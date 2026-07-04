"""Assembly of full stellar systems: star + dynamically consistent planets.

Stability rules applied to adjacent planet pairs:
- INVALID  : radial excursions overlap (a_out*(1-e_out) < a_in*(1+e_in)) ->
             orbit crossing; or Hill separation Delta < 2*sqrt(3) (Gladman
             1993 analytic instability for two planets).
- QUESTIONABLE : Delta in [2*sqrt(3), 9): Gladman-stable pair, but N-body
             studies of multi-planet chains (e.g. Pu & Wu 2015) show Gyr
             survival typically requires Delta >~ 9-12.

Geometry: a single system plane orientation is drawn isotropically
(cos i uniform), and each planet's orbit normal is tilted away from the
system-plane normal by a mutual inclination (Rayleigh, sigma_i from the
Architecture config; default 1.5 deg, Fabrycky+ 2014) at a uniformly
random nodal azimuth. The azimuth matters:
only the component of the tilt along the line of sight changes the impact
parameter, so a Rayleigh(sigma) tilt with random node projects to a
Normal(0, sigma) inclination offset. (Applying the full tilt with a random
sign, as this module did pre-v0.4, inflates the projected scatter by
sqrt(2) and distorts transit multiplicity — the observable this project
studies.)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

import numpy as np

from .architecture import DEFAULT_ARCH, Architecture
from .atmospheres import assign_atmosphere
from .constants import M_EARTH, M_SUN
from .flags import Flag, Severity
from .planets import Planet, generate_planet, sample_period
from .stars import Star, generate_star
from .stellar_noise import StellarNoise, generate_stellar_noise

GLADMAN_DELTA = 2.0 * math.sqrt(3.0)
LONGTERM_DELTA = 9.0
MAX_PLANETS = 7
MEAN_PLANETS = 2.2


@dataclass
class StellarSystem:
    name: str
    seed: int
    star: Star
    planets: List[Planet]           # sorted by period
    sys_inc_deg: float              # inclination of the system plane
    noise: StellarNoise = None      # astrophysical noise state of the star
    flags: List[Flag] = field(default_factory=list)

    def add_flag(self, severity: Severity, rule: str, message: str) -> None:
        self.flags.append(Flag(severity, rule, message))


def tilted_inclination_deg(sys_inc_deg: float, mut_deg: float,
                           node_rad: float) -> float:
    """Line-of-sight inclination of an orbit tilted mut_deg from the system
    plane at nodal azimuth node_rad about the plane normal (spherical law of
    cosines between orbit normal and line of sight). Result lies in
    [0, 180]; i > 90 is the mirror-equivalent transit chord, and clipping it
    to 90 would create an artificial b = 0 pileup."""
    i_s = math.radians(sys_inc_deg)
    mu = math.radians(mut_deg)
    cos_ip = (math.cos(i_s) * math.cos(mu)
              + math.sin(i_s) * math.sin(mu) * math.cos(node_rad))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_ip))))


def mutual_hill_delta(star: Star, p_in: Planet, p_out: Planet) -> float:
    """Separation of adjacent planets in mutual Hill radii."""
    m_sum = (p_in.mass + p_out.mass) * M_EARTH / (star.mass * M_SUN)
    r_hill = ((m_sum / 3.0) ** (1.0 / 3.0)) * 0.5 * (p_in.a + p_out.a)
    return (p_out.a - p_in.a) / r_hill


def check_pair_stability(star: Star, p_in: Planet, p_out: Planet) -> Flag | None:
    peri_out = p_out.a * (1.0 - p_out.ecc)
    apo_in = p_in.a * (1.0 + p_in.ecc)
    if peri_out <= apo_in:
        return Flag(Severity.INVALID, "stability.orbit_crossing",
                    f"Orbits of P={p_in.period:.1f} d and P={p_out.period:.1f} d "
                    "planets cross (apoapsis exceeds neighbor periastron)")
    delta = mutual_hill_delta(star, p_in, p_out)
    if delta < GLADMAN_DELTA:
        return Flag(Severity.INVALID, "stability.hill_unstable",
                    f"Adjacent pair separated by {delta:.1f} mutual Hill radii "
                    f"< 2*sqrt(3): Gladman-unstable")
    if delta < LONGTERM_DELTA:
        return Flag(Severity.QUESTIONABLE, "stability.tightly_packed",
                    f"Adjacent pair at {delta:.1f} mutual Hill radii: formally "
                    "stable but Gyr-timescale survival of packed multis "
                    "typically needs >~9 (Pu & Wu 2015)")
    return None


def _draw_planets(rng: np.random.Generator, star: Star, arch: Architecture,
                  max_planet_tries: int = 40,
                  sys_flags: List[Flag] | None = None
                  ) -> tuple[List[Planet], float]:
    """Draw the planet set and line-of-sight inclinations for one star.

    This is the exact planet-generation loop + inclination block that
    historically lived inline in generate_system, extracted so
    kepler_field.py can run it against DR25-conditioned stars while
    skipping the (expensive, detection-irrelevant) noise and atmosphere
    draws. The rng draw order is byte-identical to the pre-refactor code —
    generate_system remains bit-for-bit reproducible (guarded by
    tests/test_architecture.py::test_default_arch_bit_for_bit). In
    particular the radius_latent draw happens only when arch.sigma_r is
    not None and the f_hot coin only when arch.f_hot > 0, exactly as
    before.

    Returns (planets sorted by period, system-plane inclination in deg).
    System-level flags (multiplicity_reduced) are appended to sys_flags
    when a list is provided."""
    n_target = min(int(rng.poisson(MEAN_PLANETS)), MAX_PLANETS)
    is_single = n_target == 1

    # Latent intra-system radius scale (copula, architecture.py). Drawn only
    # when the knob is on so default populations keep the historical stream.
    radius_latent = None
    if arch.sigma_r is not None:
        radius_latent = (float(rng.standard_normal()), arch.sigma_r)

    planets: List[Planet] = []
    tries = 0
    while len(planets) < n_target and tries < max_planet_tries:
        tries += 1
        period = sample_period(rng)
        # Enforce minimum period ratio 1.2 against existing planets
        if any(max(period, q.period) / min(period, q.period) < 1.2 for q in planets):
            continue
        cand = generate_planet(rng, star, period, is_single, radius_latent)
        if cand.is_invalid:
            continue
        trial = sorted(planets + [cand], key=lambda p: p.period)
        idx = trial.index(cand)
        ok = True
        pair_flags: List[Flag] = []
        for j in (idx - 1, idx):
            if 0 <= j < len(trial) - 1:
                f = check_pair_stability(star, trial[j], trial[j + 1])
                if f is not None and f.severity == Severity.INVALID:
                    ok = False
                    break
                if f is not None:
                    pair_flags.append(f)
        if not ok:
            continue
        cand.flags.extend(pair_flags)
        planets = trial

    sys_inc = math.degrees(math.acos(rng.uniform(0.0, 1.0)))
    if arch.isotropic:
        # Exact independent-isotropic limit: no shared plane at all
        for p in planets:
            p.inc_deg = math.degrees(math.acos(rng.uniform(0.0, 1.0)))
    else:
        sigma_i = arch.sigma_i
        if arch.f_hot > 0.0 and rng.random() < arch.f_hot:
            sigma_i = arch.sigma_i_hot
        for p in planets:
            mut = rng.rayleigh(sigma_i)
            node = rng.uniform(0.0, 2.0 * math.pi)
            p.inc_deg = tilted_inclination_deg(sys_inc, mut, node)

    if sys_flags is not None and n_target > 0 and len(planets) < n_target:
        sys_flags.append(Flag(
            Severity.INFO, "system.multiplicity_reduced",
            f"Targeted {n_target} planets but only {len(planets)} "
            "satisfied stability constraints"))
    return planets, sys_inc


def generate_system(seed: int, name: str, max_planet_tries: int = 40,
                    dmax_pc: float = 300.0,
                    arch: Architecture | None = None) -> StellarSystem:
    """Generate one validated stellar system. INVALID draws are resampled.

    dmax_pc caps the host distance draw (see stars.generate_star): it only
    rescales the distance, never the random stream, so (seed, name, dmax_pc)
    is bit-for-bit reproducible and dmax_pc=300 matches historical worlds.
    arch (architecture.Architecture) sets the dichotomy-study knobs; the
    default consumes the identical rng sequence as arch=None, so default
    worlds match the historical baseline bit-for-bit. Re-generation paths
    (CLI lightcurve, web UI) must pass the dmax_pc AND arch the population
    was generated with (both stored in DB meta)."""
    arch = DEFAULT_ARCH if arch is None else arch
    rng = np.random.default_rng(seed)

    star = generate_star(rng, dmax_pc)
    while star.is_invalid:
        star = generate_star(rng, dmax_pc)

    sys_flags: List[Flag] = []
    planets, sys_inc = _draw_planets(rng, star, arch, max_planet_tries,
                                     sys_flags)

    # Astrophysical noise state and atmospheres are drawn here (not in the
    # orchestrator) so a system is bit-for-bit reproducible from (seed, name)
    noise = generate_stellar_noise(rng, star)
    for p in planets:
        p.atmosphere = assign_atmosphere(rng, star, p)

    system = StellarSystem(name=name, seed=seed, star=star, planets=planets,
                           sys_inc_deg=sys_inc, noise=noise,
                           flags=sys_flags)
    if len(planets) >= 5:
        system.add_flag(Severity.QUESTIONABLE, "system.high_multiplicity",
                        f"{len(planets)}-planet system: real analogs exist "
                        "(Kepler-90, TRAPPIST-1) but are rare; long-term "
                        "stability not verified by N-body integration")
    return system
