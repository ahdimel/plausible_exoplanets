"""Assembly of full stellar systems: star + dynamically consistent planets.

Stability rules applied to adjacent planet pairs:
- INVALID  : radial excursions overlap (a_out*(1-e_out) < a_in*(1+e_in)) ->
             orbit crossing; or Hill separation Delta < 2*sqrt(3) (Gladman
             1993 analytic instability for two planets).
- QUESTIONABLE : Delta in [2*sqrt(3), 9): Gladman-stable pair, but N-body
             studies of multi-planet chains (e.g. Pu & Wu 2015) show Gyr
             survival typically requires Delta >~ 9-12.

Geometry: a single system plane orientation is drawn isotropically
(cos i uniform), and each planet receives a small mutual inclination
(Rayleigh, sigma=1.5 deg, Fabrycky+ 2014) about that plane.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

import numpy as np

from .constants import M_EARTH, M_SUN
from .flags import Flag, Severity
from .planets import Planet, generate_planet, sample_period
from .stars import Star, generate_star

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
    flags: List[Flag] = field(default_factory=list)

    def add_flag(self, severity: Severity, rule: str, message: str) -> None:
        self.flags.append(Flag(severity, rule, message))


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


def generate_system(seed: int, name: str, max_planet_tries: int = 40) -> StellarSystem:
    """Generate one validated stellar system. INVALID draws are resampled."""
    rng = np.random.default_rng(seed)

    star = generate_star(rng)
    while star.is_invalid:
        star = generate_star(rng)

    n_target = min(int(rng.poisson(MEAN_PLANETS)), MAX_PLANETS)
    is_single = n_target == 1

    planets: List[Planet] = []
    tries = 0
    while len(planets) < n_target and tries < max_planet_tries:
        tries += 1
        period = sample_period(rng)
        # Enforce minimum period ratio 1.2 against existing planets
        if any(max(period, q.period) / min(period, q.period) < 1.2 for q in planets):
            continue
        cand = generate_planet(rng, star, period, is_single)
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
    for p in planets:
        mut = rng.rayleigh(1.5)
        sign = 1.0 if rng.random() < 0.5 else -1.0
        # i > 90 deg is geometrically equivalent (mirrored transit chord);
        # clipping at 90 would create an artificial pileup at b = 0
        p.inc_deg = float(np.clip(sys_inc + sign * mut, 0.0, 180.0))

    system = StellarSystem(name=name, seed=seed, star=star, planets=planets,
                           sys_inc_deg=sys_inc)
    if n_target > 0 and len(planets) < n_target:
        system.add_flag(Severity.INFO, "system.multiplicity_reduced",
                        f"Targeted {n_target} planets but only {len(planets)} "
                        "satisfied stability constraints")
    if len(planets) >= 5:
        system.add_flag(Severity.QUESTIONABLE, "system.high_multiplicity",
                        f"{len(planets)}-planet system: real analogs exist "
                        "(Kepler-90, TRAPPIST-1) but are rare; long-term "
                        "stability not verified by N-body integration")
    return system
