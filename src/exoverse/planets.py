"""Procedural generation of physically plausible exoplanets.

Physics & data grounding
------------------------
- Multiplicity      : truncated Poisson, mean ~2.2 planets per star with
                      P < 640 d and Rp > 0.5 Re (Kepler occurrence studies,
                      e.g. Zhu+ 2018 find ~1-3 depending on cuts).
- Periods           : occurrence rises as ~P^1.5 inside 10 d, then roughly
                      flat per log-interval out to ~640 d (Petigura+ 2013).
- Radii             : bimodal super-Earth / sub-Neptune distribution with the
                      photoevaporation "radius valley" at ~1.9 Re, whose
                      center shifts to smaller radii at longer periods
                      (Fulton+ 2017, Van Eylen+ 2018 slope d logR/d logP ~ -0.1);
                      a separate giant branch whose occurrence scales with
                      host metallicity (Fischer & Valenti 2005) and stellar mass.
- Mass from radius  : probabilistic broken power laws in the spirit of
                      Chen & Kipping (2017) / Wolfgang+ (2016), with lognormal
                      scatter. Giant-branch masses are drawn nearly
                      independently of radius (degenerate-EOS plateau).
- Eccentricity      : Rayleigh(0.05) in multis; Beta(0.867, 3.03) for singles
                      (Kipping 2013). Multi-planet systems are dynamically
                      colder, matching Kepler statistics (Xie+ 2016).
- Mutual inclination: Rayleigh(1.5 deg) around the system plane (Fabrycky+ 2014).

Hard validity rules (INVALID -> resample):
- Mass above 13 Mjup (deuterium-burning limit -> brown dwarf, not a planet).
- Bulk density above the pure-iron mass-radius curve (Zeng+ 2016): no known
  formation channel makes a planet denser than solid iron.
- Bulk density below 0.03 g/cc (well below the most extreme super-puffs).
- Periastron inside the fluid Roche limit or inside the star.
- Crossing orbits or adjacent-pair separation < 2*sqrt(3) mutual Hill radii
  (Gladman 1993 two-planet criterion).

QUESTIONABLE (kept, flagged):
- Density 0.03-0.1 g/cc ("super-puff" regime, known but poorly understood).
- Inflated radius (>1.25 Rjup) without strong irradiation (Teq < 1000 K):
  observed hot-Jupiter inflation correlates with irradiation.
- Adjacent-pair Hill separation 3.46-9: Gladman-stable but long-term
  (Gyr) stability of multi-planet chains typically needs Delta >~ 9-12.
- Ultra-short periods (< 1 d) for planets above 2 Re ("hot-Neptune desert",
  Mazeh+ 2016 - strongly depleted in the real population).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

import numpy as np

from .architecture import correlated_quantile
from .constants import (
    AU, DAY, G, M_EARTH, M_JUP, M_SUN, R_EARTH, R_JUP, R_SUN,
)
from .flags import Flag, Severity
from .stars import Star


@dataclass
class Planet:
    # Bulk
    radius: float          # Earth radii
    mass: float            # Earth masses
    density: float         # g cm^-3
    composition_class: str # rocky | sub-neptune | neptunian | giant
    # Orbit
    period: float          # days
    a: float               # AU
    ecc: float
    omega_deg: float       # argument of periastron
    inc_deg: float         # inclination to the line of sight (90 = edge-on)
    t0_frac: float         # transit-phase offset in [0,1)
    # Environment
    teq: float             # K, equilibrium temperature (A=0.3, full redistribution)
    insolation: float      # Earth flux units
    in_habitable_zone: bool
    atmosphere: object = None   # exoverse.atmospheres.Atmosphere, set at system level
    flags: List[Flag] = field(default_factory=list)

    def add_flag(self, severity: Severity, rule: str, message: str) -> None:
        self.flags.append(Flag(severity, rule, message))

    @property
    def is_invalid(self) -> bool:
        return any(f.severity == Severity.INVALID for f in self.flags)


# ----------------------------------------------------------------------------
# Distributions
# ----------------------------------------------------------------------------
def sample_period(rng: np.random.Generator, p_lo=0.5, p_hi=640.0) -> float:
    """Broken power law: dN/dlnP ~ P^1.5 for P<10 d, flat beyond (Petigura+13)."""
    while True:
        p = float(np.exp(rng.uniform(np.log(p_lo), np.log(p_hi))))
        w = min((p / 10.0) ** 1.5, 1.0)
        if rng.random() < w:
            return p


def sample_small_radius(rng: np.random.Generator, period: float) -> float:
    """Bimodal super-Earth / sub-Neptune radii with a period-dependent valley."""
    shift = (period / 10.0) ** (-0.10)  # valley moves inward with period
    if rng.random() < 0.45:
        r = float(np.exp(rng.normal(np.log(1.3), 0.20)))
    else:
        r = float(np.exp(rng.normal(np.log(2.4), 0.25)))
    return float(np.clip(r * shift ** 0.5, 0.4, 4.0))


def small_radius_cdf(r: float, period: float) -> float:
    """Exact CDF of sample_small_radius (before the [0.4, 4.0] clip atoms):
    the same two-lognormal mixture with the same period-dependent shift.
    Kept in lockstep with sample_small_radius — the sigma_r copula path
    (architecture.py) relies on this being the identical marginal."""
    f = ((period / 10.0) ** (-0.10)) ** 0.5
    x = math.log(r / f)
    phi1 = 0.5 * (1.0 + math.erf((x - math.log(1.3)) / (0.20 * math.sqrt(2.0))))
    phi2 = 0.5 * (1.0 + math.erf((x - math.log(2.4)) / (0.25 * math.sqrt(2.0))))
    return 0.45 * phi1 + 0.55 * phi2


def small_radius_ppf(u: float, period: float) -> float:
    """Inverse CDF of sample_small_radius, including the clip atoms at 0.4
    and 4.0 Re (quantile mass outside the clip maps onto the boundary, which
    is exactly what np.clip does to the sampled values)."""
    if u <= small_radius_cdf(0.4, period):
        return 0.4
    if u >= small_radius_cdf(4.0, period):
        return 4.0
    lo, hi = 0.4, 4.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if small_radius_cdf(mid, period) < u:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def giant_probability(star: Star) -> float:
    """P(planet is a giant) ~ metallicity (Fischer & Valenti 05) x stellar mass."""
    p = 0.10 * 10.0 ** (1.2 * star.feh) * min(star.mass, 1.5)
    return float(np.clip(p, 0.005, 0.30))


def mass_from_radius(rng: np.random.Generator, radius: float) -> tuple[float, str]:
    """Probabilistic mass-radius relation (Earth units). Returns (mass, class)."""
    if radius < 1.6:
        # Rocky branch: M ~ R^3.7 (Zeng-like Earth composition), 25% scatter
        m = radius ** 3.7 * float(np.exp(rng.normal(0.0, 0.25)))
        return m, "rocky"
    if radius < 4.0:
        # Volatile-rich: Wolfgang+16-ish M ~ 2.7 R^1.3, larger scatter
        m = 2.7 * radius ** 1.3 * float(np.exp(rng.normal(0.0, 0.40)))
        return m, "sub-neptune"
    if radius < 8.0:
        m = 2.7 * radius ** 1.3 * float(np.exp(rng.normal(0.0, 0.35)))
        return m, "neptunian"
    # Giant branch: degenerate EOS -> radius ~ flat in mass; draw mass in
    # 0.2-10 Mjup, lognormal centered on ~0.9 Mjup
    m_jup = float(np.exp(rng.normal(np.log(0.9), 0.75)))
    m_jup = float(np.clip(m_jup, 0.15, 15.0))  # >13 Mjup triggers INVALID below
    return m_jup * (M_JUP / M_EARTH), "giant"


def sample_giant_radius(rng: np.random.Generator, teq_guess: float) -> float:
    """Giant radii in Rjup: ~1 Rjup plateau, irradiation-driven inflation."""
    base = float(np.exp(rng.normal(np.log(1.0), 0.12)))
    if teq_guess > 1000.0:
        # empirical inflation up to ~2 Rjup for the hottest Jupiters
        base *= 1.0 + 0.6 * min((teq_guess - 1000.0) / 1500.0, 1.0) * rng.random()
    return base * (R_JUP / R_EARTH)


# ----------------------------------------------------------------------------
# Physical helper quantities
# ----------------------------------------------------------------------------
def semimajor_axis_au(period_days: float, mstar_msun: float) -> float:
    p = period_days * DAY
    a3 = G * mstar_msun * M_SUN * p * p / (4.0 * math.pi ** 2)
    return a3 ** (1.0 / 3.0) / AU


def equilibrium_temp(star: Star, a_au: float, albedo: float = 0.3) -> float:
    a_m = a_au * AU
    r_m = star.radius * R_SUN
    return star.teff * math.sqrt(r_m / (2.0 * a_m)) * (1.0 - albedo) ** 0.25


def pure_iron_radius(mass_me: float) -> float:
    """Approximate pure-iron M-R curve (Zeng+ 2016): R ~ 0.78 * M^0.27 (Earth units)."""
    return 0.78 * mass_me ** 0.27


def roche_limit_au(star: Star, planet_density_cgs: float) -> float:
    """Fluid Roche limit d = 2.44 R* (rho*/rho_p)^(1/3)."""
    rho_star_cgs = star.density / 1000.0
    return 2.44 * (star.radius * R_SUN / AU) * (rho_star_cgs / planet_density_cgs) ** (1 / 3)


def habitable_zone_au(star: Star) -> tuple[float, float]:
    """Conservative HZ (runaway greenhouse / maximum greenhouse), simplified
    Kopparapu+ 2013: S_inner=1.11, S_outer=0.36 with a mild Teff correction."""
    ts = star.teff - 5780.0
    s_inner = 1.11 + 1.4e-4 * ts
    s_outer = 0.36 + 5.0e-5 * ts
    return math.sqrt(star.luminosity / s_inner), math.sqrt(star.luminosity / s_outer)


# ----------------------------------------------------------------------------
# Single-planet generation (bulk + orbit around a given star)
# ----------------------------------------------------------------------------
def generate_planet(rng: np.random.Generator, star: Star, period: float,
                    is_single: bool,
                    radius_latent: tuple[float, float] | None = None) -> Planet:
    """radius_latent=(z_sys, sigma_r) switches the small-planet radius draw
    to the marginal-preserving intra-system copula (architecture.py); None
    keeps the historical independent draw bit-for-bit."""
    a_au = semimajor_axis_au(period, star.mass)
    teq_guess = equilibrium_temp(star, a_au)

    if rng.random() < giant_probability(star):
        radius = sample_giant_radius(rng, teq_guess)
        mass, comp = mass_from_radius(rng, radius)
    else:
        if radius_latent is None:
            radius = sample_small_radius(rng, period)
        else:
            u = correlated_quantile(rng, radius_latent[0], radius_latent[1])
            radius = small_radius_ppf(u, period)
        mass, comp = mass_from_radius(rng, radius)

    rho = (mass * M_EARTH) / (4.0 / 3.0 * math.pi * (radius * R_EARTH) ** 3) / 1000.0  # g/cc

    if is_single:
        ecc = float(rng.beta(0.867, 3.03))          # Kipping 2013
    else:
        ecc = float(min(rng.rayleigh(0.05), 0.7))   # dynamically cold multis
    omega = float(rng.uniform(0.0, 360.0))

    hz_in, hz_out = habitable_zone_au(star)
    insol = star.luminosity / a_au ** 2

    planet = Planet(
        radius=radius, mass=mass, density=rho, composition_class=comp,
        period=period, a=a_au, ecc=ecc, omega_deg=omega,
        inc_deg=90.0,  # set at system level
        t0_frac=float(rng.random()),
        teq=equilibrium_temp(star, a_au), insolation=insol,
        in_habitable_zone=hz_in <= a_au <= hz_out,
    )

    # --- Hard physics checks ------------------------------------------------
    if mass * M_EARTH > 13.0 * M_JUP:
        planet.add_flag(Severity.INVALID, "mass.deuterium_burning",
                        f"M={mass/317.8:.1f} Mjup exceeds 13 Mjup: object would "
                        "be a brown dwarf, not a planet")
    if radius < pure_iron_radius(mass) and comp != "giant":
        planet.add_flag(Severity.INVALID, "density.exceeds_pure_iron",
                        f"R={radius:.2f} Re at M={mass:.1f} Me is denser than a "
                        "pure-iron sphere (Zeng+16 limit); no known formation "
                        "channel produces this")
    if rho < 0.03:
        planet.add_flag(Severity.INVALID, "density.below_superpuff_floor",
                        f"rho={rho:.3f} g/cc is below even extreme super-puffs")
    peri_au = a_au * (1.0 - ecc)
    r_roche = roche_limit_au(star, max(rho, 0.03))
    if peri_au < r_roche:
        planet.add_flag(Severity.INVALID, "orbit.inside_roche_limit",
                        f"Periastron {peri_au:.4f} AU inside fluid Roche limit "
                        f"{r_roche:.4f} AU: planet would be tidally shredded")
    if peri_au * AU < 1.5 * star.radius * R_SUN:
        planet.add_flag(Severity.INVALID, "orbit.grazes_star",
                        "Periastron within 1.5 stellar radii")

    # --- Plausibility flags ---------------------------------------------------
    if 0.03 <= rho < 0.1:
        planet.add_flag(Severity.QUESTIONABLE, "density.super_puff",
                        f"rho={rho:.3f} g/cc: 'super-puff' regime (cf. Kepler-51); "
                        "real but poorly understood, may imply rings/high haze")
    if comp == "giant" and radius > 1.25 * (R_JUP / R_EARTH) and planet.teq < 1000.0:
        planet.add_flag(Severity.QUESTIONABLE, "radius.cold_inflated_giant",
                        f"R={radius/(R_JUP/R_EARTH):.2f} Rjup at Teq={planet.teq:.0f} K: "
                        "radius inflation without strong irradiation lacks a "
                        "known mechanism")
    if period < 1.0 and radius > 2.0:
        planet.add_flag(Severity.QUESTIONABLE, "orbit.hot_neptune_desert",
                        f"P={period:.2f} d with R={radius:.1f} Re sits in the "
                        "'hot-Neptune desert' (Mazeh+16), strongly depleted in "
                        "the observed population")
    if planet.teq > 2500.0:
        planet.add_flag(Severity.QUESTIONABLE, "environment.ultra_hot",
                        f"Teq={planet.teq:.0f} K: ultra-hot regime, atmosphere "
                        "models weakly constrained")
    if comp == "giant" and star.mass < 0.3:
        planet.add_flag(Severity.QUESTIONABLE, "formation.giant_around_late_m_dwarf",
                        f"Gas giant around a {star.mass:.2f} Msun M dwarf: real "
                        "analogs exist (GJ 3512b) but core accretion struggles "
                        "to form them; occurrence is very low")
    if comp == "rocky" and mass > 10.0:
        planet.add_flag(Severity.QUESTIONABLE, "mass.mega_earth",
                        f"Rocky planet of {mass:.0f} Me ('mega-Earth', cf. "
                        "Kepler-10c debate): rare and formation is debated")
    return planet
