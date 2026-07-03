"""Procedural generation of physically plausible main-sequence stars.

Physics grounding
-----------------
- Mass function     : Kroupa (2001) broken power-law IMF, restricted to
                      0.08-2.2 Msun (transit surveys overwhelmingly target
                      FGKM dwarfs; massive stars are rare, short-lived, and
                      poor transit hosts -> excluded, recorded as INFO flag).
- L(M), R(M)        : standard main-sequence power-law scaling relations
                      (e.g. Kippenhahn & Weigert), with lognormal scatter
                      representing age/metallicity spread.
- Teff              : from Stefan-Boltzmann, L = 4 pi R^2 sigma Teff^4.
- Lifetime          : t_MS ~ 10 Gyr * (M/Msun) / (L/Lsun).
- Metallicity       : [Fe/H] ~ N(-0.1, 0.2), clipped to [-1.0, +0.5]
                      (local thin-disk distribution, e.g. Hayden+ 2015).
- Distance          : uniform stellar density -> p(d) ~ d^2, out to 300 pc.
- Bolometric corr.  : Flower (1996) polynomials as corrected by Torres (2010).
- Colors            : coarse Teff -> (V-Ic), (V-J) lookup (Pecaut & Mamajek
                      2013 dwarf sequence, interpolated). QUESTIONABLE-free
                      but explicitly approximate (INFO flag on the star).
- Limb darkening    : quadratic (u1, u2) vs Teff, coarse fit to Claret &
                      Bloemen (2011) Kepler-band tabulations.

Known simplifications (each recorded as a flag where relevant):
- Single stars only; ~45% of FGK stars are in binaries/multiples.
- Main-sequence only; no subgiant/giant evolution. Stars drawn with an age
  in the last 10% of their MS lifetime are flagged QUESTIONABLE because
  radius inflation near turnoff is not modeled.
- Scaling-relation scatter is a proxy for real isochrone physics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .constants import (
    AU, GYR, L_SUN, M_SUN, MBOL_SUN, PC, R_SUN, RHO_SUN, SIGMA_SB, TEFF_SUN,
)
from .flags import Flag, Severity


# ----------------------------------------------------------------------------
# Empirical lookup tables (dwarf sequence, Pecaut & Mamajek 2013, coarse)
# ----------------------------------------------------------------------------
_TEFF_GRID = np.array([2800, 3000, 3500, 4000, 4500, 5000, 5500, 5772, 6000, 6500, 7000, 7500])
_V_MINUS_IC = np.array([4.20, 3.10, 2.06, 1.44, 1.06, 0.85, 0.72, 0.66, 0.62, 0.50, 0.42, 0.36])
_V_MINUS_J = np.array([5.80, 4.80, 3.72, 2.60, 2.00, 1.60, 1.32, 1.20, 1.10, 0.92, 0.78, 0.66])
# Quadratic limb darkening, Kepler/TESS-like broad optical band, ~Claret & Bloemen 2011
_LD_U1 = np.array([0.62, 0.60, 0.56, 0.55, 0.52, 0.48, 0.44, 0.41, 0.38, 0.33, 0.28, 0.24])
_LD_U2 = np.array([0.10, 0.12, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.27, 0.29, 0.31, 0.32])


def _interp_teff(teff: float, table: np.ndarray) -> float:
    return float(np.interp(teff, _TEFF_GRID, table))


def bolometric_correction_v(teff: float) -> float:
    """BC_V from Flower (1996) polynomials, coefficients per Torres (2010)."""
    lt = math.log10(teff)
    if lt < 3.70:
        c = [-0.190537291496456e5, 0.155144866764412e5, -0.421278819301717e4,
             0.381476328422343e3]
    elif lt < 3.90:
        c = [-0.370510203809015e5, 0.385672629965804e5, -0.150651486316025e5,
             0.261724637119416e4, -0.170623810323864e3]
    else:
        c = [-0.118115450538963e6, 0.137145973583929e6, -0.636233812100225e5,
             0.147412923562646e5, -0.170587278406872e4, 0.788731721804990e2]
    return float(sum(ci * lt ** i for i, ci in enumerate(c)))


# ----------------------------------------------------------------------------
# Star model
# ----------------------------------------------------------------------------
@dataclass
class Star:
    mass: float            # Msun
    radius: float          # Rsun
    luminosity: float      # Lsun
    teff: float            # K
    feh: float             # [Fe/H] dex
    age_gyr: float
    ms_lifetime_gyr: float
    distance_pc: float
    # Observables
    mag_v: float
    mag_tess: float        # ~ Cousins Ic proxy
    mag_j: float
    # Transit-modeling inputs
    u1: float              # quadratic limb-darkening coeff
    u2: float
    density: float         # kg m^-3 (mean)
    spectral_type: str
    flags: List[Flag] = field(default_factory=list)

    def add_flag(self, severity: Severity, rule: str, message: str) -> None:
        self.flags.append(Flag(severity, rule, message))

    @property
    def is_invalid(self) -> bool:
        return any(f.severity == Severity.INVALID for f in self.flags)


def spectral_type_from_teff(teff: float) -> str:
    """Coarse dwarf spectral type from Teff (Pecaut & Mamajek 2013 boundaries)."""
    bounds = [(30000, "B"), (7300, "A"), (6000, "F"), (5300, "G"), (3900, "K"), (2300, "M")]
    letters = [(7300, 30000, "A"), (6000, 7300, "F"), (5300, 6000, "G"),
               (3900, 5300, "K"), (2300, 3900, "M")]
    for lo, hi, letter in letters:
        if lo <= teff < hi:
            # subclass 0-9 linearly across the class
            frac = (hi - teff) / (hi - lo)
            sub = min(9, int(frac * 10))
            return f"{letter}{sub}V"
    return "??"


def sample_kroupa_mass(rng: np.random.Generator, m_lo: float = 0.08,
                       m_hi: float = 2.2) -> float:
    """Sample from the Kroupa (2001) IMF: dN/dm ~ m^-1.3 (m<0.5), m^-2.3 (m>=0.5)."""
    mb = 0.5
    a1, a2 = 1.3, 2.3
    # Segment integrals of m^-alpha dm (continuous at the break with weight mb^(a2-a1))
    def seg(lo, hi, a):
        p = 1.0 - a
        return (hi ** p - lo ** p) / p
    w1 = seg(m_lo, mb, a1)
    w2 = mb ** (a2 - a1) * seg(mb, m_hi, a2)
    if rng.random() < w1 / (w1 + w2):
        u, p = rng.random(), 1.0 - a1
        return float((m_lo ** p + u * (mb ** p - m_lo ** p)) ** (1 / p))
    u, p = rng.random(), 1.0 - a2
    return float((mb ** p + u * (m_hi ** p - mb ** p)) ** (1 / p))


def ms_luminosity(mass: float) -> float:
    """Main-sequence L/Lsun from piecewise mass-luminosity relation."""
    if mass < 0.43:
        return 0.23 * mass ** 2.3
    if mass < 2.0:
        return mass ** 4.0
    return 1.4 * mass ** 3.5


def ms_radius(mass: float) -> float:
    """Main-sequence R/Rsun (ZAMS-ish power laws)."""
    if mass < 1.0:
        return mass ** 0.9
    return mass ** 0.6


def generate_star(rng: np.random.Generator, dmax_pc: float = 300.0) -> Star:
    """Generate one plausible main-sequence star with validated metadata.

    dmax_pc caps the distance draw (p(d) ~ d^2 out to dmax_pc). It scales a
    single rng.random() call, so the random stream — and therefore every
    other property of the star and its system — is identical for any value;
    only the distance (and apparent magnitudes) change. Populations with
    dmax_pc < 300 model the solar neighborhood for direct-imaging studies."""
    mass = sample_kroupa_mass(rng)
    feh = float(np.clip(rng.normal(-0.1, 0.2), -1.0, 0.5))

    # Scaling relations with modest lognormal scatter (age/metallicity proxy)
    lum = ms_luminosity(mass) * float(np.exp(rng.normal(0.0, 0.10)))
    radius = ms_radius(mass) * float(np.exp(rng.normal(0.0, 0.05)))
    teff = TEFF_SUN * (lum / radius ** 2) ** 0.25

    ms_lifetime = 10.0 * mass / lum  # Gyr
    age = float(rng.uniform(0.5, min(12.0, ms_lifetime)))

    dist = dmax_pc * float(rng.random() ** (1.0 / 3.0))  # p(d) ~ d^2
    dist = max(dist, 5.0)

    mbol = MBOL_SUN - 2.5 * math.log10(lum)
    bc_v = bolometric_correction_v(teff)
    abs_v = mbol - bc_v
    mu = 5.0 * math.log10(dist / 10.0)
    mag_v = abs_v + mu
    mag_tess = mag_v - _interp_teff(teff, _V_MINUS_IC)   # TESS band ~ Cousins Ic
    mag_j = mag_v - _interp_teff(teff, _V_MINUS_J)

    density = RHO_SUN * mass / radius ** 3

    star = Star(
        mass=mass, radius=radius, luminosity=lum, teff=teff, feh=feh,
        age_gyr=age, ms_lifetime_gyr=ms_lifetime, distance_pc=dist,
        mag_v=mag_v, mag_tess=mag_tess, mag_j=mag_j,
        u1=_interp_teff(teff, _LD_U1), u2=_interp_teff(teff, _LD_U2),
        density=density, spectral_type=spectral_type_from_teff(teff),
    )

    # --- Validation & plausibility flags -----------------------------------
    if teff < 2300 or teff > 8500:
        star.add_flag(Severity.INVALID, "star.teff_out_of_range",
                      f"Teff={teff:.0f} K outside modeled dwarf range 2300-8500 K")
    if age > 0.9 * ms_lifetime:
        star.add_flag(Severity.QUESTIONABLE, "star.near_turnoff",
                      f"Age {age:.1f} Gyr is >90% of MS lifetime "
                      f"({ms_lifetime:.1f} Gyr); turnoff radius inflation is "
                      "not modeled, radius may be underestimated")
    if feh < -0.6:
        star.add_flag(Severity.QUESTIONABLE, "star.metal_poor",
                      f"[Fe/H]={feh:.2f}: scaling relations are calibrated on "
                      "near-solar-metallicity stars; R and L less reliable")
    star.add_flag(Severity.INFO, "star.single_only",
                  "Generator models single stars only (~45% of FGK stars are "
                  "actually in multiple systems)")
    star.add_flag(Severity.INFO, "star.colors_approximate",
                  "TESS/J magnitudes derived from coarse Teff-color tables "
                  "(Pecaut & Mamajek 2013); +/-0.1 mag level accuracy")
    return star
