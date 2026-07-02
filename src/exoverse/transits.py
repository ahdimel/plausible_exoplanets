"""Transit geometry and limb-darkened light-curve model.

Light-curve model: exact per-annulus occultation geometry with quadratic limb
darkening I(mu) = 1 - u1(1-mu) - u2(1-mu)^2, integrated numerically over ~400
radial annuli of the stellar disk. This reproduces the Mandel & Agol (2002)
quadratic model to ~1 ppm without external dependencies.

Geometry follows Winn (2010, "Transits and Occultations") throughout,
including the eccentricity correction factor sqrt(1-e^2)/(1+e sin w) applied
to the star-planet separation at inferior conjunction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import AU, R_EARTH, R_SUN
from .planets import Planet
from .stars import Star


def transit_flux(z: np.ndarray, k: float, u1: float, u2: float,
                 n_annuli: int = 300) -> np.ndarray:
    """Relative flux for planet/star radius ratio k at sky-projected
    separations z (in stellar radii). Quadratic limb darkening (u1, u2).

    Integrates the occulted intensity over annuli spanning only the radial
    band [z-k, z+k] actually covered by the planet, so accuracy is
    independent of k (better than ~1e-4 relative for n_annuli=300)."""
    z = np.atleast_1d(np.abs(np.asarray(z, dtype=float)))
    total = math.pi * (1.0 - u1 / 3.0 - u2 / 6.0)  # exact disk-integrated flux

    # Per-z radial grid across the covered band, clipped to the stellar disk
    lo = np.clip(z - k, 0.0, 1.0)[:, None]
    hi = np.clip(z + k, 0.0, 1.0)[:, None]
    frac = ((np.arange(n_annuli) + 0.5) / n_annuli)[None, :]
    rr = lo + (hi - lo) * frac
    dr = (hi - lo) / n_annuli

    mu = np.sqrt(np.clip(1.0 - rr * rr, 0.0, 1.0))
    inten = 1.0 - u1 * (1.0 - mu) - u2 * (1.0 - mu) ** 2

    # Half-angle of the arc of each annulus covered by the planet disk
    zz = z[:, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_phi = (zz * zz + rr * rr - k * k) / (2.0 * zz * rr)
    phi = np.arccos(np.clip(cos_phi, -1.0, 1.0))
    phi = np.where(rr <= np.maximum(k - zz, 0.0), np.pi, phi)   # annulus fully inside
    phi = np.where(zz >= rr + k, 0.0, phi)                       # no overlap

    blocked = np.sum(inten * 2.0 * phi * rr * dr, axis=1)
    return 1.0 - blocked / total


@dataclass
class TransitGeometry:
    transits: bool
    b: float                 # impact parameter (stellar radii)
    depth_ppm: float         # limb-darkened depth at mid-transit
    depth_uniform_ppm: float # geometric (Rp/R*)^2
    t14_hours: float         # total duration (1st-4th contact)
    t23_hours: float         # full duration (2nd-3rd contact); 0 if grazing
    prob_transit: float      # geometric a-priori transit probability
    k: float                 # radius ratio


def compute_geometry(star: Star, planet: Planet) -> TransitGeometry:
    a_over_rs = planet.a * AU / (star.radius * R_SUN)
    k = planet.radius * R_EARTH / (star.radius * R_SUN)
    inc = math.radians(planet.inc_deg)
    e, w = planet.ecc, math.radians(planet.omega_deg)

    # Separation scale factor at inferior conjunction (Winn 2010 eq. 7)
    ecc_fac = (1.0 - e * e) / (1.0 + e * math.sin(w))
    b = a_over_rs * math.cos(inc) * ecc_fac

    prob = min(1.0, (1.0 + k) / (a_over_rs * ecc_fac))
    transits = abs(b) < 1.0 + k
    if not transits:
        return TransitGeometry(False, b, 0.0, 0.0, 0.0, 0.0, prob, k)

    depth_uniform = k * k * 1e6
    depth_ld = float((1.0 - transit_flux(np.array([abs(b)]), k, star.u1, star.u2)[0]) * 1e6)

    # Durations (Winn 2010 eqs. 14-16 with eccentricity factor)
    vel_fac = math.sqrt(1.0 - e * e) / (1.0 + e * math.sin(w))
    def duration(contact_k: float) -> float:
        arg = (1.0 + contact_k) ** 2 - b * b
        if arg <= 0.0 or a_over_rs * math.sin(inc) == 0.0:
            return 0.0
        x = math.sqrt(arg) / (a_over_rs * math.sin(inc))
        if x >= 1.0:
            return 0.0
        return planet.period / math.pi * math.asin(x) * vel_fac * 24.0

    t14 = duration(k)
    arg23 = (1.0 - k) ** 2 - b * b
    t23 = 0.0
    if arg23 > 0.0:
        x = math.sqrt(arg23) / (a_over_rs * math.sin(inc))
        t23 = planet.period / math.pi * math.asin(min(x, 1.0)) * vel_fac * 24.0

    return TransitGeometry(True, b, depth_ld, depth_uniform, t14, t23, prob, k)


def model_light_curve(star: Star, geom: TransitGeometry,
                      n_points: int = 300, window_factor: float = 1.6):
    """Model light curve around one transit. Returns (t_hours, flux) with t=0
    at mid-transit. Uses a constant-velocity chord approximation across the
    transit window (excellent for T14 << P)."""
    if not geom.transits or geom.t14_hours <= 0.0:
        t = np.linspace(-2.0, 2.0, n_points)
        return t, np.ones_like(t)
    half = geom.t14_hours * window_factor / 2.0
    t = np.linspace(-half, half, n_points)
    # Sky-projected separation: planet crosses chord at impact parameter b
    x = t / (geom.t14_hours / 2.0) * math.sqrt(max((1.0 + geom.k) ** 2 - geom.b ** 2, 0.0))
    z = np.sqrt(x * x + geom.b ** 2)
    return t, transit_flux(z, geom.k, star.u1, star.u2)
