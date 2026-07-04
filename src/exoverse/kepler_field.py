"""Kepler-field conditioned simulation for the dichotomy study (Phase 2).

One synthetic "universe" = one pass over the DR25-cut target list: every
target hosts exactly one generator draw (system.py planet loop + inclination
block via system._draw_planets), and detection is evaluated against that
target's own empirical noise and observing window.

Noise model
-----------
Detection uses the DR25 per-target RMS CDPP (Christiansen+ 2012; the
`rrmscdpp01p5 ... rrmscdpp15p0` columns of the DR25 stellar table). CDPP is
an *empirical, total* noise measure — photon noise, instrument systematics,
and stellar variability of that actual star, combined — so the exoverse
analytic stellar-noise model (stellar_noise.py) is deliberately NOT added on
top: doing so would double-count stellar variability. Atmosphere draws are
also skipped (irrelevant to detection; large speedup).

Detection rule (pipeline-style, threshold-crossing)
---------------------------------------------------
- Geometric transit via transits.compute_geometry (Winn 2010 geometry with
  limb-darkened depth).
- n_tr = dataspan * dutycycle / period, require n_tr >= MIN_TRANSITS (3, the
  Kepler pipeline's minimum for a TCE; Jenkins+ 2010).
- sigma = CDPP interpolated to the transit duration t14 over
  CDPP_DURATIONS_HR (clamped at the ends).
- SNR = depth_ppm / sigma * sqrt(n_tr); detect if >= SNR_THRESHOLD (7.1,
  the Kepler transiting-planet-search threshold; Jenkins+ 2010).
This is a sharp threshold, not the MES-ramp detection-efficiency curve —
the ramp (and per-duration mesthres, window-function probabilities) are
pre-registered Phase 4 robustness variants (docs/robustness_plan.md).

Targets are duck-typed against kepler_data.KeplerTarget (kepid, teff, feh,
radius, mass, kepmag, dutycycle, dataspan, cdpp_ppm); kepler_data is never
imported here, so this module (and its tests) work without the DR25
snapshots on disk.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .architecture import DEFAULT_ARCH, Architecture
from .constants import RHO_SUN, TEFF_SUN
from .planets import Planet
from .stars import _LD_U1, _LD_U2, Star, _interp_teff, spectral_type_from_teff
from .system import _draw_planets
from .transits import compute_geometry

# Transit durations (hours) of the 14 DR25 rrmscdpp columns, in order. Must
# match kepler_data.KeplerTarget.CDPP_DURATIONS_HR (kept local so this module
# has no import-time dependency on the ingestion module / data snapshots).
CDPP_DURATIONS_HR: Tuple[float, ...] = (
    1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 6.0, 7.5, 9.0, 10.5, 12.0, 12.5, 15.0)

SNR_THRESHOLD = 7.1     # Kepler TPS multiple-event statistic threshold
MIN_TRANSITS = 3        # minimum transit count for a TCE

# Multiplicity histogram keys: 1..5 exact, everything above pooled as "6+"
N_K_KEYS: Tuple[object, ...] = (1, 2, 3, 4, 5, "6+")


def star_from_target(target) -> Star:
    """Build an exoverse Star from a DR25 stellar-table target.

    mass/radius/teff/feh are taken directly from the catalog;
    L/Lsun = (R/Rsun)^2 (Teff/5772 K)^4 (Stefan-Boltzmann); quadratic limb
    darkening comes from the same Teff lookup generate_star uses
    (stars._interp_teff over the Claret & Bloemen 2011-like tables); mean
    density = rho_sun * (M/Msun) / (R/Rsun)^3 in kg m^-3 as the Star
    dataclass expects (transit-duration consistency).

    distance_pc and the apparent magnitudes are filled with sane
    placeholders (100 pc; kepmag for all three bands): detection here runs
    entirely off the target's own CDPP, so neither distance nor brightness
    enters — they exist only to satisfy the dataclass. age is set to half
    the (uncomputed) MS lifetime placeholder of 5 Gyr; equally unused."""
    teff = float(target.teff)
    lum = float(target.radius) ** 2 * (teff / TEFF_SUN) ** 4
    kepmag = float(getattr(target, "kepmag", 12.0))
    return Star(
        mass=float(target.mass), radius=float(target.radius),
        luminosity=lum, teff=teff, feh=float(target.feh),
        age_gyr=5.0, ms_lifetime_gyr=10.0, distance_pc=100.0,
        mag_v=kepmag, mag_tess=kepmag, mag_j=kepmag,
        u1=_interp_teff(teff, _LD_U1), u2=_interp_teff(teff, _LD_U2),
        density=RHO_SUN * float(target.mass) / float(target.radius) ** 3,
        spectral_type=spectral_type_from_teff(teff),
    )


def interp_cdpp(cdpp_ppm: Sequence[float], t14_hours: float) -> float:
    """Effective per-transit noise (ppm) at duration t14_hours: linear
    interpolation of the 14 DR25 RMS CDPP values over CDPP_DURATIONS_HR,
    clamped to the end values outside [1.5, 15] h (np.interp semantics)."""
    return float(np.interp(t14_hours, CDPP_DURATIONS_HR, cdpp_ppm))


def detected_planets(star: Star, planets: Sequence[Planet],
                     target) -> List[Planet]:
    """Planets from `planets` that this target's Kepler observation detects.

    Applies the module-level detection rule (see module docstring):
    geometric transit, >= MIN_TRANSITS transits in dataspan * dutycycle,
    limb-darkened SNR >= SNR_THRESHOLD against the CDPP interpolated to the
    transit duration."""
    span_days = float(target.dataspan) * float(target.dutycycle)
    out: List[Planet] = []
    for p in planets:
        n_tr = span_days / p.period
        if n_tr < MIN_TRANSITS:
            continue
        geom = compute_geometry(star, p)
        if not geom.transits or geom.t14_hours <= 0.0:
            continue
        sigma = interp_cdpp(target.cdpp_ppm, geom.t14_hours)
        if geom.depth_ppm / sigma * np.sqrt(n_tr) >= SNR_THRESHOLD:
            out.append(p)
    return out


@dataclass
class UniverseResult:
    """Detected catalog of one conditioned universe.

    n_k        : detected-multiplicity histogram, keys 1..5 (int) and "6+".
    detected   : per-system lists of (period_days, radius_re) for systems
                 with >= 1 detection, each list sorted by period (adjacency
                 for pair statistics); order follows the target list.
    n_targets  : number of targets simulated (denominator for rates).
    n_detected_planets : total detected planets across all systems.
    """
    n_k: Dict[object, int]
    detected: List[List[Tuple[float, float]]]
    n_targets: int
    n_detected_planets: int


def simulate_universe(targets: Sequence, seed: int,
                      arch: Architecture | None = None) -> UniverseResult:
    """Simulate one universe over the target list. Deterministic in
    (targets order, seed, arch): each target gets its own independent rng
    from SeedSequence(seed).spawn(len(targets)) — child i belongs to target
    i, so results are reproducible and target-parallelizable.

    Skips stellar-noise and atmosphere generation entirely (see module
    docstring): only system._draw_planets + geometry + CDPP thresholding
    run per target."""
    arch = DEFAULT_ARCH if arch is None else arch
    children = np.random.SeedSequence(seed).spawn(len(targets))
    n_k: Dict[object, int] = {k: 0 for k in N_K_KEYS}
    detected: List[List[Tuple[float, float]]] = []
    n_planets = 0
    for target, child in zip(targets, children):
        rng = np.random.default_rng(child)
        star = star_from_target(target)
        planets, _ = _draw_planets(rng, star, arch)
        found = detected_planets(star, planets, target)
        if not found:
            continue
        k = len(found)
        n_k[k if k <= 5 else "6+"] += 1
        n_planets += k
        detected.append([(p.period, p.radius) for p in found])
    return UniverseResult(n_k=n_k, detected=detected,
                          n_targets=len(targets),
                          n_detected_planets=n_planets)
