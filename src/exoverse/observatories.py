"""Observatory noise models and transit detectability.

Each observatory model answers: given this star's brightness and this
transit's depth/duration, what photometric precision do we get, and is the
signal detectable?

Models (all approximate, calibrated to published performance; every result
carries the caveat that systematics are represented by simple noise floors):

- TESS          : 1-hr precision vs T magnitude fit to the Sullivan+ (2015) /
                  Stassun+ (2018) noise curves; 60 ppm/hr systematic floor;
                  27.4-day sector baseline (1 sector by default). Detection:
                  phase-folded SNR >= 7.1 and >= 2 transits observed.
- Kepler        : archival-style 4-yr baseline; 6.5-hr CDPP ~30 ppm at Kp=12
                  scaled by photon statistics with a 20 ppm floor. Kepler is
                  retired: results represent "what Kepler would have seen"
                  (INFO-level caveat in the database).
- JWST NIRISS SOSS  : single-transit spectrophotometry, white-light curve.
                  ~20 ppm/hr photon-limited precision at J=8, 10 ppm floor;
                  saturates for J < 6.5. Detection: single-transit SNR >= 5
                  (targeted follow-up, ephemeris assumed known).
- JWST NIRSpec Prism: highest sensitivity for faint hosts, saturates J < 10.5.
                  ~12 ppm/hr at J=11, 10 ppm floor. Same detection rule.
- Ground 1-m survey : 2 mmag/point + scintillation-like floor of 1 mmag/hr
                  binned; 90-night campaign, 33% usable duty cycle, transits
                  observable only if depth > ~1 mmag. Detection: folded
                  SNR >= 7 and >= 3 transits.

The per-transit noise scales as sigma_1hr / sqrt(T14 in hours); the folded
mission SNR as per-transit SNR x sqrt(N_transits). This is the standard
box-least-squares SNR scaling (e.g. Kovacs+ 2002).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List

from .planets import Planet
from .stars import Star
from .transits import TransitGeometry


@dataclass
class Observation:
    observatory: str
    usable: bool                # target observable at all (saturation, depth floor)
    note: str
    sigma_1hr_ppm: float
    snr_per_transit: float
    n_transits: float
    snr_total: float
    detectable: bool


def _tess_sigma_1hr(tmag: float) -> float:
    """1-hour combined noise in ppm vs TESS magnitude (fit to Stassun+18:
    ~150 ppm @ T=8, ~500 @ T=10, ~1750 @ T=12, ~6000 @ T=14)."""
    photon = 10.0 ** (0.267 * tmag + 0.04)
    return math.sqrt(photon ** 2 + 60.0 ** 2)


def _kepler_sigma_1hr(kp: float) -> float:
    """Kepler ~6.5h CDPP 30 ppm at Kp=12 -> convert to 1-hr, photon scaling."""
    cdpp65 = 30.0 * 10.0 ** (0.2 * (kp - 12.0))
    sigma_1hr = cdpp65 * math.sqrt(6.5)
    return math.sqrt(sigma_1hr ** 2 + 20.0 ** 2)


def _jwst_niriss_sigma_1hr(jmag: float) -> float:
    photon = 20.0 * 10.0 ** (0.2 * (jmag - 8.0))
    return math.sqrt(photon ** 2 + 10.0 ** 2)


def _jwst_nirspec_sigma_1hr(jmag: float) -> float:
    photon = 12.0 * 10.0 ** (0.2 * (jmag - 11.0))
    return math.sqrt(photon ** 2 + 10.0 ** 2)


def observe(star: Star, planet: Planet, geom: TransitGeometry) -> List[Observation]:
    """Evaluate one transiting planet against all modeled observatories."""
    out: List[Observation] = []
    if not geom.transits:
        return out
    depth = geom.depth_ppm
    t14 = max(geom.t14_hours, 0.05)

    def folded(name: str, sigma1: float, baseline_days: float, duty: float,
               min_transits: float, snr_req: float, usable: bool, note: str):
        n_tr = baseline_days * duty / planet.period
        sigma_tr = sigma1 / math.sqrt(t14)
        snr_1 = depth / sigma_tr if usable else 0.0
        snr_tot = snr_1 * math.sqrt(max(n_tr, 0.0))
        det = usable and n_tr >= min_transits and snr_tot >= snr_req
        out.append(Observation(name, usable, note, sigma1, snr_1, n_tr, snr_tot, det))

    # TESS: 1 sector
    folded("TESS (1 sector)", _tess_sigma_1hr(star.mag_tess), 27.4, 1.0, 2.0, 7.1,
           star.mag_tess > 4.0,
           "27.4-d sector; SNR>=7.1 & >=2 transits" if star.mag_tess > 4.0
           else "saturated (T<4)")

    # Kepler archival: Kp ~ V for FGK, ~V-0.3 for M; use V as proxy
    folded("Kepler (4yr, archival)", _kepler_sigma_1hr(star.mag_v), 1460.0, 0.92,
           3.0, 7.1, star.mag_v > 6.0,
           "hypothetical: Kepler retired 2018; Kp~V proxy" if star.mag_v > 6.0
           else "saturated (Kp<6)")

    # Ground 1-m
    ground_sigma = math.sqrt(2000.0 ** 2 + 1000.0 ** 2)  # ppm per 1-hr bin
    folded("Ground 1-m survey", ground_sigma, 90.0, 0.33, 3.0, 7.0,
           depth > 1000.0,
           "90 nights, 33% duty; needs depth >~1 mmag" if depth > 1000.0
           else "depth below ground-based systematics floor")

    # JWST: targeted single transit
    for name, sigma_fn, sat_limit in (
        ("JWST NIRISS SOSS", _jwst_niriss_sigma_1hr, 6.5),
        ("JWST NIRSpec Prism", _jwst_nirspec_sigma_1hr, 10.5),
    ):
        usable = star.mag_j > sat_limit
        sigma1 = sigma_fn(star.mag_j)
        sigma_tr = sigma1 / math.sqrt(t14)
        snr = depth / sigma_tr if usable else 0.0
        out.append(Observation(
            name, usable,
            f"single targeted transit; saturates J<{sat_limit}" if usable
            else f"saturated (J={star.mag_j:.1f} < {sat_limit})",
            sigma1, snr, 1.0, snr, usable and snr >= 5.0))
    return out
