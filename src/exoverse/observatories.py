"""Observatory noise models and transit/imaging detectability.

Each observatory model answers: given this star's brightness, this star's
astrophysical noise, and this transit's depth/duration, what precision do we
get and is the signal detectable?

Every noise budget is sigma_total = sqrt(sigma_instrument^2 + sigma_stellar^2)
where sigma_stellar comes from `stellar_noise.StellarNoise` (granulation,
oscillations, spots, flares) with a band-dependent spot-contrast factor.

Operating observatories (calibrated to published performance):
- TESS          : 1-hr noise fit to Stassun+ 2018 (~150 ppm @ T=8, ~6000 @
                  T=14), 60 ppm floor; 27.4-d sector. SNR>=7.1, >=2 transits.
- Kepler        : archival (retired 2018); 30 ppm CDPP@Kp=12 photon-scaled,
                  20 ppm floor, 4-yr baseline. SNR>=7.1, >=3 transits.
- JWST NIRISS SOSS  : ~20 ppm/hr @ J=8, floor 10 ppm, saturates J<6.5.
                  Targeted single transit, SNR>=5.
- JWST NIRSpec Prism: ~12 ppm/hr @ J=11, floor 10 ppm, saturates J<10.5.
- Ground 1-m    : 2 mmag + 1 mmag floors, 90-night campaign, 33% duty.

Future observatories (specs WILL evolve; see docs/OBSERVATORIES.md):
- Roman GBTDS   : Nancy Grace Roman Space Telescope (2.4 m, launch late
                  2026, science 2027) Galactic Bulge Time Domain Survey:
                  F146 (0.93-2.0 um) every ~12 min, 6 seasons x ~72 d over
                  ~5 yr, sensitive for F146 ~ 8-21; predicted to find
                  60k-200k transiting planets (Wilson+ 2023). Modeled as
                  photon noise anchored at ~700 ppm/hr @ F146=16 with a
                  200 ppm bright floor; F146 ~ J proxy. Detection SNR>=7.1,
                  >=3 transits. NOTE: GBTDS points at the bulge; applying it
                  to our local population answers "if this system were in a
                  Roman field" (INFO note attached to every result).
- HWO           : Habitable Worlds Observatory (~6-8 m, launch 2040s, specs
                  pre-Phase-A): direct imaging, not transits. Requirement:
                  raw contrast ~1e-10, best post-processed sensitivity
                  ~3e-11 on bright stars, IWA ~3 lambda/D ~ 60 mas @ V.
                  We compute the planet's reflected-light contrast at
                  quadrature (per-planet geometric albedo drawn in
                  atmospheres.py; class-mean fallback) and its angular
                  separation, and score it against a photon-limited
                  contrast floor that degrades for hosts fainter than V~7.
                  Evaluated for hosts V<11 within ~30 pc; the IWA and the
                  degraded floor - not a hard magnitude cut - decide the
                  faint M-dwarf cases. Goal: ~25 exo-Earths (Astro2020).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

from .constants import AU, R_EARTH
from .planets import Planet
from .stars import Star
from .stellar_noise import BAND_FACTORS, StellarNoise
from .transits import TransitGeometry


@dataclass
class Observation:
    observatory: str
    mode: str                   # "transit" | "imaging"
    usable: bool
    note: str
    sigma_1hr_ppm: float        # instrument-only white noise (0 for imaging)
    sigma_stellar_ppm: float    # stellar noise on this measurement
    snr_per_transit: float      # transit: per-event SNR; imaging: contrast SNR proxy
    n_transits: float
    snr_total: float
    detectable: bool
    # Imaging-specific extras (0 for transit rows)
    contrast: float = 0.0
    separation_mas: float = 0.0


# ------------------------------------------------------- instrument noise fits
def _tess_sigma_1hr(tmag: float) -> float:
    photon = 10.0 ** (0.267 * tmag + 0.04)
    return math.sqrt(photon ** 2 + 60.0 ** 2)


def _kepler_sigma_1hr(kp: float) -> float:
    cdpp65 = 30.0 * 10.0 ** (0.2 * (kp - 12.0))
    sigma_1hr = cdpp65 * math.sqrt(6.5)
    return math.sqrt(sigma_1hr ** 2 + 20.0 ** 2)


def _jwst_niriss_sigma_1hr(jmag: float) -> float:
    photon = 20.0 * 10.0 ** (0.2 * (jmag - 8.0))
    return math.sqrt(photon ** 2 + 10.0 ** 2)


def _jwst_nirspec_sigma_1hr(jmag: float) -> float:
    photon = 12.0 * 10.0 ** (0.2 * (jmag - 11.0))
    return math.sqrt(photon ** 2 + 10.0 ** 2)


def _roman_sigma_1hr(f146: float) -> float:
    """Photon-scaled around ~700 ppm/hr at F146=16 (Wilson+ 2023 regime),
    200 ppm bright-end floor (systematics, saturation handling)."""
    photon = 700.0 * 10.0 ** (0.2 * (f146 - 16.0))
    return math.sqrt(photon ** 2 + 200.0 ** 2)


JWST_INSTRUMENTS = [
    ("JWST NIRISS SOSS", _jwst_niriss_sigma_1hr, 6.5),
    ("JWST NIRSpec Prism", _jwst_nirspec_sigma_1hr, 10.5),
]


# ------------------------------------------------------------------ transits
def observe(star: Star, planet: Planet, geom: TransitGeometry,
            noise: Optional[StellarNoise] = None) -> List[Observation]:
    """Evaluate one planet against all observatories. Transit rows require a
    transiting geometry; the HWO imaging row is computed for any geometry."""
    out: List[Observation] = []

    if geom.transits:
        depth = geom.depth_ppm
        t14 = max(geom.t14_hours, 0.05)

        def stellar(band: str) -> float:
            if noise is None:
                return 0.0
            return noise.sigma_transit_ppm(t14, BAND_FACTORS[band])

        def folded(name, sigma1, band, baseline_days, duty, min_tr, snr_req,
                   usable, note):
            sig_star = stellar(band)
            n_tr = baseline_days * duty / planet.period
            sigma_tr = math.sqrt((sigma1 / math.sqrt(t14)) ** 2 + sig_star ** 2)
            snr_1 = depth / sigma_tr if usable else 0.0
            snr_tot = snr_1 * math.sqrt(max(n_tr, 0.0))
            det = usable and n_tr >= min_tr and snr_tot >= snr_req
            out.append(Observation(name, "transit", usable, note, sigma1,
                                   sig_star, snr_1, n_tr, snr_tot, det))

        folded("TESS (1 sector)", _tess_sigma_1hr(star.mag_tess), "optical_red",
               27.4, 1.0, 2.0, 7.1, star.mag_tess > 4.0,
               "27.4-d sector; SNR>=7.1 & >=2 transits"
               if star.mag_tess > 4.0 else "saturated (T<4)")

        folded("Kepler (4yr, archival)", _kepler_sigma_1hr(star.mag_v),
               "optical_blue", 1460.0, 0.92, 3.0, 7.1, star.mag_v > 6.0,
               "hypothetical: Kepler retired 2018; Kp~V proxy"
               if star.mag_v > 6.0 else "saturated (Kp<6)")

        ground_sigma = math.sqrt(2000.0 ** 2 + 1000.0 ** 2)
        folded("Ground 1-m survey", ground_sigma, "optical_blue", 90.0, 0.33,
               3.0, 7.0, depth > 1000.0,
               "90 nights, 33% duty; needs depth >~1 mmag"
               if depth > 1000.0 else "depth below ground-based systematics floor")

        # Roman GBTDS: F146 ~ J-band proxy; 6 seasons x 72 d, 12-min cadence
        f146 = star.mag_j
        roman_usable = 8.0 < f146 < 21.0
        folded("Roman GBTDS (2027+)", _roman_sigma_1hr(f146), "nir",
               6 * 72.0, 0.95, 3.0, 7.1, roman_usable,
               "FUTURE/bulge survey: 'if this system were in a Roman field'; "
               "specs will evolve" if roman_usable
               else f"outside F146 8-21 dynamic range (F146~{f146:.1f})")

        for name, sigma_fn, sat in JWST_INSTRUMENTS:
            usable = star.mag_j > sat
            sigma1 = sigma_fn(star.mag_j)
            sig_star = stellar("nir")
            sigma_tr = math.sqrt((sigma1 / math.sqrt(t14)) ** 2 + sig_star ** 2)
            snr = depth / sigma_tr if usable else 0.0
            out.append(Observation(
                name, "transit", usable,
                f"single targeted transit; saturates J<{sat}" if usable
                else f"saturated (J={star.mag_j:.1f} < {sat})",
                sigma1, sig_star, snr, 1.0, snr, usable and snr >= 5.0))

    out.append(hwo_imaging(star, planet))
    return out


# ------------------------------------------------------------------- imaging
HWO_CONTRAST_FLOOR = 3e-11    # best post-processed floor, bright hosts (req-era)
HWO_IWA_MAS = 60.0            # ~3 lambda/D at V band for ~6.5 m
HWO_OWA_MAS = 500.0
HWO_VMAG_REF = 7.0            # floor is photon-limited fainter than this
HWO_VMAG_LIMIT = 11.0         # fainter hosts: exposure times impractical
HWO_DIST_LIMIT_PC = 30.0

# Class-MEAN V-band geometric albedos — fallback only, used when a planet
# carries no per-planet draw (atmospheres.sample_geometric_albedo is the
# primary model: temperature-dependent, with measured-population scatter).
HWO_GEOMETRIC_ALBEDO = {
    "h_he": 0.50,
    "h_he_rich": 0.40,
    "steam": 0.35,
    "secondary": 0.30,
    "airless": 0.12,
}
HWO_DEFAULT_ALBEDO = 0.3


def hwo_contrast_floor(mag_v: float) -> float:
    """Post-processed 1-sigma-class contrast floor vs host brightness:
    systematics-limited at HWO_CONTRAST_FLOOR for V <= 7, photon-limited
    (10^(0.2*(V-7))) for fainter hosts."""
    return HWO_CONTRAST_FLOOR * 10.0 ** (0.2 * max(mag_v - HWO_VMAG_REF, 0.0))


def hwo_imaging(star: Star, planet: Planet) -> Observation:
    """Reflected-light direct-imaging detectability for HWO (2040s, notional).

    Contrast at quadrature: C = A_g * Phi(90deg) * (Rp / a)^2 with the
    planet's own drawn geometric albedo A_g (atmospheres.
    sample_geometric_albedo; class-mean fallback if unset) and Lambert
    phase Phi(90)=1/pi. Angular separation is the quadrature projected
    separation a/d. Detection requires the contrast to clear the
    brightness-dependent post-processed floor within the working angles."""
    atm = planet.atmosphere
    if atm is not None and getattr(atm, "geometric_albedo", 0.0) > 0.0:
        albedo = atm.geometric_albedo
    elif atm is not None:
        albedo = HWO_GEOMETRIC_ALBEDO.get(atm.atm_class, HWO_DEFAULT_ALBEDO)
    else:
        albedo = HWO_DEFAULT_ALBEDO
    a_m = planet.a * AU
    contrast = albedo * (1.0 / math.pi) * (planet.radius * R_EARTH / a_m) ** 2
    sep_mas = planet.a / star.distance_pc * 1000.0   # a[AU]/d[pc] = sep["]

    reachable = (star.mag_v < HWO_VMAG_LIMIT
                 and star.distance_pc < HWO_DIST_LIMIT_PC)
    resolved = HWO_IWA_MAS < sep_mas < HWO_OWA_MAS
    floor = hwo_contrast_floor(star.mag_v)
    bright_enough = contrast > floor
    det = reachable and resolved and bright_enough

    if not reachable:
        note = (f"host outside HWO reach (V={star.mag_v:.1f}, "
                f"d={star.distance_pc:.0f} pc; needs V<{HWO_VMAG_LIMIT:.0f}, "
                f"d<{HWO_DIST_LIMIT_PC:.0f} pc)")
    elif not resolved:
        note = (f"separation {sep_mas:.1f} mas outside working angles "
                f"[{HWO_IWA_MAS:.0f}, {HWO_OWA_MAS:.0f}] mas")
    elif not bright_enough:
        note = (f"contrast {contrast:.1e} (A_g={albedo:.2f}) below "
                f"{floor:.1e} floor at V={star.mag_v:.1f}")
    else:
        note = (f"FUTURE (2040s, pre-Phase-A specs): directly imageable "
                f"(A_g={albedo:.2f}, floor {floor:.1e} at V={star.mag_v:.1f})")
    snr_proxy = contrast / floor if reachable and resolved else 0.0
    return Observation("HWO imaging (2040s)", "imaging", reachable, note,
                       0.0, 0.0, snr_proxy, 0.0, snr_proxy, det,
                       contrast=contrast, separation_mas=sep_mas)


N_TRANSIT_OBSERVATORIES = 6   # transit rows produced per transiting planet
