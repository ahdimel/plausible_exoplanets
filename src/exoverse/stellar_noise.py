"""Astrophysical (stellar) noise model.

Real stars are not photometrically quiet. Four components are modeled, each
expressed as an effective white-noise contribution (ppm) on transit
timescales, plus a slow variability amplitude:

1. Granulation - convective surface flicker. Amplitude scales inversely
   with the peak oscillation frequency nu_max ~ g / sqrt(Teff) (Kallinger+
   2014; Bastien+ 2013 "flicker"). Sun: ~40 ppm on hours timescales.
2. p-mode oscillations - solar-like pulsations, amplitude ~ (L/M)^0.9
   (Kjeldsen & Bedding 1995). Sun: ~3 ppm. Minutes-timescale, mostly
   averages down over a transit; kept for evolved/massive dwarfs.
3. Rotational spot variability - amplitude set by activity, which declines
   with Rossby number / age (gyrochronology: Prot ~ sqrt(age), Skumanich
   1972; Barnes 2007). Only the residual over a transit duration matters
   for depth measurement: sigma ~ A_var * (T14 / Prot), plus a spot-crossing
   term for active stars. M dwarfs stay active far longer than FGK stars
   (Newton+ 2017).
4. Flares - white-light flares on active M dwarfs contaminate light curves
   stochastically; modeled as an extra white-noise term for late-M/active
   stars (Davenport 2016 statistics, crudely).

Chromaticity: spots and flares are contrast effects that weaken toward the
red/IR. Each observatory band applies a scaling factor (~1.0 in the blue
optical, ~0.75 TESS red-optical, ~0.35-0.45 in JWST/Roman NIR bands;
Rackham+ 2018).

All coefficients are order-of-magnitude calibrations to Kepler statistics
(McQuillan+ 2014 rotation; Basri+ 2013 variability fractions), adequate for
population-level detectability studies but not for fitting any single star.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .flags import Severity
from .stars import Star

NU_MAX_SUN = 3090.0   # micro-Hz
PROT_SUN = 24.5       # days
AGE_SUN = 4.57        # Gyr


@dataclass
class StellarNoise:
    prot_days: float          # rotation period
    activity_level: float     # 0 (dead quiet) .. ~1+ (very active), ~Rvar proxy
    var_amp_ppm: float        # semi-amplitude of rotational variability
    gran_1hr_ppm: float       # granulation white-noise equivalent, 1-hr bins
    osc_1hr_ppm: float        # p-mode residual, 1-hr bins
    flare_1hr_ppm: float      # flare contamination white-noise equivalent

    def sigma_transit_ppm(self, t14_hours: float, band_factor: float) -> float:
        """Effective stellar noise (ppm) on a transit-depth measurement for
        a transit of duration t14, in a band with spot/flare contrast factor
        band_factor (1.0 = Kepler-like optical)."""
        t14 = max(t14_hours, 0.1)
        # Granulation: correlated -> averages down slowly (~t^-0.35 empirically)
        gran = self.gran_1hr_ppm / t14 ** 0.35
        # Oscillations: incoherent over hours -> averages as white noise
        osc = self.osc_1hr_ppm / math.sqrt(t14)
        # Rotational variability: linear trend across the transit window that
        # survives detrending, ~ amplitude * fraction of rotation elapsed
        spot = band_factor * self.var_amp_ppm * min(t14 / (self.prot_days * 24.0), 0.5)
        # Spot-crossing bumps for active stars: ~5% of variability amplitude
        crossing = band_factor * 0.05 * self.var_amp_ppm
        flare = band_factor * self.flare_1hr_ppm / math.sqrt(t14)
        return math.sqrt(gran ** 2 + osc ** 2 + spot ** 2 + crossing ** 2 + flare ** 2)


def generate_stellar_noise(rng: np.random.Generator, star: Star) -> StellarNoise:
    """Draw a noise state consistent with the star's mass, age, and structure."""
    # --- Rotation from gyrochronology --------------------------------------
    # Prot ~ sqrt(age), normalized to the Sun, steeper mass dependence for
    # low-mass stars (which spin down to ~100 d at late ages, Newton+ 2017)
    mass_fac = (max(star.mass, 0.1)) ** (-0.6)
    prot = PROT_SUN * math.sqrt(max(star.age_gyr, 0.1) / AGE_SUN) * mass_fac
    prot *= float(np.exp(rng.normal(0.0, 0.25)))
    prot = float(np.clip(prot, 0.3, 180.0))

    # --- Activity: declines with age; M dwarfs stay active longer ----------
    if star.mass < 0.35:
        spin_down_gyr = 5.0        # fully convective stars stay active ~Gyrs
    elif star.mass < 0.6:
        spin_down_gyr = 2.5
    else:
        spin_down_gyr = 1.0
    activity = math.exp(-star.age_gyr / spin_down_gyr)
    activity = float(np.clip(activity + rng.normal(0.0, 0.08), 0.005, 1.5))

    # Rotational variability semi-amplitude: quiet Sun ~200-400 ppm; young /
    # active stars 5,000-20,000 ppm (McQuillan+ 2014 span)
    var_amp = 300.0 + 15000.0 * activity ** 1.5
    var_amp *= float(np.exp(rng.normal(0.0, 0.5)))
    var_amp = float(np.clip(var_amp, 50.0, 50000.0))

    # --- Granulation via nu_max ---------------------------------------------
    g_ratio = star.mass / star.radius ** 2                  # g/g_sun
    nu_max = NU_MAX_SUN * g_ratio / math.sqrt(star.teff / 5772.0)
    gran = 40.0 * (NU_MAX_SUN / nu_max) ** 0.6              # Kallinger+ 2014
    gran = float(np.clip(gran * np.exp(rng.normal(0.0, 0.15)), 5.0, 500.0))

    # --- Oscillations ---------------------------------------------------------
    osc = 3.0 * (star.luminosity / star.mass) ** 0.9
    osc = float(np.clip(osc, 0.2, 100.0))

    # --- Flares (active M dwarfs) ---------------------------------------------
    flare = 0.0
    if star.mass < 0.5:
        flare = 200.0 * activity * float(np.exp(rng.normal(0.0, 0.4)))

    noise = StellarNoise(prot_days=prot, activity_level=activity,
                         var_amp_ppm=var_amp, gran_1hr_ppm=gran,
                         osc_1hr_ppm=osc, flare_1hr_ppm=flare)

    # --- Flags on the star ----------------------------------------------------
    if activity > 0.5:
        star.add_flag(Severity.QUESTIONABLE, "star.high_activity",
                      f"Young/active star (activity index {activity:.2f}, "
                      f"Prot={prot:.1f} d, variability ~{var_amp:.0f} ppm): "
                      "transit depths biased by spots (Rackham+ 2018 transit "
                      "light source effect); atmosphere retrievals unreliable")
    if star.mass < 0.5 and flare > 150.0:
        star.add_flag(Severity.INFO, "star.flaring_m_dwarf",
                      f"Flaring M dwarf (~{flare:.0f} ppm/hr flare noise): "
                      "individual transits may need flare masking")
    return noise


# Band-dependent spot/flare contrast factors (Rackham+ 2018-like)
BAND_FACTORS = {
    "optical_blue": 1.0,      # Kepler, ground V
    "optical_red": 0.75,      # TESS
    "nir": 0.40,              # JWST NIRISS/NIRSpec, Roman F146
}
