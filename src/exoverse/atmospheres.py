"""Exoplanet atmospheres: composition class, structure, and detectability.

Atmosphere assignment
---------------------
Class is assigned from bulk properties plus an escape criterion:

- giant / neptunian      -> primary H/He envelope (mu ~ 2.3)
- sub-neptune, rho < 2   -> H/He-rich envelope (mu ~ 3.5)
- sub-neptune, rho >= 2  -> steam / volatile "water world" (mu ~ 18)
- rocky                  -> "cosmic shoreline" test (Zahnle & Catling 2017):
                            a planet retains an atmosphere when its escape
                            velocity is high enough for its irradiation,
                            I_crit ~ (v_esc / v_esc_earth)^4 in Earth
                            insolation units. Retained -> secondary CO2/N2
                            (mu ~ 33); stripped -> airless.
  Rocky planets around M dwarfs get a QUESTIONABLE flag either way: whether
  they keep atmospheres at all is one of the biggest open questions JWST is
  currently working on (e.g. TRAPPIST-1 b/c results).

Structure: scale height H = k*Teq / (mu * m_H * g).

Transmission signal: the canonical estimate for one scale height is
delta_1H = 2 * H * Rp / R*^2; a cloud-free feature spans ~5 H. A random
cloud/haze suppression factor in [0.15, 1.0] is drawn (clouds are ubiquitous
and unpredictable - flagged as the dominant uncertainty).

Metrics
-------
- TSM: Transmission Spectroscopy Metric of Kempton+ (2018), eq. 1, with the
  published per-radius-bin scale factors and J-band magnitude.
- ESM: Emission Spectroscopy Metric of Kempton+ (2018), eq. 4: 7.5-micron
  blackbody flux ratio scaled to K magnitude ~ approximated here from J
  (mK ~ mJ - 0.5, INFO-level approximation).

Detectability: for each JWST instrument, the per-transit white-light noise
is converted to a per-spectral-bin noise (~sqrt(15) penalty for R~30-ish
binning across the band), compared against the expected feature amplitude
(~2 scale heights after cloud suppression), and the number of transits
needed for a 5-sigma feature detection is reported (capped evaluation at 25;
>25 is scored "impractical"). Stellar activity adds a spot-contamination
noise floor (transit light source effect, Rackham+ 2018).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List

import numpy as np

from .constants import G, M_EARTH, R_EARTH
from .flags import Flag, Severity
from .planets import Planet
from .stars import Star
from .stellar_noise import StellarNoise

K_B = 1.380649e-23
M_H = 1.6735575e-27
V_ESC_EARTH = 11186.0  # m/s

# White-light -> effective feature-detection noise penalty. Features are
# fit across the full band (not read off one spectral bin), so the penalty
# is mild; sqrt(5) reproduces published JWST results (e.g. K2-18b: features
# detected at ~5 sigma in ~2-3 transits of a J~9 M dwarf).
SPECTRAL_BIN_PENALTY = math.sqrt(5.0)
FEATURE_SCALE_HEIGHTS = 2.0              # typical single-feature amplitude
MAX_PRACTICAL_TRANSITS = 25


@dataclass
class Atmosphere:
    atm_class: str            # h_he | h_he_rich | steam | secondary | airless
    mu: float                 # mean molecular weight (amu); 0 for airless
    scale_height_km: float
    delta_1h_ppm: float       # transmission signal of one scale height
    feature_ppm: float        # expected feature after cloud suppression
    cloud_factor: float       # 0.15-1.0 multiplier applied to features
    tsm: float
    esm: float
    tsm_priority: bool        # above Kempton+18 threshold for its size bin
    geometric_albedo: float = 0.0   # V-band A_g, per-planet draw (0 = unset)
    flags: List[Flag] = field(default_factory=list)

    def add_flag(self, severity: Severity, rule: str, message: str) -> None:
        self.flags.append(Flag(severity, rule, message))


@dataclass
class AtmosphereObservation:
    observatory: str
    feature_snr_per_transit: float
    n_transits_5sigma: float   # inf -> impractical
    practical: bool            # <= MAX_PRACTICAL_TRANSITS and instrument usable
    note: str


# V-band geometric-albedo model: mean A_g vs Teq for H/He envelopes,
# following the Sudarsky/Marley cloud sequence bracketed by measurements:
# cold ammonia-cloud giants (Jupiter 0.52, Saturn 0.47, Neptune 0.44 at
# Teq ~ 50-130 K); predicted bright water-cloud regime near ~250 K; clear
# alkali-absorbing atmospheres 700-1500 K are DARK (Kepler hot-Jupiter
# population A_g ~ 0.03-0.11, TrES-2b < 0.04); ultra-hot silicate-cloud
# objects recover moderate brightness (Kepler-7b ~ 0.35).
_GIANT_ALBEDO_TEQ = [
    (50.0, 0.47), (110.0, 0.52), (250.0, 0.65), (400.0, 0.30),
    (700.0, 0.12), (1000.0, 0.06), (1500.0, 0.10), (2200.0, 0.30),
]


def _giant_albedo_mean(teq: float) -> float:
    ts, alb = zip(*_GIANT_ALBEDO_TEQ)
    return float(np.interp(teq, ts, alb))


def sample_geometric_albedo(rng: np.random.Generator, atm_class: str,
                            teq: float) -> float:
    """Draw a per-planet V-band geometric albedo for reflected-light imaging.

    Real albedos span far more than the solar-system anchors: measured hot
    Jupiters are nearly black (A_g 0.03-0.11) while Venus reaches 0.70 and
    icy surfaces (Europa 0.67, Enceladus > 1) are brighter still. Rocky
    anchors: Kepler super-Earth statistics 0.16-0.30 (Demory 2014), Earth
    0.24 (2026 phase-curve re-analysis; yield studies conventionally use
    0.2), Mars 0.17, Mercury 0.14, Moon 0.07.

    Exactly two rng draws are consumed for every planet regardless of class
    (a mixture selector and a scatter deviate), keeping the random stream
    layout uniform."""
    u = float(rng.random())     # mixture / branch selector
    z = float(rng.normal())     # scatter deviate

    if atm_class in ("h_he", "h_he_rich"):
        mean = _giant_albedo_mean(teq)
        if atm_class == "h_he_rich":
            mean *= 0.9   # metal-enriched envelopes run slightly darker
        # Wide lognormal scatter: measured hot-Jupiter spread is ~x3
        return float(np.clip(mean * math.exp(0.35 * z), 0.02, 0.85))
    if atm_class == "steam":
        # No measurements exist; water clouds plausibly bright, hazes dark
        return float(np.clip(0.30 * math.exp(0.30 * z), 0.08, 0.60))
    if atm_class == "secondary":
        if u < 0.22:   # Venus-like global cloud deck
            return float(np.clip(0.65 + 0.07 * z, 0.45, 0.80))
        return float(np.clip(0.24 + 0.08 * z, 0.08, 0.45))
    # airless: dark regolith, with an icy-bright branch for cold worlds
    if teq < 200.0 and u < 0.15:
        return float(np.clip(0.55 + 0.15 * z, 0.25, 0.90))
    return float(np.clip(0.12 + 0.05 * z, 0.04, 0.30))


def surface_gravity(planet: Planet) -> float:
    return G * planet.mass * M_EARTH / (planet.radius * R_EARTH) ** 2


def escape_velocity(planet: Planet) -> float:
    return math.sqrt(2.0 * G * planet.mass * M_EARTH / (planet.radius * R_EARTH))


def assign_atmosphere(rng: np.random.Generator, star: Star,
                      planet: Planet) -> Atmosphere:
    g = surface_gravity(planet)
    comp = planet.composition_class

    if comp in ("giant", "neptunian"):
        atm_class, mu = "h_he", 2.3
    elif comp == "sub-neptune":
        if planet.density < 2.0:
            atm_class, mu = "h_he_rich", 3.5
        else:
            atm_class, mu = "steam", 18.0
    else:  # rocky: cosmic shoreline, I_crit ~ v_esc^4 (Zahnle & Catling 2017)
        # Normalized through Mars (v_esc=5.03 km/s at S=0.43), the solar-system
        # body sitting on the empirical line: I_crit = 10.5 (v/v_earth)^4.
        # Earth (I_crit=10.5 vs S=1) retains; Moon/Mercury strip. Normalizing
        # through Earth instead would put Earth itself on the knife's edge.
        i_crit = 10.5 * (escape_velocity(planet) / V_ESC_EARTH) ** 4
        if planet.insolation > i_crit:
            atm_class, mu = "airless", 0.0
        else:
            atm_class, mu = "secondary", 33.0

    if mu > 0.0:
        h_m = K_B * planet.teq / (mu * M_H * g)
    else:
        h_m = 0.0
    rp_m = planet.radius * R_EARTH
    rs_m = star.radius * 6.957e8
    delta_1h = 2.0 * h_m * rp_m / rs_m ** 2 * 1e6 if mu > 0 else 0.0
    cloud = float(rng.uniform(0.15, 1.0)) if mu > 0 else 1.0
    feature = delta_1h * FEATURE_SCALE_HEIGHTS * cloud
    albedo = sample_geometric_albedo(rng, atm_class, planet.teq)

    atm = Atmosphere(atm_class=atm_class, mu=mu,
                     scale_height_km=h_m / 1000.0,
                     delta_1h_ppm=delta_1h, feature_ppm=feature,
                     cloud_factor=cloud,
                     tsm=compute_tsm(star, planet),
                     esm=compute_esm(star, planet),
                     tsm_priority=False,
                     geometric_albedo=albedo)
    atm.tsm_priority = is_tsm_priority(planet, atm.tsm)

    # ---- flags -------------------------------------------------------------
    if comp == "rocky" and star.teff < 3900:
        atm.add_flag(Severity.QUESTIONABLE, "atmosphere.m_dwarf_rocky_retention",
                     "Whether rocky planets around M dwarfs retain atmospheres "
                     "at all is an open question (XUV history; cf. JWST "
                     "TRAPPIST-1 b/c null results); class assignment here is "
                     "a coin-toss informed by the cosmic shoreline")
    if atm_class == "steam":
        atm.add_flag(Severity.QUESTIONABLE, "atmosphere.water_world_degeneracy",
                     "Bulk density is degenerate between volatile-rich "
                     "interiors and rock+H/He; steam classification is one of "
                     "several allowed solutions")
    if mu > 0 and cloud < 0.35:
        atm.add_flag(Severity.INFO, "atmosphere.high_cloud_suppression",
                     f"Drawn cloud/haze factor {cloud:.2f}: features strongly "
                     "muted (cf. GJ 1214b); this is common and unpredictable")
    if albedo >= 0.45 and atm_class in ("secondary", "airless"):
        atm.add_flag(Severity.INFO, "atmosphere.high_albedo_draw",
                     f"Drawn A_g={albedo:.2f}: Venus-like cloud deck or icy "
                     "surface branch; strongly boosts reflected-light "
                     "detectability and is unconstrained for exoplanets")
    if planet.teq > 2000 and atm_class in ("h_he", "h_he_rich"):
        atm.add_flag(Severity.INFO, "atmosphere.thermal_dissociation",
                     "Ultra-hot: molecular features partially dissociated; "
                     "simple scale-height signal estimate is optimistic")
    return atm


def compute_tsm(star: Star, planet: Planet) -> float:
    """Kempton+ 2018 eq. 1 (0 if outside table or missing prerequisites)."""
    r = planet.radius
    if r < 1.5:
        sf = 0.19
    elif r < 2.75:
        sf = 1.26
    elif r < 4.0:
        sf = 1.28
    elif r <= 10.0:
        sf = 1.19
    else:
        return 0.0
    return (sf * r ** 3 * planet.teq /
            (planet.mass * star.radius ** 2) * 10.0 ** (-star.mag_j / 5.0))


def is_tsm_priority(planet: Planet, tsm: float) -> bool:
    """Kempton+ 2018 first-tier thresholds: 10 for R<1.5, 90 above."""
    return tsm > (10.0 if planet.radius < 1.5 else 90.0)


def _planck_ratio_7p5um(t_planet: float, t_star: float) -> float:
    """B(7.5um, Tp) / B(7.5um, Ts)."""
    if t_planet <= 0:
        return 0.0
    x = 1.98645e-25 / (7.5e-6 * K_B)  # h*c/(lambda*k)
    def b(t):
        return 1.0 / (math.exp(x / t) - 1.0)
    return b(t_planet) / b(t_star)


def compute_esm(star: Star, planet: Planet) -> float:
    """Kempton+ 2018 eq. 4; K magnitude approximated as J - 0.5."""
    t_day = 1.10 * planet.teq
    k = planet.radius * R_EARTH / (star.radius * 6.957e8)
    mag_k = star.mag_j - 0.5
    return 4.29e6 * _planck_ratio_7p5um(t_day, star.teff) * k * k \
        * 10.0 ** (-mag_k / 5.0)


def score_atmosphere_observability(
        atm: Atmosphere, t14_hours: float, noise: StellarNoise,
        jwst_instruments: List[tuple[str, float, bool]],
) -> List[AtmosphereObservation]:
    """Score spectroscopic detectability per JWST instrument.

    jwst_instruments: list of (name, sigma_1hr_white_ppm, usable).
    """
    out: List[AtmosphereObservation] = []
    t14 = max(t14_hours, 0.3)
    # Stellar contamination floor (transit light source effect): scales with
    # activity; NIR band factor already partially mitigates -> use 0.4
    contamination = 0.4 * 0.02 * noise.var_amp_ppm  # ~2% of variability amp
    for name, sigma_1hr, usable in jwst_instruments:
        if not usable or atm.mu <= 0.0 or atm.feature_ppm <= 0.0:
            out.append(AtmosphereObservation(
                name, 0.0, float("inf"), False,
                "airless / no feature" if atm.mu <= 0 else "instrument unusable"))
            continue
        sigma_bin = sigma_1hr / math.sqrt(t14) * SPECTRAL_BIN_PENALTY
        sigma_bin = math.sqrt(sigma_bin ** 2 + contamination ** 2)
        snr1 = atm.feature_ppm / sigma_bin
        n_needed = (5.0 / snr1) ** 2 if snr1 > 0 else float("inf")
        practical = n_needed <= MAX_PRACTICAL_TRANSITS
        note = (f"feature {atm.feature_ppm:.0f} ppm vs {sigma_bin:.0f} ppm/bin"
                f"/transit; stellar contamination floor "
                f"{contamination:.0f} ppm")
        out.append(AtmosphereObservation(name, snr1, n_needed, practical, note))
    return out
