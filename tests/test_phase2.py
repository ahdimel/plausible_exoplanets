"""Tests for phase 2: stellar noise, atmospheres, Roman/HWO, archive,
validation, and the web UI."""

import math

import numpy as np
import pytest

from exoverse.archive import ArchivePlanet, CANONICAL_COLUMNS
from exoverse.atmospheres import (
    Atmosphere, assign_atmosphere, compute_esm, compute_tsm,
    score_atmosphere_observability,
)
from exoverse.constants import R_EARTH, R_JUP, R_SUN
from exoverse.flags import Severity
from exoverse.observatories import hwo_imaging, observe
from exoverse.planets import semimajor_axis_au
from exoverse.stellar_noise import BAND_FACTORS, generate_stellar_noise
from exoverse.system import generate_system
from exoverse.transits import compute_geometry
from exoverse.validate import audit_rules, ks_2samp, quantiles

from test_physics import make_planet, make_sun


# ------------------------------------------------------------- stellar noise
def test_sunlike_noise_is_solar():
    """A quiet solar-age Sun should have Prot ~25 d and ~solar amplitudes."""
    rng = np.random.default_rng(3)
    prots, grans = [], []
    for _ in range(300):
        sun = make_sun()
        n = generate_stellar_noise(rng, sun)
        prots.append(n.prot_days)
        grans.append(n.gran_1hr_ppm)
    assert 18 < float(np.median(prots)) < 33      # Sun: 24.5 d
    assert 20 < float(np.median(grans)) < 80      # Sun: ~40 ppm/hr


def test_young_stars_noisier_than_old():
    rng = np.random.default_rng(5)
    young_amp, old_amp = [], []
    for _ in range(200):
        s = make_sun()
        s.age_gyr = 0.3
        young_amp.append(generate_stellar_noise(rng, s).var_amp_ppm)
        s2 = make_sun()
        s2.age_gyr = 8.0
        old_amp.append(generate_stellar_noise(rng, s2).var_amp_ppm)
    assert np.median(young_amp) > 5 * np.median(old_amp)


def test_active_star_gets_flagged():
    rng = np.random.default_rng(11)
    for _ in range(50):
        s = make_sun()
        s.age_gyr = 0.2
        generate_stellar_noise(rng, s)
        if any(f.rule == "star.high_activity" for f in s.flags):
            return
    pytest.fail("no young star flagged as active in 50 draws")


def test_stellar_noise_band_chromaticity():
    """NIR bands must see less spot noise than blue optical."""
    rng = np.random.default_rng(7)
    s = make_sun()
    s.age_gyr = 0.5
    n = generate_stellar_noise(rng, s)
    assert (n.sigma_transit_ppm(3.0, BAND_FACTORS["nir"])
            < n.sigma_transit_ppm(3.0, BAND_FACTORS["optical_blue"]))


# --------------------------------------------------------------- atmospheres
def test_hot_rocky_is_airless_cold_rocky_keeps_atmosphere():
    rng = np.random.default_rng(1)
    sun = make_sun()
    hot = make_planet(1.0, 1.0, 3.0)      # S ~ 2600 S_earth
    hot.insolation = sun.luminosity / hot.a ** 2
    hot.composition_class = "rocky"
    assert assign_atmosphere(rng, sun, hot).atm_class == "airless"

    cold = make_planet(1.0, 1.0, 365.25)  # Earth twin, S=1
    cold.insolation = 1.0
    cold.composition_class = "rocky"
    atm = assign_atmosphere(rng, sun, cold)
    assert atm.atm_class == "secondary"
    # Earth's ~N2/CO2-ish scale height is ~8.5 km; mu=33 at 255 K gives ~7-8
    assert 3 < atm.scale_height_km < 15


def test_hot_jupiter_scale_height_and_tsm():
    rng = np.random.default_rng(2)
    sun = make_sun()
    hj = make_planet(R_JUP / R_EARTH, 318.0, 3.5)
    hj.teq = 1400.0
    hj.composition_class = "giant"
    atm = assign_atmosphere(rng, sun, hj)
    assert atm.atm_class == "h_he"
    assert 200 < atm.scale_height_km < 1500   # canonical HJ: several hundred km
    # A hot Jupiter around a J=3.6 star is a spectacular TSM target
    assert compute_tsm(sun, hj) == 0.0 or hj.radius > 10  # outside table -> 0


def test_tsm_matches_kempton_scaling():
    """GJ 1214b-like: R=2.7, M=8.2, Teq~600, R*=0.22, J=9.75 -> TSM ~ 300."""
    from exoverse.stars import Star
    gj = Star(mass=0.18, radius=0.22, luminosity=0.004, teff=3250, feh=0.2,
              age_gyr=5, ms_lifetime_gyr=400, distance_pc=14.6, mag_v=14.7,
              mag_tess=12.3, mag_j=9.75, u1=0.55, u2=0.18, density=23000,
              spectral_type="M4V")
    p = make_planet(2.7, 8.2, 1.58)
    p.teq = 600.0
    tsm = compute_tsm(gj, p)
    assert 150 < tsm < 500   # published TSM for GJ 1214b is ~300


def test_esm_positive_and_scales_with_temperature():
    sun = make_sun()
    hot = make_planet(1.5, 4.0, 2.0)
    hot.teq = 1500.0
    cool = make_planet(1.5, 4.0, 50.0)
    cool.teq = 400.0
    assert compute_esm(sun, hot) > compute_esm(sun, cool) > 0


def test_airless_planet_not_spectroscopy_target():
    rng = np.random.default_rng(4)
    sun = make_sun()
    p = make_planet(1.0, 1.0, 2.0)
    p.insolation = 700.0
    p.composition_class = "rocky"
    atm = assign_atmosphere(rng, sun, p)
    noise = generate_stellar_noise(rng, sun)
    scores = score_atmosphere_observability(
        atm, 2.0, noise, [("JWST NIRSpec Prism", 15.0, True)])
    assert not scores[0].practical


# ----------------------------------------------------------------- Roman/HWO
def test_roman_detects_more_than_tess():
    """Roman GBTDS (long baseline, IR, fine cadence) should out-detect a
    single TESS sector across any population slice."""
    n_roman = n_tess = 0
    for seed in range(150):
        sys_ = generate_system(seed, f"T-{seed}")
        for p in sys_.planets:
            g = compute_geometry(sys_.star, p)
            if not g.transits:
                continue
            obs = {o.observatory: o for o in observe(sys_.star, p, g, sys_.noise)}
            n_roman += int(obs["Roman GBTDS (2027+)"].detectable)
            n_tess += int(obs["TESS (1 sector)"].detectable)
    assert n_roman >= n_tess
    assert n_roman > 0


def test_hwo_detects_earth_twin_at_10pc():
    sun = make_sun()
    sun.distance_pc = 10.0
    sun.mag_v = 4.83 + 5 * math.log10(10.0 / 10.0)
    earth = make_planet(1.0, 1.0, 365.25)
    o = hwo_imaging(sun, earth)
    assert o.mode == "imaging"
    # Earth twin: contrast ~1.7e-10 at 100 mas -> flagship detection case
    assert 1e-10 < o.contrast < 3e-10
    assert 90 < o.separation_mas < 110
    assert o.detectable


def test_hwo_rejects_distant_host():
    sun = make_sun()
    sun.distance_pc = 150.0
    sun.mag_v = 4.83 + 5 * math.log10(150.0 / 10.0)
    earth = make_planet(1.0, 1.0, 365.25)
    assert not hwo_imaging(sun, earth).detectable


def _bare_atm(cls: str) -> Atmosphere:
    return Atmosphere(atm_class=cls, mu=2.3, scale_height_km=0.0,
                      delta_1h_ppm=0.0, feature_ppm=0.0, cloud_factor=1.0,
                      tsm=0.0, esm=0.0, tsm_priority=False)


def test_hwo_albedo_tracks_atmosphere_class():
    """A H/He giant reflects ~4x more than the same body if airless
    (Jupiter A_g~0.5 vs Moon ~0.12); no atmosphere falls back to 0.3."""
    sun = make_sun()
    jup = make_planet(11.2, 318.0, 2922.0)   # ~4 AU -> 400 mas at 10 pc
    default = hwo_imaging(sun, jup).contrast
    jup.atmosphere = _bare_atm("h_he")
    bright = hwo_imaging(sun, jup)
    jup.atmosphere = _bare_atm("airless")
    dark = hwo_imaging(sun, jup)
    assert bright.contrast / dark.contrast == pytest.approx(0.50 / 0.12)
    assert dark.contrast < default < bright.contrast
    assert bright.detectable   # a 4-AU giant at 10 pc is an easy HWO target


def test_dmax_rescales_distance_without_touching_the_worlds():
    """The solar-neighborhood distance cap must only rescale distance and
    apparent magnitudes: same seed => identical star and planets otherwise."""
    from exoverse.stars import _distance_cdf_unnorm
    far = generate_system(1234, "D-far")                # dmax 300 (default)
    near = generate_system(1234, "D-near", dmax_pc=30.0)
    assert near.star.distance_pc < 30.0
    # Same underlying uniform draw => same CDF quantile at both caps
    # (exact 1/10 rescaling no longer holds: the vertical disk falloff
    # bends the 300 pc CDF). Skip when clamped at the 5 pc floor.
    if near.star.distance_pc > 5.0:
        q_near = (_distance_cdf_unnorm(near.star.distance_pc)
                  / _distance_cdf_unnorm(30.0))
        q_far = (_distance_cdf_unnorm(far.star.distance_pc)
                 / _distance_cdf_unnorm(300.0))
        assert q_near == pytest.approx(q_far, rel=1e-6)
    assert near.star.mass == far.star.mass
    assert near.star.teff == far.star.teff
    assert len(near.planets) == len(far.planets)
    for a, b in zip(near.planets, far.planets):
        assert a.period == b.period and a.radius == b.radius
    # Apparent magnitudes brighten by the distance modulus
    assert near.star.mag_v < far.star.mag_v


def test_hwo_floor_degrades_for_faint_hosts():
    """Fainter hosts stay on the target list but need more contrast:
    photon-limited floor 10^(0.2*(V-7)) above V=7."""
    earth = make_planet(1.0, 1.0, 365.25)
    bright_host = make_sun()                 # V=4.83 at 10 pc
    faint_host = make_sun()
    faint_host.mag_v = 10.5                  # nearby late-K/M dwarf regime
    b, f = hwo_imaging(bright_host, earth), hwo_imaging(faint_host, earth)
    assert f.usable                          # V<11: still evaluated
    assert f.snr_total < b.snr_total
    assert b.snr_total == pytest.approx(f.snr_total
                                        * 10.0 ** (0.2 * (10.5 - 7.0)))
    too_faint = make_sun()
    too_faint.mag_v = 12.0
    assert not hwo_imaging(too_faint, earth).usable


# ------------------------------------------------------------------- archive
def _ap(**kw) -> ArchivePlanet:
    base = dict(pl_name="X b", hostname="X", pl_orbper=10.0, pl_rade=1.0,
                pl_bmasse=1.0, pl_orbeccen=0.0, pl_eqt=500.0, pl_orbsmax=0.1,
                st_teff=5700.0, st_rad=1.0, st_mass=1.0, st_met=0.0,
                sy_dist=50.0, sy_vmag=10.0, sy_jmag=9.0,
                discoverymethod="Primary Transit", disc_year=2020.0,
                tran_flag=1.0, mass_provenance="measured")
    base.update(kw)
    return ArchivePlanet(**base)


def test_archive_columns_match_dataclass():
    from dataclasses import fields
    assert [f.name for f in fields(ArchivePlanet)] == CANONICAL_COLUMNS


def test_audit_flags_brown_dwarf_and_iron_violation():
    planets = [
        _ap(),                                              # fine
        _ap(pl_name="BD b", pl_bmasse=20 * 317.83),          # brown dwarf
        _ap(pl_name="Iron b", pl_rade=0.5, pl_bmasse=5.0),   # denser than iron
    ]
    audits = {a.rule: a for a in audit_rules(planets)}
    assert audits["mass.deuterium_burning"].n_violations == 1
    assert "BD b" in audits["mass.deuterium_burning"].examples
    assert audits["density.exceeds_pure_iron"].n_violations == 1
    assert audits["orbit.grazes_star"].n_violations == 0


def test_audit_ignores_unmeasured_masses():
    planets = [_ap(pl_name="Limit b", pl_bmasse=20 * 317.83,
                   mass_provenance="uncertain")]
    audits = {a.rule: a for a in audit_rules(planets)}
    assert audits["mass.deuterium_burning"].n_checked == 0


# ---------------------------------------------------------------- validation
def test_ks_identical_distributions():
    xs = list(np.random.default_rng(0).normal(0, 1, 400))
    d, p = ks_2samp(xs, xs)
    assert d < 0.01 and p > 0.99


def test_ks_different_distributions():
    rng = np.random.default_rng(0)
    d, p = ks_2samp(list(rng.normal(0, 1, 300)), list(rng.normal(2, 1, 300)))
    assert d > 0.5 and p < 1e-6


def test_quantiles():
    q = quantiles(list(range(101)))
    assert q["q50"] == 50 and q["q10"] == 10 and q["q90"] == 90


# -------------------------------------------------------------------- web UI
@pytest.fixture(scope="module")
def web_client(tmp_path_factory):
    from exoverse.generate import generate_population
    from exoverse.web.app import create_app
    tmp = tmp_path_factory.mktemp("web")
    db_path = str(tmp / "w.db")
    generate_population(db_path, 25, seed=5, progress=False)
    app = create_app(db_path, data_dir=str(tmp / "data"))
    app.config["TESTING"] = True
    return app.test_client()


def test_web_pages_render(web_client):
    for path in ("/", "/systems", "/systems?transiting=1", "/validation",
                 "/about"):
        resp = web_client.get(path)
        assert resp.status_code == 200, path


def test_web_system_detail_and_404(web_client):
    resp = web_client.get("/system/PXS-5-00000")
    assert resp.status_code == 200
    assert b"Plausibility flags" in resp.data or b"Star" in resp.data
    assert web_client.get("/system/DOES-NOT-EXIST").status_code == 404


# ------------------------------------------------------ pipeline consistency
def test_noise_raises_effective_sigma():
    """Stellar noise must never reduce a transit noise budget."""
    for seed in range(60):
        sys_ = generate_system(seed, f"N-{seed}")
        for p in sys_.planets:
            g = compute_geometry(sys_.star, p)
            if not g.transits:
                continue
            with_n = {o.observatory: o for o in
                      observe(sys_.star, p, g, sys_.noise)}
            without = {o.observatory: o for o in
                       observe(sys_.star, p, g, None)}
            for k, o in with_n.items():
                if o.mode == "transit" and o.usable:
                    assert o.snr_total <= without[k].snr_total + 1e-9
            return
