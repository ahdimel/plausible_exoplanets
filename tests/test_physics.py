"""Physics sanity tests: the generator must reproduce known real-world values."""

import math

import numpy as np
import pytest

from exoverse.constants import AU, M_SUN, R_EARTH, R_JUP, R_SUN
from exoverse.flags import Severity
from exoverse.planets import (
    Planet, equilibrium_temp, habitable_zone_au, pure_iron_radius,
    semimajor_axis_au,
)
from exoverse.stars import (
    Star, bolometric_correction_v, generate_star, ms_luminosity, ms_radius,
    sample_kroupa_mass, spectral_type_from_teff,
)
from exoverse.system import generate_system, mutual_hill_delta
from exoverse.transits import compute_geometry, transit_flux


def make_sun() -> Star:
    return Star(mass=1.0, radius=1.0, luminosity=1.0, teff=5772.0, feh=0.0,
                age_gyr=4.6, ms_lifetime_gyr=10.0, distance_pc=10.0,
                mag_v=4.83, mag_tess=4.17, mag_j=3.63, u1=0.41, u2=0.26,
                density=1408.0, spectral_type="G2V")


def make_planet(radius_re, mass_me, period_d, mstar=1.0, ecc=0.0, inc=90.0):
    a = semimajor_axis_au(period_d, mstar)
    return Planet(radius=radius_re, mass=mass_me, density=1.0,
                  composition_class="test", period=period_d, a=a, ecc=ecc,
                  omega_deg=90.0, inc_deg=inc, t0_frac=0.0, teq=300.0,
                  insolation=1.0, in_habitable_zone=False)


# ---------------------------------------------------------------- stars
def test_sun_like_scaling_relations():
    assert ms_luminosity(1.0) == pytest.approx(1.0)
    assert ms_radius(1.0) == pytest.approx(1.0)


def test_bolometric_correction_sun():
    # Flower/Torres BC_V for the Sun should be ~ -0.07 mag
    assert bolometric_correction_v(5772.0) == pytest.approx(-0.07, abs=0.05)


def test_spectral_types():
    assert spectral_type_from_teff(5772.0).startswith("G")
    assert spectral_type_from_teff(3200.0).startswith("M")
    assert spectral_type_from_teff(6500.0).startswith("F")


def test_kroupa_masses_in_range():
    rng = np.random.default_rng(1)
    m = np.array([sample_kroupa_mass(rng) for _ in range(2000)])
    assert m.min() >= 0.08 and m.max() <= 2.2
    # IMF is bottom-heavy: median well below 0.5 Msun
    assert np.median(m) < 0.5


def test_generated_stars_valid():
    rng = np.random.default_rng(7)
    for _ in range(200):
        s = generate_star(rng)
        if s.is_invalid:
            continue
        assert 2300 < s.teff < 8500
        assert s.ms_lifetime_gyr > 0
        assert s.mag_v > s.mag_j  # stars are brighter in J than V


# ---------------------------------------------------------------- orbits
def test_kepler_third_law_earth():
    assert semimajor_axis_au(365.25, 1.0) == pytest.approx(1.0, rel=1e-3)


def test_equilibrium_temp_earth():
    sun = make_sun()
    # Earth: Teq ~ 255 K with A=0.3
    assert equilibrium_temp(sun, 1.0) == pytest.approx(255.0, abs=8.0)


def test_habitable_zone_contains_earth():
    hz_in, hz_out = habitable_zone_au(make_sun())
    assert hz_in < 1.0 < hz_out
    assert 0.8 < hz_in < 1.0   # conservative inner edge ~0.95 AU
    assert 1.3 < hz_out < 2.0  # conservative outer edge ~1.67 AU


def test_pure_iron_limit():
    # Earth (M=1, R=1) must be less dense than pure iron sphere of same mass
    assert pure_iron_radius(1.0) < 1.0


# ---------------------------------------------------------------- transits
def test_uniform_depth_jupiter():
    # No limb darkening: depth = (Rp/R*)^2; Jupiter/Sun ~ 1.06%
    k = R_JUP / R_SUN
    flux = transit_flux(np.array([0.0]), k, 0.0, 0.0)
    assert 1.0 - flux[0] == pytest.approx(k * k, rel=2e-3)


def test_earth_transit_depth_and_duration():
    sun, earth = make_sun(), make_planet(1.0, 1.0, 365.25)
    g = compute_geometry(sun, earth)
    assert g.transits
    # Geometric depth ~84 ppm; limb darkening boosts central depth ~10-30%
    assert g.depth_uniform_ppm == pytest.approx(84.0, abs=3.0)
    assert 84.0 < g.depth_ppm < 120.0
    # Central Earth transit lasts ~13 hours
    assert g.t14_hours == pytest.approx(13.0, abs=0.7)


def test_hot_jupiter_depth():
    sun = make_sun()
    hj = make_planet(R_JUP / R_EARTH, 318.0, 3.5)
    g = compute_geometry(sun, hj)
    assert g.transits
    assert 9000 < g.depth_uniform_ppm < 12000  # ~1%


def test_limb_darkening_shape():
    # LD makes central transits deeper and grazing transits shallower than uniform
    k = 0.1
    z = np.array([0.0, 0.95])
    ld = transit_flux(z, k, 0.5, 0.2)
    uni = transit_flux(z, k, 0.0, 0.0)
    assert (1 - ld[0]) > (1 - uni[0])   # deeper at center
    assert (1 - ld[1]) < (1 - uni[1])   # shallower at limb


def test_no_transit_when_face_on():
    sun = make_sun()
    p = make_planet(1.0, 1.0, 365.25, inc=45.0)
    g = compute_geometry(sun, p)
    assert not g.transits
    assert g.prob_transit < 0.01  # Earth-Sun geometric prob ~0.47%


def test_flux_out_of_transit_is_unity():
    f = transit_flux(np.array([2.0, 5.0]), 0.1, 0.4, 0.25)
    assert np.allclose(f, 1.0)


# ---------------------------------------------------------------- stability
def test_mutual_hill_solar_system():
    sun = make_sun()
    earth = make_planet(1.0, 1.0, 365.25)
    venus = make_planet(0.95, 0.815, 224.7)
    delta = mutual_hill_delta(sun, venus, earth)
    assert delta > 25  # Venus-Earth are very widely separated dynamically


def test_generated_systems_are_stable():
    for seed in range(30):
        sys_ = generate_system(seed, f"T-{seed}")
        # No INVALID flags may survive generation
        for p in sys_.planets:
            assert not p.is_invalid
        # Adjacent pairs obey the Gladman criterion
        for a, b in zip(sys_.planets, sys_.planets[1:]):
            assert mutual_hill_delta(sys_.star, a, b) > 2 * math.sqrt(3) - 1e-9
            assert b.a * (1 - b.ecc) > a.a * (1 + a.ecc)


def test_generation_deterministic():
    s1 = generate_system(123, "A")
    s2 = generate_system(123, "A")
    assert s1.star.mass == s2.star.mass
    assert len(s1.planets) == len(s2.planets)
    for p, q in zip(s1.planets, s2.planets):
        assert p.period == q.period and p.radius == q.radius
