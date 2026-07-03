"""Tests for phase 3: realistic albedo distributions, disk-structure
distance sampling with absolute normalization, and system-architecture
(structure) comparisons against the real catalog."""

import math

import numpy as np
import pytest

from exoverse.atmospheres import sample_geometric_albedo
from exoverse.observatories import hwo_imaging
from exoverse.stars import (
    DISK_SCALE_HEIGHT_PC, LOCAL_MS_SYSTEM_DENSITY_PC3,
    _distance_cdf_unnorm, expected_systems_within, sample_distance,
)
from exoverse.validate import (
    _adjacent_pairs, _desert_fraction, _multiplicity_hist, compare_structure,
)

from test_phase2 import _ap
from test_physics import make_planet, make_sun


# ------------------------------------------------------------------- albedos
def _draws(atm_class, teq, n=800, seed=0):
    rng = np.random.default_rng(seed)
    return np.array([sample_geometric_albedo(rng, atm_class, teq)
                     for _ in range(n)])


def test_hot_jupiters_are_dark():
    """Kepler hot-Jupiter population: A_g ~ 0.03-0.11 median."""
    a = _draws("h_he", 1100.0)
    assert 0.03 < np.median(a) < 0.12
    assert a.min() >= 0.02


def test_cold_giants_are_bright():
    """Jupiter/Saturn regime: A_g ~ 0.5."""
    a = _draws("h_he", 110.0)
    assert 0.40 < np.median(a) < 0.62


def test_ultra_hot_silicate_clouds_brighten():
    """Kepler-7b class: ultra-hot giants recover moderate albedo."""
    assert np.median(_draws("h_he", 2200.0)) > np.median(_draws("h_he", 1000.0))


def test_secondary_has_venus_tail():
    """~22% of secondary atmospheres draw a Venus-like bright cloud deck."""
    a = _draws("secondary", 300.0)
    frac_bright = float(np.mean(a > 0.45))
    assert 0.12 < frac_bright < 0.32
    # the dark branch centers near Earth's measured 0.24
    assert 0.15 < np.median(a[a <= 0.45]) < 0.33


def test_airless_mostly_dark_with_icy_branch_when_cold():
    hot = _draws("airless", 800.0)
    assert hot.max() <= 0.30                    # no ice on hot rocks
    cold = _draws("airless", 150.0)
    assert float(np.mean(cold > 0.30)) > 0.05   # icy-bright branch exists
    assert np.median(cold) < 0.20               # but regolith-dark dominates


def test_all_albedos_physical():
    for cls in ("h_he", "h_he_rich", "steam", "secondary", "airless"):
        for teq in (100.0, 700.0, 2000.0):
            a = _draws(cls, teq, n=300, seed=7)
            assert a.min() > 0.0 and a.max() <= 0.9


def test_draw_count_is_uniform_across_classes():
    """Every albedo draw must consume exactly two rng values regardless of
    class/branch, or downstream worlds lose reproducibility."""
    for cls in ("h_he", "steam", "secondary", "airless"):
        r1 = np.random.default_rng(99)
        sample_geometric_albedo(r1, cls, 150.0)
        r2 = np.random.default_rng(99)
        r2.random(); r2.normal()
        assert r1.random() == r2.random()


def test_hwo_uses_stored_per_planet_albedo():
    from exoverse.atmospheres import Atmosphere
    sun = make_sun()
    earth = make_planet(1.0, 1.0, 365.25)
    earth.atmosphere = Atmosphere(
        atm_class="secondary", mu=33.0, scale_height_km=8.0, delta_1h_ppm=0.0,
        feature_ppm=0.0, cloud_factor=1.0, tsm=0.0, esm=0.0,
        tsm_priority=False, geometric_albedo=0.65)
    bright = hwo_imaging(sun, earth).contrast
    earth.atmosphere.geometric_albedo = 0.13
    dark = hwo_imaging(sun, earth).contrast
    assert bright / dark == pytest.approx(0.65 / 0.13)


# ------------------------------------------- disk structure / normalization
def test_expected_systems_match_census():
    """~7,100 systems within 30 pc (RECONS-scaled), ~5.2M within 300 pc."""
    assert 6300 < expected_systems_within(30.0) < 7900
    assert 4.4e6 < expected_systems_within(300.0) < 6.0e6


def test_vertical_falloff_reduces_far_counts():
    """At d = h the sphere-averaged density is ~71% of midplane."""
    uniform = (4.0 * math.pi / 3.0) * LOCAL_MS_SYSTEM_DENSITY_PC3 * 300.0 ** 3
    ratio = expected_systems_within(300.0) / uniform
    assert 0.68 < ratio < 0.73
    # ...but a 30 pc sphere is essentially uniform (< 5% correction)
    uniform30 = (4.0 * math.pi / 3.0) * LOCAL_MS_SYSTEM_DENSITY_PC3 * 30.0 ** 3
    assert expected_systems_within(30.0) / uniform30 > 0.95


def test_sample_distance_inverts_cdf():
    assert sample_distance(1.0, 300.0) == pytest.approx(300.0, abs=1e-3)
    assert sample_distance(0.0, 300.0) == pytest.approx(0.0, abs=1e-3)
    # monotone, and biased nearer than the uniform-density draw at 300 pc
    ds = [sample_distance(u, 300.0) for u in (0.2, 0.5, 0.8)]
    assert ds == sorted(ds)
    assert sample_distance(0.5, 300.0) < 300.0 * 0.5 ** (1.0 / 3.0)
    # quantile consistency with the unnormalized CDF
    d = sample_distance(0.37, 300.0)
    assert (_distance_cdf_unnorm(d) / _distance_cdf_unnorm(300.0)
            == pytest.approx(0.37, rel=1e-4))
    assert DISK_SCALE_HEIGHT_PC == 300.0


# --------------------------------------------------- structure comparisons
def test_adjacent_pairs_and_multiplicity():
    groups = {
        "a": [(10.0, 1.0), (20.0, 2.0), (45.0, 1.0)],
        "b": [(5.0, 3.0)],
    }
    pr, rr = _adjacent_pairs(groups)
    assert sorted(pr) == [2.0, 2.25]
    assert rr == [abs(math.log10(2.0)), abs(math.log10(0.5))]
    assert _multiplicity_hist(groups) == {"1": 1, "2": 0, "3+": 1}


def test_desert_fraction():
    hits, tot = _desert_fraction([(2.0, 5.0), (2.0, 1.0), (10.0, 5.0),
                                  (None, 5.0)])
    assert (hits, tot) == (1, 3)


def test_compare_structure_runs_on_small_population(tmp_path):
    from exoverse.generate import generate_population
    db_path = str(tmp_path / "s.db")
    generate_population(db_path, 150, seed=9, progress=False)
    real = [
        _ap(pl_name="K-1 b", hostname="K-1", pl_orbper=5.0, pl_rade=1.5),
        _ap(pl_name="K-1 c", hostname="K-1", pl_orbper=10.0, pl_rade=1.6),
        _ap(pl_name="K-2 b", hostname="K-2", pl_orbper=3.0, pl_rade=5.0),
    ]
    out = {c.metric: c for c in compare_structure(db_path, real)}
    assert set(out) == {"peas_in_a_pod_radius_uniformity",
                        "adjacent_period_ratios", "transit_multiplicity",
                        "hot_neptune_desert_occupancy"}
    assert out["transit_multiplicity"].real == {"1": 1, "2": 1, "3+": 0,
                                                "singles_per_multi": 1.0}
    assert out["hot_neptune_desert_occupancy"].real["in_desert"] == 1
