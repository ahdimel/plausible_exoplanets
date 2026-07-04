"""Tests for the dichotomy observables + distance metric (dichotomy.py).

All catalogs are built in code — no network, no data/ snapshots."""

import math
import random

import pytest

from exoverse.dichotomy import (
    Observables, combined_distance, multiplicity_distance,
    observables_from_koi, observables_from_universe, size_distance, spearman,
)
from exoverse.kepler_field import UniverseResult


def make_obs(n_k, dlogr=(), n_targets=None):
    counts = {k: 0 for k in (1, 2, 3, 4, 5, "6+")}
    counts.update(n_k)
    return Observables(n_k=counts, dlogr=sorted(dlogr), monotonicity=None,
                       n_systems=sum(counts.values()),
                       n_planets=sum(k * v for k, v in counts.items()
                                     if k != "6+"),
                       n_targets=n_targets)


# ----------------------------------------------------------------- builders
GROUPS = [
    [(3.0, 1.0)],                                   # single
    [(5.0, 2.0), (11.0, 2.5)],                      # pair
    [(4.0, 1.0), (9.0, 1.5), (20.0, 2.0)],          # monotone triple
    [(2.0, 1.0), (6.0, 35.0), (15.0, 1.2)],         # one radius >= 30 Re
]


def test_builders_share_one_code_path():
    """UniverseResult and KOI-dict builders must agree exactly on the same
    grouped catalog — detection bias cancels only if both sides are
    computed identically."""
    res = UniverseResult(n_k={1: 2, 2: 1, 3: 2, 4: 0, 5: 0, "6+": 0},
                         detected=[list(g) for g in GROUPS],
                         n_targets=100,
                         n_detected_planets=9)
    obs_syn = observables_from_universe(res)
    obs_real = observables_from_koi({i: g for i, g in enumerate(GROUPS)},
                                    n_targets=100)
    assert obs_syn == obs_real


def test_observables_statistics():
    obs = observables_from_koi({i: g for i, g in enumerate(GROUPS)})
    assert obs.n_k[1] == 1 and obs.n_k[2] == 1 and obs.n_k[3] == 2
    assert obs.n_systems == 4 and obs.n_planets == 9
    # dlogr: pair (2.0->2.5), triple (1->1.5, 1.5->2); the 35 Re planet's two
    # adjacencies are excluded (both radii must be < 30 Re)
    expected = sorted([abs(math.log10(2.5 / 2.0)), abs(math.log10(1.5)),
                       abs(math.log10(2.0 / 1.5))])
    assert obs.dlogr == pytest.approx(expected)
    # monotonicity: mean Spearman over the two >= 3-planet systems:
    # perfectly increasing (+1) and 1.0 -> 35 -> 1.2 i.e. ranks 1,3,2 (+0.5)
    assert obs.monotonicity == pytest.approx((1.0 + 0.5) / 2.0)


def test_monotonicity_none_when_no_triples():
    obs = observables_from_koi({0: [(3.0, 1.0)], 1: [(5.0, 2.0), (9.0, 2.2)]})
    assert obs.monotonicity is None


def test_koi_builder_accepts_objects_with_attrs():
    class KOI:
        def __init__(self, period, prad):
            self.period, self.prad = period, prad
    obs_t = observables_from_koi({7: [(3.0, 1.0), (8.0, 2.0)]})
    obs_o = observables_from_koi({7: [KOI(3.0, 1.0), KOI(8.0, 2.0)]})
    assert obs_t == obs_o


# ----------------------------------------------------------------- spearman
def test_spearman_monotone_and_ties():
    assert spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)
    assert spearman([1, 2, 3, 4], [40, 30, 20, 10]) == pytest.approx(-1.0)
    # Tied values get average ranks: y ranks (1.5, 1.5, 3.5, 3.5)
    assert spearman([1, 2, 3, 4], [5, 5, 9, 9]) == pytest.approx(
        4.0 / math.sqrt(5.0 * 4.0))
    assert math.isnan(spearman([1, 2, 3], [7, 7, 7]))   # constant variable


# ------------------------------------------------------------- multiplicity
def test_multinomial_zero_for_identical_shapes():
    real = make_obs({1: 100, 2: 30, 3: 8})
    same_shape_bigger = make_obs({1: 1000, 2: 300, 3: 80})
    d = multiplicity_distance(same_shape_bigger, real, mode="multinomial")
    assert 0.0 <= d < 0.02   # ~0 up to the Laplace-smoothing residual


def test_multinomial_grows_with_shape_difference():
    real = make_obs({1: 100, 2: 30, 3: 8})
    near = make_obs({1: 95, 2: 35, 3: 8})
    far = make_obs({1: 135, 2: 3, 3: 0})    # dichotomy-like excess of singles
    d_near = multiplicity_distance(near, real)
    d_far = multiplicity_distance(far, real)
    assert d_near < d_far
    assert d_far > 0.2


def test_multinomial_is_shape_only():
    """Scaling the synthetic counts must not change the distance: the fixed
    planets-per-star rate cannot contaminate the fit. (Only asymptotically
    exact: the +0.5 Laplace smoothing fades with count size.)"""
    real = make_obs({1: 100, 2: 30, 3: 8})
    syn = make_obs({1: 600, 2: 250, 3: 40})
    syn10 = make_obs({1: 6000, 2: 2500, 3: 400})
    assert multiplicity_distance(syn, real) == pytest.approx(
        multiplicity_distance(syn10, real), abs=0.01)


def test_poisson_mode_uses_absolute_rates():
    real = make_obs({1: 100, 2: 30}, n_targets=10000)
    syn_match = make_obs({1: 100, 2: 30}, n_targets=10000)
    syn_low = make_obs({1: 50, 2: 15}, n_targets=10000)   # same shape, half rate
    d_match = multiplicity_distance(syn_match, real, mode="poisson")
    d_low = multiplicity_distance(syn_low, real, mode="poisson")
    assert d_match < d_low   # multinomial mode would call these equal
    with pytest.raises(ValueError):
        multiplicity_distance(make_obs({1: 5}), real, mode="poisson")
    with pytest.raises(ValueError):
        multiplicity_distance(syn_match, real, mode="nope")


# --------------------------------------------------------------------- size
def _gauss(n, mu, seed):
    rng = random.Random(seed)
    return [rng.gauss(mu, 1.0) for _ in range(n)]


def test_size_distance_ks_sanity():
    real = make_obs({1: 1}, dlogr=_gauss(400, 0.0, 1))
    same = make_obs({1: 1}, dlogr=_gauss(400, 0.0, 2))
    shifted = make_obs({1: 1}, dlogr=_gauss(400, 1.0, 3))
    d_same = size_distance(same, real, mode="ks")
    d_shift = size_distance(shifted, real, mode="ks")
    assert 0.0 <= d_same < 0.12       # ~KS null fluctuation at n=400
    assert d_shift > 0.3
    assert d_shift > d_same


def test_size_distance_ad_sanity():
    real = make_obs({1: 1}, dlogr=_gauss(400, 0.0, 4))
    same = make_obs({1: 1}, dlogr=_gauss(400, 0.0, 5))
    shifted = make_obs({1: 1}, dlogr=_gauss(400, 1.0, 6))
    t_same = size_distance(same, real, mode="ad")
    t_shift = size_distance(shifted, real, mode="ad")
    assert abs(t_same) < 3.0          # standardized: ~0 under the null
    assert t_shift > 10.0             # gross shift blows past the null scale
    with pytest.raises(ValueError):
        size_distance(same, real, mode="nope")


# ----------------------------------------------------------------- combined
def test_combined_distance_components_and_sum():
    real = make_obs({1: 100, 2: 30, 3: 8}, dlogr=_gauss(200, 0.0, 7))
    syn = make_obs({1: 140, 2: 10, 3: 1}, dlogr=_gauss(200, 0.5, 8))
    out = combined_distance(syn, real, w_mult=2.0, w_size=0.5)
    assert out["multiplicity"] == pytest.approx(
        multiplicity_distance(syn, real))
    assert out["size"] == pytest.approx(size_distance(syn, real))
    assert out["total"] == pytest.approx(
        2.0 * out["multiplicity"] + 0.5 * out["size"])
    # Perfect match: both components (hence the sum) are ~0
    self_out = combined_distance(real, real)
    assert self_out["total"] == pytest.approx(0.0, abs=0.02)
