"""Architecture-knob tests (sigma_r radius copula, sigma_i dispersion):
bit-for-bit default neutrality, exact marginal preservation, inclination
limits, and determinism. See docs/sigma_r_note.md for the copula math."""

import math

import numpy as np
import pytest

from exoverse.architecture import (
    Architecture, DEFAULT_ARCH, correlated_quantile, normal_cdf,
)
from exoverse.planets import (
    sample_small_radius, small_radius_cdf, small_radius_ppf,
)
from exoverse.system import generate_system


def systems_equal(s1, s2) -> bool:
    if (s1.star.mass != s2.star.mass or s1.sys_inc_deg != s2.sys_inc_deg
            or len(s1.planets) != len(s2.planets)):
        return False
    return all(p.period == q.period and p.radius == q.radius
               and p.mass == q.mass and p.ecc == q.ecc
               and p.inc_deg == q.inc_deg and p.t0_frac == q.t0_frac
               for p, q in zip(s1.planets, s2.planets))


# ------------------------------------------------------------ defaults
def test_default_arch_bit_for_bit():
    """Architecture() must consume the identical rng sequence as arch=None:
    default populations reproduce the post-geometry-fix baseline exactly."""
    for seed in range(40):
        assert systems_equal(generate_system(seed, "A"),
                             generate_system(seed, "A", arch=Architecture()))


def test_arch_determinism():
    arch = Architecture(sigma_r=0.3, sigma_i=5.0, f_hot=0.4, sigma_i_hot=30.0)
    for seed in (7, 99):
        assert systems_equal(generate_system(seed, "A", arch=arch),
                             generate_system(seed, "A", arch=arch))


def test_arch_meta_roundtrip():
    arch = Architecture(sigma_r=0.25, sigma_i=3.0, f_hot=0.1,
                        sigma_i_hot=45.0, isotropic=False)
    meta = arch.meta_items()
    assert Architecture.from_meta(lambda k, d=None: meta.get(k, d)) == arch
    # Old DBs without arch keys must yield the default
    assert Architecture.from_meta(lambda k, d=None: d) == DEFAULT_ARCH


# ------------------------------------------------- radius marginal (sigma_r)
def test_small_radius_ppf_cdf_roundtrip():
    for period in (1.0, 10.0, 200.0):
        for u in np.linspace(0.001, 0.999, 41):
            r = small_radius_ppf(float(u), period)
            if 0.4 < r < 4.0:   # off the clip atoms the inverse is exact
                assert small_radius_cdf(r, period) == pytest.approx(
                    float(u), abs=1e-9)


def test_small_radius_marginal_matches_sampler():
    """ppf(U) must reproduce sample_small_radius's distribution exactly:
    KS between the two at fixed period, plus matching clip-atom masses."""
    rng = np.random.default_rng(11)
    period = 10.0
    direct = np.sort([sample_small_radius(rng, period) for _ in range(4000)])
    via_ppf = np.sort([small_radius_ppf(float(u), period)
                       for u in rng.uniform(0.0, 1.0, 4000)])
    # Two-sample KS by hand (repo has no scipy)
    grid = np.linspace(0.4, 4.0, 200)
    d = max(abs(np.searchsorted(direct, g, "right") / 4000
                - np.searchsorted(via_ppf, g, "right") / 4000) for g in grid)
    assert d < 0.04   # ~2x the n=4000 KS 1% critical value


def test_marginal_radius_invariant_under_sigma_r():
    """The one-point radius distribution must not move as sigma_r varies:
    correlated quantiles are still marginally uniform."""
    rng = np.random.default_rng(5)
    period = 10.0
    pooled = {}
    for sigma_r in (0.05, 0.5, 2.0):
        radii = []
        for _ in range(1500):
            z_sys = float(rng.standard_normal())
            radii.append(small_radius_ppf(
                correlated_quantile(rng, z_sys, sigma_r), period))
        pooled[sigma_r] = np.sort(radii)
    grid = np.linspace(0.4, 4.0, 200)
    for a in (0.05, 0.5):
        d = max(abs(np.searchsorted(pooled[a], g, "right") / 1500
                    - np.searchsorted(pooled[2.0], g, "right") / 1500)
                for g in grid)
        assert d < 0.07


def test_sigma_r_uniformity_monotonic():
    """Sibling |dlogR| dispersion must shrink as sigma_r -> 0."""
    rng = np.random.default_rng(17)
    med = {}
    for sigma_r in (0.05, 0.5, 5.0):
        deltas = []
        for _ in range(2000):
            z_sys = float(rng.standard_normal())
            r1 = small_radius_ppf(correlated_quantile(rng, z_sys, sigma_r), 8.0)
            r2 = small_radius_ppf(correlated_quantile(rng, z_sys, sigma_r), 20.0)
            deltas.append(abs(math.log10(r2 / r1)))
        med[sigma_r] = float(np.median(deltas))
    assert med[0.05] < med[0.5] < med[5.0]
    assert med[0.05] < 0.3 * med[5.0]   # near-peas limit is much tighter


def test_correlated_quantile_marginally_uniform():
    rng = np.random.default_rng(23)
    us = np.sort([correlated_quantile(rng, float(rng.standard_normal()), 0.3)
                  for _ in range(4000)])
    d = max(abs(np.searchsorted(us, g, "right") / 4000 - g)
            for g in np.linspace(0.01, 0.99, 99))
    assert d < 0.03
    assert normal_cdf(0.0) == pytest.approx(0.5)


# ------------------------------------------------- inclinations (sigma_i)
def test_sigma_i_coplanar_limit():
    """sigma_i=0 must put every planet exactly in the system plane."""
    arch = Architecture(sigma_i=0.0)
    found_multi = False
    for seed in range(60):
        sys_ = generate_system(seed, "C", arch=arch)
        for p in sys_.planets:
            assert p.inc_deg == pytest.approx(sys_.sys_inc_deg, abs=1e-9)
        found_multi = found_multi or len(sys_.planets) >= 2
    assert found_multi


def test_isotropic_limit():
    """isotropic=True: cos(i) uniform on [0, 1] and siblings independent."""
    arch = Architecture(isotropic=True)
    cosi, pair_diffs = [], []
    for seed in range(700):
        sys_ = generate_system(seed, "ISO", arch=arch)
        for p in sys_.planets:
            cosi.append(math.cos(math.radians(p.inc_deg)))
        for a, b in zip(sys_.planets, sys_.planets[1:]):
            pair_diffs.append(abs(a.inc_deg - b.inc_deg))
    cosi = np.sort(cosi)
    d = max(abs(np.searchsorted(cosi, g, "right") / len(cosi) - g)
            for g in np.linspace(0.02, 0.98, 49))
    assert d < 0.05
    # Independent isotropic siblings differ by ~30 deg on median; a shared
    # plane with Rayleigh(1.5) would put this at ~2 deg
    assert float(np.median(pair_diffs)) > 10.0


def test_sigma_i_scales_planet_scatter():
    """Pooled |inc - sys_inc| grows with sigma_i."""
    def pooled_offsets(sigma_i):
        arch = Architecture(sigma_i=sigma_i)
        out = []
        for seed in range(150):
            sys_ = generate_system(seed, "S", arch=arch)
            out.extend(abs(p.inc_deg - sys_.sys_inc_deg)
                       for p in sys_.planets)
        return float(np.median(out))
    m1, m5, m20 = (pooled_offsets(s) for s in (1.0, 5.0, 20.0))
    assert m1 < m5 < m20
    assert m5 == pytest.approx(5.0 * m1, rel=0.25)   # Rayleigh medians scale


def test_hot_mixture_two_populations():
    """f_hot systems draw the hot Rayleigh: per-system max offset separates
    the two components cleanly for sigma_i_hot >> sigma_i."""
    arch = Architecture(sigma_i=1.0, f_hot=0.5, sigma_i_hot=40.0)
    n_hot, n_multi = 0, 0
    for seed in range(400):
        sys_ = generate_system(seed, "H", arch=arch)
        if len(sys_.planets) < 1:
            continue
        n_multi += 1
        if max(abs(p.inc_deg - sys_.sys_inc_deg)
               for p in sys_.planets) > 8.0:
            n_hot += 1
    frac = n_hot / n_multi
    assert 0.3 < frac < 0.7   # ~f_hot (binomial scatter + Rayleigh tails)