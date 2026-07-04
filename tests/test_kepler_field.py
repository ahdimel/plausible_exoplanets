"""Tests for the DR25-conditioned Kepler-field simulation (kepler_field.py).

Targets are fake, built in code (duck-typed against
kepler_data.KeplerTarget) — no network, no data/ snapshots needed."""

import math
from dataclasses import dataclass, field
from typing import Tuple

import numpy as np
import pytest

from exoverse.architecture import Architecture
from exoverse.constants import RHO_SUN
from exoverse.kepler_field import (
    CDPP_DURATIONS_HR, MIN_TRANSITS, SNR_THRESHOLD, UniverseResult,
    detected_planets, interp_cdpp, simulate_universe, star_from_target,
)
from exoverse.planets import Planet


@dataclass
class FakeTarget:
    """Duck-typed stand-in for kepler_data.KeplerTarget."""
    kepid: int
    teff: float = 5772.0
    logg: float = 4.44
    feh: float = 0.0
    radius: float = 1.0
    mass: float = 1.0
    kepmag: float = 12.0
    dutycycle: float = 0.9
    dataspan: float = 1400.0
    cdpp_ppm: Tuple[float, ...] = field(default=tuple([30.0] * 14))


def sun_like_targets(n: int) -> list:
    return [FakeTarget(kepid=i) for i in range(n)]


def make_planet(radius: float, period: float, inc_deg: float = 90.0,
                mass_star: float = 1.0) -> Planet:
    """Minimal circular edge-on-capable planet for detection tests."""
    a = (period / 365.25) ** (2.0 / 3.0) * mass_star ** (1.0 / 3.0)
    return Planet(radius=radius, mass=radius ** 2.06, density=5.0,
                  composition_class="test", period=period, a=a, ecc=0.0,
                  omega_deg=90.0, inc_deg=inc_deg, t0_frac=0.0,
                  teq=500.0, insolation=1.0, in_habitable_zone=False)


# ------------------------------------------------------- star construction
def test_star_from_target_physics():
    star = star_from_target(FakeTarget(kepid=1, teff=5772.0, radius=1.0,
                                       mass=1.0, feh=0.1))
    assert star.luminosity == pytest.approx(1.0, rel=1e-6)
    assert star.density == pytest.approx(RHO_SUN, rel=1e-9)
    assert star.feh == 0.1
    # Limb darkening must come from the same Teff tables generate_star uses
    assert 0.3 < star.u1 < 0.5 and 0.2 < star.u2 < 0.3
    # Cooler, smaller star: L = R^2 (T/5772)^4
    star_k = star_from_target(FakeTarget(kepid=2, teff=4500.0, radius=0.7,
                                         mass=0.75))
    assert star_k.luminosity == pytest.approx(
        0.7 ** 2 * (4500.0 / 5772.0) ** 4, rel=1e-6)
    assert star_k.density == pytest.approx(RHO_SUN * 0.75 / 0.7 ** 3)


# --------------------------------------------------------- CDPP interpolation
def test_interp_cdpp_endpoints_clamp():
    cdpp = tuple(float(i) for i in range(14))   # 0..13 over 1.5..15 h
    assert interp_cdpp(cdpp, 0.5) == 1.5 * 0.0 + cdpp[0]     # below -> first
    assert interp_cdpp(cdpp, CDPP_DURATIONS_HR[0]) == cdpp[0]
    assert interp_cdpp(cdpp, 30.0) == cdpp[-1]               # above -> last
    assert interp_cdpp(cdpp, CDPP_DURATIONS_HR[-1]) == cdpp[-1]


def test_interp_cdpp_interior_linear():
    cdpp = tuple(float(i) for i in range(14))
    # Halfway between the 1.5 h and 2.0 h columns
    assert interp_cdpp(cdpp, 1.75) == pytest.approx(0.5)
    # Exactly on a grid point
    assert interp_cdpp(cdpp, 6.0) == pytest.approx(7.0)


# ----------------------------------------------------------- detection rule
def test_deep_planet_on_quiet_star_detected():
    target = FakeTarget(kepid=1, cdpp_ppm=tuple([30.0] * 14))
    star = star_from_target(target)
    planet = make_planet(radius=10.0, period=3.0)   # ~8400 ppm, 420 transits
    assert detected_planets(star, [planet], target) == [planet]


def test_shallow_planet_on_noisy_star_missed():
    target = FakeTarget(kepid=2, cdpp_ppm=tuple([500.0] * 14))
    star = star_from_target(target)
    planet = make_planet(radius=1.0, period=100.0)  # ~84 ppm, ~13 transits
    assert detected_planets(star, [planet], target) == []


def test_min_transit_count_enforced():
    # 1400 d * 0.9 / 600 d = 2.1 transits < 3: rejected however deep
    target = FakeTarget(kepid=3, cdpp_ppm=tuple([1.0] * 14))
    star = star_from_target(target)
    planet = make_planet(radius=10.0, period=600.0)
    assert detected_planets(star, [planet], target) == []
    # Same planet at 400 d gives 3.15 transits and a huge SNR: detected
    ok = make_planet(radius=10.0, period=400.0)
    assert detected_planets(star, [ok], target) == [ok]


def test_non_transiting_planet_missed():
    target = FakeTarget(kepid=4)
    star = star_from_target(target)
    planet = make_planet(radius=10.0, period=3.0, inc_deg=45.0)
    assert detected_planets(star, [planet], target) == []


# ---------------------------------------------------------- universe results
def test_simulate_universe_deterministic():
    targets = sun_like_targets(200)
    r1 = simulate_universe(targets, seed=42)
    r2 = simulate_universe(targets, seed=42)
    assert r1.n_k == r2.n_k
    assert r1.detected == r2.detected
    assert r1.n_detected_planets == r2.n_detected_planets
    assert r1.n_targets == 200
    # A different seed draws different universes
    r3 = simulate_universe(targets, seed=43)
    assert r1.detected != r3.detected


def test_universe_result_shape():
    res = simulate_universe(sun_like_targets(400), seed=7)
    assert isinstance(res, UniverseResult)
    assert set(res.n_k) == {1, 2, 3, 4, 5, "6+"}
    assert sum(res.n_k.values()) == len(res.detected)
    assert res.n_detected_planets == sum(len(g) for g in res.detected)
    for group in res.detected:
        assert len(group) >= 1
        periods = [p for p, _ in group]
        assert periods == sorted(periods)   # adjacency for pair statistics
        assert all(r > 0 for _, r in group)


def test_hot_universe_boosts_singles_per_multi():
    """The dichotomy axis: large mutual inclinations destroy transit multis
    on the same targets/seed, raising the singles-per-multi ratio."""
    targets = sun_like_targets(3000)
    cold = simulate_universe(targets, seed=11, arch=Architecture(sigma_i=1.0))
    hot = simulate_universe(targets, seed=11, arch=Architecture(sigma_i=30.0))

    def singles_per_multi(res):
        multis = sum(v for k, v in res.n_k.items() if k != 1)
        assert multis > 0 or res is hot   # cold must yield real multis
        return res.n_k[1] / max(1, multis)

    assert singles_per_multi(hot) > 2.0 * singles_per_multi(cold)


# ------------------------------------------------------ detection variants
def test_detection_default_bit_for_bit():
    """An explicitly constructed default Detection must reproduce the
    det=None path exactly (guards the Phase 4 parameterization)."""
    from exoverse.kepler_field import DEFAULT_DET, Detection
    targets = sun_like_targets(600)
    base = simulate_universe(targets, seed=5)
    for det in (DEFAULT_DET, Detection(), Detection(snr_threshold=7.1)):
        res = simulate_universe(targets, seed=5, det=det)
        assert res.n_k == base.n_k and res.detected == base.detected
    assert not DEFAULT_DET.needs_rng


def noisy_targets(n: int) -> list:
    """Targets spanning quiet to very noisy, so planet SNRs straddle the
    7.1 threshold (the flat 30 ppm default leaves nothing near it)."""
    cdpps = np.geomspace(20.0, 3000.0, n)
    return [FakeTarget(kepid=i, cdpp_ppm=tuple([float(c)] * 14))
            for i, c in enumerate(cdpps)]


def test_detection_threshold_monotone():
    """Lowering the SNR threshold can only add detections."""
    from exoverse.kepler_field import Detection
    targets = noisy_targets(1500)
    lo = simulate_universe(targets, seed=9, det=Detection(snr_threshold=6.5))
    hi = simulate_universe(targets, seed=9, det=Detection(snr_threshold=8.0))
    assert lo.n_detected_planets > hi.n_detected_planets


def test_detection_mes_ramp_reproducible_and_close_to_step():
    """The logistic MES ramp is seed-reproducible, differs from the hard
    step, and stays within ~25% of its detection count (the ramp only
    reshuffles near-threshold candidates)."""
    from exoverse.kepler_field import Detection
    targets = noisy_targets(1500)
    det = Detection(mes_ramp_width=1.0)
    a = simulate_universe(targets, seed=13, det=det)
    b = simulate_universe(targets, seed=13, det=det)
    step = simulate_universe(targets, seed=13)
    assert a.n_k == b.n_k and a.detected == b.detected
    assert a.detected != step.detected
    assert abs(a.n_detected_planets - step.n_detected_planets) \
        < 0.25 * step.n_detected_planets


def test_detection_binomial_window_reproducible():
    """The binomial window variant draws from its own stream: results are
    seed-reproducible and near the deterministic expected-count variant."""
    from exoverse.kepler_field import Detection
    targets = sun_like_targets(1500)
    det = Detection(window="binomial")
    a = simulate_universe(targets, seed=17, det=det)
    b = simulate_universe(targets, seed=17, det=det)
    step = simulate_universe(targets, seed=17)
    assert a.n_k == b.n_k
    assert abs(a.n_detected_planets - step.n_detected_planets) \
        < 0.25 * step.n_detected_planets
