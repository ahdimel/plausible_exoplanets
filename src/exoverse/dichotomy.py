"""Observables and distance metric for the Kepler-dichotomy study (Phase 2).

The comparison philosophy (Zhu 2020; docs/phase0_audit.md): compute the SAME
summary statistics on the synthetic detected catalog (kepler_field.
UniverseResult) and on the real detected catalog (DR25 KOIs grouped by host),
so the Kepler detection bias enters both sides identically and cancels by
construction. Both builders therefore funnel through one shared code path
(_observables_from_groups) operating on per-system lists of
(period_days, radius_re).

Statistics
----------
- n_k          : detected-multiplicity histogram {1..5, "6+"}.
- dlogr        : |log10(R_out / R_in)| over period-adjacent detected pairs
                 with both radii < 30 Re (peas-in-a-pod size uniformity;
                 Weiss+ 2018).
- monotonicity : mean Spearman rank correlation of radius against period
                 order within systems of >= 3 detections (intra-system size
                 ordering; Weiss+ 2018, Millholland+ 2017).

Distances
---------
- multiplicity_distance "multinomial" (fiducial): shape-only — the real N_k
  counts scored against the synthetic multiplicity *proportions*, so our
  fixed planets-per-star rate (MEAN_PLANETS) cannot contaminate the
  (sigma_r, sigma_i) fit. Reported as the per-real-system negative log
  multinomial likelihood minus its minimum (attained when the synthetic
  proportions equal the real ones), i.e. the KL divergence
  D(real || synthetic): >= 0, and ~0 when the shapes match. Synthetic
  proportions are Laplace-smoothed (+0.5 per bin) so empty synthetic bins
  stay finite.
- multiplicity_distance "poisson" (variant): negative log Poisson
  likelihood of the real per-bin counts under synthetic counts scaled to
  equal target numbers — sensitive to the absolute detected-system rate,
  not just the shape.
- size_distance "ks" (fiducial): two-sample Kolmogorov-Smirnov D on the
  dlogr samples (implemented by hand as in validate.ks_2samp; no scipy).
- size_distance "ad": two-sample Anderson-Darling, the standardized
  k-sample statistic of Scholz & Stephens (1987) with k=2 — more weight in
  the tails than KS.
- combined_distance: components plus a weighted sum; weights are explicit
  and revisited in Phase 3 (keep modular).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from .kepler_field import N_K_KEYS, UniverseResult

MAX_PAIR_RADIUS_RE = 30.0   # dlogr pairs require both radii below this
MIN_MONOTONICITY_N = 3      # planets per system for the Spearman statistic


# ---------------------------------------------------------------- Spearman
def _average_ranks(xs: Sequence[float]) -> List[float]:
    """Ranks 1..n with ties assigned the average of their positions."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for m in range(i, j + 1):
            ranks[order[m]] = avg
        i = j + 1
    return ranks


def spearman(x: Sequence[float], y: Sequence[float]) -> float:
    """Spearman rank correlation: Pearson correlation of the (tie-averaged)
    ranks. Returns nan for n < 2 or when either variable is constant."""
    n = len(x)
    if n < 2 or n != len(y):
        return float("nan")
    rx, ry = _average_ranks(x), _average_ranks(y)
    mx, my = sum(rx) / n, sum(ry) / n
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0.0 or syy == 0.0:
        return float("nan")
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    return sxy / math.sqrt(sxx * syy)


# ------------------------------------------------------------- Observables
@dataclass
class Observables:
    """Detection-biased summary statistics of one detected catalog (real or
    synthetic — built by the identical code path, so biases cancel in the
    comparison). n_targets is the number of stars searched (available for
    the synthetic side and for the DR25 stellar-cut count; None when
    unknown) — used only by the "poisson" multiplicity mode."""
    n_k: Dict[object, int]
    dlogr: List[float]                  # sorted |log10(R_out/R_in)|
    monotonicity: Optional[float]       # mean Spearman; None if too few
    n_systems: int
    n_planets: int
    n_targets: Optional[int] = field(default=None)


def _observables_from_groups(groups: Sequence[Sequence[Tuple[float, float]]],
                             n_targets: Optional[int] = None) -> Observables:
    """Shared statistics path. groups = per-system lists of
    (period_days, radius_re) detected planets; each group is sorted by
    period here so real and synthetic catalogs are treated identically."""
    n_k: Dict[object, int] = {k: 0 for k in N_K_KEYS}
    dlogr: List[float] = []
    rhos: List[float] = []
    n_systems = 0
    n_planets = 0
    for group in groups:
        planets = sorted(group, key=lambda pr: pr[0])
        if not planets:
            continue
        n_systems += 1
        k = len(planets)
        n_planets += k
        n_k[k if k <= 5 else "6+"] += 1
        for (_, r_in), (_, r_out) in zip(planets, planets[1:]):
            if r_in < MAX_PAIR_RADIUS_RE and r_out < MAX_PAIR_RADIUS_RE:
                dlogr.append(abs(math.log10(r_out / r_in)))
        if k >= MIN_MONOTONICITY_N:
            rho = spearman([p for p, _ in planets], [r for _, r in planets])
            if not math.isnan(rho):
                rhos.append(rho)
    mono = sum(rhos) / len(rhos) if rhos else None
    return Observables(n_k=n_k, dlogr=sorted(dlogr), monotonicity=mono,
                       n_systems=n_systems, n_planets=n_planets,
                       n_targets=n_targets)


def observables_from_universe(result: UniverseResult) -> Observables:
    """Observables of a synthetic kepler_field.UniverseResult."""
    return _observables_from_groups(result.detected,
                                    n_targets=result.n_targets)


def observables_from_koi(koi_systems: Dict[int, Sequence],
                         n_targets: Optional[int] = None) -> Observables:
    """Observables of the real KOI catalog: {kepid: [(period, prad), ...]}.

    Entries may be (period, prad) tuples or KOI-like objects exposing
    .period and .prad (kepler_data.load_koi_systems output). Pass
    n_targets = len(load_stellar_targets(...)) to enable the "poisson"
    multiplicity mode on absolute rates."""
    groups = []
    for kois in koi_systems.values():
        group = []
        for k in kois:
            if hasattr(k, "period"):
                group.append((float(k.period), float(k.prad)))
            else:
                period, prad = k
                group.append((float(period), float(prad)))
        groups.append(group)
    return _observables_from_groups(groups, n_targets=n_targets)


# --------------------------------------------------------------- distances
def multiplicity_distance(syn: Observables, real: Observables,
                          mode: str = "multinomial") -> float:
    """Distance between detected-multiplicity histograms (see module
    docstring for the two modes). Fiducial mode is "multinomial"."""
    if mode == "multinomial":
        n_real = sum(real.n_k.get(k, 0) for k in N_K_KEYS)
        if n_real == 0:
            return float("nan")
        # Laplace-smoothed synthetic proportions (add 0.5 per bin)
        n_syn = sum(syn.n_k.get(k, 0) for k in N_K_KEYS)
        denom = n_syn + 0.5 * len(N_K_KEYS)
        d = 0.0
        for k in N_K_KEYS:
            q = real.n_k.get(k, 0) / n_real
            if q > 0.0:
                p = (syn.n_k.get(k, 0) + 0.5) / denom
                d += q * math.log(q / p)
        return d
    if mode == "poisson":
        if not syn.n_targets or not real.n_targets:
            raise ValueError("poisson mode needs n_targets on both sides")
        n_real = sum(real.n_k.get(k, 0) for k in N_K_KEYS)
        if n_real == 0:
            return float("nan")
        scale = real.n_targets / syn.n_targets
        nll = 0.0
        for k in N_K_KEYS:
            lam = (syn.n_k.get(k, 0) + 0.5) * scale   # smoothed, lam > 0
            n = real.n_k.get(k, 0)
            nll += lam - n * math.log(lam) + math.lgamma(n + 1.0)
        return nll / n_real
    raise ValueError(f"unknown multiplicity mode: {mode!r}")


def _ks_2samp_stat(x1: Sequence[float], x2: Sequence[float]) -> float:
    """Two-sample KS statistic D (merge walk over the sorted samples, same
    construction as validate.ks_2samp; no scipy)."""
    a, b = sorted(x1), sorted(x2)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return float("nan")
    i = j = 0
    d = 0.0
    while i < n1 and j < n2:
        if a[i] <= b[j]:
            i += 1
        else:
            j += 1
        d = max(d, abs(i / n1 - j / n2))
    return d


def _ad_2samp_stat(x1: Sequence[float], x2: Sequence[float]) -> float:
    """Standardized two-sample Anderson-Darling statistic
    T = (A^2_kN - (k-1)) / sigma_N with k=2 (Scholz & Stephens 1987,
    eqs. 3-7; continuous version, ties broken by sort order — dlogr samples
    are continuous so exact ties are measure-zero). T ~ 0 for samples from
    one distribution and grows without bound as they separate."""
    n1, n2 = len(x1), len(x2)
    n = n1 + n2
    if n1 == 0 or n2 == 0 or n < 4:   # sigma_N needs (N-1)(N-2)(N-3) > 0
        return float("nan")
    combined = sorted([(v, 0) for v in x1] + [(v, 1) for v in x2])
    a2 = 0.0
    m1 = 0   # running count of sample-1 members among the first j order stats
    for j in range(1, n):            # j = 1 .. N-1
        if combined[j - 1][1] == 0:
            m1 += 1
        m2 = j - m1
        denom = j * (n - j)
        a2 += (n * m1 - j * n1) ** 2 / denom / n1
        a2 += (n * m2 - j * n2) ** 2 / denom / n2
    a2 /= n
    # Variance of A^2_kN under H0 (Scholz & Stephens 1987, eq. 4), k=2
    k = 2
    h_cap = 1.0 / n1 + 1.0 / n2
    h = sum(1.0 / i for i in range(1, n))
    g = sum((1.0 / (n - i)) * sum(1.0 / j for j in range(i + 1, n))
            for i in range(1, n - 1))
    a = (4.0 * g - 6.0) * (k - 1) + (10.0 - 6.0 * g) * h_cap
    b = ((2.0 * g - 4.0) * k * k + 8.0 * h * k
         + (2.0 * g - 14.0 * h - 4.0) * h_cap - 8.0 * h + 4.0 * g - 6.0)
    c = ((6.0 * h + 2.0 * g - 2.0) * k * k + (4.0 * h - 4.0 * g + 6.0) * k
         + (2.0 * h - 6.0) * h_cap + 4.0 * h)
    d = (2.0 * h + 6.0) * k * k - 4.0 * h * k
    var = ((a * n ** 3 + b * n ** 2 + c * n + d)
           / ((n - 1.0) * (n - 2.0) * (n - 3.0)))
    return (a2 - (k - 1)) / math.sqrt(var)


def size_distance(syn: Observables, real: Observables,
                  mode: str = "ks") -> float:
    """Distance between the |dlogR| (size-uniformity) samples: two-sample
    KS D in [0, 1] (fiducial) or the standardized two-sample
    Anderson-Darling statistic ("ad")."""
    if mode == "ks":
        return _ks_2samp_stat(syn.dlogr, real.dlogr)
    if mode == "ad":
        return _ad_2samp_stat(syn.dlogr, real.dlogr)
    raise ValueError(f"unknown size mode: {mode!r}")


def combined_distance(syn: Observables, real: Observables,
                      w_mult: float = 1.0, w_size: float = 1.0,
                      mult_mode: str = "multinomial",
                      size_mode: str = "ks") -> Dict[str, float]:
    """Component distances and their weighted sum. The components live on
    different scales (KL nats vs KS D), so the weights are explicit knobs —
    fiducial 1.0/1.0, revisited once (sigma_r, sigma_i) grids are run."""
    d_mult = multiplicity_distance(syn, real, mode=mult_mode)
    d_size = size_distance(syn, real, mode=size_mode)
    return {"multiplicity": d_mult, "size": d_size,
            "w_mult": w_mult, "w_size": w_size,
            "total": w_mult * d_mult + w_size * d_size}
