# Checkpoint 4 — results summary and outcome decision (2026-07-04)

*All Phase 3 + Phase 4 computation is complete. Per the study brief, no
paper text is drafted until we agree which pre-registered Outcome this
is. This document is the decision input.*

## The question (verbatim from the brief)

> How much mutual-inclination dispersion does the Kepler multiplicity
> distribution actually require once intra-system radius correlation is
> modeled jointly — and how degenerate are the two effects?

## Headline numbers (fine grid, M=16 in the valley, smoothed likelihood)

| Quantity | Value |
|---|---|
| σ_R (radius-correlation) | 1.55; 68% [1.50, 1.64]; 95% [1.36, 1.86] |
| implied sibling correlation ρ = 1/(1+σ_R²) | ≈ 0.29 |
| σ_i (Rayleigh, cold population) | upper limit: ≤ 0.26° (68%), ≤ 0.44° (95%), σ_R free |
| σ_i with σ_R fixed uncorrelated | ≤ 0.21° (68%), ≤ 0.24° (95%) |
| correlated vs uncorrelated radii | Δln L = 10.3 (uncorrelated firmly disfavored) |
| posterior correlation (ln σ_R, ln σ_i) | **0.07 — essentially orthogonal** |
| single-population multiplicity GoF | **rejected**: G² = 32.6 / 3 dof (p ≈ 4×10⁻⁷) |
| best mixture | f_hot = 0.10, σ_i,hot = 10° (30° similar), cold σ_i = 0.2° |
| mixture vs single | Δln L = +9.2, ΔAIC = −14.4 (seed noise ~4.4) |
| mixture residual GoF | G² = 17.2 (improved from 25.8, still > 7.8 at p=0.05) |
| SPM (singles per multi) | real 3.40; single-pop best 2.90; mixture 3.29 |

Decomposition (the direct answer to "how degenerate"): σ_R is
constrained almost entirely by |ΔlogR| (~1084 ln L of leverage across
the row; the multiplicity term moves ~5, within noise) and σ_i almost
entirely by the multiplicity vector (~461 ln L; the size term is flat).
Under DR25 conditioning the two effects are **not degenerate**: each
statistic pins its own parameter.

## Robustness table (pre-registered, docs/robustness_plan.md)

Operationalized pass rule: |Δln σ_R(best)| < 0.10 (the fiducial 68%
width), |Δln UL95(σ_i)| < 0.31, and 95%-region overlap ≥ 0.5.

| variant | best σ_R (Δln) | σ_i UL95 (Δln) | overlap95 | verdict |
|---|---|---|---|---|
| score > 0.0 | 1.71 (+0.00) | 0.42 (+0.00) | 0.99 | pass |
| score > 0.9 | 1.84 (+0.07) | 0.54 (+0.27) | 0.91 | pass |
| Teff 4700–6500 (GK) | 1.75 (+0.02) | 0.39 (−0.06) | 0.89 | pass |
| logg ≥ 4.2 | 1.62 (−0.06) | 0.43 (+0.05) | 1.00 | pass |
| no dataspan cut | 1.72 (+0.01) | 0.37 (−0.11) | 0.66 | pass |
| SNR 6.5 | 1.71 (+0.00) | 0.42 (+0.02) | 0.93 | pass |
| SNR 8.0 | 1.72 (+0.01) | 0.37 (−0.10) | 0.78 | pass |
| logistic MES ramp | 1.68 (−0.02) | 0.43 (+0.04) | 0.94 | pass |
| binomial window | 1.71 (+0.00) | 0.41 (+0.00) | 1.00 | pass |
| N_k truncated 4+ | 1.72 (+0.01) | 0.44 (+0.07) | 0.98 | pass |
| KDE bandwidth ×0.5 | 1.80 (+0.05) | 0.41 (−0.01) | 0.87 | pass |
| **radius-error injection** | 2.26 (+0.28) | 0.44 (+0.06) | 0.42 | **FAIL** |
| **KDE bandwidth ×2** | 1.42 (−0.18) | 0.45 (+0.07) | 0.57 | **FAIL** |
| **monotonicity as size obs.** | 3.00 (+0.56) | 0.72 (+0.55) | 0.03 | **FAIL** |
| **Poisson occurrence** | 1.58 (−0.08) | 0.22 (−0.63) | 0.08 | **FAIL** |

Reading of the failures (all four localize away from the headline):

1. **radius_err**: the synthetic radii carry no measurement noise while
   the real |ΔlogR| does, so σ_R's *central value* absorbs the
   difference (+0.28 ln when real radii are perturbed by catalog
   errors). σ_i is untouched (+0.06). Proposed fix: convolve synthetic
   detected radii with the DR25 fractional-error distribution
   (scoring-time change; ~1 h of re-simulation of the valley).
2. **bw ×2**: σ_R central value carries ~±0.2 ln metric systematic. The
   ×0.5 direction passes; only the heavy-smoothing direction fails, and
   only on σ_R.
3. **monotonicity**: pre-documented generator limitation (the real
   Spearman monotonicity +0.47 is unreachable anywhere in the model
   family) — the variant chases misspecification noise; unusable as the
   size observable rather than evidence against the result.
4. **poisson**: quantifies exactly what it was designed to quantify —
   the fixed intrinsic occurrence (Poisson 2.2/star) leaks into
   occurrence-sensitive fits (σ_i UL *tightens* to 0.22°). Conclusions
   are explicitly conditional on shape-only comparison; note the failure
   direction makes the headline upper limit conservative.

All nine sample/detection/completeness/window variants pass — the
brief's Outcome-D trigger ("the detection model's geometry or
completeness simplifications dominate the error budget") does **not**
fire. The robustness plan's stricter "any fail" clause does, on the four
metric/real-side variants above.

## Pre-registered outcomes (verbatim from the brief)

- **Outcome A:** σ_i constraint weakens substantially (>50% wider, or
  best-fit σ_i drops materially) when σ_R is free → headline result.
- **Outcome B:** σ_i constraint barely moves → still publishable as a
  null: the dichotomy is robust to the size-correlation degeneracy;
  quantify the ridge anyway.
- **Outcome C:** two-component mixture remains required at all σ_R →
  the dichotomy survives joint modeling; report the hot fraction f_hot
  as a function of σ_R.
- **Outcome D (failure mode):** the detection model's geometry or
  completeness simplifications dominate the error budget → stop, report
  which simplification, propose the fix before drawing conclusions.

## Where the results sit

- **A's trigger fires numerically**: the 95% σ_i upper limit widens from
  0.24° to 0.44° (+83%, > the 50% threshold) when σ_R is freed. But the
  mechanism is *not* a degeneracy (posterior correlation 0.07); the
  correlated model simply has a flatter σ_i profile. And qualitatively
  nothing changes: the cold population is near-coplanar either way.
- **C's trigger fires**: the mixture is still statistically demanded
  (ΔAIC = −14.4; single-population G² rejected at p ≈ 4×10⁻⁷), with a
  small hot fraction f_hot ≈ 0.05–0.2 — though the mixture itself still
  leaves G² = 17.2, so it is an improvement, not a complete resolution.
  (Caveat: f_hot was swept at the best-fit σ_R only, justified by the
  measured orthogonality, not "at all σ_R" as the outcome text asks.)
- **D per the brief's definition does not fire** (all detection /
  completeness / sample variants pass); the four metric-side failures
  bound the σ_R central value to ~±0.2–0.3 ln and flag the shape-only
  conditioning — caveats, or one cheap fix (radius-error convolution),
  not a stop, unless we choose to honor the robustness plan's stricter
  letter.

## Proposed reading (for discussion, not yet adopted)

A joint **A + C** paper: "Once intra-system radius correlation is
modeled jointly, the Kepler multiplicity distribution requires almost no
mutual-inclination dispersion in the dominant population (σ_i ≤ 0.44°
at 95%), and the two effects — often assumed degenerate — decouple
cleanly under DR25 conditioning. A small misaligned fraction
(f_hot ≈ 10%) is still demanded, so the dichotomy survives joint
modeling in attenuated form." With: the σ_R systematic honestly boxed
(radius errors + bandwidth: central value 1.4–2.3), the Poisson
conditionality stated, and the radius-error convolution either run
before drafting (~1 h compute + small code change) or listed as the
named follow-up.

Alternative stricter reading: honor the robustness plan's "any fail"
letter → run the radius-error fix first, re-derive the σ_R region, then
draft with the remaining fails (bandwidth, monotonicity, poisson) as
documented conditionalities.

## Decision needed at this checkpoint

1. Which Outcome do we declare (A+C as proposed / B+C / strict-D pause)?
2. Run the radius-error convolution fix before drafting, or after as
   follow-up?
3. Any changes to the headline figure set before the paper section?
