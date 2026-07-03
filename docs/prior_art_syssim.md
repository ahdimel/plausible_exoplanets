# Prior-art memo: SysSim vs. the planned (σ_R, σ_i) decomposition

*Phase 0 deliverable for the Kepler-dichotomy decomposition study. Written
2026-07-03 from the arXiv abstracts and ar5iv methods sections of the
"Architectures of Exoplanetary Systems" (AE) series; exact fitted values
that were not retrievable online are marked TO-VERIFY (human: PDFs of AE I
Table 2 / AE III results tables would close these).*

## Research question we are positioning

> How much mutual-inclination dispersion does the Kepler multiplicity
> distribution actually require once intra-system radius correlation is
> modeled jointly — and how degenerate are the two effects?

## What the SysSim papers did

**AE I — He, Ford & Ragozzine 2019, MNRAS 490, 4575 (arXiv:1907.07773).**
Clustered Poisson point-process forward model for Kepler's FGK dwarfs,
pushed through a model of the Kepler DR25 detection pipeline. Verified
methods facts:

- **Radii are clustered**: each cluster draws a radius scale `R_p,c` from a
  broken power law; each planet draws `R_p,i ~ Lognormal(R_p,c, σ_R)`.
  `σ_R` is a free parameter (fitted value TO-VERIFY).
- **Mutual inclinations are a two-component Rayleigh mixture**: a fraction
  `f_σi,high` of systems draw from Rayleigh(σ_i,high), the rest from
  Rayleigh(σ_i,low). Three free parameters (fitted values TO-VERIFY).
- **Distance**: weighted sum of KS/AD distances over ~9 observed summary
  statistics (planet rate, transit multiplicity, periods, period ratios,
  depths, depth ratios, durations, period-normalized duration ratios),
  each normalized by its Monte Carlo scatter.
- Headline: clustered periods+sizes fit Kepler far better than independent
  draws; the Kepler dichotomy is read as evidence for a high-inclination
  population; 0.56 (+0.18/−0.15) of FGK stars host a planet >0.5 R⊕ within
  3–300 d.

**AE II — He et al. 2020, AJ (arXiv:2003.04348).** Adds a spectral-type
(Teff) gradient in occurrence. Orthogonal to our question.

**AE III — He, Ford, Ragozzine & Carrera 2020, AJ 161, 16
(arXiv:2007.14473).** Replaces the two-population mixture with mutual
inclinations set by the angular-momentum-deficit (AMD) stability limit:
median mutual inclination scales as `μ̃_i,n ∝ n^(−1.73 ± 0.09)` with
intrinsic multiplicity n. Conclusion: a *single* population with
multiplicity-dependent excitation "can also match the observed Kepler
population" — i.e., the dichotomy does not demand two populations. They
also note size orderings/spacings are "more extreme than what can be
produced by the detection biases of the Kepler mission alone" (their
answer to Zhu 2020), and sizes remain clustered in this model.

## What SysSim did NOT do (the gap we occupy)

1. **No explicit (σ_R, σ_i) degeneracy map.** Both papers fit ~10-plus
   parameters jointly and report marginal posteriors; neither presents the
   2D likelihood surface over size-correlation strength vs. inclination
   dispersion, quantifies the covariance between them, or asks how the
   σ_i constraint degrades when σ_R is freed vs. pinned.
2. **No decomposition of the excess-singles signal.** Neither paper reports
   what fraction of the dichotomy signal (observed singles-per-multi excess
   over a coplanar/uncorrelated null) is absorbed by radius clustering
   *alone* at fixed σ_i. Our exoverse null-model result — an uncorrelated
   universe *overproduces* singles (5.52 ± 0.10 synthetic vs 4.06 ± 0.17
   real, current heterogeneous-census definition) — is the effect they
   never isolate: size correlation changes which systems yield multiple
   detections, so multiplicity and size uniformity are entangled
   observables.
3. **The two-population question is unsettled between their own papers.**
   AE I invokes a dichotomous mixture; AE III says a continuum suffices.
   Our secondary question (is the mixture statistically demanded *at each
   σ_R*, tested along the degeneracy ridge) directly addresses that
   tension with a controlled, two-knob experiment.

## Honest scoping (what they do better)

- SysSim fits intrinsic multiplicity, period structure, eccentricities, and
  (AE III) an AMD-physics prior jointly; we hold the intrinsic multiplicity
  function (truncated Poisson(2.2)) and period occurrence fixed. Our σ_i
  numbers are therefore *conditional* on the exoverse generative baseline,
  and the paper section must say so — our contribution is the decomposition
  and degeneracy geometry, not a superseding measurement of σ_i.
- Their Kepler detection model (window function, per-target MES ramp,
  vetting efficiency) is more faithful than exoverse's current
  folded-SNR threshold; Phase 0 audit lists the gaps and Phase 4 perturbs
  them. If contour motion under those perturbations exceeds the contour
  width, we stop (pre-registered Outcome D) rather than publish.

## Verdict

**Not redundant — proceed.** The planned result (2D degeneracy surface,
marginal σ_i with/without free σ_R, mixture test along the ridge, and the
"fraction of excess singles absorbed by size correlation" number) is
contained in neither AE I nor AE III, and speaks directly to the AE I vs
AE III disagreement. Position the work as a decomposition study built on a
deliberately minimal null-model generator, citing He+ 2019/2020 as the
joint-modeling prior art it complements.

## References

- He, Ford & Ragozzine 2019, MNRAS 490, 4575 — arXiv:1907.07773
- He et al. 2020, AJ 160, 276 — arXiv:2003.04348 (AE II)
- He, Ford, Ragozzine & Carrera 2020, AJ 161, 16 — arXiv:2007.14473 (AE III)
- Lissauer et al. 2011; Johansen et al. 2012 — the dichotomy
- Ballard & Johnson 2016 — dichotomy for M dwarfs, two-population framing
- Weiss et al. 2018; Millholland et al. 2017 — peas in a pod
- Gilbert & Fabrycky 2020 — architecture similarity metrics
- Zhu 2020; Weiss & Petigura 2020 — detection-bias debate
- Zhu et al. 2018 — independent multiplicity-dependent excitation claim
