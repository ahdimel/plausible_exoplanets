# Pre-registered robustness plan (Phase 4)

*Written at the end of Phase 2, 2026-07-03, BEFORE any (σ_R, σ_i) grid
results were computed or examined. Each variant re-derives the credible
region on the (σ_R, σ_i) plane; the pass criterion is that the region
moves by less than its own width (else → pre-registered Outcome D
handling: identify the dominating simplification, propose the fix, do not
interpret the science).*

Fiducial configuration: cuts and metric per docs/phase2_design.md
(stellar FGK-dwarf cuts on DR25; KOI CONFIRMED+CANDIDATE with
koi_score > 0.5, P ∈ [0.5, 640] d, R_p ≤ 30 R⊕; multinomial N_k shape
distance + KS |Δlog R| distance, equal weights).

## Variants (minimum set, brief §Phase 4 + Phase 0/2 additions)

1. **Robovetter score cut**: koi_score > 0.0 and > 0.9 (two alternatives
   to the fiducial 0.5).
2. **Stellar sample cuts**: (a) Teff window narrowed to 4700–6500 K
   (GK); (b) logg ≥ 4.2; (c) dataspan cut dropped.
3. **Completeness perturbation**: detection SNR threshold 7.1 → 6.5 and
   8.0 (≈ ±20% efficiency near threshold under the SNR^~steep local
   slope); plus a smooth MES-ramp variant (logistic with width 1.0
   centered at 7.1) vs the hard step.
4. **Window function**: probabilistic transit-count variant — draw
   observed transit epochs as Binomial(round(dataspan/P), dutycycle) and
   require ≥ 3 — vs the deterministic expected-count fiducial.
5. **Radius uncertainty injection**: perturb catalog koi_prad by its
   reported errors (Gaussian, symmetrized) and re-derive the real
   |Δlog R| distribution; repeat 25×; fold the induced spread into the
   size-distance comparison.
6. **Multiplicity vector truncation**: N_k truncated at 4+ (merge
   4, 5, 6+) vs the full k = 1..6+ fiducial.
7. **Uniformity metric**: Gilbert & Fabrycky (2020)-style monotonicity
   as the size observable instead of |Δlog R| KS; also AD instead of KS
   on |Δlog R|.
8. **Occurrence normalization**: Poisson likelihood on absolute N_k
   counts (occurrence-sensitive) vs the fiducial shape-only multinomial —
   quantifies how much the fixed planets-per-star rate (Poisson 2.2)
   leaks into the constraint.
9. **Stellar properties source**: Berger et al. 2020 (Gaia–Kepler,
   Vizier J/AJ/159/280) radii/Teff joined on kepid, if network allows —
   affects both the stellar cuts and koi_prad (∝ R_star). If the join is
   unavailable, document as a known limitation instead.

## Reporting

One table: variant | Δ(best-fit σ_i) | Δ(best-fit σ_R) | credible-region
overlap fraction with fiducial | pass/fail. Plus one sentence per failed
variant. The table is a paper deliverable regardless of outcome.

## Decision rule (pre-registered)

- All pass → interpret against Outcomes A/B/C as defined in the study
  brief.
- Any fail → Outcome D path: report the dominating simplification and a
  proposed fix; no headline claims until resolved.
