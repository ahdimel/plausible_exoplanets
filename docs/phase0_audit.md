# Phase 0 audit — Kepler dichotomy decomposition study

*2026-07-03. Baseline audit for the study:*

> How much mutual-inclination dispersion does the Kepler multiplicity
> distribution actually require once intra-system radius correlation is
> modeled jointly — and how degenerate are the two effects?

Companion memo: [prior_art_syssim.md](prior_art_syssim.md) (verdict:
not redundant with SysSim — proceed).

## 1. Pipeline map

| Stage | Module | Notes |
|---|---|---|
| Orchestrator | `generate.py::generate_population` | `SeedSequence(seed)` → per-system child seeds; population is a bit-exact prefix of any larger run. ~2,080 sys/s with full persistence (1M ≈ 8 min, 2.2 GB). |
| System assembly | `system.py::generate_system` | One `default_rng(child)` per system; strict draw order star → planet count (truncated Poisson(2.2), max 7) → per-planet period/radius/mass (independent draws; resample on INVALID) → pair stability → **inclinations** → stellar noise → atmospheres. |
| Transit geometry | `transits.py::compute_geometry` | Winn 2010 incl. eccentricity factors; per-planet `inc_deg`. |
| Detection | `observatories.py::observe` | Kepler row: folded SNR ≥ 7.1, n_transits ≥ 3, baseline 1460 d × 0.92 duty, CDPP fit 30 ppm@Kp12 (√6.5-scaled) + 20 ppm floor + stellar noise, Kp≈V proxy, usable V > 6. |
| Real catalog | `archive.py` | Confirmed-planet snapshot (NASA `ps` TAP primary, exoplanet.eu CSV fallback). Current snapshot: **eu, 2026-07-02, 8,240 planets** (NASA TAP was down). |
| §5.5 experiment | `validate.py::compare_structure` | Synthetic Kepler-detectable systems vs confirmed transit-discovered planets grouped by hostname. |

## 2. Baseline reproduction (pre-fix, worlds.db v0.3, seed 42, N=1M)

Bootstrap: 2,000 resamples over systems/hosts, 68% CI, seed 20260703.

| Quantity | Value |
|---|---|
| Synthetic singles-per-multi | **5.52** (5.42–5.63); N_k = {1: 18019, 2: 2776, 3: 435, 4: 44, 5: 8} |
| Real (confirmed transit-discovered, all surveys) | **4.06** (3.90–4.24); N_k = {1: 3024, 2: 488, 3: 171, 4: 56, 5: 21, 6+: 9} |

Paper §5.5's "5.5 vs 4.1" reproduced exactly.

## 3. Geometry audit → fix (commit `5db5410`)

**Found:** mutual inclinations were applied as `i = i_sys ± mut` (full
Rayleigh tilt, random sign) — no nodal azimuth. A tilt `mut` at uniform
node Ω projects onto the line of sight as ≈ `mut·cos Ω`; Rayleigh(σ)·cos Ω
is exactly Normal(0, σ), while ±Rayleigh(σ) has dispersion σ√2 (and the
wrong shape). The old model therefore inflated the effective inclination
dispersion by √2 — a direct bias on the parameter this study constrains.

**Fix:** `system.py::tilted_inclination_deg` (spherical law of cosines,
node ~ U[0, 2π)). Same RNG draw count, so **stars, planets, noise, and
atmospheres are bit-identical pre/post fix; only inclinations changed**
(planets 2,197,155 unchanged; transiting 32,096 → 32,155). Unit tests:
coplanar limit, node extremes, projected-scatter σ (not σ√2), generator
clustering. 63/63 tests pass.

**Post-fix baseline (new reference for Phase 1 bit-for-bit checks):**

| Quantity | Pre-fix | Post-fix |
|---|---|---|
| Synthetic singles-per-multi | 5.52 (5.42–5.63) | **4.43 (4.36–4.52)** |
| Synthetic N_k | {18019, 2776, 435, 44, 8} | {16787, 3135, 560, 81, 12} |
| Kepler-detectable planets | 25,092 | 25,121 |
| Peas-in-a-pod null median \|ΔlogR\| | 0.175 | 0.179 |
| Real comparator (unchanged) | 4.06 (3.90–4.24) | 4.06 (3.90–4.24) |

**Scientific consequence:** roughly half of the paper's §5.5 excess
(5.52 → 4.06) was the √2 geometry artifact. The residual null-over-real
excess survives (4.43 vs 4.06 heterogeneous census; vs **3.63** on DR25
KOIs, §5) and the size-uniformity entanglement argument stands, but
`docs/paper.html` §5.5 numbers are now stale and must be re-derived after
Phase 2 fixes the comparison catalog (don't churn it twice).

## 4. Detection-model simplifications (audit register)

Ordered by expected impact on the multiplicity vector; Phase 4 perturbs
each unless noted.

1. **No window function**: `n_tr = 1460·0.92/P` is deterministic; real
   per-target `dutycycle`/`dataspan` vary and gate long-period detections
   probabilistically. *DR25 stellar table provides both per target — fix
   available (§5).*
2. **Hard SNR step at 7.1**: real DR25 pipeline efficiency is a gradual
   MES ramp (Christiansen-style gamma CDF) plus vetting efficiency.
3. **Kp ≈ V proxy** and analytic CDPP fit vs real per-target
   `rrmscdpp*` (available per target, 14 durations).
4. **No vetting/reliability layer** (robovetter completeness ~0.9-ish near
   threshold, score-dependent).
5. Depth uses limb-darkened *central* depth rather than duration-averaged
   depth (few-% level; note only).
6. Detection treats each planet independently given geometry (no
   multiplicity-dependent vetting boost as in the real pipeline; note in
   limitations).

## 5. Catalog audit (real side)

**Current pipeline is NOT a Kepler-dichotomy measurement.** It groups all
confirmed transit-discovered planets (Kepler+K2+TESS+ground) by hostname;
TESS's 27-d baseline floods the singles bin. Measured on DR25 KOIs
(CONFIRMED+CANDIDATE, no score/stellar cuts yet):

- N_k = {1: 2423, 2: 440, 3: 153, 4: 52, 5: 20, 6: 2, 7: 1}, 3,091 hosts
- singles-per-multi = **3.63** — vs 4.06 for the heterogeneous census.
  The proper real-side number *widens* the gap to the null.

**TAP access verified 2026-07-03** (exoplanetarchive.ipac.caltech.edu):

- `q1_q17_dr25_koi`: 8,054 rows; dispositions CONFIRMED 2,730 /
  CANDIDATE 1,359 / FALSE POSITIVE 3,965; `koi_score` present on all
  4,089 C+C rows. Query strings to be recorded in `data/PROVENANCE.md`
  at Phase 2 ingestion.
- DR25 stellar: `keplerstellar WHERE st_delivname='q1_q17_dr25_stellar'`
  (200,038 targets; the literal table name `q1_q17_dr25_stellar` is not
  exposed on TAP). Columns include per-target `teff/logg/radius/kepmag`,
  **`dutycycle`, `dataspan`, `rrmscdpp01p5..15p0`, `mesthres*`** — i.e.,
  everything needed to drive detection with real per-target completeness
  inputs.
- **Berger et al. 2020 is NOT on this TAP.** Proposed join: Vizier
  J/AJ/159/280 (Gaia–Kepler stellar properties) ← `kepid` → DR25 stellar.
  Fallback stellar cuts can use DR25 columns alone; Berger radii become a
  Phase 4 robustness variant if the join is deferred.

## 6. Stellar-sample mismatch (decision for Phase 2)

Synthetic stars: volume-limited (d ≤ 300 pc), IMF-weighted (M-dwarf
heavy). Kepler targets: magnitude-limited, FGK-skewed, 0.5–3 kpc. The
multiplicity vector is Teff- and noise-dependent, so a quantitative
(σ_R, σ_i) fit must condition both sides identically. Options:

- **(a) Cut both sides to FGK dwarfs** (Teff 3900–7300 K, logg > 4 or
  Berger evolutionary state) and keep exoverse's analytic noise model.
  Cheap; leaves the Kp-distribution mismatch in place.
- **(b) Condition on the real DR25 target list** (recommended): for each
  synthetic universe, draw host stars from the 200,038-row DR25 stellar
  table (with its `rrmscdpp`, `dutycycle`, `dataspan`), attach exoverse
  planetary systems, and detect with per-target thresholds. Kills the
  target-selection, noise-model, and window-function mismatches in one
  move; simulation observable becomes directly comparable to the KOI
  N_k. Costs: a new `kepler_field` generation mode (Phase 1 flag), and
  planet occurrence priors were calibrated on FGK-ish samples anyway.

## 7. Runtime & grid sizing

Stripped path (generate_system + geometry + Kepler detection, no DB):
**~12,000 sys/s** single-core (20k-system benchmark: 1.68 s; detection
+geometry cost is <5% — generation dominates). Sizing for Phase 3:

- 200k-target universe ≈ 17 s → 20×20 grid × M=10 ≈ 19 h single-core.
  Feasible but tight; planned reductions: (i) reuse physical systems
  along the σ_i axis and redraw only orientations/nodes (common random
  numbers — σ_i enters after all physical draws), (ii) coarse 12×12 grid
  + local refinement, (iii) `multiprocessing` over cells if needed.
  Convergence test (M vs cell-to-cell variation) pre-registered for the
  first coarse grid.

## 8. RNG / reproducibility invariants (verified)

- One `default_rng(child_seed)` per system; (seed, name, dmax_pc)
  regenerates bit-for-bit (`test_generation_deterministic`).
- Draw *count* at the inclination step is unchanged by the geometry fix →
  post-fix worlds share all physical draws with v0.3 (only `inc_deg`
  moved). Post-fix `validate` output is the new frozen baseline.
- New Phase 1 parameters must draw from clearly separated stream
  positions (after existing draws, gated by config) so defaults remain
  bit-for-bit vs the post-fix baseline.

## 9. Actions carried into Phase 1/2

1. σ_R via marginal-preserving construction (Gaussian copula on radius
   quantiles — literal lognormal hierarchy would distort the
   period-dependent bimodal marginal; renormalization note to follow).
2. σ_i as free parameter replacing hardcoded Rayleigh(1.5°); two-component
   mixture flag for the secondary question.
3. Kepler-field conditioning mode (§6b) + full N_k (k = 1..6+) observable.
4. DR25 KOI + stellar ingestion with provenance capture; Berger join
   proposal (§5).
5. Paper §5.5 refresh deferred until the Phase 2 catalog lands.
