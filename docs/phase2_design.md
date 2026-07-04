# Phase 2 design: DR25 conditioning, observables, distance metric

*Authoritative spec for the dichotomy study's Phase 2 (see
docs/phase0_audit.md for the audit, docs/sigma_r_note.md for sigma_r).
Sub-agents implement exactly these APIs.*

## Fiducial cuts (Phase 4 varies them; keep as named constants)

Stellar (DR25 stellar delivery, `keplerstellar` TAP table filtered to
`st_delivname='q1_q17_dr25_stellar'`):
- 3900 <= teff <= 7300 K, logg >= 4.0 (FGK dwarfs)
- non-null: radius, mass, dataspan, dutycycle, all 14 rrmscdpp columns
- dataspan > 365 d (excludes barely-observed targets)

KOI (`q1_q17_dr25_koi`): koi_disposition IN (CONFIRMED, CANDIDATE),
koi_score > 0.5, 0.5 <= koi_period <= 640 d (generator support),
koi_prad <= 30 Re, host kepid passes the stellar cuts.

## Module 1: `src/exoverse/kepler_data.py` (ingestion)

- `fetch_dr25(data_dir="data")` — TAP queries (exact strings recorded in
  `data/PROVENANCE.md` with dates + row counts; PROVENANCE.md is
  force-added to git, snapshots stay gitignored). Keep only needed
  columns; snapshots `data/dr25_stellar.csv`, `data/dr25_koi.csv`.
- `KeplerTarget` dataclass: kepid, teff, logg, feh, radius, mass, kepmag,
  dutycycle, dataspan, cdpp_ppm (tuple of 14 floats),
  CDPP_DURATIONS_HR = (1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 6.0, 7.5, 9.0,
  10.5, 12.0, 12.5, 15.0).
- `load_stellar_targets(data_dir, cuts=FIDUCIAL) -> list[KeplerTarget]`
- `load_koi_systems(data_dir, cuts=FIDUCIAL) -> dict[int, list[KOI]]`
  (KOI: kepoi_name, period, prad, score, disposition), hosts filtered by
  the stellar cuts.
- feh may be null in DR25 → default 0.0 (flag in docstring).

## Module 2: `src/exoverse/kepler_field.py` (conditioned simulation)

One universe = one pass over the (cut) target list: each target hosts one
generator draw. Empirical per-target noise replaces the analytic model —
CDPP already contains stellar+instrument noise, so the exoverse
stellar-noise model is NOT added on top.

- Star construction from a target: mass/radius/teff/feh direct;
  luminosity = radius^2 (teff/5772)^4; limb darkening via the existing
  helper in stars.py; distance/mags irrelevant to detection (fill
  something sane).
- Planet+inclination draws: `system.py` is refactored to expose
  `_draw_planets(rng, star, arch) -> (planets, sys_inc_deg)` — the
  existing planet loop + inclination block, byte-identical draw order —
  and `generate_system` calls it; kepler_field calls it directly and
  SKIPS noise/atmospheres (not needed, big speedup).
- Detection per planet: geometric transit via compute_geometry; then
  n_tr = dataspan * dutycycle / period, require n_tr >= 3;
  sigma = interp(t14 in CDPP_DURATIONS_HR -> cdpp_ppm) (clamp ends);
  snr_total = depth_ppm / sigma * sqrt(n_tr); detect if >= 7.1.
- `simulate_universe(targets, seed, arch) -> UniverseResult`:
  n_k dict {1..5, "6+"}, detected (period, radius) lists grouped per
  system (for pair statistics), n_targets, n_detected_planets. One
  rng per system from SeedSequence(seed) children, deterministic.

## Module 3: `src/exoverse/dichotomy.py` (observables + distance)

Computed IDENTICALLY on real (KOI dict) and synthetic (UniverseResult)
detected catalogs — detection bias cancels by construction (the Zhu 2020
answer).

- `Observables`: n_k (dict), dlogr (sorted list of |log10(R_out/R_in)|
  over adjacent detected pairs, both R < 30), monotonicity (mean Spearman
  rank corr of radius vs period order over systems with >= 3 detected;
  None if too few), n_systems, n_planets.
- `multiplicity_distance(syn, real, mode)`:
  - "multinomial" (fiducial): -log multinomial likelihood of real N_k
    under synthetic N_k proportions (shape-only — our fixed
    planets-per-star rate must not contaminate the fit), normalized per
    real system. Laplace-smooth synthetic proportions (add 0.5).
  - "poisson" (variant): Poisson LL on absolute counts scaled to equal
    target numbers.
- `size_distance(syn, real, mode="ks"|"ad")`: two-sample KS D (fiducial)
  or Anderson-Darling on the |dlogR| samples.
- `combined_distance(syn, real, w_mult=1.0, w_size=1.0) -> dict` with
  components AND the weighted sum (keep modular; weights revisited).

## Testing rules

- No network in tests. kepler_data tests run against a tiny checked-in
  fixture CSV (tests/fixtures/, ~20 rows, same columns). kepler_field +
  dichotomy tests build fake targets in code.
- Existing 75 tests must stay green; default-arch bit-for-bit test guards
  the system.py refactor.
- Keep the suite < 1 min.

## Out of scope for Phase 2 (pre-registered for Phase 4)

Berger 2020 stellar join (Vizier); probabilistic window function;
MES-ramp detection efficiency; per-duration mesthres; score-cut variants
(0.0, 0.9); completeness +-20% near threshold; radius-error injection;
N_k truncation at 4+; AD vs KS. Listed in docs/robustness_plan.md before
any (sigma_r, sigma_i) results are examined.
