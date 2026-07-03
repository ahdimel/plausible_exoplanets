# exoverse — project context for Claude

Procedurally generated, physics-validated exoplanetary systems with transit
and direct-imaging observability modeling, validated against real exoplanet
catalogs, browsable via CLI and a Flask web UI.

## Purpose (the research question)

Transit surveys only find the planets their geometry, noise, and baselines
allow. This project inverts that: generate the full spectrum of *plausible*
stellar systems bottom-up from physics + occurrence statistics, observe them
with modeled instruments (real and future), and study the selection effects.
The database of generated worlds is the product; each world is individually
inspectable with full metadata and plausibility flags.

## Environment

- Python 3.13 venv at `.venv/` (NOT the system python; NOT another project's
  venv). All commands: `.venv/bin/...`
- Deps: numpy, flask (runtime); pytest, matplotlib (dev). Editable install:
  `.venv/bin/pip install -e ".[dev]"`
- Tests: `.venv/bin/python -m pytest` — keep them green; add tests alongside
  any new feature (user's explicit standing preference).
- The user's other workspace folders (visual_dna, worlds, REI...) are
  unrelated. Everything for this project lives in this repo.

## Commands

```bash
.venv/bin/exoverse --db worlds.db generate --n 100000 --seed 42
.venv/bin/exoverse --db neighborhood.db generate --n 20000 --seed 137 --dmax 30   # solar neighborhood (HWO)
.venv/bin/exoverse --db worlds.db stats | list | inspect <name> | lightcurve <name> <letter>
.venv/bin/exoverse archive --refresh          # snapshot real-planet catalog -> data/
.venv/bin/exoverse --db worlds.db validate    # audit vs real planets -> data/validation_report.json
.venv/bin/exoverse --db worlds.db serve       # web UI on http://127.0.0.1:8321
```

Two working datasets: `worlds.db` (100k systems, d<=300 pc — transit-survey
statistics) and `neighborhood.db` (20k systems, d<=30 pc — direct-imaging /
HWO statistics; a 300 pc sample starves HWO of targets, ~0.001% detectable).
`--dmax` only rescales the distance draw (same rng stream), and is stored in
DB meta so (seed, name) re-generation stays exact — pass `db.dmax_pc` to any
new re-generation call site.

## Architecture (src/exoverse/)

Generation flows: `generate.py` (orchestrator) → `system.py::generate_system`
(one seed = one system, **bit-for-bit deterministic**: star → planets →
stability → inclinations → stellar noise → atmospheres, all from one
`np.random.default_rng(seed)` — never reorder these draws or add draws
mid-sequence, it breaks reproducibility of every stored world) →
`transits.py` geometry → `observatories.py` detectability →
`database.py` persistence.

| Module | Role | Key decisions |
|---|---|---|
| `constants.py` | CODATA/IAU constants, SI | |
| `flags.py` | Flag(severity, rule, message); severities INFO / QUESTIONABLE / INVALID | INVALID ⇒ resample; never persisted. QUESTIONABLE ⇒ kept + recorded. This is the core "plausibility" mechanism. |
| `stars.py` | Kroupa IMF (0.08–2.2 M☉), MS scaling relations, Teff via Stefan-Boltzmann, Flower/Torres BC_V, Pecaut-Mamajek colors, Claret-like limb darkening | Single MS stars only (no binaries/evolution — flagged). Distances ≤300 pc, p(d)∝d². |
| `stellar_noise.py` | granulation (ν_max scaling), p-modes, spot variability (gyrochronology Prot~√age), M-dwarf flares | Chromatic: BAND_FACTORS blue 1.0 / TESS 0.75 / NIR 0.40 (Rackham+18). Coefficients are population-level Kepler calibrations. |
| `planets.py` | occurrence-shaped periods (Petigura+13) & radii (Fulton+17 valley), metallicity-scaled giants, probabilistic M-R, hard validity rules | INVALID rules: >13 Mjup, denser than pure iron (Zeng+16 R=0.78·M^0.27), ρ<0.03 g/cc, inside Roche/1.5 R*. |
| `atmospheres.py` | class (h_he / h_he_rich / steam / secondary / airless), scale height, TSM/ESM (Kempton+18), JWST spectroscopy scoring | Cosmic shoreline is **Mars-normalized**: I_crit = 10.5·(v_esc/v_earth)⁴ — normalizing through Earth puts Earth on the knife's edge (bug we fixed). SPECTRAL_BIN_PENALTY=√5 calibrated to real K2-18b-class results. |
| `system.py` | assembly; pair stability: orbit crossing ⇒ INVALID, Δ<2√3 mutual Hill radii ⇒ INVALID (Gladman), Δ<9 ⇒ QUESTIONABLE (Pu & Wu) | Inclinations may exceed 90° (mirror-equivalent); clipping at 90 creates a fake b=0 pileup (bug we fixed). |
| `transits.py` | Winn 2010 geometry incl. eccentricity factors; limb-darkened light curve by per-annulus integration over [z−k, z+k] only | Grid scaled to planet size ⇒ accuracy independent of k (~1e-13 vs analytic). A fixed whole-disk grid gave 20% errors for Earth-size planets (bug we fixed). |
| `observatories.py` | TESS, Kepler(archival), JWST NIRISS/NIRSpec, ground 1-m, **Roman GBTDS**, **HWO imaging**. σ_total² = σ_instr² + σ_stellar² | Roman/HWO specs are requirement-era, WILL change → docs/OBSERVATORIES.md is the source of truth w/ dates+links. Roman rows mean "if this system were in a Roman bulge field". |
| `archive.py` | real-catalog client: NASA Exoplanet Archive TAP primary, exoplanet.eu CSV fallback; canonical snapshot in data/ | NASA TAP served ORA errors AND HTML-with-HTTP-200 maintenance pages (2026-07-02) → response must start with `pl_name` and have >3000 rows. exoplanet.eu `mass` silently contains upper limits → provenance "measured" only if errors imply >2σ. |
| `validate.py` | (1) run real planets through our INVALID rules; (2) selection-matched KS: synthetic Kepler-DETECTABLE vs real transit-discovered | Never compare raw synthetic vs raw archive — forward-model the selection. Expected residuals: ~4% deuterium violations = real transiting brown dwarfs; ~3% iron = low-SNR masses; host-Teff mismatch is BY DESIGN (IMF vs FGK targeting). |
| `database.py` | SQLite schema v2: systems/planets/observations/atmospheres/atm_observations/flags + meta(schema) | Schema version checked on open; on change, bump SCHEMA_VERSION and regenerate (no migrations — DBs are cheap to rebuild). n_transits_5sigma = -1 encodes infinity. |
| `web/` | Flask app (`create_app(db_path, data_dir)`), Jinja templates, server-rendered SVG charts in `charts.py` | Read-only over the DB (deployable later). Light curves re-simulated from (seed, name) on request. Charts follow the dataviz method: CSS-variable palette, light+dark, fixed class→color slots in `CLASS_SLOTS`, hover via data-tip. |

## Non-obvious invariants

1. **Determinism**: a stored system re-generates exactly from (seed, name);
   `cli.py lightcurve` and `web/app.py` rely on this instead of storing arrays.
2. **INVALID never reaches the DB** — `test_pipeline` asserts it.
3. Transiting planets get 6 transit observation rows + 1 HWO imaging row = 7;
   non-transiting get just the HWO row.
4. Every approximation gets a flag or a docstring citation; "plausible but
   questionable" is a first-class output, not an error.
5. `data/` (catalog snapshots, validation reports) and `*.db` are gitignored —
   regenerable artifacts stay out of the repo.

## Where the last session left off / natural next steps

- Working datasets: worlds.db (100k, seed 42) + neighborhood.db (20k,
  seed 137, dmax 30 pc). HWO model upgraded: per-atmosphere-class geometric
  albedo, photon-limited contrast floor 3e-11*10^(0.2*(V-7)), V<11
  evaluation envelope; stats now report hwo_exo_earth_candidates (HZ,
  0.8-1.4 Re). Neighborhood run: 522 HWO detections, 44 EECs — scaled to
  the real ~10k stars within 30 pc that's ~22 exo-Earths, matching the
  Astro2020 ">=25" goal. Dominant EEC loss: M-dwarf hosts (79% too faint,
  V>=11) then IWA (20%); EEC hosts come out K>G>F, no M.
- Validation at 100k: audit residuals unchanged (3.8% deuterium = brown
  dwarfs, 2.7% iron); Kepler-detectable radius median refined to 0.316 dex
  (0.338 on the 10k subset — sample size, not a regression).
- Project paper: `docs/paper.html` — self-contained HTML (motivation, goals,
  validation, findings, researched future directions). Its figures/numbers
  are HARDCODED snapshots of the 2026-07-02 populations; if you regenerate
  either DB with different physics, the paper's numbers must be re-derived
  and updated by hand (see SKILLS.md).
- Known model deficiency (found at 10k scale, still present): the INTRINSIC
  radius distribution has no dip at 1.7-2.0 Re — the Fulton valley is
  smeared flat by the two wide lognormal peaks in
  `planets.py::sample_small_radius`. Top candidate physics fix. (The
  DETECTED population still matches real surveys well enough for KS.)
- Candidate research directions (paper section 5, with literature refs):
  (1) survey trade studies / yield forecasting (PLATO arXiv:2407.15917,
  Roman arXiv:2305.16204); (2) HWO IWA-vs-contrast-floor design trades +
  exozodi + multi-visit scheduling (TSS25 target list arXiv:2509.20544);
  (3) closed-loop occurrence-rate inference validation (inject synthetic
  universe -> modeled survey -> eta-Earth pipeline -> check recovery);
  (4) exporting labeled light-curve corpora for ML vetting
  (arXiv:2507.19520); (5) public web deploy of the read-only UI.
- Open scope: N-body verification of packed multis, binary hosts, TTVs,
  radiative-transfer spectra instead of scale-height heuristics, refreshing
  Roman/HWO specs as the missions evolve (docs/OBSERVATORIES.md has "last
  researched" dates).
- NASA TAP may recover: `exoverse archive --refresh --prefer nasa` will then
  produce a cleaner audit than exoplanet.eu (true mass provenance flags).
