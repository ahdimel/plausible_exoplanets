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
.venv/bin/exoverse --db worlds.db generate --n 10000 --seed 42
.venv/bin/exoverse --db worlds.db stats | list | inspect <name> | lightcurve <name> <letter>
.venv/bin/exoverse archive --refresh          # snapshot real-planet catalog -> data/
.venv/bin/exoverse --db worlds.db validate    # audit vs real planets -> data/validation_report.json
.venv/bin/exoverse --db worlds.db serve       # web UI on http://127.0.0.1:8321
```

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

- 10,000-system population is the working dataset (seed 42).
- Open scope: N-body verification of packed multis, binary hosts, TTVs,
  radiative-transfer spectra instead of scale-height heuristics, deploying
  the web UI publicly, refreshing Roman/HWO specs as the missions evolve
  (docs/OBSERVATORIES.md has "last researched" dates).
- NASA TAP may recover: `exoverse archive --refresh --prefer nasa` will then
  produce a cleaner audit than exoplanet.eu (true mass provenance flags).
