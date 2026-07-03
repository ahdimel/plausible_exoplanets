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
.venv/bin/exoverse --db worlds.db generate --n 1000000 --seed 42
.venv/bin/exoverse --db neighborhood.db generate --n 7083 --seed 137 --dmax 30   # 1:1 solar neighborhood (HWO)
.venv/bin/exoverse --db worlds.db stats | list | inspect <name> | lightcurve <name> <letter>
.venv/bin/exoverse archive --refresh          # snapshot real-planet catalog -> data/
.venv/bin/exoverse --db worlds.db validate    # audit vs real planets -> data/validation_report.json
.venv/bin/exoverse --db worlds.db serve       # web UI on http://127.0.0.1:8321
```

Two working datasets: `worlds.db` (1M systems, d<=300 pc — a 19% random
sample of the ~5.2M real systems in that volume; transit-survey statistics)
and `neighborhood.db` (7,083 systems, d<=30 pc — a 1:1 analog of the real
solar neighborhood: 7083 = round(stars.expected_systems_within(30)), so HWO
counts are absolute yield estimates). Distances follow a disk-structure
density (exponential z falloff, h=300 pc, n0=0.065 MS systems/pc^3 from
RECONS/GCNS). `--dmax` only reshapes the one distance draw (same rng
stream), is stored in DB meta so (seed, name) re-generation stays exact —
pass `db.dmax_pc` to any new re-generation call site — but is NOT a linear
rescale (CDF inversion in stars.sample_distance).

## Workflow (token discipline)

Dispatch a **sub-agent** for bulky, self-contained work so its intermediate
context (screenshots, full-file reads, Chrome runs) never accumulates in the
main thread — the main agent only sees the sub-agent's final report.

- **Updating the findings paper (`docs/paper.html`)**: hand the whole job to a
  `general-purpose` sub-agent. Instruct it to re-derive the numbers (see the
  "update the project paper" skill in SKILLS.md), edit the prose + figure
  arrays, verify light AND dark mode in headless Chrome, and report only what
  it changed. Do NOT read full-page PNGs into the main context.
- **Any graphics rendering or visual checking** (headless-Chrome screenshots
  of the web UI or paper, image cropping, layout/figure verification):
  delegate to a sub-agent (`Explore` if read-only "look and report",
  `general-purpose` if it should also fix what it finds). It returns a
  paragraph; the images stay in its throwaway context.
- Sub-agents start COLD (they re-read this file + SKILLS.md), so use them for
  genuinely separable chunks, not two-line tweaks. Keep CLAUDE.md/SKILLS.md
  current so a cold sub-agent can orient from them alone.

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
| `stars.py` | Kroupa IMF (0.08–2.2 M☉), MS scaling relations, Teff via Stefan-Boltzmann, Flower/Torres BC_V, Pecaut-Mamajek colors, Claret-like limb darkening, disk-structure distances + absolute density normalization | Single MS stars only (no binaries/evolution — flagged). Distance: exponential z falloff (h=300 pc), one rng draw inverted by bisection (sample_distance). expected_systems_within() anchors reality (0.065 sys/pc³): ~5.2M within 300 pc, 7,083 within 30 pc. |
| `stellar_noise.py` | granulation (ν_max scaling), p-modes, spot variability (gyrochronology Prot~√age), M-dwarf flares | Chromatic: BAND_FACTORS blue 1.0 / TESS 0.75 / NIR 0.40 (Rackham+18). Coefficients are population-level Kepler calibrations. |
| `planets.py` | occurrence-shaped periods (Petigura+13) & radii (Fulton+17 valley), metallicity-scaled giants, probabilistic M-R, hard validity rules | INVALID rules: >13 Mjup, denser than pure iron (Zeng+16 R=0.78·M^0.27), ρ<0.03 g/cc, inside Roche/1.5 R*. |
| `atmospheres.py` | class (h_he / h_he_rich / steam / secondary / airless), scale height, TSM/ESM (Kempton+18), JWST spectroscopy scoring, per-planet V-band geometric albedo | Cosmic shoreline is **Mars-normalized**: I_crit = 10.5·(v_esc/v_earth)⁴ — normalizing through Earth puts Earth on the knife's edge (bug we fixed). SPECTRAL_BIN_PENALTY=√5 calibrated to real K2-18b-class results. Albedo draws (sample_geometric_albedo) are Teq-dependent w/ measured anchors (hot Jupiters DARK 0.03–0.11; 22% Venus branch on secondary; icy branch on cold airless) and consume exactly 2 rng draws for every class — keep that invariant. |
| `system.py` | assembly; pair stability: orbit crossing ⇒ INVALID, Δ<2√3 mutual Hill radii ⇒ INVALID (Gladman), Δ<9 ⇒ QUESTIONABLE (Pu & Wu) | Inclinations may exceed 90° (mirror-equivalent); clipping at 90 creates a fake b=0 pileup (bug we fixed). |
| `transits.py` | Winn 2010 geometry incl. eccentricity factors; limb-darkened light curve by per-annulus integration over [z−k, z+k] only | Grid scaled to planet size ⇒ accuracy independent of k (~1e-13 vs analytic). A fixed whole-disk grid gave 20% errors for Earth-size planets (bug we fixed). |
| `observatories.py` | TESS, Kepler(archival), JWST NIRISS/NIRSpec, ground 1-m, **Roman GBTDS**, **HWO imaging**. σ_total² = σ_instr² + σ_stellar² | Roman/HWO specs are requirement-era, WILL change → docs/OBSERVATORIES.md is the source of truth w/ dates+links. Roman rows mean "if this system were in a Roman bulge field". HWO uses the planet's stored geometric_albedo (class-mean fallback); its detection rule has NO integration-time budget, so yields are geometry/floor-limited — albedo moves contrast margins, not counts. |
| `archive.py` | real-catalog client: NASA Exoplanet Archive TAP primary, exoplanet.eu CSV fallback; canonical snapshot in data/ | NASA TAP served ORA errors AND HTML-with-HTTP-200 maintenance pages (2026-07-02) → response must start with `pl_name` and have >3000 rows. exoplanet.eu `mass` silently contains upper limits → provenance "measured" only if errors imply >2σ. |
| `validate.py` | (1) run real planets through our INVALID rules; (2) selection-matched KS: synthetic Kepler-DETECTABLE vs real transit-discovered; (3) architecture null-model comparisons (compare_structure) | Never compare raw synthetic vs raw archive — forward-model the selection. Expected residuals: ~4% deuterium violations = real transiting brown dwarfs; ~3% iron = low-SNR masses; host-Teff mismatch is BY DESIGN (IMF vs FGK targeting). Structure metrics (peas-in-a-pod, period-ratio resonances, Kepler dichotomy, hot-Neptune desert) are SUPPOSED to disagree: synthetic = no-formation-physics null, the residual is the science. |
| `database.py` | SQLite schema v3: systems/planets/observations/atmospheres(+geometric_albedo)/atm_observations/flags + meta(schema) | Schema version checked on open; on change, bump SCHEMA_VERSION and regenerate (no migrations — DBs are cheap to rebuild). n_transits_5sigma = -1 encodes infinity. stats() reports expected_real_systems_within_dmax + fraction_of_reality_sampled. |
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

- v0.3 (2026-07-02): three upgrades landed together, ALL worlds regenerated
  (schema v3). (1) Per-planet V-band geometric albedo
  (atmospheres.sample_geometric_albedo): Teq-dependent Sudarsky-like curve
  for H/He (hot Jupiters dark 0.03-0.11, cold giants ~0.5), 22% Venus
  branch for secondary, icy branch for cold airless; stored in DB, used by
  hwo_imaging. (2) Disk-structure distances + absolute normalization:
  n0=0.065 MS systems/pc^3, h=300 pc → 5.2M real systems within 300 pc,
  7,083 within 30 pc (stars.expected_systems_within). (3) Architecture
  null-model comparisons in validate (compare_structure).
- Working datasets: worlds.db (1M, seed 42; 19% sample of reality, ~2.1 GB)
  + neighborhood.db (7,083, seed 137, dmax 30 — exact 1:1 solar
  neighborhood). Headline results: HWO images 209 planets, 19 EECs
  (absolute yield, vs Astro2020 >=25); EEC loss budget 79% host-too-faint /
  19% inside IWA; yield is IWA-DOMINATED: 60→45 mas triples to 48, 30 mas
  gives 99, while 3x deeper floor or V<12.5 envelope adds ZERO (same stars).
  Flat A_g=0.2 vs drawn albedos changes counts by ~1% (no exposure-time
  budget in the detection rule — that's WHY; adding one is the top HWO-model
  next step, then exozodi + multi-visit).
- Structure findings at 1M (validate → §4.5 of paper): real peas-in-a-pod
  median |dlogR| 0.121 vs null 0.175; ~2x real pileup wide of 3:2/2:1;
  singles-per-multi 4.1 real vs 5.5 null (dichotomy metric entangled with
  size uniformity — null overproduces singles); desert occupancy 2.1% real
  vs 3.7% null. KS p-values collapse at n~25k: watch D, not p.
- Validation at 1M: audit residuals unchanged (3.8% deuterium, 2.7% iron);
  Kepler-detectable radius median 0.313 dex (n=25,092); log-period D=0.074.
- Project paper `docs/paper.html` updated to v0.3 (IWA trade + fig4,
  architecture table); numbers are HARDCODED snapshots — see SKILLS.md for
  the re-derivation recipe. Verified light+dark via headless Chrome.
- Paper restructured for scientific review (2026-07-02, later session):
  new §3 "The generative model" (full methods: IMF/scalings, distance CDF
  equation, occurrence/M-R/eccentricity parameters, INVALID-rule table,
  noise model, shoreline equation, albedo table, transit accuracy,
  instrument-noise-model table, determinism/tests), validation now §4
  (adds selection-matched KS table incl. the honest radius-D=0.21 row),
  findings §5, NEW §6 Limitations, future §7, refs §8 = ~40-entry
  alphabetical bibliography with author-year inline citations (no more
  numbered superscripts). Section cross-refs all renumbered — grep for
  "§5.4"-style refs before renumbering again.
- Known model deficiency (still present): the INTRINSIC radius distribution
  has no dip at 1.7-2.0 Re — Fulton valley smeared flat in
  `planets.py::sample_small_radius`. Top candidate physics fix. Related new
  target: give multis correlated radii (peas-in-a-pod knob) so the §4.5
  null can be dialed between "no memory" and "full memory".
- Density-model caveats (documented in stars.py): n0 uncertain ±20%;
  single-h thin disk (no thick disk, Sun's z-offset ignored). A future
  full-reality 300 pc run is `generate --n 5199422` (~45 min, ~11 GB) and
  keeps the current 1M worlds as an exact prefix.
- NASA TAP may recover: `exoverse archive --refresh --prefer nasa` will then
  produce a cleaner audit than exoplanet.eu (true mass provenance flags).
