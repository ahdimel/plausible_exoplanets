# exoverse — operational skills

Task-shaped recipes for continuing this work. Read CLAUDE.md first for
architecture; this file is "how do I do X" with the traps called out.

## Skill: regenerate the population

```bash
rm -f worlds.db   # schema version mismatch raises; DBs are disposable
.venv/bin/exoverse --db worlds.db generate --n 1000000 --seed 42
.venv/bin/exoverse --db neighborhood.db generate --n 7083 --seed 137 --dmax 30
```
- Generation runs ~2000 systems/s (1M ≈ 9 min, ~2.1 GB). SeedSequence
  children are a stream: growing N keeps all earlier (seed, name) worlds
  identical — 1M is a bit-exact prefix of a future full-reality run
  (`--n 5199422` would be 1:1 with the ~5.2M real systems within 300 pc,
  ~45 min / ~11 GB; stats reports the sampled fraction either way).
- The neighborhood N is not arbitrary: 7,083 =
  `round(stars.expected_systems_within(30.0))`, making that population a
  1:1 solar-neighborhood analog whose HWO counts are absolute yield
  estimates. Recompute it if the density model changes.
- `--dmax` (pc) makes a solar-neighborhood population for direct-imaging
  studies. It consumes ONE rng.random() call, so it never shifts the random
  stream; it is stored in DB meta and must be passed to re-generation
  (`db.dmax_pc` — cli lightcurve and web app already do). Note the draw is
  a disk-structure CDF inversion (stars.sample_distance), NOT a linear
  rescale: distances no longer scale exactly by dmax ratios.
- Changing ANY rng draw inside `generate_system` (order, count, distribution)
  changes every world downstream of that draw. That's allowed — worlds are
  versionless — but re-run `generate`, `validate`, and the tests afterwards,
  and expect all system names to map to different worlds.
- After regenerating, refresh the validation report:
  `.venv/bin/exoverse --db worlds.db validate`

## Skill: add a new observatory

1. Add a noise fit `_<name>_sigma_1hr(mag)` in `observatories.py` and a
   `folded(...)` call inside `observe()` (transit) or a new imaging function
   modeled on `hwo_imaging()`.
2. Pick the band factor for stellar noise: `optical_blue` / `optical_red` /
   `nir` (spot contrast chromaticity).
3. Update `N_TRANSIT_OBSERVATORIES` and the `== 7` row-count assertion in
   `tests/test_pipeline.py`.
4. Document assumptions + sources + volatility in `docs/OBSERVATORIES.md`
   (future missions: include a "last researched" date; their specs drift).
5. Add a detection sanity test in `tests/test_phase2.py` (pattern:
   `test_roman_detects_more_than_tess`, `test_hwo_detects_earth_twin_at_10pc`).

## Skill: add a plausibility rule

- Hard physics bound → `Severity.INVALID` in `planets.py::generate_planet`
  (or `system.py::check_pair_stability` for dynamics). The generator
  resamples; INVALID must never appear in the DB.
- Weak-model regime → `Severity.QUESTIONABLE` with a message that cites the
  literature and says WHY it's questionable. These are product output.
- Then extend `validate.py::audit_rules` so the same rule is audited against
  real planets — a rule that rejects many real planets is a wrong rule.
- Rule ids are dotted paths (`density.exceeds_pure_iron`); flags table is
  queryable by rule.

## Skill: refresh / debug the real-planet catalog

```bash
.venv/bin/exoverse archive --refresh                # NASA first, eu fallback
.venv/bin/exoverse archive --refresh --prefer eu    # force fallback
```
Traps encountered in the wild (all handled in `archive.py`, keep the guards):
- NASA TAP returns ORA-xxxx errors and sometimes a full HTML maintenance page
  **with HTTP 200**. Validation requires the body to start with `pl_name`
  and contain >3000 rows.
- exoplanet.eu CSV: `mass` column silently mixes real measurements with
  upper limits; radius/mass are in Jupiter units; `planet_status` must be
  "Confirmed". Mass counts as provenance "measured" only when its error
  columns imply a >2σ detection.
- Snapshot header line records source + fetch date; `validate` reports carry
  it forward. Never mix conclusions across sources without noting it.

## Skill: interpret validation output

- `mass.deuterium_burning` violations ≈ 4%: these are REAL transiting brown
  dwarfs (CoRoT-3b class) that catalogs file with planets. Not a bug.
- `density.exceeds_pure_iron` ≈ 3%: low-SNR catalog masses. Rising well above
  that after a catalog refresh ⇒ suspect the provenance filter broke.
- KS host_teff_K will ALWAYS fail (D~0.5): real surveys target FGK, our
  sample is IMF-weighted. That row exists to prove we forward-model
  selection, not to pass.
- Meaningful regressions to watch: log-period KS p collapsing (period
  occurrence broken) or synthetic radius median drifting off ~0.36 dex
  (radius valley / M-R relations broken).
- The "System architecture" section compares metrics the generator
  deliberately does NOT model (peas-in-a-pod uniformity, resonance pileups,
  Kepler dichotomy, desert depletion). These are EXPECTED to disagree —
  the synthetic side is the no-formation-physics null; the residual is the
  science. A synthetic-vs-real match here would mean the metric has no
  discriminating power, not that the generator is right.

## Skill: run and extend the web UI

> **Do visual verification in a sub-agent.** Any headless-Chrome screenshot,
> image crop, or layout/chart check should be delegated (`Explore` to look and
> report, `general-purpose` to also fix) so the PNGs never enter the main
> context. See "Workflow (token discipline)" in CLAUDE.md.

```bash
.venv/bin/exoverse --db worlds.db serve --port 8321
```
- Templates in `src/exoverse/web/templates/`, charts in `web/charts.py`
  (server-rendered SVG strings; no JS chart libs).
- Chart conventions: colors ONLY via CSS vars (`--series-N`, defined
  light+dark in base.html); composition classes map to fixed slots in
  `charts.CLASS_SLOTS` — never reassign colors by rank; tooltips via
  `data-tip` attribute (base.html has the global handler).
- Verify visually: headless Chrome —
  `"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless
  --disable-gpu --window-size=1200,1400 --screenshot=out.png <url>`
- The app is read-only; anything that writes belongs in the CLI, not the web.

## Skill: update the project paper (docs/paper.html)

> **Run this in a sub-agent.** Hand the whole job to a `general-purpose`
> sub-agent (re-derive numbers, edit prose + figure arrays, verify light AND
> dark in headless Chrome, report only the diff). The screenshots and
> full-file reads stay in its throwaway context — never pull full-page PNGs
> into the main thread. See "Workflow (token discipline)" in CLAUDE.md.

- Fully self-contained: inline CSS/JS, no external requests (fonts are
  system stacks; figures are inline SVG built by the `barChart()` helper at
  the bottom of the file). Keep it that way — it must render offline.
- All figures and prose numbers are HARDCODED snapshots (2026-07-02 v0.3
  populations, seeds 42/137, N=1,000,000 + 7,083). After any regeneration
  that changes physics, re-derive them: `stats` for both DBs, the validate
  report (KS + structure sections feed §3 and §4.5), plus the analysis
  queries (HZ-rocky funnel, rare-regime flag counts, HWO loss modes for HZ
  0.8-1.4 Re planets, EEC hosts by spectral type, the IWA×floor sweep of
  §4.4 over stored contrast/separation_mas, Kepler-detectable radius
  median) — pattern: JOIN planets/observations, `observatory LIKE 'HWO%'`.
- Chart data lives in the arrays `fig1/fig2/fig3a/fig3b/fig4` in the script
  tag; each figure has a matching table fallback (`fillTable`).
- Light AND dark mode exist (CSS vars + prefers-color-scheme). Verify both
  with headless Chrome before committing; to force light, strip the dark
  @media block into a temp copy (Chrome headless follows OS dark mode).

## Skill: query the database directly

```python
import sqlite3; c = sqlite3.connect("worlds.db"); c.row_factory = sqlite3.Row
```
- `observations.mode` = 'transit' | 'imaging'; imaging rows use
  contrast/separation_mas, transit rows use sigma/SNR columns.
- `atm_observations.n_transits_5sigma = -1` means infinite/impractical.
- Flags: `SELECT rule, COUNT(*) FROM flags GROUP BY rule ORDER BY 2 DESC`.
- Joins: planets.system_id → systems.id; observations/atmospheres/
  atm_observations.planet_id → planets.id.

## Skill: performance & testing discipline

- Tests must stay fast (<1 min): population tests use n=25..150 systems.
- Physics tests assert against KNOWN real values (Sun BC_V ≈ −0.07, Earth
  depth 84 ppm / 13 h, Jupiter ~1%, Venus–Earth Hill Δ>25, GJ 1214b TSM
  ~300, Earth-twin HWO contrast ~1.7e-10 @ 100 mas). When adding physics,
  find the real-world anchor first, then write the test.
- `pytest` runs from repo root; test imports use helpers in
  `tests/test_physics.py` (`make_sun`, `make_planet`).
