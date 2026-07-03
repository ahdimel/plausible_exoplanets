# exoverse — operational skills

Task-shaped recipes for continuing this work. Read CLAUDE.md first for
architecture; this file is "how do I do X" with the traps called out.

## Skill: regenerate the population

```bash
rm -f worlds.db   # schema version mismatch raises; DBs are disposable
.venv/bin/exoverse --db worlds.db generate --n 10000 --seed 42
```
- ~10k systems takes a couple of minutes; run in background for larger N.
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

## Skill: run and extend the web UI

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
