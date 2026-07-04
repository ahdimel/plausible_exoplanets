# plausible_exoplanets (`exoverse`)

Procedurally generated, physics-validated exoplanetary systems — a bottom-up
research sandbox for transit detectability.

Transit surveys only find the planets their geometry, noise floors, and
baselines allow. This project inverts the problem: it generates the **full
spectrum of plausible stellar systems** consistent with the laws of physics
and current occurrence statistics, inserts validated planets, simulates the
resulting light signatures, and asks *which observatories could actually see
each world*. The output is a queryable SQLite database of worlds, each with
rich validated metadata and explicit plausibility flags.

**Project paper**: [docs/paper.html](docs/paper.html) — motivation, methods,
validation, findings from the 1,007,083-system populations (a 19% sample of
the ~5.2M real systems within 300 pc plus a 1:1 solar-neighborhood analog),
and future directions (fully self-contained; open it in any browser).

## Quick start

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"

.venv/bin/exoverse --db worlds.db generate --n 500 --seed 42
.venv/bin/exoverse --db worlds.db stats
.venv/bin/exoverse --db worlds.db list --transiting
.venv/bin/exoverse --db worlds.db inspect PXS-42-00481      # full metadata + flags
.venv/bin/exoverse --db worlds.db lightcurve PXS-42-00481 d # PNG light curve

.venv/bin/exoverse archive --refresh        # fetch real-planet catalog snapshot
.venv/bin/exoverse --db worlds.db validate  # audit our physics vs real planets
.venv/bin/exoverse --db worlds.db serve     # browser UI at localhost:8321

.venv/bin/python -m pytest      # physics sanity tests vs known real values
```

Generation is fully deterministic per seed: any stored system can be
re-generated bit-for-bit from its `(seed, name)`, which is how `lightcurve`
re-simulates without storing arrays.

## Architecture knobs (Kepler-dichotomy study)

Two optional generator parameters decompose the Kepler dichotomy into
size-uniformity and inclination-dispersion effects (study plan and audit:
`docs/phase0_audit.md`; math: `docs/sigma_r_note.md`):

```bash
.venv/bin/exoverse --db hot.db generate --n 10000 --sigma-r 0.3 --sigma-i 5
```

- `--sigma-r` — intra-system radius correlation via a marginal-preserving
  Gaussian copula: small values give peas-in-a-pod uniformity (Weiss+
  2018), unset keeps the historical independent draws. The population's
  one-point radius distribution is identical for every value.
- `--sigma-i` — Rayleigh scale (deg) of mutual inclinations about the
  system plane (default 1.5). `--f-hot` / `--sigma-i-hot` add a second,
  dynamically hot component (the classic dichotomy model); `--isotropic`
  removes the shared plane entirely.

Defaults are bit-for-bit neutral (enforced by `tests/test_architecture.py`),
and the knobs are stored in DB meta so `(seed, name)` re-generation stays
exact.

## What gets generated

**Stars** — Kroupa (2001) IMF over 0.08–2.2 M☉; main-sequence L(M), R(M)
scaling relations with scatter; T_eff from Stefan–Boltzmann; thin-disk
metallicity distribution; distances from a disk-structure density model
(uniform in the plane, exponential vertical falloff with 300 pc scale
height) with an absolute normalization of 0.065 main-sequence systems/pc³
(RECONS / Gaia GCNS) — so a population is an explicit random sample of the
~5.2 million real systems within 300 pc (~7,100 within 30 pc), and `stats`
reports the sampled fraction of reality;
V/TESS/J magnitudes via Flower (1996)/Torres (2010) bolometric corrections and
Pecaut & Mamajek (2013) colors; quadratic limb-darkening coefficients vs
T_eff (Claret & Bloemen 2011-like).

**Planets** — multiplicity ~Poisson(2.2); period distribution rising to ~10 d
then flat in log P (Petigura+ 2013); bimodal radius distribution with the
period-dependent photoevaporation valley at ~1.9 R⊕ (Fulton+ 2017); giant
occurrence scaling with host metallicity (Fischer & Valenti 2005) and stellar
mass; probabilistic mass–radius relations (Chen & Kipping 2017 / Wolfgang+
2016 style); Kipping (2013) eccentricities for singles, dynamically cold
Rayleigh eccentricities for multis; Rayleigh(1.5°) mutual inclinations
(Fabrycky+ 2014); equilibrium temperature, insolation, and conservative
habitable-zone membership (simplified Kopparapu+ 2013).

**Stellar noise** — every star carries an astrophysical noise state:
granulation via ν_max scaling (Kallinger+ 2014), p-mode oscillations
(Kjeldsen & Bedding 1995), rotational spot variability with gyrochronology
(Prot ~ √age; McQuillan+ 2014 amplitude span; M dwarfs stay active longer),
and M-dwarf flares. Spot/flare terms are chromatic (Rackham+ 2018): full
strength in the blue optical, ~75% in the TESS band, ~40% in JWST/Roman NIR.
Stellar noise is added in quadrature to every observatory's budget.

**Atmospheres** — each planet gets an atmosphere class (H/He, H/He-rich,
steam, secondary CO₂/N₂, or airless via the Mars-normalized cosmic shoreline
of Zahnle & Catling 2017), scale height, expected transmission feature with a
random cloud/haze suppression factor, TSM and ESM metrics (Kempton+ 2018),
per-JWST-instrument spectroscopy scoring (transits needed for a 5σ
feature detection, with a stellar-contamination noise floor), and a
**per-planet V-band geometric albedo** drawn from class- and
temperature-dependent distributions spanning the measured range — dark
alkali-absorbing hot Jupiters (A_g 0.03–0.11) through bright ammonia/water
cloud decks, a 22% Venus-like branch for secondary atmospheres, and an
icy-bright branch for cold airless worlds (see docs/OBSERVATORIES.md).

**Transits** — Winn (2010) geometry with full eccentricity corrections
(impact parameter, T14/T23 durations, a-priori probability) and a
quadratic-limb-darkened light-curve model integrated per-annulus over exactly
the radial band the planet covers (matches the analytic uniform-source depth
to ~1e-13 relative; equivalent to Mandel & Agol 2002).

**Observability** — per-transit and phase-folded SNR for:

| Observatory | Model | Detection rule |
|---|---|---|
| TESS (1 sector) | noise curve fit to Stassun+ 2018, 60 ppm floor | SNR ≥ 7.1, ≥ 2 transits |
| Kepler (4 yr, archival) | 30 ppm CDPP@Kp=12 photon-scaled, 20 ppm floor | SNR ≥ 7.1, ≥ 3 transits |
| JWST NIRISS SOSS | ~20 ppm/hr @ J=8, saturates J<6.5 | single targeted transit, SNR ≥ 5 |
| JWST NIRSpec Prism | ~12 ppm/hr @ J=11, saturates J<10.5 | single targeted transit, SNR ≥ 5 |
| Ground 1-m survey | 2 mmag + scintillation floor, 90 nights | SNR ≥ 7, ≥ 3 transits, depth > 1 mmag |
| **Roman GBTDS** (2027+) | ~700 ppm/hr @ F146=16, 6×72 d seasons | SNR ≥ 7.1, ≥ 3 transits |
| **HWO imaging** (2040s) | reflected light: C = A_g·Φ·(Rp/a)² with per-planet drawn A_g, floor 3e-11·10^(0.2·max(V−7,0)), 60 mas IWA | host V<11, d<30 pc, resolved, above floor |

The two future facilities are documented — with sources and explicit
"specs will evolve" caveats — in [docs/OBSERVATORIES.md](docs/OBSERVATORIES.md).

## Validation against real exoplanets

`exoverse archive --refresh` snapshots the real confirmed-planet catalog
(NASA Exoplanet Archive TAP when its service is healthy; The Extrasolar
Planets Encyclopaedia as automatic fallback — provenance recorded in the
snapshot). `exoverse validate` then runs two checks:

1. **Physics-rule audit**: every real planet with sufficient measured data is
   run through our generator's hard INVALID rules. Violation rates are
   1–4% and fully attributable to catalog artifacts (transiting brown dwarfs
   cataloged as planets; low-SNR mass measurements) — i.e. our validity
   rules correctly bound reality.
2. **Selection-matched population comparison**: synthetic planets detectable
   by our modeled Kepler are compared to real transit-discovered planets via
   two-sample KS tests. Periods agree (p ≈ 0.9); host-star temperatures
   deliberately disagree (real surveys target FGK stars, our sample is
   IMF-weighted) — the report distinguishes matches from designed mismatches.
3. **System-architecture comparison**: metrics that probe physics the
   generator deliberately leaves out, so the synthetic population acts as a
   *null model*: intra-system radius uniformity ("peas in a pod"), adjacent
   period ratios (resonance pileups), transit multiplicity (the Kepler
   dichotomy), and hot-Neptune-desert occupancy. Real-vs-synthetic residuals
   on these metrics isolate formation/migration signatures that survey
   selection alone cannot produce.

## Browser UI

`exoverse --db worlds.db serve` starts a Flask app (dashboard with population
charts, filterable system list, per-system detail pages with light curves,
observability and atmosphere tables, plausibility flags, and the validation
report). Read-only over the database, ready to be deployed for collaboration.

## The plausibility framework

Every object carries typed flags at three severities:

- **INVALID** — violates physics or robust empirical limits; the generator
  rejects and resamples, so these *never appear in the database*. Rules:
  mass above the 13 M_Jup deuterium-burning limit; bulk density above the
  pure-iron mass–radius curve (Zeng+ 2016) or below 0.03 g/cc; periastron
  inside the fluid Roche limit or grazing the star; crossing orbits; adjacent
  pairs closer than 2√3 mutual Hill radii (Gladman 1993).
- **QUESTIONABLE** — allowed by physics but in regimes where models or data
  are weak; kept, with the doubt recorded. Rules include: super-puff
  densities (0.03–0.1 g/cc); inflated giants without strong irradiation;
  planets in the hot-Neptune desert (Mazeh+ 2016); Hill spacing 3.46–9
  (Gladman-stable but questionable over Gyr, Pu & Wu 2015); giants around
  late M dwarfs (GJ 3512b regime); mega-Earths; ultra-hot atmospheres;
  near-turnoff host stars; metal-poor hosts.
- **INFO** — honest bookkeeping of modeling simplifications (single stars
  only, approximate colors, reduced multiplicity after stability rejection).

`inspect` prints all flags for a system; the `flags` table is queryable
(e.g. `SELECT rule, COUNT(*) FROM flags GROUP BY rule`).

## Findings from the default 500-system population

- ~2.1 planets/star drawn; ~2–3% transit geometrically — matching the naive
  a/R★ expectation for a population dominated by P < 100 d planets.
- **Single-sector TESS detected 0 of 20 transiting planets** in this
  (unbiased, volume-limited-ish) sample: the deepest signals sit on faint
  M dwarfs, and the best bright-star candidate transits only once per
  27-day sector. This reproduces the real reason TESS's yield is dominated
  by short-period planets around bright stars — and is exactly the
  selection effect this project exists to study.
- JWST, pointed at the same transits as targeted single-visit observations,
  recovers ~80% of them; hypothetical-Kepler recovers ~95% given its 4-yr
  baseline.

## Known limitations (deliberate scope)

- Single main-sequence stars only: no binaries, no evolved hosts.
- Occurrence distributions are literature-shaped parametric fits, not fits to
  the actual Kepler/TESS posterior samples.
- No transit-timing variations, no mean-motion resonances, no N-body
  verification of multi-planet stability (Hill-spacing criteria only).
- Observatory systematics are noise floors; no red noise or window functions.
- Stellar noise coefficients are order-of-magnitude Kepler calibrations,
  fine for populations, not for fitting an individual star.
- Atmosphere classes are probabilistic assignments over genuinely degenerate
  bulk compositions (flagged per planet); spectroscopy scoring is a
  scale-height heuristic, not radiative transfer.
- Roman/HWO models are requirement-era and WILL drift as the missions evolve
  (see docs/OBSERVATORIES.md, "Last researched" date).

## Layout

```
src/exoverse/
  constants.py       physical constants (CODATA/IAU)
  flags.py           severity-typed plausibility flag framework
  stars.py           IMF, scaling relations, magnitudes, limb darkening
  stellar_noise.py   granulation / oscillations / spots / flares model
  planets.py         occurrence distributions, M-R relations, validity rules
  atmospheres.py     atmosphere classes, scale heights, TSM/ESM, spectroscopy
  system.py          assembly + Hill/orbit-crossing stability
  transits.py        Winn-2010 geometry + limb-darkened light curves
  observatories.py   TESS / Kepler / JWST / ground / Roman / HWO
  archive.py         real-planet catalog client (NASA TAP + exoplanet.eu)
  validate.py        physics-rule audit + KS population comparison
  database.py        SQLite schema and queries
  generate.py        population orchestration
  cli.py             generate/list/inspect/stats/lightcurve/archive/validate/serve
  web/               Flask browser UI (templates + server-rendered SVG charts)
docs/OBSERVATORIES.md  observatory assumptions, sources, volatility notes
tests/               sanity tests against known real-world values
```
