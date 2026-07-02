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

## Quick start

```bash
python3 -m venv .venv && .venv/bin/pip install -e ".[dev]"

.venv/bin/exoverse --db worlds.db generate --n 500 --seed 42
.venv/bin/exoverse --db worlds.db stats
.venv/bin/exoverse --db worlds.db list --transiting
.venv/bin/exoverse --db worlds.db inspect PXS-42-00481      # full metadata + flags
.venv/bin/exoverse --db worlds.db lightcurve PXS-42-00481 d # PNG light curve

.venv/bin/python -m pytest      # physics sanity tests vs known real values
```

Generation is fully deterministic per seed: any stored system can be
re-generated bit-for-bit from its `(seed, name)`, which is how `lightcurve`
re-simulates without storing arrays.

## What gets generated

**Stars** — Kroupa (2001) IMF over 0.08–2.2 M☉; main-sequence L(M), R(M)
scaling relations with scatter; T_eff from Stefan–Boltzmann; thin-disk
metallicity distribution; distances from uniform local density (≤300 pc);
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

## Known limitations (deliberate v1 scope)

- Single main-sequence stars only: no binaries, no evolved hosts, no spots or
  stellar activity noise.
- Occurrence distributions are literature-shaped parametric fits, not fits to
  the actual Kepler/TESS posterior samples.
- No transit-timing variations, no mean-motion resonances, no N-body
  verification of multi-planet stability (Hill-spacing criteria only).
- Observatory systematics are single noise floors; no red noise, no
  window-function/duty-cycle modeling beyond a scalar.
- No transmission spectroscopy yet — JWST is modeled as a white-light
  transit detector; atmosphere characterization SNR is the natural next step.

## Layout

```
src/exoverse/
  constants.py       physical constants (CODATA/IAU)
  flags.py           severity-typed plausibility flag framework
  stars.py           IMF, scaling relations, magnitudes, limb darkening
  planets.py         occurrence distributions, M-R relations, validity rules
  system.py          assembly + Hill/orbit-crossing stability
  transits.py        Winn-2010 geometry + limb-darkened light curves
  observatories.py   TESS / Kepler / JWST / ground noise + detection
  database.py        SQLite schema and queries
  generate.py        population orchestration
  cli.py             generate / list / inspect / stats / lightcurve
tests/               sanity tests against known real-world values
```
