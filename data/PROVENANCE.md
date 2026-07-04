# DR25 snapshot provenance

Data source: NASA Exoplanet Archive TAP sync service
(`https://exoplanetarchive.ipac.caltech.edu/TAP/sync`, `format=csv`),
fetched by `exoverse.kepler_data.fetch_dr25()`.

- **Retrieved:** 2026-07-04 01:46:12 UTC (single un-paged sync fetch;
  row counts cross-checked against `count(*)` queries)
- **Stellar reference:** Mathur et al. 2017, ApJS 229, 30 (DR25 stellar
  properties, delivery `q1_q17_dr25_stellar`)
- **KOI reference:** Thompson et al. 2018, ApJS 235, 38 (DR25 uniform KOI
  catalog with Robovetter disposition scores)

## Exact queries

`data/dr25_stellar.csv` (200,038 rows; count(*) verified 200,038):

```sql
select kepid,teff,logg,feh,radius,mass,kepmag,dutycycle,dataspan,
       rrmscdpp01p5,rrmscdpp02p0,rrmscdpp02p5,rrmscdpp03p0,rrmscdpp03p5,
       rrmscdpp04p5,rrmscdpp05p0,rrmscdpp06p0,rrmscdpp07p5,rrmscdpp09p0,
       rrmscdpp10p5,rrmscdpp12p0,rrmscdpp12p5,rrmscdpp15p0
from keplerstellar where st_delivname='q1_q17_dr25_stellar'
```

(submitted as a single line, columns comma-joined without whitespace)

`data/dr25_koi.csv` (8,054 rows; count(*) verified 8,054):

```sql
select kepoi_name,kepid,koi_disposition,koi_score,koi_period,koi_prad
from q1_q17_dr25_koi
```

## Fiducial cuts (docs/phase2_design.md; constants in `kepler_data.Cuts`)

Stellar (FGK dwarfs with complete detection metadata):

1. `3900 <= teff <= 7300` K
2. `logg >= 4.0`
3. non-null `radius`, `mass`, `dutycycle`, `dataspan`, and all 14
   `rrmscdpp*` columns
4. `dataspan > 365` d (strict)

KOI: `koi_disposition in (CONFIRMED, CANDIDATE)`, `koi_score > 0.5`
(strict), `0.5 <= koi_period <= 640` d, `koi_prad <= 30` Re (non-null),
host kepid passes the stellar cuts.

`feh` null defaults to 0.0 (solar) — note the 2026-07-04 snapshot has zero
null feh values (DR25 fills feh from fits/priors), so the default is
defensive. `kepmag` is not cut on (null -> NaN).

## Row counts through the cuts (2026-07-04 snapshot)

Stellar (cuts applied cumulatively, in order):

| stage | rows |
|---|---|
| raw (`q1_q17_dr25_stellar`) | 200,038 |
| teff in [3900, 7300] | 190,688 |
| + logg >= 4.0 | 152,370 |
| + non-null radius/mass/dutycycle/dataspan/cdpp | 148,725 |
| + dataspan > 365 d | **137,493** |

All 200,038 kepids are unique. The 3,645 rows dropped for nulls are 2,855
mass-only nulls plus 790 with no long-cadence metadata (dutycycle,
dataspan, and all cdpp null together, some also missing mass).

KOI (cuts applied cumulatively, in order):

| stage | rows |
|---|---|
| raw (`q1_q17_dr25_koi`) | 8,054 |
| disposition CONFIRMED/CANDIDATE | 4,089 |
| + koi_score > 0.5 | 3,914 |
| + 0.5 <= period <= 640 d | 3,904 |
| + prad <= 30 Re (non-null) | 3,842 |
| + host passes stellar cuts | **3,400** |

The 3,400 surviving KOIs group into **2,547 systems**; observed
multiplicity N_k for k = 1..5, 6+: **1968, 389, 127, 45, 15, 3**
(singles-per-multi = 1968/579 = 3.40).

Score-cut variants (all other cuts fiducial):

| score cut | KOIs | systems | N_k (1..5, 6+) | singles/multi |
|---|---|---|---|---|
| > 0.5 (fiducial) | 3,400 | 2,547 | 1968, 389, 127, 45, 15, 3 | 3.40 |
| > 0.0 | 3,506 | 2,631 | 2041, 390, 136, 46, 15, 3 | 3.46 |
| > 0.9 | 3,021 | 2,303 | 1803, 340, 115, 34, 9, 2 | 3.61 |

Adjacent-pair size uniformity under the fiducial cuts (pairs adjacent in
period within a system, both radii < 30 Re): **853 pairs**, median
|log10(R_out/R_in)| = **0.122**.

Snapshots (`dr25_stellar.csv`, `dr25_koi.csv`) are gitignored; only this
provenance file is committed (`git add -f data/PROVENANCE.md`). Re-running
`fetch_dr25()` must be followed by an update of this file.

---

`archive_snapshot.csv` (real confirmed-planet catalog for validation) is
documented in its own header line and `src/exoverse/archive.py`.
