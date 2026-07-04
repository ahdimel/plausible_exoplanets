"""Kepler DR25 stellar + KOI catalog ingestion for the dichotomy study.

Sources (NASA Exoplanet Archive TAP sync service, same endpoint as
archive.py):

- Stellar : table `keplerstellar` filtered to
            st_delivname='q1_q17_dr25_stellar' — the DR25 stellar
            properties delivery (Mathur et al. 2017, ApJS 229, 30), one row
            per observed Kepler target with per-target duty cycle, data
            span, and the 14-duration rrmscdpp* robust-RMS CDPP values used
            for empirical detection modeling.
- KOI     : table `q1_q17_dr25_koi` — the final uniform DR25 KOI catalog
            (Thompson et al. 2018, ApJS 235, 38), including the Robovetter
            disposition score `koi_score`.

Hardening (learned in archive.py, 2026-07-02): the TAP service can return
whole HTML maintenance pages with HTTP 200, XML error documents, or ORA-
Oracle errors. Every response is validated by its header line and row
count; the stellar fetch is cross-checked against a count(*) query and
falls back to paging by kepid ranges if the sync endpoint ever truncates.

Snapshots are written to data/dr25_stellar.csv and data/dr25_koi.csv
(gitignored); exact query strings, UTC timestamps, and row counts are
recorded in data/PROVENANCE.md (force-added to git).

Fiducial cuts (docs/phase2_design.md; Phase 4 varies them via `Cuts`):

- Stellar: 3900 <= teff <= 7300 K, logg >= 4.0 (FGK dwarfs); non-null
  radius, mass, dataspan, dutycycle and all 14 rrmscdpp columns;
  dataspan > 365 d.
- KOI: koi_disposition in (CONFIRMED, CANDIDATE), koi_score > 0.5,
  0.5 <= koi_period <= 640 d (generator support), koi_prad <= 30 Re,
  and the host kepid passes the stellar cuts.

feh may be null in DR25 (no spectroscopic/photometric constraint) — it
defaults to 0.0 (solar) in KeplerTarget, matching the DR25 pipeline's own
solar-neighborhood prior. kepmag is not cut on; a null kepmag becomes NaN.
"""

from __future__ import annotations

import csv
import io
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from .archive import NASA_TAP, _f

STELLAR_TABLE = "keplerstellar"
STELLAR_DELIVERY = "q1_q17_dr25_stellar"
KOI_TABLE = "q1_q17_dr25_koi"

STELLAR_SNAPSHOT = "dr25_stellar.csv"
KOI_SNAPSHOT = "dr25_koi.csv"

# Transit-duration grid (hours) of the DR25 rrmscdpp* columns, in column
# order: rrmscdpp01p5 ... rrmscdpp15p0 (Mathur et al. 2017, Table 2).
CDPP_DURATIONS_HR = (1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 6.0,
                     7.5, 9.0, 10.5, 12.0, 12.5, 15.0)
CDPP_COLUMNS = tuple(
    f"rrmscdpp{int(d):02d}p{int(round(d * 10)) % 10}" for d in CDPP_DURATIONS_HR
)  # rrmscdpp01p5, rrmscdpp02p0, ..., rrmscdpp15p0

STELLAR_COLUMNS = ("kepid", "teff", "logg", "feh", "radius", "mass",
                   "kepmag", "dutycycle", "dataspan") + CDPP_COLUMNS
KOI_COLUMNS = ("kepoi_name", "kepid", "koi_disposition", "koi_score",
               "koi_period", "koi_prad")


# ------------------------------------------------------------------- cuts
@dataclass(frozen=True)
class Cuts:
    """Sample-selection cuts. FIDUCIAL is the pre-registered default;
    Phase 4 robustness variants use dataclasses.replace on it."""
    teff_min: float = 3900.0        # K, FGK lower edge
    teff_max: float = 7300.0        # K, FGK upper edge
    logg_min: float = 4.0           # dwarfs only
    dataspan_min: float = 365.0     # d, exclude barely-observed targets
    koi_dispositions: tuple = ("CONFIRMED", "CANDIDATE")
    koi_score_min: float = 0.5      # Robovetter disposition score, strict >
    koi_period_min: float = 0.5     # d, generator support
    koi_period_max: float = 640.0   # d, generator support
    koi_prad_max: float = 30.0      # Re


FIDUCIAL = Cuts()


# ------------------------------------------------------------- dataclasses
@dataclass
class KeplerTarget:
    """One DR25 stellar target that passed the cuts.

    cdpp_ppm holds the 14 rrmscdpp values (ppm) in CDPP_DURATIONS_HR
    order. feh is 0.0 when DR25 reports no metallicity (null -> solar);
    kepmag is NaN when missing (never cut on).
    """
    kepid: int
    teff: float
    logg: float
    feh: float
    radius: float      # Rsun
    mass: float        # Msun
    kepmag: float
    dutycycle: float   # fraction of dataspan with valid data
    dataspan: float    # d, first-to-last cadence
    cdpp_ppm: tuple    # 14 floats, CDPP_DURATIONS_HR order


@dataclass
class KOI:
    """One DR25 KOI that passed the cuts (host kepid is the dict key in
    load_koi_systems)."""
    kepoi_name: str
    period: float      # d
    prad: float        # Re
    score: float       # Robovetter disposition score
    disposition: str   # CONFIRMED | CANDIDATE


# ------------------------------------------------------------------- fetch
def _tap_csv(query: str, first_col: str, min_rows: int,
             timeout: float = 300.0, retries: int = 2) -> str:
    """Run a TAP sync query, returning validated CSV text.

    Rejects HTML maintenance pages / XML error documents served with
    HTTP 200 (header check), ORA- Oracle errors, and truncated results
    (row-count check). Retries with backoff like archive.fetch_nasa.
    """
    url = NASA_TAP + "?" + urllib.parse.urlencode(
        {"query": query, "format": "csv"})
    last_err: Exception | None = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            head = text.lstrip()[:500]
            first_line = head.splitlines()[0] if head else ""
            if "ORA-" in head:
                raise RuntimeError(f"TAP returned Oracle error: {head[:200]!r}")
            if not first_line.startswith(first_col):
                raise RuntimeError(f"TAP returned non-CSV response: {head[:200]!r}")
            n_rows = text.count("\n") - 1  # data rows (header excluded)
            if n_rows < min_rows:
                raise RuntimeError(
                    f"TAP returned only {n_rows} rows (expected >= {min_rows}); "
                    "treating as truncated/partial")
            return text
        except Exception as e:  # noqa: BLE001 - report last failure to caller
            last_err = e
            if attempt < retries:
                time.sleep(5.0 * (attempt + 1))
    raise RuntimeError(f"NASA Exoplanet Archive TAP unavailable: {last_err}")


def _tap_scalar(query: str, timeout: float = 120.0) -> float:
    """Run a single-value TAP query (count/min/max) and return the number."""
    text = _tap_csv(query, first_col="n", min_rows=1, timeout=timeout)
    return float(text.strip().splitlines()[1])


def _fetch_stellar_paged(base_query: str, where: str, n_expected: int,
                         n_chunks: int = 10) -> str:
    """Fallback: page the stellar query by kepid ranges if the sync
    endpoint ever caps result size (it did not on 2026-07-03)."""
    lo = int(_tap_scalar(f"select min(kepid) as n from {STELLAR_TABLE} where {where}"))
    hi = int(_tap_scalar(f"select max(kepid) as n from {STELLAR_TABLE} where {where}"))
    edges = [lo + (hi + 1 - lo) * i // n_chunks for i in range(n_chunks + 1)]
    parts: List[str] = []
    for a, b in zip(edges[:-1], edges[1:]):
        text = _tap_csv(
            f"{base_query} and kepid >= {a} and kepid < {b}",
            first_col="kepid", min_rows=0)
        lines = text.splitlines()
        if not parts:
            parts.append(lines[0])          # keep one header
        parts.extend(lines[1:])
    joined = "\n".join(parts) + "\n"
    n_got = joined.count("\n") - 1
    if n_got != n_expected:
        raise RuntimeError(
            f"paged stellar fetch returned {n_got} rows, expected {n_expected}")
    return joined


def stellar_query() -> str:
    """Exact ADQL for the stellar snapshot (recorded in PROVENANCE.md)."""
    return (f"select {','.join(STELLAR_COLUMNS)} from {STELLAR_TABLE} "
            f"where st_delivname='{STELLAR_DELIVERY}'")


def koi_query() -> str:
    """Exact ADQL for the KOI snapshot (recorded in PROVENANCE.md)."""
    return f"select {','.join(KOI_COLUMNS)} from {KOI_TABLE}"


def fetch_dr25(data_dir: str | Path = "data",
               timeout: float = 300.0, retries: int = 2) -> dict:
    """Fetch the DR25 stellar and KOI tables and write CSV snapshots.

    Returns {"stellar_path", "stellar_rows", "koi_path", "koi_rows",
    "fetched_utc"}. Row counts are cross-checked against count(*)
    queries; a mismatch on the (wide) stellar table triggers paging by
    kepid ranges. Update data/PROVENANCE.md whenever this is re-run.
    """
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    fetched = time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime())

    where = f"st_delivname='{STELLAR_DELIVERY}'"
    n_stellar = int(_tap_scalar(
        f"select count(*) as n from {STELLAR_TABLE} where {where}"))
    try:
        stellar_text = _tap_csv(stellar_query(), first_col="kepid",
                                min_rows=n_stellar, timeout=timeout,
                                retries=retries)
    except RuntimeError:
        stellar_text = _fetch_stellar_paged(stellar_query(), where, n_stellar)

    n_koi = int(_tap_scalar(f"select count(*) as n from {KOI_TABLE}"))
    koi_text = _tap_csv(koi_query(), first_col="kepoi_name",
                        min_rows=n_koi, timeout=timeout, retries=retries)

    stellar_path = data_dir / STELLAR_SNAPSHOT
    koi_path = data_dir / KOI_SNAPSHOT
    stellar_path.write_text(
        f"# table={STELLAR_TABLE} delivery={STELLAR_DELIVERY} "
        f"fetched={fetched.replace(' ', 'T')} rows={n_stellar}\n" + stellar_text)
    koi_path.write_text(
        f"# table={KOI_TABLE} fetched={fetched.replace(' ', 'T')} "
        f"rows={n_koi}\n" + koi_text)
    return {"stellar_path": stellar_path, "stellar_rows": n_stellar,
            "koi_path": koi_path, "koi_rows": n_koi, "fetched_utc": fetched}


# -------------------------------------------------------------------- load
def _read_snapshot(path: Path) -> csv.DictReader:
    if not path.exists():
        raise FileNotFoundError(
            f"No DR25 snapshot at {path}. Run exoverse.kepler_data.fetch_dr25().")
    lines = path.read_text().splitlines()
    while lines and lines[0].startswith("#"):
        lines = lines[1:]
    return csv.DictReader(io.StringIO("\n".join(lines)))


def _passes_stellar(teff: Optional[float], logg: Optional[float],
                    radius: Optional[float], mass: Optional[float],
                    dutycycle: Optional[float], dataspan: Optional[float],
                    cdpp: tuple, cuts: Cuts) -> bool:
    """Fiducial stellar cuts: FGK dwarfs with complete noise metadata."""
    if teff is None or not (cuts.teff_min <= teff <= cuts.teff_max):
        return False
    if logg is None or logg < cuts.logg_min:
        return False
    if radius is None or mass is None or dutycycle is None:
        return False
    if dataspan is None or dataspan <= cuts.dataspan_min:
        return False
    if any(c is None for c in cdpp):
        return False
    return True


def load_stellar_targets(data_dir: str | Path = "data",
                         cuts: Cuts = FIDUCIAL) -> List[KeplerTarget]:
    """Load the stellar snapshot and apply the stellar cuts.

    feh null -> 0.0 (solar; DR25 lacks metallicity for many targets).
    Returned in snapshot (kepid) order.
    """
    targets: List[KeplerTarget] = []
    for row in _read_snapshot(Path(data_dir) / STELLAR_SNAPSHOT):
        cdpp = tuple(_f(row[c]) for c in CDPP_COLUMNS)
        teff, logg = _f(row["teff"]), _f(row["logg"])
        radius, mass = _f(row["radius"]), _f(row["mass"])
        dutycycle, dataspan = _f(row["dutycycle"]), _f(row["dataspan"])
        if not _passes_stellar(teff, logg, radius, mass, dutycycle,
                               dataspan, cdpp, cuts):
            continue
        feh = _f(row["feh"])
        kepmag = _f(row["kepmag"])
        targets.append(KeplerTarget(
            kepid=int(row["kepid"]),
            teff=teff, logg=logg,
            feh=0.0 if feh is None else feh,
            radius=radius, mass=mass,
            kepmag=float("nan") if kepmag is None else kepmag,
            dutycycle=dutycycle, dataspan=dataspan,
            cdpp_ppm=cdpp))
    return targets


def load_koi_systems(data_dir: str | Path = "data",
                     cuts: Cuts = FIDUCIAL) -> Dict[int, List[KOI]]:
    """Load the KOI snapshot, apply the KOI cuts, and keep only hosts
    that pass the stellar cuts. Returns {kepid: [KOI, ...]} with each
    system's KOIs sorted by period."""
    good_hosts = {t.kepid for t in load_stellar_targets(data_dir, cuts)}
    systems: Dict[int, List[KOI]] = {}
    for row in _read_snapshot(Path(data_dir) / KOI_SNAPSHOT):
        disposition = (row["koi_disposition"] or "").strip()
        if disposition not in cuts.koi_dispositions:
            continue
        score = _f(row["koi_score"])
        if score is None or score <= cuts.koi_score_min:
            continue
        period = _f(row["koi_period"])
        if period is None or not (cuts.koi_period_min <= period
                                  <= cuts.koi_period_max):
            continue
        prad = _f(row["koi_prad"])
        if prad is None or prad > cuts.koi_prad_max:
            continue
        kepid = int(row["kepid"])
        if kepid not in good_hosts:
            continue
        systems.setdefault(kepid, []).append(KOI(
            kepoi_name=(row["kepoi_name"] or "").strip(),
            period=period, prad=prad, score=score,
            disposition=disposition))
    for kois in systems.values():
        kois.sort(key=lambda k: k.period)
    return systems


__all__ = [
    "CDPP_COLUMNS", "CDPP_DURATIONS_HR", "Cuts", "FIDUCIAL", "KOI",
    "KOI_SNAPSHOT", "KeplerTarget", "STELLAR_SNAPSHOT", "fetch_dr25",
    "koi_query", "load_koi_systems", "load_stellar_targets",
    "stellar_query",
]
