"""Client for real exoplanet catalogs, used to validate the generator.

Primary source : NASA Exoplanet Archive TAP service (table `pscomppars`,
                 one composite row per confirmed planet).
                 https://exoplanetarchive.ipac.caltech.edu/TAP
Fallback source: The Extrasolar Planets Encyclopaedia (exoplanet.eu) full
                 catalog CSV, filtered to status "Confirmed", with columns
                 mapped (and Jupiter units converted) to the NASA schema.

The NASA TAP service returned server-side errors (ORA-04063, broken view)
throughout initial development (2026-07-02); the client therefore tries NASA
with retries and falls back automatically. `source` is recorded in the
snapshot header so the provenance is always known.

Snapshots are cached at data/archive_snapshot.csv with canonical columns:

    pl_name, hostname, pl_orbper (d), pl_rade (Re), pl_bmasse (Me),
    pl_orbeccen, pl_eqt (K), pl_orbsmax (AU), st_teff (K), st_rad (Rsun),
    st_mass (Msun), st_met (dex), sy_dist (pc), sy_vmag, sy_jmag,
    discoverymethod, disc_year, tran_flag (0/1)
"""

from __future__ import annotations

import csv
import io
import math
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, fields
from pathlib import Path
from typing import List, Optional

M_JUP_IN_EARTH = 317.83
R_JUP_IN_EARTH = 11.209

NASA_TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
EU_CSV = "https://exoplanet.eu/catalog/csv/"

CANONICAL_COLUMNS = [
    "pl_name", "hostname", "pl_orbper", "pl_rade", "pl_bmasse", "pl_orbeccen",
    "pl_eqt", "pl_orbsmax", "st_teff", "st_rad", "st_mass", "st_met",
    "sy_dist", "sy_vmag", "sy_jmag", "discoverymethod", "disc_year", "tran_flag",
    "mass_provenance",   # "measured" | "msini" | "" (unknown)
]


@dataclass
class ArchivePlanet:
    pl_name: str
    hostname: str
    pl_orbper: Optional[float]
    pl_rade: Optional[float]
    pl_bmasse: Optional[float]
    pl_orbeccen: Optional[float]
    pl_eqt: Optional[float]
    pl_orbsmax: Optional[float]
    st_teff: Optional[float]
    st_rad: Optional[float]
    st_mass: Optional[float]
    st_met: Optional[float]
    sy_dist: Optional[float]
    sy_vmag: Optional[float]
    sy_jmag: Optional[float]
    discoverymethod: str
    disc_year: Optional[float]
    tran_flag: Optional[float]
    mass_provenance: str = ""


def _f(value: str) -> Optional[float]:
    value = (value or "").strip()
    if not value:
        return None
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except ValueError:
        return None


# ---------------------------------------------------------------------- NASA
def fetch_nasa(timeout: float = 120.0, retries: int = 2) -> str:
    """Fetch pscomppars from NASA TAP; returns canonical CSV text."""
    # pl_bmassprov distinguishes true masses from Msini lower limits
    cols = ("pl_name,hostname,pl_orbper,pl_rade,pl_bmasse,pl_orbeccen,pl_eqt,"
            "pl_orbsmax,st_teff,st_rad,st_mass,st_met,sy_dist,sy_vmag,sy_jmag,"
            "discoverymethod,disc_year,tran_flag,pl_bmassprov as mass_provenance")
    query = f"select {cols} from pscomppars"
    url = NASA_TAP + "?" + urllib.parse.urlencode({"query": query, "format": "csv"})
    last_err: Exception | None = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                text = resp.read().decode("utf-8", errors="replace")
            # The service can return XML error documents or whole HTML
            # maintenance pages with HTTP 200 - accept only real CSV.
            first_line = text.lstrip().splitlines()[0] if text.strip() else ""
            if not first_line.startswith("pl_name"):
                raise RuntimeError(f"TAP returned non-CSV response: {text[:200]!r}")
            n_rows = text.count("\n")
            if n_rows < 3000:
                raise RuntimeError(
                    f"TAP returned only ~{n_rows} rows (expected >5000); "
                    "treating as truncated/partial")
            return text
        except Exception as e:  # noqa: BLE001 - report last failure to caller
            last_err = e
            if attempt < retries:
                time.sleep(5.0 * (attempt + 1))
    raise RuntimeError(f"NASA Exoplanet Archive TAP unavailable: {last_err}")


# ---------------------------------------------------------------- exoplanet.eu
def fetch_eu(timeout: float = 180.0) -> str:
    """Fetch exoplanet.eu catalog and convert to canonical CSV text."""
    req = urllib.request.Request(EU_CSV, headers={"User-Agent": "exoverse/0.2"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read().decode("utf-8", errors="replace")

    reader = csv.DictReader(io.StringIO(raw))
    out = io.StringIO()
    writer = csv.writer(out)
    writer.writerow(CANONICAL_COLUMNS)
    for row in reader:
        if (row.get("planet_status") or "").strip().lower() != "confirmed":
            continue
        mass_meas = _f(row.get("mass", ""))
        mass_mj = mass_meas if mass_meas is not None else _f(row.get("mass_sini", ""))
        # exoplanet.eu's CSV stores upper limits indistinguishably in `mass`;
        # a mass only counts as "measured" if its errors exist and imply a
        # >2-sigma detection (error < 50% of the value)
        if mass_meas is not None:
            err = max(abs(_f(row.get("mass_error_min", "")) or 0.0),
                      abs(_f(row.get("mass_error_max", "")) or 0.0))
            provenance = "measured" if 0.0 < err < 0.5 * mass_meas else "uncertain"
        elif mass_mj is not None:
            provenance = "msini"
        else:
            provenance = ""
        radius_rj = _f(row.get("radius", ""))
        teq = _f(row.get("temp_calculated", "")) or _f(row.get("temp_measured", ""))
        method = (row.get("detection_type") or "").strip()
        tran = 1.0 if "transit" in method.lower() else 0.0
        writer.writerow([
            (row.get("name") or "").strip(),
            (row.get("star_name") or "").strip(),
            _f(row.get("orbital_period", "")),
            radius_rj * R_JUP_IN_EARTH if radius_rj is not None else None,
            mass_mj * M_JUP_IN_EARTH if mass_mj is not None else None,
            _f(row.get("eccentricity", "")),
            teq,
            _f(row.get("semi_major_axis", "")),
            _f(row.get("star_teff", "")),
            _f(row.get("star_radius", "")),
            _f(row.get("star_mass", "")),
            _f(row.get("star_metallicity", "")),
            _f(row.get("star_distance", "")),
            _f(row.get("mag_v", "")),
            _f(row.get("mag_j", "")),
            method,
            _f(row.get("discovered", "")),
            tran,
            provenance,
        ])
    return out.getvalue()


# ------------------------------------------------------------------ interface
def refresh_snapshot(data_dir: str | Path = "data",
                     prefer: str = "nasa") -> tuple[Path, str]:
    """Fetch a fresh catalog snapshot. Returns (path, source_used)."""
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    path = data_dir / "archive_snapshot.csv"

    source_used = None
    text = None
    errors = []
    order = ["nasa", "eu"] if prefer == "nasa" else ["eu", "nasa"]
    for source in order:
        try:
            text = fetch_nasa() if source == "nasa" else fetch_eu()
            source_used = source
            break
        except Exception as e:  # noqa: BLE001
            errors.append(f"{source}: {e}")
    if text is None:
        raise RuntimeError("All catalog sources failed: " + " | ".join(errors))

    header = (f"# source={source_used} fetched={time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    path.write_text(header + text)
    return path, source_used


def load_snapshot(data_dir: str | Path = "data") -> tuple[List[ArchivePlanet], dict]:
    """Load the cached snapshot. Returns (planets, provenance dict)."""
    path = Path(data_dir) / "archive_snapshot.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"No archive snapshot at {path}. Run: exoverse archive --refresh")
    lines = path.read_text().splitlines()
    meta = {"source": "unknown", "fetched": "unknown"}
    if lines and lines[0].startswith("#"):
        for part in lines[0].lstrip("# ").split():
            if "=" in part:
                k, v = part.split("=", 1)
                meta[k] = v
        lines = lines[1:]
    reader = csv.DictReader(io.StringIO("\n".join(lines)))
    planets: List[ArchivePlanet] = []
    float_fields = {f.name for f in fields(ArchivePlanet)} - {
        "pl_name", "hostname", "discoverymethod", "mass_provenance"}
    for row in reader:
        kwargs = {}
        for f_ in fields(ArchivePlanet):
            v = row.get(f_.name, "")
            kwargs[f_.name] = _f(v) if f_.name in float_fields else (v or "")
        planets.append(ArchivePlanet(**kwargs))
    meta["count"] = len(planets)
    return planets, meta
