"""Validation of the generator against real, confirmed exoplanets.

Two independent checks:

1. Physics-rule audit ("do real planets pass our INVALID rules?")
   Every confirmed planet with sufficient measured data is run through the
   same hard validity rules the generator enforces (deuterium-burning mass
   limit, pure-iron density bound, super-puff floor, Roche limit). If our
   rules rejected a meaningful fraction of *real* planets, the rules would
   be wrong. A ~1% violation rate is expected and healthy: it reflects
   provisional/erroneous catalog measurements (and the known brown-dwarf
   boundary ambiguity), not generator physics.

2. Population comparison ("does our synthetic detected population look like
   the real detected population?")
   The raw archive is heavily selection-biased, so raw-vs-raw comparisons
   are meaningless. Instead we forward-model the selection: our synthetic
   planets are filtered to those DETECTABLE by the modeled Kepler survey,
   and compared against the real *transit-discovered* archive planets on
   log period, log radius, and host Teff using two-sample KS statistics and
   quantile tables. Agreement will never be perfect (Kepler's real target
   selection favored FGK stars; ours is IMF-weighted and volume-limited) -
   the report says what matches and what deliberately does not.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import List

from .archive import ArchivePlanet, load_snapshot
from .constants import AU, M_EARTH, M_JUP, R_SUN
from .database import WorldDB
from .planets import pure_iron_radius

M_JUP_IN_EARTH = M_JUP / M_EARTH


# --------------------------------------------------------------- KS utilities
def ks_2samp(x1: List[float], x2: List[float]) -> tuple[float, float]:
    """Two-sample Kolmogorov-Smirnov statistic and asymptotic p-value."""
    a, b = sorted(x1), sorted(x2)
    n1, n2 = len(a), len(b)
    if n1 == 0 or n2 == 0:
        return float("nan"), float("nan")
    i = j = 0
    d = 0.0
    while i < n1 and j < n2:
        if a[i] <= b[j]:
            i += 1
        else:
            j += 1
        d = max(d, abs(i / n1 - j / n2))
    n_eff = n1 * n2 / (n1 + n2)
    lam = (math.sqrt(n_eff) + 0.12 + 0.11 / math.sqrt(n_eff)) * d
    p = 2.0 * sum((-1) ** (k - 1) * math.exp(-2.0 * k * k * lam * lam)
                  for k in range(1, 101))
    return d, max(0.0, min(1.0, p))


def quantiles(xs: List[float], qs=(0.1, 0.25, 0.5, 0.75, 0.9)) -> dict:
    if not xs:
        return {}
    s = sorted(xs)
    out = {}
    for q in qs:
        idx = q * (len(s) - 1)
        lo, hi = int(idx), min(int(idx) + 1, len(s) - 1)
        out[f"q{int(q*100)}"] = round(s[lo] + (s[hi] - s[lo]) * (idx - lo), 4)
    return out


# ----------------------------------------------------------- rule audit
@dataclass
class RuleAudit:
    rule: str
    n_checked: int
    n_violations: int
    examples: List[str] = field(default_factory=list)
    interpretation: str = ""


def audit_rules(planets: List[ArchivePlanet]) -> List[RuleAudit]:
    audits: List[RuleAudit] = []

    # Deuterium-burning mass limit. Measured masses only: Msini lower limits
    # cannot prove a violation, and exoplanet.eu deliberately catalogs
    # objects up to ~60 Mjup (brown-dwarf candidates NASA would exclude).
    checked = [p for p in planets
               if p.pl_bmasse is not None and p.mass_provenance == "measured"
               and p.tran_flag == 1.0]
    viol = [p for p in checked if p.pl_bmasse > 13.0 * M_JUP_IN_EARTH]
    audits.append(RuleAudit(
        "mass.deuterium_burning", len(checked), len(viol),
        [p.pl_name for p in viol[:8]],
        "Restricted to transiting, >2-sigma masses. Violators are the known "
        "transiting brown dwarfs cataloged alongside planets (KELT-1b, "
        "CoRoT-3b class); our generator excludes this regime by design. "
        "(exoplanet.eu additionally lists imaged companions to ~60 Mjup, "
        "excluded here.)"))

    # Denser than pure iron. Measured-mass transiting planets only: Msini
    # substitutions and provisional RV masses dominate false violations.
    checked = [p for p in planets
               if p.pl_bmasse is not None and p.pl_rade is not None
               and p.pl_rade < 5.0 and p.pl_bmasse < 13.0 * M_JUP_IN_EARTH
               and p.mass_provenance == "measured" and p.tran_flag == 1.0]
    viol = [p for p in checked if p.pl_rade < pure_iron_radius(p.pl_bmasse)]
    audits.append(RuleAudit(
        "density.exceeds_pure_iron", len(checked), len(viol),
        [p.pl_name for p in viol[:8]],
        "Remaining violations are planets with low-SNR mass measurements; "
        "a ~1% rate reflects catalog measurement error, not physics."))

    # Super-puff density floor
    viol = []
    n_checked = 0
    for p in checked:
        rho = (p.pl_bmasse * M_EARTH) / (
            4.0 / 3.0 * math.pi * (p.pl_rade * 6.3781e6) ** 3) / 1000.0
        n_checked += 1
        if rho < 0.03:
            viol.append(p)
    audits.append(RuleAudit(
        "density.below_superpuff_floor", n_checked, len(viol),
        [p.pl_name for p in viol[:8]],
        "Densities below 0.03 g/cc in the catalog are usually inflated-radius "
        "young planets with weak mass limits."))

    # Periastron inside ~1.5 stellar radii (grazes star)
    checked = [p for p in planets
               if None not in (p.pl_orbsmax, p.pl_orbeccen, p.st_rad)]
    viol = []
    for p in checked:
        peri_m = p.pl_orbsmax * AU * (1.0 - p.pl_orbeccen)
        if peri_m < 1.5 * p.st_rad * R_SUN:
            viol.append(p)
    audits.append(RuleAudit(
        "orbit.grazes_star", len(checked), len(viol),
        [p.pl_name for p in viol[:8]],
        "Real planets this close would be actively disrupted; violations "
        "indicate inconsistent catalog (a, e, R*) combinations."))
    return audits


# ------------------------------------------------- population comparison
@dataclass
class DistComparison:
    quantity: str
    n_synthetic: int
    n_real: int
    ks_stat: float
    ks_pvalue: float
    synthetic_quantiles: dict
    real_quantiles: dict
    note: str = ""


def compare_populations(db_path: str, planets: List[ArchivePlanet]
                        ) -> List[DistComparison]:
    db = WorldDB(db_path)
    # Synthetic planets detectable by the modeled Kepler survey (the least
    # noise-limited long-baseline instrument we model = best analogy to the
    # real transit-discovered census)
    rows = db.conn.execute(
        """SELECT p.period_d, p.radius_re, s.st_teff FROM planets p
           JOIN systems s ON s.id = p.system_id
           JOIN observations o ON o.planet_id = p.id
           WHERE o.observatory LIKE 'Kepler%' AND o.detectable = 1""").fetchall()
    db.close()
    syn_p = [math.log10(r[0]) for r in rows]
    syn_r = [math.log10(r[1]) for r in rows]
    syn_t = [r[2] for r in rows]

    real = [p for p in planets if p.tran_flag == 1.0
            and "transit" in p.discoverymethod.lower()]
    real_p = [math.log10(p.pl_orbper) for p in real
              if p.pl_orbper and 0.3 < p.pl_orbper < 1000]
    real_r = [math.log10(p.pl_rade) for p in real
              if p.pl_rade and 0.3 < p.pl_rade < 25]
    real_t = [p.st_teff for p in real if p.st_teff and 2300 < p.st_teff < 8500]

    out = []
    for name, syn, rl, note in (
        ("log10(period_days)", syn_p, real_p,
         "Kepler-detectable synthetic vs transit-discovered real planets"),
        ("log10(radius_Re)", syn_r, real_r,
         "Radius valley position and giant fraction are the key features"),
        ("host_teff_K", syn_t, real_t,
         "Expected mismatch: real surveys target FGK; our sample is "
         "IMF-weighted (M-dwarf heavy) and volume-limited"),
    ):
        d, pval = ks_2samp(syn, rl)
        out.append(DistComparison(name, len(syn), len(rl), round(d, 4),
                                  round(pval, 6), quantiles(syn),
                                  quantiles(rl), note))
    return out


# ------------------------------------------- system-architecture comparison
@dataclass
class StructureComparison:
    metric: str
    synthetic: dict
    real: dict
    ks_stat: float
    ks_pvalue: float
    note: str = ""


def _adjacent_pairs(groups: dict) -> tuple[List[float], List[float]]:
    """Given {system: [(period, radius), ...]}, return adjacent-pair period
    ratios (outer/inner) and |log10 radius ratios|, sorted by period."""
    pratios: List[float] = []
    rratios: List[float] = []
    for members in groups.values():
        members = sorted(m for m in members if m[0] and m[0] > 0)
        for inner, outer in zip(members, members[1:]):
            pratios.append(outer[0] / inner[0])
            if inner[1] and outer[1] and inner[1] > 0 and outer[1] > 0:
                rratios.append(abs(math.log10(outer[1] / inner[1])))
    return pratios, rratios


def _multiplicity_hist(groups: dict) -> dict:
    hist = {"1": 0, "2": 0, "3+": 0}
    for members in groups.values():
        n = len(members)
        hist["1" if n == 1 else "2" if n == 2 else "3+"] += 1
    return hist


def _desert_fraction(planets_pr: List[tuple]) -> tuple[int, int]:
    """(count in hot-Neptune desert, total) for [(period, radius), ...].
    Desert box: P < 4 d and 3 <= R < 8 Re (Mazeh+ 2016 core)."""
    with_data = [(p, r) for p, r in planets_pr if p and r]
    hits = sum(1 for p, r in with_data if p < 4.0 and 3.0 <= r < 8.0)
    return hits, len(with_data)


def compare_structure(db_path: str, planets: List[ArchivePlanet]
                      ) -> List[StructureComparison]:
    """Selection-matched comparison of system ARCHITECTURE, not just 1-D
    occurrence: these metrics probe physics the generator deliberately does
    not put in (formation memory, migration), so real-vs-synthetic residuals
    isolate signatures that survey selection alone cannot explain.

    - peas-in-a-pod: real Kepler multis have correlated planet sizes
      (Weiss+ 2018); our planets are drawn independently, so the synthetic
      dispersion of adjacent radius ratios is the *no-formation-memory null*.
    - period ratios: real multis pile up just wide of 2:1 and 3:2 resonances
      (Fabrycky+ 2014, migration + tidal damping); our spacing is smooth
      (only Hill-stability carving), so any real excess is dynamics.
    - multiplicity (the 'Kepler dichotomy'): tests the single Rayleigh(1.5
      deg) mutual-inclination population against reality's apparent excess
      of single-transiting systems.
    - hot-Neptune desert: we only flag the desert, we don't deplete it;
      real photoevaporation/tidal stripping does. The occupancy gap measures
      that missing physics in the DETECTED population."""
    db = WorldDB(db_path)
    rows = db.conn.execute(
        """SELECT p.system_id, p.period_d, p.radius_re FROM planets p
           JOIN observations o ON o.planet_id = p.id
           WHERE o.observatory LIKE 'Kepler%' AND o.detectable = 1
           ORDER BY p.system_id, p.period_d""").fetchall()
    db.close()
    syn_groups: dict = {}
    for sid, per, rad in rows:
        syn_groups.setdefault(sid, []).append((per, rad))

    real_groups: dict = {}
    for p in planets:
        if p.tran_flag == 1.0 and "transit" in p.discoverymethod.lower():
            real_groups.setdefault(p.hostname or p.pl_name, []).append(
                (p.pl_orbper, p.pl_rade))

    syn_pr, syn_rr = _adjacent_pairs(syn_groups)
    real_pr, real_rr = _adjacent_pairs(real_groups)
    syn_pr_in = [x for x in syn_pr if 1.0 < x < 5.0]
    real_pr_in = [x for x in real_pr if 1.0 < x < 5.0]

    out: List[StructureComparison] = []

    d, pval = ks_2samp(syn_rr, real_rr)
    out.append(StructureComparison(
        "peas_in_a_pod_radius_uniformity",
        {"n_pairs": len(syn_rr), **quantiles(syn_rr)},
        {"n_pairs": len(real_rr), **quantiles(real_rr)},
        round(d, 4), round(pval, 6),
        "|log10(R_out/R_in)| of adjacent detected pairs. Synthetic = "
        "independent-draw null; if real medians are smaller, real systems "
        "have intra-system size uniformity (formation memory, Weiss+ 2018) "
        "that selection effects alone cannot produce."))

    def near(vals, lo, hi):
        return sum(1 for x in vals if lo <= x < hi)

    def res_stats(vals):
        n = len(vals)
        return {
            "n_pairs": n,
            "frac_near_3to2": round(near(vals, 1.5, 1.6) / n, 4) if n else 0.0,
            "frac_near_2to1": round(near(vals, 2.0, 2.1) / n, 4) if n else 0.0,
            **quantiles(vals),
        }

    d, pval = ks_2samp(syn_pr_in, real_pr_in)
    out.append(StructureComparison(
        "adjacent_period_ratios",
        res_stats(syn_pr_in), res_stats(real_pr_in),
        round(d, 4), round(pval, 6),
        "Period ratios of adjacent detected pairs (1-5). Synthetic spacing "
        "is smooth (Hill stability only); a real excess just wide of 3:2 "
        "and 2:1 is the fingerprint of migration + resonance capture."))

    syn_mult = _multiplicity_hist(syn_groups)
    real_mult = _multiplicity_hist(real_groups)

    def dichotomy(h):
        multi = h["2"] + h["3+"]
        return round(h["1"] / multi, 2) if multi else float("inf")

    out.append(StructureComparison(
        "transit_multiplicity",
        {**syn_mult, "singles_per_multi": dichotomy(syn_mult)},
        {**real_mult, "singles_per_multi": dichotomy(real_mult)},
        float("nan"), float("nan"),
        "Detected-planets-per-system histogram (the 'Kepler dichotomy' "
        "test). Our single Rayleigh(1.5 deg) mutual-inclination population "
        "predicts fewer observed singles than reality if real systems have "
        "a second, dynamically hot component."))

    syn_hits, syn_tot = _desert_fraction(
        [(p, r) for members in syn_groups.values() for p, r in members])
    real_hits, real_tot = _desert_fraction(
        [(p, r) for members in real_groups.values() for p, r in members])
    out.append(StructureComparison(
        "hot_neptune_desert_occupancy",
        {"in_desert": syn_hits, "total": syn_tot,
         "fraction": round(syn_hits / syn_tot, 4) if syn_tot else 0.0},
        {"in_desert": real_hits, "total": real_tot,
         "fraction": round(real_hits / real_tot, 4) if real_tot else 0.0},
        float("nan"), float("nan"),
        "Fraction of detected planets with P < 4 d and 3 <= R < 8 Re. The "
        "generator flags but does not deplete the desert; the synthetic "
        "excess over reality measures the strength of the real "
        "photoevaporation/tidal-stripping sink (Mazeh+ 2016)."))
    return out


# --------------------------------------------------------------- reporting
def run_validation(db_path: str, data_dir: str | Path = "data") -> dict:
    planets, meta = load_snapshot(data_dir)
    report = {
        "archive": meta,
        "rule_audits": [asdict(a) for a in audit_rules(planets)],
        "distributions": [asdict(c) for c in compare_populations(db_path, planets)],
        "structure": [asdict(c) for c in compare_structure(db_path, planets)],
    }
    out_path = Path(data_dir) / "validation_report.json"
    out_path.write_text(json.dumps(report, indent=2))
    report["report_path"] = str(out_path)
    return report


def format_report(report: dict) -> str:
    lines = []
    a = report["archive"]
    lines.append(f"Catalog: {a.get('count')} confirmed planets "
                 f"(source={a.get('source')}, fetched={a.get('fetched')})")
    lines.append("\n-- Physics-rule audit (our INVALID rules vs real planets) --")
    for r in report["rule_audits"]:
        frac = r["n_violations"] / r["n_checked"] * 100 if r["n_checked"] else 0.0
        lines.append(f"  {r['rule']:<32} {r['n_violations']:>4}/{r['n_checked']:<5} "
                     f"({frac:.2f}%) violations")
        if r["examples"]:
            lines.append(f"      e.g. {', '.join(r['examples'][:5])}")
        lines.append(f"      -> {r['interpretation']}")
    lines.append("\n-- Population comparison (selection-matched) --")
    for c in report["distributions"]:
        lines.append(f"  {c['quantity']:<22} KS D={c['ks_stat']:.3f} "
                     f"p={c['ks_pvalue']:.4f}  "
                     f"(n_syn={c['n_synthetic']}, n_real={c['n_real']})")
        lines.append(f"      synthetic median={c['synthetic_quantiles'].get('q50')}"
                     f"  real median={c['real_quantiles'].get('q50')}")
        if c["note"]:
            lines.append(f"      note: {c['note']}")
    lines.append("\n-- System architecture (physics the generator leaves out"
                 " on purpose) --")
    for c in report.get("structure", []):
        ks = ""
        if c["ks_stat"] == c["ks_stat"]:   # not NaN
            ks = f"  KS D={c['ks_stat']:.3f} p={c['ks_pvalue']:.4f}"
        lines.append(f"  {c['metric']}{ks}")
        lines.append(f"      synthetic: {c['synthetic']}")
        lines.append(f"      real     : {c['real']}")
        if c["note"]:
            lines.append(f"      note: {c['note']}")
    return "\n".join(lines)
