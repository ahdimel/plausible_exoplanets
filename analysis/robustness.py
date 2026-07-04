"""Phase 4 pre-registered robustness variants (docs/robustness_plan.md).

Each variant re-derives the (sigma_r, sigma_i) credible region with the
same likelihood + cubic-smoothing machinery as grid_inference.py, over
the fiducial valley cells (Delta 2lnL < 100 in the fiducial fit). Three
kinds:

  metric : rescore the cached fiducial simulations (V6 truncation,
           V7 size-observable swaps, V8 Poisson occurrence)
  real   : rescore cached simulations against a re-cut real catalog
           (V1 robovetter score cuts, V5 radius-error injection)
  resim  : re-simulate the valley cells with changed stellar cuts or
           detection rule (V2 a-c, V3 a-c, V4), paired seeds with the
           fiducial grid (same cell_seed(i_sr, i_si, m) stream)

The fiducial comparison uses m=0..7 only, matching the variants' M=8,
so all comparisons are seed-paired and sample-size-matched.

Note on V7: the plan named AD-vs-KS on the *distance* metric; the
inference layer replaced KS with a fixed-bandwidth KDE likelihood, so
the corresponding sensitivity variants here are KDE bandwidth x0.5/x2
plus the pre-registered monotonicity swap (scored as a Gaussian
likelihood of the real Spearman monotonicity under the per-seed
synthetic distribution). This substitution is documented, not silent.

Pass rule (operationalizing "region moves by less than its own width"):
  |Delta ln sigma_r(best)| < full 68% width in ln sigma_r  AND
  |Delta ln UL95(sigma_i)| < ln(UL95/UL68) width of the fiducial
  sigma_i marginal AND 95%-region overlap fraction >= 0.5.

Usage:
  .venv/bin/python analysis/robustness.py --run teff_gk --workers 4
  .venv/bin/python analysis/robustness.py --run all --workers 4
  .venv/bin/python analysis/robustness.py --analyze
"""
import argparse
import dataclasses
import json
import multiprocessing as mp
import time
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from grid_inference import (  # noqa: E402
    KCATS, LVL_1SIG, LVL_2SIG, kde_loglike, mult_loglike, poly_design,
    silverman_bw, upper_quantiles,
)
from grid_sweep import FINE_SR, FINE_SI, cell_seed  # noqa: E402
from exoverse.architecture import Architecture
from exoverse.dichotomy import observables_from_koi, observables_from_universe
from exoverse.kepler_data import (
    FIDUCIAL, Cuts, load_koi_systems, load_stellar_targets,
)
from exoverse.kepler_field import Detection, simulate_universe

ROOT = Path(__file__).resolve().parent.parent
GRID = ROOT / "results" / "grid"
ROB = ROOT / "results" / "robustness"
VALLEY_D2 = 100.0
M = 8

VARIANTS = {
    # V1 robovetter score
    "score0":     {"kind": "real", "cuts": {"koi_score_min": 0.0}},
    "score9":     {"kind": "real", "cuts": {"koi_score_min": 0.9}},
    # V2 stellar cuts
    "teff_gk":    {"kind": "resim", "cuts": {"teff_min": 4700.0,
                                             "teff_max": 6500.0}},
    "logg42":     {"kind": "resim", "cuts": {"logg_min": 4.2}},
    "nodataspan": {"kind": "resim", "cuts": {"dataspan_min": 0.0}},
    # V3 completeness
    "snr65":      {"kind": "resim", "det": {"snr_threshold": 6.5}},
    "snr80":      {"kind": "resim", "det": {"snr_threshold": 8.0}},
    "mesramp":    {"kind": "resim", "det": {"mes_ramp_width": 1.0}},
    # V4 window function
    "binwindow":  {"kind": "resim", "det": {"window": "binomial"}},
    # V5 radius uncertainty (real-side; needs data/dr25_koi_prad_err.csv)
    "radius_err": {"kind": "real", "metric": "radius_perturb"},
    # V6 multiplicity truncation
    "trunc4":     {"kind": "metric", "metric": "trunc4"},
    # V7 size-observable sensitivity
    "bw_half":    {"kind": "metric", "metric": "bw0.5"},
    "bw_double":  {"kind": "metric", "metric": "bw2.0"},
    "monotonicity": {"kind": "metric", "metric": "monotonicity"},
    # V8 occurrence normalization
    "poisson":    {"kind": "metric", "metric": "poisson"},
}

_TARGETS = None
_DET = None


def variant_cuts(cfg: dict) -> Cuts:
    return dataclasses.replace(FIDUCIAL, **cfg.get("cuts", {}))


def variant_det(cfg: dict) -> Detection:
    return Detection(**cfg.get("det", {}))


def fiducial_valley() -> list:
    """(sr, si) cells with Delta 2lnL < VALLEY_D2 in the fiducial fit."""
    inf = json.loads((ROOT / "results" / "grid_inference.json").read_text())
    free = {k: v["total"] for k, v in inf["logL"].items()
            if not k.startswith("None|")}
    mx = max(free.values())
    out = []
    for k, v in free.items():
        if 2 * (mx - v) < VALLEY_D2:
            sr_s, si_s = k.split("|")
            out.append((float(sr_s), float(si_s)))
    return sorted(out)


# ------------------------------------------------------------ simulation
def _init_worker(cuts_cfg: dict, det_cfg: dict) -> None:
    global _TARGETS, _DET
    _TARGETS = load_stellar_targets(
        "data", dataclasses.replace(FIDUCIAL, **cuts_cfg))
    _DET = Detection(**det_cfg)


def _run_cell(job: tuple) -> str:
    sr, si, m, seed, out_dir = job
    obs = observables_from_universe(simulate_universe(
        _TARGETS, seed, Architecture(sigma_r=sr, sigma_i=si), det=_DET))
    tag = f"sr{sr}_si{si}_m{m}"
    (Path(out_dir) / f"cell_{tag}.json").write_text(json.dumps({
        "sigma_r": sr, "sigma_i": si, "m": m, "seed": seed,
        "n_k": {str(k): v for k, v in obs.n_k.items()},
        "dlogr": [round(x, 6) for x in obs.dlogr],
        "monotonicity": obs.monotonicity,
        "n_systems": obs.n_systems, "n_planets": obs.n_planets,
        "n_targets": obs.n_targets}))
    return tag


def run_variant(name: str, workers: int) -> None:
    cfg = VARIANTS[name]
    assert cfg["kind"] == "resim", f"{name} needs no simulation"
    out_dir = ROB / name
    out_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for sr, si in fiducial_valley():
        i_sr, i_si = FINE_SR.index(sr), FINE_SI.index(si)
        for m in range(M):
            if (out_dir / f"cell_sr{sr}_si{si}_m{m}.json").exists():
                continue
            jobs.append((sr, si, m, cell_seed(i_sr, i_si, m), str(out_dir)))
    print(f"[{name}] {len(jobs)} universes, {workers} workers", flush=True)
    t0 = time.time()
    with mp.get_context("spawn").Pool(
            workers, initializer=_init_worker,
            initargs=(cfg.get("cuts", {}), cfg.get("det", {}))) as pool:
        for done, tag in enumerate(
                pool.imap_unordered(_run_cell, jobs, chunksize=1), 1):
            rate = done / (time.time() - t0)
            print(f"  [{name}] {done}/{len(jobs)} {tag} "
                  f"(eta {(len(jobs)-done)/rate/60:.0f} min)", flush=True)
    print(f"[{name}] complete in {(time.time()-t0)/60:.1f} min", flush=True)


# ------------------------------------------------------------- analysis
def load_cells(folder: Path, sr, si, ms=range(M)) -> list:
    return [json.loads((folder / f"cell_sr{sr}_si{si}_m{m}.json").read_text())
            for m in ms]


def perturbed_real_dlogr(n_draws: int = 25, seed: int = 20260704):
    """Real |dlogR| samples with koi_prad perturbed by its catalog errors
    (symmetrized Gaussian). Requires data/dr25_koi_prad_err.csv with
    kepoi_name, koi_prad_err1, koi_prad_err2."""
    import csv
    err_file = ROOT / "data" / "dr25_koi_prad_err.csv"
    if not err_file.exists():
        return None
    errs = {}
    with err_file.open() as f:
        for row in csv.DictReader(f):
            try:
                e1 = abs(float(row["koi_prad_err1"]))
                e2 = abs(float(row["koi_prad_err2"]))
                errs[row["kepoi_name"].strip()] = 0.5 * (e1 + e2)
            except (TypeError, ValueError):
                continue
    systems = load_koi_systems("data")
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_draws):
        dlogr = []
        for kois in systems.values():
            if len(kois) < 2:
                continue
            rads = [max(0.05, k.prad + rng.normal(0.0, errs.get(
                k.kepoi_name, 0.0))) for k in kois]
            for a, b in zip(rads[:-1], rads[1:]):
                dlogr.append(abs(np.log10(b) - np.log10(a)))
        draws.append(np.array(dlogr))
    return draws


def score_cells(cells, real, bw, metric: str, perturbed=None) -> float:
    counts = np.zeros(len(KCATS))
    dlogr, monos = [], []
    for d in cells:
        for i, k in enumerate(KCATS):
            counts[i] += d["n_k"].get(k, 0)
        dlogr.extend(d["dlogr"])
        if d["monotonicity"] is not None:
            monos.append(d["monotonicity"])
    syn = np.array(dlogr)
    n_obs = np.array([real.n_k.get("6+" if k == "6+" else int(k), 0)
                      for k in KCATS], dtype=float)

    if metric == "trunc4":
        c4 = np.array([counts[0], counts[1], counts[2], counts[3:].sum()])
        n4 = np.array([n_obs[0], n_obs[1], n_obs[2], n_obs[3:].sum()])
        p = (c4 + 0.5) / (c4.sum() + 2.0)
        lm = float((n4 * np.log(p)).sum())
        return lm + kde_loglike(np.array(real.dlogr), syn, bw)
    if metric == "poisson":
        lam = np.maximum(counts / len(cells), 0.5)
        lm = float((n_obs * np.log(lam) - lam).sum())
        return lm + kde_loglike(np.array(real.dlogr), syn, bw)
    if metric in ("bw0.5", "bw2.0"):
        f = 0.5 if metric == "bw0.5" else 2.0
        return mult_loglike(real.n_k, counts) + kde_loglike(
            np.array(real.dlogr), syn, bw * f)
    if metric == "monotonicity":
        mu, sd = float(np.mean(monos)), float(np.std(monos, ddof=1))
        sd = max(sd, 1e-3)
        ls = -0.5 * ((real.monotonicity - mu) / sd) ** 2 - np.log(sd)
        return mult_loglike(real.n_k, counts) + float(ls)
    if metric == "radius_perturb":
        ls = float(np.mean([kde_loglike(rd, syn, bw) for rd in perturbed]))
        return mult_loglike(real.n_k, counts) + ls
    # fiducial scoring
    return mult_loglike(real.n_k, counts) + kde_loglike(
        np.array(real.dlogr), syn, bw)


def fit_region(points: list) -> dict:
    """Cubic smooth over (ln sr, ln si, lnL) points -> best fit, 68/95
    region masks on a fixed mesh, and sigma_i/sigma_r marginals."""
    vx = np.log([p[0] for p in points])
    vy = np.log([p[1] for p in points])
    vz = np.array([p[2] for p in points])
    beta, *_ = np.linalg.lstsq(poly_design(vx, vy), vz, rcond=None)
    gx = np.linspace(vx.min(), vx.max(), 160)
    gy = np.linspace(vy.min(), vy.max(), 160)
    GX, GY = np.meshgrid(gx, gy, indexing="ij")
    s = (poly_design(GX.ravel(), GY.ravel()) @ beta).reshape(GX.shape)
    i_b, j_b = np.unravel_index(np.argmax(s), s.shape)
    d2 = 2 * (s.max() - s)
    w = np.exp(s - s.max())
    post_si = w.sum(axis=0)
    post_sr = w.sum(axis=1)
    resid = vz - poly_design(vx, vy) @ beta
    return {"best_sr": float(np.exp(gx[i_b])),
            "best_si": float(np.exp(gy[j_b])),
            "mask95": d2 < LVL_2SIG, "mask68": d2 < LVL_1SIG,
            "si": upper_quantiles(np.exp(gy), post_si),
            "sr": upper_quantiles(np.exp(gx), post_sr),
            "edge": bool(i_b in (0, 159) or j_b in (0, 159)),
            "resid_rms": float(np.std(resid))}


def analyze() -> None:
    valley = fiducial_valley()
    targets = load_stellar_targets("data")
    real_fid = observables_from_koi(load_koi_systems("data"),
                                    n_targets=len(targets))
    bw = silverman_bw(np.array(real_fid.dlogr))
    perturbed = perturbed_real_dlogr()

    def region_for(name: str) -> dict | None:
        cfg = VARIANTS.get(name, {"kind": "fiducial"})
        kind = cfg.get("kind", "fiducial")
        metric = cfg.get("metric", "fiducial")
        folder = ROB / name if kind == "resim" else GRID
        if kind == "resim" and not folder.exists():
            return None
        if metric == "radius_perturb" and perturbed is None:
            return None
        cuts = variant_cuts(cfg) if kind in ("real", "resim") else FIDUCIAL
        tg = (load_stellar_targets("data", cuts)
              if kind in ("real", "resim") else targets)
        real = (observables_from_koi(load_koi_systems("data", cuts),
                                     n_targets=len(tg))
                if kind in ("real", "resim") else real_fid)
        pts = []
        for sr, si in valley:
            try:
                cells = load_cells(folder, sr, si)
            except FileNotFoundError:
                return None
            pts.append((sr, si, score_cells(cells, real, bw,
                                            metric, perturbed)))
        r = fit_region(pts)
        r["n_real_systems"] = real.n_systems
        return r

    fid = region_for("fiducial")
    rows, results = [], {"fiducial": {k: v for k, v in fid.items()
                                      if not k.startswith("mask")}}
    w68_sr = np.log(fid["sr"]["q68"] / fid["sr"]["q32"])
    w_si = np.log(fid["si"]["q95"] / fid["si"]["q68"])
    for name in VARIANTS:
        r = region_for(name)
        if r is None:
            rows.append((name, "SKIPPED (no data)", "", "", "", ""))
            results[name] = "skipped"
            continue
        dsr = np.log(r["best_sr"] / fid["best_sr"])
        dsi = np.log(r["si"]["q95"] / fid["si"]["q95"])
        overlap = float((r["mask95"] & fid["mask95"]).sum()
                        / max(1, fid["mask95"].sum()))
        ok = (abs(dsr) < w68_sr and abs(dsi) < w_si and overlap >= 0.5)
        rows.append((name, f"{r['best_sr']:.2f}", f"{dsr:+.2f}",
                     f"{r['si']['q95']:.2f} ({dsi:+.2f})",
                     f"{overlap:.2f}", "pass" if ok else "FAIL"))
        results[name] = {k: v for k, v in r.items()
                         if not k.startswith("mask")}
        results[name].update({"d_ln_best_sr": float(dsr),
                              "d_ln_ul95_si": float(dsi),
                              "overlap95": overlap, "pass": bool(ok)})

    (ROOT / "results" / "robustness.json").write_text(json.dumps({
        "valley_cells": len(valley), "pass_widths":
            {"ln_sr_68_width": float(w68_sr), "ln_si_ul_width": float(w_si)},
        "fiducial": results["fiducial"],
        "variants": {k: v for k, v in results.items() if k != "fiducial"},
    }, indent=1))

    print(f"fiducial (m=0..7, {len(valley)} valley cells): "
          f"best sr={fid['best_sr']:.2f}, si UL95={fid['si']['q95']:.2f}, "
          f"sr 68% [{fid['sr']['q32']:.2f},{fid['sr']['q68']:.2f}]")
    print(f"pass widths: |d ln sr|<{w68_sr:.2f}, "
          f"|d ln UL95(si)|<{w_si:.2f}, overlap95>=0.5\n")
    hdr = ("variant", "best_sr", "d_ln_sr", "si_UL95 (d_ln)",
           "overlap95", "verdict")
    wid = [13, 8, 8, 16, 10, 8]
    print(" | ".join(h.ljust(w) for h, w in zip(hdr, wid)))
    print("-" * 74)
    for row in rows:
        print(" | ".join(str(c).ljust(w) for c, w in zip(row, wid)))
    print("\nwritten: results/robustness.json")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", help="resim variant name or 'all'")
    ap.add_argument("--analyze", action="store_true")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    if a.run:
        names = ([n for n, c in VARIANTS.items() if c["kind"] == "resim"]
                 if a.run == "all" else [a.run])
        for n in names:
            run_variant(n, a.workers)
    if a.analyze:
        analyze()


if __name__ == "__main__":
    main()
