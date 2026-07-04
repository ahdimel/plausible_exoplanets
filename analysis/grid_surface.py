"""Assemble the (sigma_r, sigma_i) distance surface from cached grid
observables (results/grid/) against the real DR25 observables. Writes
results/grid_surface.json and prints the raw table. Distances recomputed
here so Phase 4 metric variants never require re-simulation.

Run: .venv/bin/python analysis/grid_surface.py [--mult-mode multinomial]
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from exoverse.dichotomy import (
    Observables, combined_distance, observables_from_koi,
)
from exoverse.kepler_data import load_koi_systems, load_stellar_targets

GRID = Path(__file__).resolve().parent.parent / "results" / "grid"


def obs_from_cell(d: dict) -> Observables:
    n_k = {(int(k) if k != "6+" else "6+"): v for k, v in d["n_k"].items()}
    return Observables(n_k=n_k, dlogr=d["dlogr"],
                       monotonicity=d["monotonicity"],
                       n_systems=d["n_systems"], n_planets=d["n_planets"],
                       n_targets=d["n_targets"])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mult-mode", default="multinomial")
    ap.add_argument("--size-mode", default="ks")
    ap.add_argument("--w-mult", type=float, default=1.0)
    ap.add_argument("--w-size", type=float, default=1.0)
    a = ap.parse_args()

    axes = json.loads((GRID / "axes.json").read_text())
    sr_axis = [None if s == "None" else float(s) for s in axes["sigma_r"]]
    si_axis = [float(s) for s in axes["sigma_i"]]
    n_seeds = axes["n_seeds"]

    targets = load_stellar_targets("data")
    real = observables_from_koi(load_koi_systems("data"),
                                n_targets=len(targets))

    surface, cells = {}, {}
    for sr in sr_axis:
        for si in si_axis:
            comps, spms = [], []
            for m in range(n_seeds):
                tag = f"sr{'None' if sr is None else sr}_si{si}_m{m}"
                d = json.loads((GRID / f"cell_{tag}.json").read_text())
                obs = obs_from_cell(d)
                comps.append(combined_distance(
                    obs, real, w_mult=a.w_mult, w_size=a.w_size,
                    mult_mode=a.mult_mode, size_mode=a.size_mode))
                multis = sum(v for k, v in obs.n_k.items() if k != 1)
                spms.append(obs.n_k.get(1, 0) / multis if multis else np.inf)
            key = f"{'None' if sr is None else sr}|{si}"
            totals = [c["total"] for c in comps]
            cells[key] = {
                "sigma_r": sr, "sigma_i": si,
                "total_mean": float(np.mean(totals)),
                "total_std": float(np.std(totals)),
                "mult_mean": float(np.mean([c["multiplicity"] for c in comps])),
                "size_mean": float(np.mean([c["size"] for c in comps])),
                "spm_mean": float(np.mean(spms)),
            }
            surface[key] = cells[key]["total_mean"]

    out = {"axes": axes, "metric": vars(a),
           "real_n_k": {str(k): v for k, v in real.n_k.items()},
           "cells": cells}
    (GRID.parent / "grid_surface.json").write_text(json.dumps(out, indent=1))

    # Raw table: rows sigma_r, cols sigma_i, entries total distance (mean)
    best = min(cells.values(), key=lambda c: c["total_mean"])
    mc_noise = float(np.median([c["total_std"] for c in cells.values()]))
    print("total distance (mean over seeds); rows sigma_r, cols sigma_i")
    print("sig_r\\i " + " ".join(f"{si:>7.2f}" for si in si_axis))
    for sr in sr_axis:
        row = [cells[f"{'None' if sr is None else sr}|{si}"]["total_mean"]
               for si in si_axis]
        lab = " None" if sr is None else f"{sr:5.3f}"
        print(f"{lab:>7} " + " ".join(f"{v:7.4f}" for v in row))
    print(f"\nbest cell: sigma_r={best['sigma_r']} sigma_i={best['sigma_i']}"
          f"  total={best['total_mean']:.4f} (mult {best['mult_mean']:.4f}"
          f" + size {best['size_mean']:.4f})  spm={best['spm_mean']:.2f}")
    print(f"median MC noise per cell (std over {n_seeds} seeds): "
          f"{mc_noise:.4f}")
    print(f"real spm: 3.399   surface written: results/grid_surface.json")


if __name__ == "__main__":
    main()
