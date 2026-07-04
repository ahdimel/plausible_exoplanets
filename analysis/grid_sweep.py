"""(sigma_r, sigma_i) grid sweep on the DR25-conditioned Kepler field.

Caches per-(cell, seed) OBSERVABLES (not distances) as JSON in
results/grid/, so the distance metric can be recomputed under any Phase 4
variant without re-simulating. Resumable: existing cache files are
skipped. Deterministic: universe seed = f(base_seed, i_sr, i_si, m).

Usage:
  .venv/bin/python analysis/grid_sweep.py                        # coarse 8x8, M=3
  .venv/bin/python analysis/grid_sweep.py --fine --workers 4     # 14x16, M=8
"""
import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from exoverse.architecture import Architecture
from exoverse.dichotomy import observables_from_universe
from exoverse.kepler_data import load_stellar_targets
from exoverse.kepler_field import simulate_universe

BASE_SEED = 20260704
OUT = Path(__file__).resolve().parent.parent / "results" / "grid"

# sigma_r axis: index 0 = None (correlation off, the exact baseline null),
# then log-spaced. sigma_i axis: log-spaced degrees.
COARSE_SR = [None, 0.05, 0.093, 0.17, 0.32, 0.59, 1.08, 2.0]
COARSE_SI = [0.5, 0.87, 1.5, 2.6, 4.5, 7.9, 13.7, 24.0]
# Fine axes fixed at checkpoint 3 (user-approved): the coarse surface
# excludes sigma_r < 0.2 and sigma_i > 8 deg firmly; the low-sigma_i edge
# must extend to 0.2 deg because the coarse minimum touched 0.5 deg.
FINE_SR = [None] + list(np.round(np.geomspace(0.2, 3.0, 13), 4))
FINE_SI = list(np.round(np.geomspace(0.2, 8.0, 16), 4))

_TARGETS = None


def _init_worker() -> None:
    global _TARGETS
    _TARGETS = load_stellar_targets("data")


def _run_cell(job: tuple) -> str:
    sr, si, m, seed, tag = job
    res = simulate_universe(_TARGETS, seed, Architecture(sigma_r=sr, sigma_i=si))
    obs = observables_from_universe(res)
    (OUT / f"cell_{tag}.json").write_text(json.dumps({
        "sigma_r": sr, "sigma_i": si, "m": m, "seed": seed,
        "n_k": {str(k): v for k, v in obs.n_k.items()},
        "dlogr": [round(x, 6) for x in obs.dlogr],
        "monotonicity": obs.monotonicity,
        "n_systems": obs.n_systems,
        "n_planets": obs.n_planets,
        "n_targets": obs.n_targets}))
    return tag


def cell_seed(i_sr: int, i_si: int, m: int) -> int:
    return int(np.random.SeedSequence(
        [BASE_SEED, i_sr, i_si, m]).generate_state(1)[0])


def run(sr_axis, si_axis, n_seeds: int, workers: int = 1) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "axes.json").write_text(json.dumps(
        {"sigma_r": ["None" if s is None else s for s in sr_axis],
         "sigma_i": si_axis, "n_seeds": n_seeds, "base_seed": BASE_SEED},
        indent=2))
    jobs = []
    for i_sr, sr in enumerate(sr_axis):
        for i_si, si in enumerate(si_axis):
            for m in range(n_seeds):
                tag = f"sr{'None' if sr is None else sr}_si{si}_m{m}"
                if not (OUT / f"cell_{tag}.json").exists():
                    jobs.append((sr, si, m, cell_seed(i_sr, i_si, m), tag))
    total = len(sr_axis) * len(si_axis) * n_seeds
    print(f"{total} cells, {len(jobs)} to run, {workers} workers",
          flush=True)
    t0 = time.time()
    if workers <= 1:
        _init_worker()
        for done, job in enumerate(jobs, 1):
            _run_cell(job)
            rate = done / (time.time() - t0)
            print(f"  {done}/{len(jobs)} {job[4]} "
                  f"(eta {(len(jobs)-done)/rate/60:.0f} min)", flush=True)
    else:
        with mp.get_context("spawn").Pool(workers,
                                          initializer=_init_worker) as pool:
            for done, tag in enumerate(
                    pool.imap_unordered(_run_cell, jobs, chunksize=1), 1):
                rate = done / (time.time() - t0)
                print(f"  {done}/{len(jobs)} {tag} "
                      f"(eta {(len(jobs)-done)/rate/60:.0f} min)",
                      flush=True)
    print(f"grid complete: {total} cells in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fine", action="store_true")
    ap.add_argument("--workers", type=int, default=1)
    a = ap.parse_args()
    if a.fine:
        run(FINE_SR, FINE_SI, 8, workers=a.workers)
    else:
        run(COARSE_SR, COARSE_SI, 3, workers=a.workers)
