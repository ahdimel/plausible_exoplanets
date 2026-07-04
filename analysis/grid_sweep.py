"""(sigma_r, sigma_i) grid sweep on the DR25-conditioned Kepler field.

Caches per-(cell, seed) OBSERVABLES (not distances) as JSON in
results/grid/, so the distance metric can be recomputed under any Phase 4
variant without re-simulating. Resumable: existing cache files are
skipped. Deterministic: universe seed = f(base_seed, i_sr, i_si, m).

Usage:
  .venv/bin/python analysis/grid_sweep.py            # coarse 8x8, M=3
  .venv/bin/python analysis/grid_sweep.py --fine     # 20x20, M=5
"""
import argparse
import json
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
FINE_SR = [None] + list(np.round(np.geomspace(0.03, 3.0, 19), 4))
FINE_SI = list(np.round(np.geomspace(0.3, 30.0, 20), 4))


def cell_seed(i_sr: int, i_si: int, m: int) -> int:
    return int(np.random.SeedSequence(
        [BASE_SEED, i_sr, i_si, m]).generate_state(1)[0])


def run(sr_axis, si_axis, n_seeds: int) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    targets = load_stellar_targets("data")
    (OUT / "axes.json").write_text(json.dumps(
        {"sigma_r": ["None" if s is None else s for s in sr_axis],
         "sigma_i": si_axis, "n_seeds": n_seeds, "base_seed": BASE_SEED,
         "n_targets": len(targets)}, indent=2))
    total = len(sr_axis) * len(si_axis) * n_seeds
    done = 0
    t0 = time.time()
    for i_sr, sr in enumerate(sr_axis):
        for i_si, si in enumerate(si_axis):
            for m in range(n_seeds):
                done += 1
                tag = f"sr{'None' if sr is None else sr}_si{si}_m{m}"
                path = OUT / f"cell_{tag}.json"
                if path.exists():
                    continue
                arch = Architecture(sigma_r=sr, sigma_i=si)
                res = simulate_universe(targets, cell_seed(i_sr, i_si, m),
                                        arch)
                obs = observables_from_universe(res)
                path.write_text(json.dumps({
                    "sigma_r": sr, "sigma_i": si, "m": m,
                    "seed": cell_seed(i_sr, i_si, m),
                    "n_k": {str(k): v for k, v in obs.n_k.items()},
                    "dlogr": [round(x, 6) for x in obs.dlogr],
                    "monotonicity": obs.monotonicity,
                    "n_systems": obs.n_systems,
                    "n_planets": obs.n_planets,
                    "n_targets": obs.n_targets}))
                rate = done / (time.time() - t0)
                print(f"  {done}/{total} {tag} "
                      f"(eta {(total-done)/rate/60:.0f} min)", flush=True)
    print(f"grid complete: {total} cells in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fine", action="store_true")
    a = ap.parse_args()
    if a.fine:
        run(FINE_SR, FINE_SI, 5)
    else:
        run(COARSE_SR, COARSE_SI, 3)
