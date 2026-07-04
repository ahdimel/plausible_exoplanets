"""Two-population inclination mixture test (the secondary question):
is an isotropic-ish 'hot' fraction still statistically demanded once
radius correlation is modeled jointly?

Design: sigma_r is FIXED at the joint best fit (1.5244). Justified by
the fine-grid decomposition (results/grid_inference.json): sigma_r is
constrained almost entirely by the |dlogR| statistic and is orthogonal
to the inclination structure, so the mixture axes are swept at fixed
radius correlation. Sweep:
  sigma_i (cold)  in {0.2, 0.33, 0.53, 0.87, 1.43} deg
  f_hot           in {0.05, 0.1, 0.15, 0.2, 0.3, 0.4}
  sigma_i_hot     in {10, 30} deg
M=8 seeds per cell -> 480 universes. f_hot=0 baselines are the existing
fine-grid cells at sr=1.5244. Observables cached per cell in
results/grid_mixture/ in the grid_sweep cell format.

Seeds: SeedSequence([BASE_SEED, 7, i_si, i_f, i_h, m]) — the leading 7
namespaces the mixture away from the fine-grid stream.

Run: .venv/bin/python analysis/mixture_sweep.py --workers 4
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
SIGMA_R = 1.5244
SI_COLD = [0.2, 0.3271, 0.5349, 0.8747, 1.4304]   # fine-grid values
F_HOT = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4]
SI_HOT = [10.0, 30.0]
M = 8
OUT = Path(__file__).resolve().parent.parent / "results" / "grid_mixture"

_TARGETS = None


def _init_worker() -> None:
    global _TARGETS
    _TARGETS = load_stellar_targets("data")


def _run_cell(job: tuple) -> str:
    si, f, sih, m, seed, tag = job
    arch = Architecture(sigma_r=SIGMA_R, sigma_i=si, f_hot=f,
                        sigma_i_hot=sih)
    obs = observables_from_universe(simulate_universe(_TARGETS, seed, arch))
    (OUT / f"cell_{tag}.json").write_text(json.dumps({
        "sigma_r": SIGMA_R, "sigma_i": si, "f_hot": f, "sigma_i_hot": sih,
        "m": m, "seed": seed,
        "n_k": {str(k): v for k, v in obs.n_k.items()},
        "dlogr": [round(x, 6) for x in obs.dlogr],
        "monotonicity": obs.monotonicity,
        "n_systems": obs.n_systems,
        "n_planets": obs.n_planets,
        "n_targets": obs.n_targets}))
    return tag


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "axes.json").write_text(json.dumps(
        {"sigma_r": SIGMA_R, "sigma_i": SI_COLD, "f_hot": F_HOT,
         "sigma_i_hot": SI_HOT, "n_seeds": M, "base_seed": BASE_SEED},
        indent=2))
    jobs = []
    for i_si, si in enumerate(SI_COLD):
        for i_f, f in enumerate(F_HOT):
            for i_h, sih in enumerate(SI_HOT):
                for m in range(M):
                    tag = f"si{si}_f{f}_sih{sih}_m{m}"
                    if (OUT / f"cell_{tag}.json").exists():
                        continue
                    seed = int(np.random.SeedSequence(
                        [BASE_SEED, 7, i_si, i_f, i_h, m]
                    ).generate_state(1)[0])
                    jobs.append((si, f, sih, m, seed, tag))
    print(f"{len(jobs)} mixture universes, {a.workers} workers", flush=True)
    t0 = time.time()
    with mp.get_context("spawn").Pool(a.workers,
                                      initializer=_init_worker) as pool:
        for done, tag in enumerate(
                pool.imap_unordered(_run_cell, jobs, chunksize=1), 1):
            rate = done / (time.time() - t0)
            print(f"  {done}/{len(jobs)} {tag} "
                  f"(eta {(len(jobs)-done)/rate/60:.0f} min)", flush=True)
    print(f"mixture sweep complete in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
