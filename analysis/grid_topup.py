"""Top up Monte-Carlo seeds (m=8..15) for fine-grid cells in the
likelihood valley, where the Wilks 1-sigma level (Delta 2lnL = 2.3) is
below the M=8 seed-noise floor. Cells are selected from
results/grid_inference.json: sigma_r>0 cells within DELTA_MAX of the
maximum, plus the sigma_r=None row within DELTA_MAX of its own row
maximum (it is the fixed-uncorrelated comparison model).

Seeds stay on the same deterministic stream as grid_sweep.py:
cell_seed(i_sr, i_si, m) with the fine-axis indices, m in 8..15.

Run: .venv/bin/python analysis/grid_topup.py --workers 4
"""
import argparse
import json
import multiprocessing as mp
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grid_sweep import (  # noqa: E402
    FINE_SR, FINE_SI, OUT, _init_worker, _run_cell, cell_seed,
)

DELTA_MAX = 200.0
M_RANGE = range(8, 16)


def valley_cells() -> list:
    inf = json.loads(
        (OUT.parent / "grid_inference.json").read_text())
    ll = {k: v["total"] for k, v in inf["logL"].items()}
    free = {k: v for k, v in ll.items() if not k.startswith("None|")}
    none = {k: v for k, v in ll.items() if k.startswith("None|")}
    keep = [k for k, v in free.items()
            if 2 * (max(free.values()) - v) < DELTA_MAX]
    keep += [k for k, v in none.items()
             if 2 * (max(none.values()) - v) < DELTA_MAX]
    return keep


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()

    jobs = []
    for key in valley_cells():
        sr_s, si_s = key.split("|")
        sr = None if sr_s == "None" else float(sr_s)
        si = float(si_s)
        i_sr = FINE_SR.index(sr)
        i_si = FINE_SI.index(si)
        for m in M_RANGE:
            tag = f"sr{sr_s}_si{si}_m{m}"
            if not (OUT / f"cell_{tag}.json").exists():
                jobs.append((sr, si, m, cell_seed(i_sr, i_si, m), tag))
    print(f"{len(jobs)} top-up universes, {a.workers} workers", flush=True)
    t0 = time.time()
    with mp.get_context("spawn").Pool(a.workers,
                                      initializer=_init_worker) as pool:
        for done, tag in enumerate(
                pool.imap_unordered(_run_cell, jobs, chunksize=1), 1):
            rate = done / (time.time() - t0)
            print(f"  {done}/{len(jobs)} {tag} "
                  f"(eta {(len(jobs)-done)/rate/60:.0f} min)", flush=True)
    print(f"top-up complete in {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
