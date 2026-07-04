"""Phase 2 integration check: one default-architecture universe on the real
DR25 target sample vs the real KOI observables, end to end.

This is NOT the grid — a single smoke universe proving the conditioned
pipeline runs and the distance metric returns finite, sensible components.
Run: .venv/bin/python analysis/phase2_fiducial_check.py
"""
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from exoverse.architecture import Architecture
from exoverse.dichotomy import (
    combined_distance, observables_from_koi, observables_from_universe,
)
from exoverse.kepler_data import load_koi_systems, load_stellar_targets
from exoverse.kepler_field import simulate_universe

SEED = 20260703


def spm(n_k: dict) -> float:
    multis = sum(v for k, v in n_k.items() if k != 1)
    return n_k.get(1, 0) / multis if multis else float("inf")


def main() -> None:
    targets = load_stellar_targets("data")
    koi = load_koi_systems("data")
    real = observables_from_koi(koi, n_targets=len(targets))
    print(f"targets after cuts : {len(targets)}")
    print(f"real: N_k={real.n_k}  spm={spm(real.n_k):.3f}  "
          f"pairs={len(real.dlogr)}  monotonicity={real.monotonicity}")

    for label, arch in (
        ("default (sigma_i=1.5, no sigma_r)", None),
        ("cold+peas (sigma_i=1.5, sigma_r=0.3)", Architecture(sigma_r=0.3)),
        ("hot (sigma_i=8)", Architecture(sigma_i=8.0)),
    ):
        t0 = time.perf_counter()
        syn = observables_from_universe(simulate_universe(targets, SEED, arch))
        dt = time.perf_counter() - t0
        d = combined_distance(syn, real)
        print(f"\n[{label}]  ({dt:.1f} s)")
        print(f"  syn: N_k={syn.n_k}  spm={spm(syn.n_k):.3f}  "
              f"pairs={len(syn.dlogr)}  monotonicity={syn.monotonicity}")
        print(f"  distance: {json.dumps({k: round(v, 5) for k, v in d.items() if isinstance(v, float)})}")


if __name__ == "__main__":
    main()
