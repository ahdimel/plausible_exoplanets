"""Population generation: build N validated systems and persist them."""

from __future__ import annotations

import sys
import time
from typing import Optional

import numpy as np

from .database import WorldDB
from .observatories import observe
from .system import generate_system
from .transits import compute_geometry


def generate_population(db_path: str, n_systems: int, seed: int = 42,
                        progress: bool = True) -> dict:
    """Generate n_systems validated systems and store them in db_path."""
    db = WorldDB(db_path)
    master = np.random.SeedSequence(seed)
    child_seeds = master.generate_state(n_systems)

    t0 = time.time()
    for i, child in enumerate(child_seeds):
        name = f"PXS-{seed}-{i:05d}"
        system = generate_system(int(child), name)
        geoms = [compute_geometry(system.star, p) for p in system.planets]
        obs = [observe(system.star, p, g) if g.transits else []
               for p, g in zip(system.planets, geoms)]
        db.save_system(system, geoms, obs)
        if progress and (i + 1) % 100 == 0:
            rate = (i + 1) / (time.time() - t0)
            print(f"  {i+1}/{n_systems} systems ({rate:.0f}/s)", file=sys.stderr)

    stats = db.stats()
    db.close()
    return stats
