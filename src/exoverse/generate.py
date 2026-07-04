"""Population generation: build N validated systems and persist them."""

from __future__ import annotations

import sys
import time

import numpy as np

from .architecture import Architecture
from .atmospheres import score_atmosphere_observability
from .database import WorldDB
from .observatories import JWST_INSTRUMENTS, observe
from .system import generate_system
from .transits import compute_geometry


def generate_population(db_path: str, n_systems: int, seed: int = 42,
                        progress: bool = True,
                        dmax_pc: float = 300.0,
                        arch: Architecture | None = None) -> dict:
    """Generate n_systems validated systems and store them in db_path.

    dmax_pc caps host distances (solar-neighborhood mode when small) and
    arch sets the architecture knobs (sigma_r/sigma_i, architecture.py);
    both are recorded in DB meta so re-generation from (seed, name) stays
    exact."""
    db = WorldDB(db_path)
    db.set_meta("dmax_pc", str(dmax_pc))
    if arch is not None:
        for key, value in arch.meta_items().items():
            db.set_meta(key, value)
    master = np.random.SeedSequence(seed)
    child_seeds = master.generate_state(n_systems)

    t0 = time.time()
    for i, child in enumerate(child_seeds):
        name = f"PXS-{seed}-{i:05d}"
        system = generate_system(int(child), name, dmax_pc=dmax_pc, arch=arch)
        star = system.star
        geoms = [compute_geometry(star, p) for p in system.planets]
        obs = [observe(star, p, g, system.noise)
               for p, g in zip(system.planets, geoms)]

        # Atmosphere spectroscopy scoring (transiting planets only)
        jwst = [(nm, fn(star.mag_j), star.mag_j > sat)
                for nm, fn, sat in JWST_INSTRUMENTS]
        atm_obs = []
        for p, g in zip(system.planets, geoms):
            if g.transits and p.atmosphere is not None:
                atm_obs.append(score_atmosphere_observability(
                    p.atmosphere, g.t14_hours, system.noise, jwst))
            else:
                atm_obs.append([])

        db.save_system(system, geoms, obs, atm_obs)
        if progress and (i + 1) % 100 == 0:
            rate = (i + 1) / (time.time() - t0)
            print(f"  {i+1}/{n_systems} systems ({rate:.0f}/s)", file=sys.stderr)

    stats = db.stats()
    db.close()
    return stats
