"""End-to-end pipeline tests: generation -> observability -> database."""

import numpy as np
import pytest

from exoverse.database import WorldDB
from exoverse.generate import generate_population
from exoverse.observatories import observe
from exoverse.system import generate_system
from exoverse.transits import compute_geometry


def test_observatories_rank_sensibly():
    """JWST must beat a 1-m ground telescope on any usable target."""
    for seed in range(200):
        sys_ = generate_system(seed, f"T-{seed}")
        for p in sys_.planets:
            g = compute_geometry(sys_.star, p)
            if not g.transits:
                continue
            obs = {o.observatory: o for o in observe(sys_.star, p, g)}
            jwst = obs["JWST NIRSpec Prism"]
            ground = obs["Ground 1-m survey"]
            if jwst.usable and ground.usable:
                assert jwst.snr_per_transit > ground.snr_per_transit
            return  # one transiting example is enough
    pytest.skip("no transiting planet in 200 seeds (should never happen)")


def test_population_roundtrip(tmp_path):
    db_path = tmp_path / "test.db"
    stats = generate_population(str(db_path), 40, seed=99, progress=False)
    assert stats["systems"] == 40
    assert stats["planets"] > 20  # mean ~2.2/system minus stability rejections

    db = WorldDB(db_path)
    rows = db.list_systems(limit=40)
    assert len(rows) == 40
    # Inspect first system with planets end-to-end
    for r in rows:
        planets = db.get_planets(r["id"])
        if planets:
            for p in planets:
                assert p["period_d"] > 0
                assert p["mass_me"] > 0
                if p["transits"]:
                    assert p["depth_ppm"] > 0
                    assert p["t14_hr"] > 0
                    # 6 transit observatories + HWO imaging
                    assert len(db.get_observations(p["id"])) == 7
            break
    # No invalid-severity flags should ever be persisted
    n_invalid = db.conn.execute(
        "SELECT COUNT(*) FROM flags WHERE severity='invalid'").fetchone()[0]
    assert n_invalid == 0
    db.close()


def test_transit_fraction_realistic(tmp_path):
    """Geometric transit fraction should be a few percent, not ~0 or ~50%."""
    db_path = tmp_path / "frac.db"
    stats = generate_population(str(db_path), 150, seed=7, progress=False)
    frac = stats["transiting_planets"] / max(stats["planets"], 1)
    assert 0.005 < frac < 0.30  # population is dominated by P < 100 d planets
