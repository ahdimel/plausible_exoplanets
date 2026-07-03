"""SQLite persistence for generated systems.

Schema v2: normalized tables (systems / planets / observations / atmospheres /
atm_observations / flags). Full metadata is retained column-by-column so each
world can be inspected individually without any re-computation.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import List, Optional

from .observatories import Observation
from .system import StellarSystem
from .transits import TransitGeometry

SCHEMA_VERSION = 2

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS systems (
    id INTEGER PRIMARY KEY,
    name TEXT UNIQUE NOT NULL,
    seed INTEGER NOT NULL,
    n_planets INTEGER NOT NULL,
    n_transiting INTEGER NOT NULL,
    sys_inc_deg REAL,
    -- star
    st_mass REAL, st_radius REAL, st_lum REAL, st_teff REAL, st_feh REAL,
    st_age_gyr REAL, st_ms_lifetime_gyr REAL, st_distance_pc REAL,
    st_mag_v REAL, st_mag_tess REAL, st_mag_j REAL,
    st_u1 REAL, st_u2 REAL, st_density REAL, st_sptype TEXT,
    -- stellar noise state
    st_prot_days REAL, st_activity REAL, st_var_amp_ppm REAL,
    st_gran_1hr_ppm REAL, st_osc_1hr_ppm REAL, st_flare_1hr_ppm REAL
);
CREATE TABLE IF NOT EXISTS planets (
    id INTEGER PRIMARY KEY,
    system_id INTEGER NOT NULL REFERENCES systems(id),
    letter TEXT NOT NULL,
    radius_re REAL, mass_me REAL, density_gcc REAL, comp_class TEXT,
    period_d REAL, a_au REAL, ecc REAL, omega_deg REAL, inc_deg REAL,
    teq_k REAL, insolation_se REAL, in_hz INTEGER,
    transits INTEGER, impact_b REAL, depth_ppm REAL, depth_uniform_ppm REAL,
    t14_hr REAL, t23_hr REAL, prob_transit REAL, radius_ratio REAL
);
CREATE TABLE IF NOT EXISTS observations (
    id INTEGER PRIMARY KEY,
    planet_id INTEGER NOT NULL REFERENCES planets(id),
    observatory TEXT NOT NULL,
    mode TEXT NOT NULL DEFAULT 'transit',
    usable INTEGER, note TEXT,
    sigma_1hr_ppm REAL, sigma_stellar_ppm REAL,
    snr_per_transit REAL, n_transits REAL,
    snr_total REAL, detectable INTEGER,
    contrast REAL, separation_mas REAL
);
CREATE TABLE IF NOT EXISTS atmospheres (
    id INTEGER PRIMARY KEY,
    planet_id INTEGER UNIQUE NOT NULL REFERENCES planets(id),
    atm_class TEXT, mu REAL, scale_height_km REAL,
    delta_1h_ppm REAL, feature_ppm REAL, cloud_factor REAL,
    tsm REAL, esm REAL, tsm_priority INTEGER
);
CREATE TABLE IF NOT EXISTS atm_observations (
    id INTEGER PRIMARY KEY,
    planet_id INTEGER NOT NULL REFERENCES planets(id),
    observatory TEXT NOT NULL,
    feature_snr_per_transit REAL,
    n_transits_5sigma REAL,           -- -1 encodes 'impractical/infinite'
    practical INTEGER,
    note TEXT
);
CREATE TABLE IF NOT EXISTS flags (
    id INTEGER PRIMARY KEY,
    scope TEXT NOT NULL,               -- 'system' | 'star' | 'planet' | 'atmosphere'
    system_id INTEGER NOT NULL REFERENCES systems(id),
    planet_id INTEGER REFERENCES planets(id),
    severity TEXT NOT NULL,
    rule TEXT NOT NULL,
    message TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_planets_system ON planets(system_id);
CREATE INDEX IF NOT EXISTS idx_obs_planet ON observations(planet_id);
CREATE INDEX IF NOT EXISTS idx_atm_planet ON atmospheres(planet_id);
CREATE INDEX IF NOT EXISTS idx_atmobs_planet ON atm_observations(planet_id);
CREATE INDEX IF NOT EXISTS idx_flags_system ON flags(system_id);
"""

PLANET_LETTERS = "bcdefghijklmn"


class WorldDB:
    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.conn = sqlite3.connect(self.path)
        self.conn.row_factory = sqlite3.Row
        self.conn.executescript(SCHEMA)
        row = self.conn.execute("SELECT value FROM meta WHERE key='schema'").fetchone()
        if row is None:
            self.conn.execute("INSERT INTO meta VALUES ('schema', ?)",
                              (str(SCHEMA_VERSION),))
            self.conn.commit()
        elif int(row[0]) != SCHEMA_VERSION:
            raise RuntimeError(
                f"{self.path} has schema v{row[0]}, this code needs "
                f"v{SCHEMA_VERSION}. Regenerate the database.")

    def close(self) -> None:
        self.conn.close()

    # ------------------------------------------------------------------ meta
    def set_meta(self, key: str, value: str) -> None:
        self.conn.execute("INSERT OR REPLACE INTO meta VALUES (?, ?)",
                          (key, str(value)))
        self.conn.commit()

    def get_meta(self, key: str, default: Optional[str] = None) -> Optional[str]:
        row = self.conn.execute("SELECT value FROM meta WHERE key=?",
                                (key,)).fetchone()
        return row[0] if row is not None else default

    @property
    def dmax_pc(self) -> float:
        """Distance cap the population was generated with (pre-existing DBs
        without the key were all generated at 300 pc)."""
        return float(self.get_meta("dmax_pc", "300.0"))

    # ------------------------------------------------------------------ save
    def save_system(self, system: StellarSystem,
                    geometries: List[TransitGeometry],
                    observations: List[List[Observation]],
                    atm_observations: Optional[List[list]] = None) -> int:
        st = system.star
        nz = system.noise
        n_transiting = sum(1 for g in geometries if g.transits)
        cur = self.conn.execute(
            """INSERT INTO systems (name, seed, n_planets, n_transiting,
               sys_inc_deg, st_mass, st_radius, st_lum, st_teff, st_feh,
               st_age_gyr, st_ms_lifetime_gyr, st_distance_pc, st_mag_v,
               st_mag_tess, st_mag_j, st_u1, st_u2, st_density, st_sptype,
               st_prot_days, st_activity, st_var_amp_ppm, st_gran_1hr_ppm,
               st_osc_1hr_ppm, st_flare_1hr_ppm)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (system.name, system.seed, len(system.planets), n_transiting,
             system.sys_inc_deg, st.mass, st.radius, st.luminosity, st.teff,
             st.feh, st.age_gyr, st.ms_lifetime_gyr, st.distance_pc,
             st.mag_v, st.mag_tess, st.mag_j, st.u1, st.u2, st.density,
             st.spectral_type,
             nz.prot_days if nz else None,
             nz.activity_level if nz else None,
             nz.var_amp_ppm if nz else None,
             nz.gran_1hr_ppm if nz else None,
             nz.osc_1hr_ppm if nz else None,
             nz.flare_1hr_ppm if nz else None))
        system_id = cur.lastrowid

        for f in system.flags:
            self._save_flag("system", system_id, None, f)
        for f in st.flags:
            self._save_flag("star", system_id, None, f)

        for i, (p, g, obs) in enumerate(zip(system.planets, geometries, observations)):
            cur = self.conn.execute(
                """INSERT INTO planets (system_id, letter, radius_re, mass_me,
                   density_gcc, comp_class, period_d, a_au, ecc, omega_deg,
                   inc_deg, teq_k, insolation_se, in_hz, transits, impact_b,
                   depth_ppm, depth_uniform_ppm, t14_hr, t23_hr, prob_transit,
                   radius_ratio)
                   VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                (system_id, PLANET_LETTERS[i], p.radius, p.mass, p.density,
                 p.composition_class, p.period, p.a, p.ecc, p.omega_deg,
                 p.inc_deg, p.teq, p.insolation, int(p.in_habitable_zone),
                 int(g.transits), g.b, g.depth_ppm, g.depth_uniform_ppm,
                 g.t14_hours, g.t23_hours, g.prob_transit, g.k))
            planet_id = cur.lastrowid
            for f in p.flags:
                self._save_flag("planet", system_id, planet_id, f)
            for o in obs:
                self.conn.execute(
                    """INSERT INTO observations (planet_id, observatory, mode,
                       usable, note, sigma_1hr_ppm, sigma_stellar_ppm,
                       snr_per_transit, n_transits, snr_total, detectable,
                       contrast, separation_mas)
                       VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)""",
                    (planet_id, o.observatory, o.mode, int(o.usable), o.note,
                     o.sigma_1hr_ppm, o.sigma_stellar_ppm, o.snr_per_transit,
                     o.n_transits, o.snr_total, int(o.detectable),
                     o.contrast, o.separation_mas))

            atm = p.atmosphere
            if atm is not None:
                self.conn.execute(
                    """INSERT INTO atmospheres (planet_id, atm_class, mu,
                       scale_height_km, delta_1h_ppm, feature_ppm,
                       cloud_factor, tsm, esm, tsm_priority)
                       VALUES (?,?,?,?,?,?,?,?,?,?)""",
                    (planet_id, atm.atm_class, atm.mu, atm.scale_height_km,
                     atm.delta_1h_ppm, atm.feature_ppm, atm.cloud_factor,
                     atm.tsm, atm.esm, int(atm.tsm_priority)))
                for f in atm.flags:
                    self._save_flag("atmosphere", system_id, planet_id, f)
                if atm_observations is not None:
                    for ao in atm_observations[i]:
                        n5 = (-1.0 if ao.n_transits_5sigma == float("inf")
                              else ao.n_transits_5sigma)
                        self.conn.execute(
                            """INSERT INTO atm_observations (planet_id,
                               observatory, feature_snr_per_transit,
                               n_transits_5sigma, practical, note)
                               VALUES (?,?,?,?,?,?)""",
                            (planet_id, ao.observatory,
                             ao.feature_snr_per_transit, n5,
                             int(ao.practical), ao.note))
        self.conn.commit()
        return system_id

    def _save_flag(self, scope: str, system_id: int, planet_id: Optional[int], f) -> None:
        self.conn.execute(
            "INSERT INTO flags (scope, system_id, planet_id, severity, rule, message)"
            " VALUES (?,?,?,?,?,?)",
            (scope, system_id, planet_id, f.severity.value, f.rule, f.message))

    # ----------------------------------------------------------------- query
    def get_system(self, name_or_id: str) -> Optional[sqlite3.Row]:
        if str(name_or_id).isdigit():
            row = self.conn.execute("SELECT * FROM systems WHERE id=?",
                                    (int(name_or_id),)).fetchone()
            if row:
                return row
        return self.conn.execute("SELECT * FROM systems WHERE name=?",
                                 (str(name_or_id),)).fetchone()

    def get_planets(self, system_id: int) -> List[sqlite3.Row]:
        return self.conn.execute(
            "SELECT * FROM planets WHERE system_id=? ORDER BY period_d",
            (system_id,)).fetchall()

    def get_observations(self, planet_id: int) -> List[sqlite3.Row]:
        return self.conn.execute(
            "SELECT * FROM observations WHERE planet_id=?", (planet_id,)).fetchall()

    def get_atmosphere(self, planet_id: int) -> Optional[sqlite3.Row]:
        return self.conn.execute(
            "SELECT * FROM atmospheres WHERE planet_id=?", (planet_id,)).fetchone()

    def get_atm_observations(self, planet_id: int) -> List[sqlite3.Row]:
        return self.conn.execute(
            "SELECT * FROM atm_observations WHERE planet_id=?",
            (planet_id,)).fetchall()

    def get_flags(self, system_id: int) -> List[sqlite3.Row]:
        return self.conn.execute(
            "SELECT * FROM flags WHERE system_id=? ORDER BY scope, planet_id",
            (system_id,)).fetchall()

    def list_systems(self, limit: int = 50, offset: int = 0,
                     transiting_only: bool = False,
                     detectable_only: bool = False) -> List[sqlite3.Row]:
        q = "SELECT DISTINCT s.* FROM systems s"
        where = []
        if detectable_only:
            q += (" JOIN planets p ON p.system_id = s.id"
                  " JOIN observations o ON o.planet_id = p.id")
            where.append("o.detectable = 1")
        if transiting_only:
            where.append("s.n_transiting > 0")
        if where:
            q += " WHERE " + " AND ".join(where)
        q += " ORDER BY s.id LIMIT ? OFFSET ?"
        return self.conn.execute(q, (limit, offset)).fetchall()

    def stats(self) -> dict:
        c = self.conn
        row = lambda q: c.execute(q).fetchone()[0]
        return {
            "systems": row("SELECT COUNT(*) FROM systems"),
            "planets": row("SELECT COUNT(*) FROM planets"),
            "transiting_planets": row("SELECT COUNT(*) FROM planets WHERE transits=1"),
            "hz_planets": row("SELECT COUNT(*) FROM planets WHERE in_hz=1"),
            "questionable_flags": row(
                "SELECT COUNT(*) FROM flags WHERE severity='questionable'"),
            "detections_by_observatory": dict(c.execute(
                "SELECT observatory, SUM(detectable) FROM observations"
                " GROUP BY observatory").fetchall()),
            "planet_classes": dict(c.execute(
                "SELECT comp_class, COUNT(*) FROM planets GROUP BY comp_class"
            ).fetchall()),
            "atmosphere_classes": dict(c.execute(
                "SELECT atm_class, COUNT(*) FROM atmospheres GROUP BY atm_class"
            ).fetchall()),
            "tsm_priority_targets": row(
                "SELECT COUNT(*) FROM atmospheres WHERE tsm_priority=1"),
            "practical_atm_targets": row(
                "SELECT COUNT(DISTINCT planet_id) FROM atm_observations"
                " WHERE practical=1"),
            "hwo_detections_by_class": dict(c.execute(
                "SELECT p.comp_class, COUNT(*) FROM planets p"
                " JOIN observations o ON o.planet_id = p.id"
                " WHERE o.observatory LIKE 'HWO%' AND o.detectable=1"
                " GROUP BY p.comp_class").fetchall()),
            # Exo-Earth candidates: HZ planets of 0.8-1.4 R_earth directly
            # imageable by HWO (the mission's headline ">=25 exo-Earths" metric)
            "hwo_exo_earth_candidates": row(
                "SELECT COUNT(*) FROM planets p"
                " JOIN observations o ON o.planet_id = p.id"
                " WHERE o.observatory LIKE 'HWO%' AND o.detectable=1"
                " AND p.in_hz=1 AND p.radius_re BETWEEN 0.8 AND 1.4"),
        }
