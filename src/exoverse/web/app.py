"""Flask web UI for browsing the generated worlds database.

Read-only over the SQLite database plus the validation report JSON; light
curves are re-simulated on demand from each system's deterministic seed.
Designed to be deployable later for collaboration (no local-only paths in
templates; DB path injected at create_app time).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from flask import Flask, abort, render_template, request

from ..constants import R_EARTH, R_JUP
from ..database import WorldDB
from . import charts


def create_app(db_path: str, data_dir: str = "data") -> Flask:
    app = Flask(__name__)
    app.config["DB_PATH"] = db_path
    app.config["DATA_DIR"] = data_dir

    def db() -> WorldDB:
        return WorldDB(app.config["DB_PATH"])

    @app.template_filter("num")
    def num_filter(v, digits=1):
        if v is None:
            return "—"
        return f"{v:,.{digits}f}"

    @app.template_filter("sizefmt")
    def sizefmt(radius_re):
        rj = radius_re / (R_JUP / R_EARTH)
        return f"{rj:.2f} R♃" if rj >= 0.5 else f"{radius_re:.2f} R⊕"

    # ------------------------------------------------------------ dashboard
    @app.route("/")
    def dashboard():
        d = db()
        stats = d.stats()
        det = stats["detections_by_observatory"]
        det_rows = [(k, v or 0, f"{v or 0} detections") for k, v in
                    sorted(det.items(), key=lambda kv: -(kv[1] or 0))]
        pts = []
        for r in d.conn.execute(
                """SELECT p.period_d, p.radius_re, p.comp_class, p.transits,
                          p.letter, s.name FROM planets p
                   JOIN systems s ON s.id=p.system_id""").fetchall():
            pts.append({
                "x": r["period_d"], "y": r["radius_re"], "cls": r["comp_class"],
                "ring": bool(r["transits"]),
                "tip": (f"{r['name']} {r['letter']} — {r['comp_class']}, "
                        f"{r['radius_re']:.2f} R⊕, P={r['period_d']:.2f} d"
                        + (" — transits" if r["transits"] else "")),
                "href": f"/system/{r['name']}",
            })
        scatter = charts.scatter_log_log(pts)
        legend = [(cls, charts.CLASS_SLOTS[cls]) for cls in
                  ("rocky", "sub-neptune", "neptunian", "giant")]
        d.close()
        return render_template("dashboard.html", stats=stats,
                               det_chart=charts.hbar(det_rows),
                               scatter=scatter, legend=legend)

    # ------------------------------------------------------------- systems
    @app.route("/systems")
    def systems():
        d = db()
        page = max(int(request.args.get("page", 1)), 1)
        per = 50
        transiting = request.args.get("transiting") == "1"
        detectable = request.args.get("detectable") == "1"
        rows = d.list_systems(limit=per, offset=(page - 1) * per,
                              transiting_only=transiting,
                              detectable_only=detectable)
        d.close()
        return render_template("systems.html", rows=rows, page=page,
                               transiting=transiting, detectable=detectable,
                               has_next=len(rows) == per)

    # ------------------------------------------------------------- detail
    @app.route("/system/<name>")
    def system(name):
        d = db()
        s = d.get_system(name)
        if s is None:
            abort(404)
        planets = []
        for p in d.get_planets(s["id"]):
            planets.append({
                "row": p,
                "atm": d.get_atmosphere(p["id"]),
                "obs": d.get_observations(p["id"]),
                "atm_obs": d.get_atm_observations(p["id"]),
            })
        flags = d.get_flags(s["id"])
        letters = {p["row"]["id"]: p["row"]["letter"] for p in planets}
        d.close()

        # Light curves re-simulated deterministically for transiting planets
        curves = []
        from ..system import generate_system
        from ..transits import compute_geometry, model_light_curve
        gen = generate_system(s["seed"], s["name"])
        for i, gp in enumerate(gen.planets):
            geom = compute_geometry(gen.star, gp)
            if geom.transits and geom.t14_hours > 0:
                t, flux = model_light_curve(gen.star, geom, n_points=260)
                curves.append({
                    "letter": "bcdefghijklmn"[i],
                    "svg": charts.light_curve_svg(list(t), list(flux),
                                                  geom.depth_ppm),
                })
        return render_template("system.html", s=s, planets=planets,
                               flags=flags, letters=letters, curves=curves)

    # ---------------------------------------------------------- validation
    @app.route("/validation")
    def validation():
        path = Path(app.config["DATA_DIR"]) / "validation_report.json"
        if not path.exists():
            return render_template("validation.html", report=None,
                                   hist_p=None, hist_r=None)
        report = json.loads(path.read_text())

        d = db()
        syn = d.conn.execute(
            """SELECT p.period_d, p.radius_re FROM planets p
               JOIN observations o ON o.planet_id=p.id
               WHERE o.observatory LIKE 'Kepler%' AND o.detectable=1""").fetchall()
        d.close()
        real_p, real_r = [], []
        try:
            from ..archive import load_snapshot
            planets, _ = load_snapshot(app.config["DATA_DIR"])
            for p in planets:
                if p.tran_flag == 1.0 and "transit" in p.discoverymethod.lower():
                    if p.pl_orbper and 0.3 < p.pl_orbper < 1000:
                        real_p.append(math.log10(p.pl_orbper))
                    if p.pl_rade and 0.3 < p.pl_rade < 25:
                        real_r.append(math.log10(p.pl_rade))
        except FileNotFoundError:
            pass

        hist_p = charts.overlay_step_hist(
            [{"name": "synthetic (Kepler-detectable)", "slot": 1,
              "values": [math.log10(r[0]) for r in syn]},
             {"name": "real (transit-discovered)", "slot": 3,
              "values": real_p}],
            lo=-0.5, hi=3.0, x_label="log10 orbital period (days)")
        hist_r = charts.overlay_step_hist(
            [{"name": "synthetic (Kepler-detectable)", "slot": 1,
              "values": [math.log10(r[1]) for r in syn]},
             {"name": "real (transit-discovered)", "slot": 3,
              "values": real_r}],
            lo=-0.4, hi=1.4, x_label="log10 planet radius (R⊕)")
        return render_template("validation.html", report=report,
                               hist_p=hist_p, hist_r=hist_r)

    @app.route("/about")
    def about():
        return render_template("about.html")

    return app
