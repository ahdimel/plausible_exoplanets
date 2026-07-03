"""Command-line interface: generate, list, inspect, stats, lightcurve,
archive, validate, serve."""

from __future__ import annotations

import argparse
import json
import sys

from .constants import R_EARTH, R_JUP
from .database import WorldDB
from .generate import generate_population


def cmd_generate(args) -> None:
    print(f"Generating {args.n} systems (seed={args.seed}, "
          f"dmax={args.dmax:.0f} pc) -> {args.db}")
    stats = generate_population(args.db, args.n, args.seed, dmax_pc=args.dmax)
    print(json.dumps(stats, indent=2, default=str))


def cmd_list(args) -> None:
    db = WorldDB(args.db)
    rows = db.list_systems(limit=args.limit, transiting_only=args.transiting,
                           detectable_only=args.detectable)
    print(f"{'id':>5} {'name':<16} {'sptype':<5} {'Teff':>6} {'d(pc)':>6} "
          f"{'V':>5} {'Npl':>3} {'Ntr':>3}")
    for r in rows:
        print(f"{r['id']:>5} {r['name']:<16} {r['st_sptype']:<5} "
              f"{r['st_teff']:>6.0f} {r['st_distance_pc']:>6.1f} "
              f"{r['st_mag_v']:>5.1f} {r['n_planets']:>3} {r['n_transiting']:>3}")
    db.close()


def _severity_marker(sev: str) -> str:
    return {"info": "i", "questionable": "?", "invalid": "X"}.get(sev, " ")


def cmd_inspect(args) -> None:
    db = WorldDB(args.db)
    s = db.get_system(args.system)
    if s is None:
        print(f"System '{args.system}' not found", file=sys.stderr)
        sys.exit(1)

    print(f"\n=== {s['name']} (id {s['id']}, seed {s['seed']}) ===")
    print(f"Star : {s['st_sptype']}  M={s['st_mass']:.2f} Msun  "
          f"R={s['st_radius']:.2f} Rsun  Teff={s['st_teff']:.0f} K  "
          f"L={s['st_lum']:.3f} Lsun")
    print(f"       [Fe/H]={s['st_feh']:+.2f}  age={s['st_age_gyr']:.1f} Gyr "
          f"(MS lifetime {s['st_ms_lifetime_gyr']:.1f} Gyr)  "
          f"d={s['st_distance_pc']:.1f} pc")
    print(f"       V={s['st_mag_v']:.2f}  T={s['st_mag_tess']:.2f}  "
          f"J={s['st_mag_j']:.2f}  limb dark u1={s['st_u1']:.2f} u2={s['st_u2']:.2f}")
    if s["st_prot_days"] is not None:
        print(f"Noise: Prot={s['st_prot_days']:.1f} d  "
              f"activity={s['st_activity']:.2f}  "
              f"spot variability~{s['st_var_amp_ppm']:.0f} ppm  "
              f"granulation={s['st_gran_1hr_ppm']:.0f} ppm/hr  "
              f"flares={s['st_flare_1hr_ppm']:.0f} ppm/hr")
    print(f"       system plane inclination: {s['sys_inc_deg']:.2f} deg "
          f"(90 = edge-on)")

    planets = db.get_planets(s["id"])
    print(f"\nPlanets ({len(planets)}):")
    for p in planets:
        rj = p["radius_re"] / (R_JUP / R_EARTH)
        size = (f"{p['radius_re']:.2f} Re" if rj < 0.5 else f"{rj:.2f} Rjup")
        hz = "  [HABITABLE ZONE]" if p["in_hz"] else ""
        print(f"\n  {s['name']} {p['letter']}  ({p['comp_class']}){hz}")
        print(f"    R={size}  M={p['mass_me']:.1f} Me  "
              f"rho={p['density_gcc']:.2f} g/cc")
        print(f"    P={p['period_d']:.3f} d  a={p['a_au']:.4f} AU  "
              f"e={p['ecc']:.3f}  i={p['inc_deg']:.2f} deg")
        print(f"    Teq={p['teq_k']:.0f} K  S={p['insolation_se']:.1f} S_earth")

        atm = db.get_atmosphere(p["id"])
        if atm:
            if atm["mu"] > 0:
                print(f"    Atmosphere: {atm['atm_class']}  mu={atm['mu']:.1f}  "
                      f"H={atm['scale_height_km']:.0f} km  "
                      f"feature~{atm['feature_ppm']:.0f} ppm "
                      f"(clouds x{atm['cloud_factor']:.2f})  "
                      f"TSM={atm['tsm']:.0f}{' *priority*' if atm['tsm_priority'] else ''}  "
                      f"ESM={atm['esm']:.0f}")
            else:
                print(f"    Atmosphere: {atm['atm_class']} (likely stripped; "
                      "cosmic shoreline)")

        if p["transits"]:
            print(f"    TRANSITS: depth={p['depth_ppm']:.0f} ppm "
                  f"(uniform {p['depth_uniform_ppm']:.0f})  b={p['impact_b']:.2f}  "
                  f"T14={p['t14_hr']:.2f} h  T23={p['t23_hr']:.2f} h  "
                  f"a-priori prob={p['prob_transit']*100:.1f}%")
        else:
            print(f"    no transit (b={p['impact_b']:.2f}, a-priori prob "
                  f"{p['prob_transit']*100:.1f}%)")
        for o in db.get_observations(p["id"]):
            if o["mode"] == "imaging":
                status = ("DETECTABLE" if o["detectable"] else "not detectable")
                print(f"      {o['observatory']:<22} contrast={o['contrast']:.1e}  "
                      f"sep={o['separation_mas']:.1f} mas  {status}  [{o['note']}]")
            else:
                status = "DETECTABLE" if o["detectable"] else (
                    "usable, below threshold" if o["usable"] else "unusable")
                print(f"      {o['observatory']:<22} "
                      f"sigma_1hr={o['sigma_1hr_ppm']:>7.0f} ppm  "
                      f"star={o['sigma_stellar_ppm']:>5.0f} ppm  "
                      f"SNR/tr={o['snr_per_transit']:>6.1f}  "
                      f"Ntr={o['n_transits']:>6.1f}  "
                      f"SNR={o['snr_total']:>7.1f}  {status}  [{o['note']}]")
        for ao in db.get_atm_observations(p["id"]):
            n5 = ao["n_transits_5sigma"]
            n5s = "impractical" if n5 < 0 or n5 > 1e5 else f"{n5:.1f} transits"
            print(f"      [spectroscopy] {ao['observatory']:<22} "
                  f"feature SNR/tr={ao['feature_snr_per_transit']:.2f}  "
                  f"5-sigma in {n5s}"
                  f"{'  PRACTICAL' if ao['practical'] else ''}")

    flags = db.get_flags(s["id"])
    if flags:
        print("\nPlausibility flags:")
        for f in flags:
            scope = f["scope"]
            if f["planet_id"]:
                letter = next((p["letter"] for p in planets
                               if p["id"] == f["planet_id"]), "?")
                scope = f"{f['scope']} {letter}"
            print(f"  [{_severity_marker(f['severity'])}] ({scope}) "
                  f"{f['rule']}: {f['message']}")
    print()
    db.close()


def cmd_stats(args) -> None:
    db = WorldDB(args.db)
    print(json.dumps(db.stats(), indent=2, default=str))
    db.close()


def cmd_lightcurve(args) -> None:
    from .system import generate_system
    from .transits import compute_geometry, model_light_curve

    db = WorldDB(args.db)
    s = db.get_system(args.system)
    if s is None:
        print(f"System '{args.system}' not found", file=sys.stderr)
        sys.exit(1)
    # deterministic re-generation (dmax_pc from meta keeps distances exact)
    system = generate_system(s["seed"], s["name"], dmax_pc=db.dmax_pc)
    letters = "bcdefghijklmn"
    idx = letters.index(args.planet)
    if idx >= len(system.planets):
        print(f"No planet '{args.planet}' in {s['name']}", file=sys.stderr)
        sys.exit(1)
    planet = system.planets[idx]
    geom = compute_geometry(system.star, planet)
    if not geom.transits:
        print(f"{s['name']} {args.planet} does not transit", file=sys.stderr)
        sys.exit(1)
    t, flux = model_light_curve(system.star, geom, n_points=500)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, (flux - 1.0) * 1e6, color="#1f4e8c", lw=1.8)
    ax.set_xlabel("Time from mid-transit (hours)")
    ax.set_ylabel("Relative flux change (ppm)")
    ax.set_title(f"{s['name']} {args.planet}: depth {geom.depth_ppm:.0f} ppm, "
                 f"T14 {geom.t14_hours:.2f} h, b={geom.b:.2f}")
    ax.grid(alpha=0.3)
    out = args.out or f"{s['name']}_{args.planet}_lightcurve.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"Saved {out}")
    db.close()


def cmd_archive(args) -> None:
    from .archive import load_snapshot, refresh_snapshot
    if args.refresh:
        try:
            path, source = refresh_snapshot(args.data_dir, prefer=args.prefer)
            print(f"Snapshot refreshed from source '{source}' -> {path}")
        except RuntimeError as e:
            print(f"Refresh failed: {e}", file=sys.stderr)
            sys.exit(1)
    try:
        planets, meta = load_snapshot(args.data_dir)
        n_tran = sum(1 for p in planets if p.tran_flag == 1.0)
        print(f"Snapshot: {meta['count']} confirmed planets "
              f"({n_tran} transiting), source={meta['source']}, "
              f"fetched={meta['fetched']}")
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)


def cmd_validate(args) -> None:
    from .validate import format_report, run_validation
    report = run_validation(args.db, args.data_dir)
    print(format_report(report))
    print(f"\nFull report: {report['report_path']}")


def cmd_serve(args) -> None:
    from .web.app import create_app
    app = create_app(args.db, data_dir=args.data_dir)
    print(f"exoverse web UI: http://127.0.0.1:{args.port}")
    app.run(host=args.host, port=args.port, debug=args.debug)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        prog="exoverse",
        description="Procedurally generated, physics-validated exoplanet worlds")
    ap.add_argument("--db", default="worlds.db", help="SQLite database path")
    ap.add_argument("--data-dir", default="data",
                    help="directory for archive snapshots and reports")
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("generate", help="generate a population")
    g.add_argument("--n", type=int, default=500)
    g.add_argument("--seed", type=int, default=42)
    g.add_argument("--dmax", type=float, default=300.0,
                   help="max host distance in pc (30 = solar neighborhood "
                        "for direct-imaging studies; default 300)")
    g.set_defaults(func=cmd_generate)

    l = sub.add_parser("list", help="list systems")
    l.add_argument("--limit", type=int, default=50)
    l.add_argument("--transiting", action="store_true")
    l.add_argument("--detectable", action="store_true")
    l.set_defaults(func=cmd_list)

    i = sub.add_parser("inspect", help="inspect one system in full detail")
    i.add_argument("system", help="system name or id")
    i.set_defaults(func=cmd_inspect)

    st = sub.add_parser("stats", help="population statistics")
    st.set_defaults(func=cmd_stats)

    lc = sub.add_parser("lightcurve", help="plot a planet's transit light curve")
    lc.add_argument("system")
    lc.add_argument("planet", help="planet letter, e.g. b")
    lc.add_argument("--out", default=None)
    lc.set_defaults(func=cmd_lightcurve)

    ar = sub.add_parser("archive", help="fetch/inspect the real-planet catalog")
    ar.add_argument("--refresh", action="store_true")
    ar.add_argument("--prefer", choices=["nasa", "eu"], default="nasa")
    ar.set_defaults(func=cmd_archive)

    va = sub.add_parser("validate",
                        help="validate generator against real exoplanets")
    va.set_defaults(func=cmd_validate)

    se = sub.add_parser("serve", help="run the browser UI")
    se.add_argument("--host", default="127.0.0.1")
    se.add_argument("--port", type=int, default=8321)
    se.add_argument("--debug", action="store_true")
    se.set_defaults(func=cmd_serve)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
