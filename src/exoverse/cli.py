"""Command-line interface: generate, list, inspect, stats, lightcurve."""

from __future__ import annotations

import argparse
import json
import sys

from .constants import R_EARTH, R_JUP
from .database import WorldDB
from .generate import generate_population


def cmd_generate(args) -> None:
    print(f"Generating {args.n} systems (seed={args.seed}) -> {args.db}")
    stats = generate_population(args.db, args.n, args.seed)
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
        if p["transits"]:
            print(f"    TRANSITS: depth={p['depth_ppm']:.0f} ppm "
                  f"(uniform {p['depth_uniform_ppm']:.0f})  b={p['impact_b']:.2f}  "
                  f"T14={p['t14_hr']:.2f} h  T23={p['t23_hr']:.2f} h  "
                  f"a-priori prob={p['prob_transit']*100:.1f}%")
            for o in db.get_observations(p["id"]):
                status = "DETECTABLE" if o["detectable"] else (
                    "usable, below threshold" if o["usable"] else "unusable")
                print(f"      {o['observatory']:<22} sigma_1hr={o['sigma_1hr_ppm']:>7.0f} ppm  "
                      f"SNR/tr={o['snr_per_transit']:>6.1f}  Ntr={o['n_transits']:>6.1f}  "
                      f"SNR={o['snr_total']:>7.1f}  {status}  [{o['note']}]")
        else:
            print(f"    no transit (b={p['impact_b']:.2f}, a-priori prob "
                  f"{p['prob_transit']*100:.1f}%)")

    flags = db.get_flags(s["id"])
    if flags:
        print("\nPlausibility flags:")
        for f in flags:
            scope = f["scope"]
            if f["planet_id"]:
                letter = next((p["letter"] for p in planets if p["id"] == f["planet_id"]), "?")
                scope = f"planet {letter}"
            print(f"  [{_severity_marker(f['severity'])}] ({scope}) "
                  f"{f['rule']}: {f['message']}")
    print()
    db.close()


def cmd_stats(args) -> None:
    db = WorldDB(args.db)
    print(json.dumps(db.stats(), indent=2, default=str))
    db.close()


def cmd_lightcurve(args) -> None:
    """Re-simulate and plot the light curve for one planet (requires matplotlib)."""
    from .system import generate_system
    from .transits import compute_geometry, model_light_curve

    db = WorldDB(args.db)
    s = db.get_system(args.system)
    if s is None:
        print(f"System '{args.system}' not found", file=sys.stderr)
        sys.exit(1)
    system = generate_system(s["seed"], s["name"])  # deterministic re-generation
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


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        prog="exoverse",
        description="Procedurally generated, physics-validated exoplanet worlds")
    ap.add_argument("--db", default="worlds.db", help="SQLite database path")
    sub = ap.add_subparsers(dest="cmd", required=True)

    g = sub.add_parser("generate", help="generate a population")
    g.add_argument("--n", type=int, default=500)
    g.add_argument("--seed", type=int, default=42)
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

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
