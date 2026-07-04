"""Phase 1 validation plots for the dichotomy-decomposition architecture knobs.

Generates four figures + stats.json into results/phase1_validation/:

  fig1_marginal_radius.png   — marginal small-planet radius CDF is invariant
                               under sigma_r (Gaussian-copula guarantee,
                               docs/sigma_r_note.md); KS D vs the None case.
  fig2_size_uniformity.png   — |log10(R_out/R_in)| over adjacent small-planet
                               pairs tightens monotonically as sigma_r drops.
  fig3_inclination_scatter.png — projected inclination scatter (inc - sys_inc)
                               has std ~ sigma_i; isotropic limit has uniform
                               cos(inc).
  fig4_dichotomy_spm.png     — geometric singles-per-multi (transit geometry
                               only) rises monotonically with sigma_i;
                               isotropic limit as reference.

Deterministic: every cell owns a np.random.SeedSequence([BASE_SEED, fig, cell])
stream (same style as generate.py). Run end-to-end with

    .venv/bin/python analysis/phase1_validation_plots.py
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from exoverse.architecture import Architecture  # noqa: E402
from exoverse.system import generate_system  # noqa: E402
from exoverse.transits import compute_geometry  # noqa: E402

BASE_SEED = 20260703
SMALL_R_MAX = 4.05          # R_earth: the small-planet branch (giants excluded)
OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "phase1_validation"

N_FIG12 = 30_000            # systems per sigma_r cell (figures 1 & 2)
N_FIG3 = 20_000             # systems per sigma_i cell (figure 3)
N_FIG4 = 100_000            # systems per sigma_i cell (figure 4)

# ---------------------------------------------------------------------------
# Style: dataviz-skill reference palette (light mode), roles not raw intent.
# ---------------------------------------------------------------------------
SERIES = ["#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRIDLINE = "#e1e0d9"
BASELINE = "#c3c2b7"

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "savefig.dpi": 150,
    "font.family": "sans-serif",
    "font.size": 10,
    "text.color": INK,
    "axes.edgecolor": BASELINE,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "axes.titlesize": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": GRIDLINE,
    "grid.linewidth": 0.8,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "legend.frameon": False,
    "lines.linewidth": 2.0,
})


def cell_seeds(fig: int, cell: int, n: int) -> np.ndarray:
    """Distinct deterministic seed stream per (figure, cell)."""
    return np.random.SeedSequence([BASE_SEED, fig, cell]).generate_state(n)


def ks_two_sample(a: np.ndarray, b: np.ndarray) -> float:
    a = np.sort(a)
    b = np.sort(b)
    grid = np.concatenate([a, b])
    grid.sort(kind="mergesort")
    cdf_a = np.searchsorted(a, grid, side="right") / a.size
    cdf_b = np.searchsorted(b, grid, side="right") / b.size
    return float(np.max(np.abs(cdf_a - cdf_b)))


def ecdf_curve(x: np.ndarray, n_points: int = 1500):
    """Decimated ECDF for plotting (quantile-spaced, exact at both ends)."""
    x = np.sort(x)
    q = np.linspace(0.0, 1.0, min(n_points, x.size))
    idx = np.minimum((q * (x.size - 1)).astype(int), x.size - 1)
    return x[idx], (idx + 1) / x.size


def label_sigma_r(v) -> str:
    return "None (independent)" if v is None else f"{v:g}"


# ---------------------------------------------------------------------------
# Generation cells
# ---------------------------------------------------------------------------

def gen_radius_cells():
    """One pass serves figures 1 and 2: pooled small radii + adjacent-pair
    |dlog10 R| per sigma_r value (union of both figures' grids)."""
    sigma_r_values = [None, 1.0, 0.5, 0.3, 0.1]
    radii = {}
    dlogr = {}
    for cell, sr in enumerate(sigma_r_values):
        t0 = time.time()
        arch = Architecture() if sr is None else Architecture(sigma_r=sr)
        pooled = []
        pairs = []
        for s in cell_seeds(fig=1, cell=cell, n=N_FIG12):
            system = generate_system(int(s), "VAL1", arch=arch)
            planets = sorted(system.planets, key=lambda p: p.period)
            rs = [p.radius for p in planets]
            pooled.extend(r for r in rs if r < SMALL_R_MAX)
            for r_in, r_out in zip(rs, rs[1:]):
                if r_in < SMALL_R_MAX and r_out < SMALL_R_MAX:
                    pairs.append(abs(math.log10(r_out / r_in)))
        radii[sr] = np.asarray(pooled)
        dlogr[sr] = np.asarray(pairs)
        print(f"  sigma_r={label_sigma_r(sr):>18}: {N_FIG12} systems, "
              f"{radii[sr].size} small radii, {dlogr[sr].size} pairs "
              f"({time.time() - t0:.1f} s)", flush=True)
    return radii, dlogr


def gen_inclination_cells():
    """Figure 3: projected inclination offsets per sigma_i + isotropic cos i."""
    sigma_i_values = [0.5, 1.5, 5.0, 30.0]
    offsets = {}
    for cell, si in enumerate(sigma_i_values):
        t0 = time.time()
        arch = Architecture(sigma_i=si)
        vals = []
        for s in cell_seeds(fig=3, cell=cell, n=N_FIG3):
            system = generate_system(int(s), "VAL3", arch=arch)
            vals.extend(p.inc_deg - system.sys_inc_deg for p in system.planets)
        offsets[si] = np.asarray(vals)
        print(f"  sigma_i={si:>5}: {offsets[si].size} planets "
              f"({time.time() - t0:.1f} s)", flush=True)

    t0 = time.time()
    cos_i = []
    for s in cell_seeds(fig=3, cell=len(sigma_i_values), n=N_FIG3):
        system = generate_system(int(s), "VAL3", arch=Architecture(isotropic=True))
        cos_i.extend(math.cos(math.radians(p.inc_deg)) for p in system.planets)
    cos_i = np.asarray(cos_i)
    print(f"  isotropic: {cos_i.size} planets ({time.time() - t0:.1f} s)",
          flush=True)
    return offsets, cos_i


def gen_dichotomy_cells():
    """Figure 4: geometric transit multiplicity per sigma_i (+ isotropic)."""
    sigma_i_values = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0]
    cells = list(sigma_i_values) + ["isotropic"]
    spm = {}
    multi_frac = {}
    for cell, key in enumerate(cells):
        t0 = time.time()
        arch = (Architecture(isotropic=True) if key == "isotropic"
                else Architecture(sigma_i=key))
        n_single = 0
        n_multi = 0
        for s in cell_seeds(fig=4, cell=cell, n=N_FIG4):
            system = generate_system(int(s), "VAL4", arch=arch)
            n_tr = sum(compute_geometry(system.star, p).transits
                       for p in system.planets)
            if n_tr == 1:
                n_single += 1
            elif n_tr >= 2:
                n_multi += 1
        spm[key] = n_single / n_multi
        multi_frac[key] = n_multi / (n_single + n_multi)
        print(f"  sigma_i={key!s:>9}: singles={n_single} multis={n_multi} "
              f"SPM={spm[key]:.2f} ({time.time() - t0:.1f} s)", flush=True)
    return spm, multi_frac


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def fig1_marginal(radii: dict, stats: dict) -> None:
    show = [None, 0.5, 0.1]
    base = radii[None]
    ks = {f"{sr:g}": ks_two_sample(radii[sr], base) for sr in show[1:]}
    stats["marginal_ks"] = ks
    stats["marginal_n"] = {label_sigma_r(sr): int(radii[sr].size) for sr in show}

    worst = max(ks.values())
    if worst > 0.02:
        # Load-bearing check failed: report, do not prettify.
        stats["marginal_ks_FAILED"] = True
        (OUT_DIR / "stats.json").write_text(json.dumps(stats, indent=2))
        sys.exit(f"FAIL: marginal radius distribution moves under sigma_r "
                 f"(max KS D = {worst:.4f} > 0.02). The copula marginal-"
                 f"preservation guarantee is broken; see stats.json.")

    fig, (ax, axd) = plt.subplots(
        2, 1, figsize=(7.2, 6.4), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1.0], "hspace": 0.12})

    dashes = [(), (5, 2), (1.5, 1.5)]
    for i, sr in enumerate(show):
        x, y = ecdf_curve(radii[sr])
        ax.plot(x, y, color=SERIES[i], dashes=dashes[i],
                label=f"sigma_r = {label_sigma_r(sr)}")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Marginal small-planet radius CDF is invariant under sigma_r")
    ax.legend(loc="lower right", title=f"{N_FIG12:,} systems each")
    ax.set_ylim(0, 1)

    grid = np.linspace(0.4, SMALL_R_MAX, 800)
    base_sorted = np.sort(base)
    cdf_base = np.searchsorted(base_sorted, grid, side="right") / base.size
    for i, sr in enumerate(show[1:], start=1):
        a = np.sort(radii[sr])
        cdf = np.searchsorted(a, grid, side="right") / a.size
        axd.plot(grid, (cdf - cdf_base) * 1e3, color=SERIES[i],
                 dashes=dashes[i],
                 label=f"sigma_r = {sr:g}  (KS D = {ks[f'{sr:g}']:.4f})")
    axd.axhline(0.0, color=BASELINE, lw=1.0)
    axd.set_xlabel("Planet radius (Earth radii)")
    axd.set_ylabel("Delta CDF vs None (x 10^-3)")
    axd.set_ylim(-8.5, 8.5)
    axd.legend(loc="upper right")
    fig.savefig(OUT_DIR / "fig1_marginal_radius.png", bbox_inches="tight")
    plt.close(fig)


def fig2_uniformity(dlogr: dict, stats: dict) -> None:
    show = [None, 1.0, 0.3, 0.1]
    medians = {label_sigma_r(sr) if sr is None else f"{sr:g}":
               float(np.median(dlogr[sr])) for sr in show}
    stats["dlogr_median"] = {("None" if sr is None else f"{sr:g}"):
                             float(np.median(dlogr[sr])) for sr in show}
    stats["dlogr_n_pairs"] = {("None" if sr is None else f"{sr:g}"):
                              int(dlogr[sr].size) for sr in show}

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for i, sr in enumerate(show):
        x, y = ecdf_curve(dlogr[sr])
        med = np.median(dlogr[sr])
        ax.plot(x, y, color=SERIES[i],
                label=f"sigma_r = {label_sigma_r(sr)}   median {med:.3f}")
    ax.set_xlabel("|log10(R_out / R_in)| of adjacent small-planet pairs (dex)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Intra-system size uniformity tightens as sigma_r drops")
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 1)
    ax.legend(loc="lower right", title=f"{N_FIG12:,} systems per curve")
    fig.savefig(OUT_DIR / "fig2_size_uniformity.png", bbox_inches="tight")
    plt.close(fig)
    _ = medians


def fig3_inclinations(offsets: dict, cos_i: np.ndarray, stats: dict) -> None:
    stats["inc_proj_std"] = {f"{k:g}": float(np.std(v))
                             for k, v in offsets.items()}
    stats["inc_proj_n"] = {f"{k:g}": int(v.size) for k, v in offsets.items()}
    u = np.sort(cos_i)
    d_uni = float(np.max(np.abs((np.arange(1, u.size + 1) / u.size) - u)))
    stats["isotropic_cos_inc_ks_uniform"] = d_uni

    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(9.6, 4.4), gridspec_kw={"width_ratios": [1.9, 1.0],
                                               "wspace": 0.28})
    for i, (si, vals) in enumerate(offsets.items()):
        lim = 60.0
        bins = np.linspace(-lim, lim, 241)
        hist, edges = np.histogram(vals, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        keep = hist > 0
        ax.plot(centers[keep], hist[keep], color=SERIES[i],
                label=f"{si:g} deg   {np.std(vals):.2f} deg")
    ax.set_yscale("log")
    ax.set_xlim(-60, 60)
    ax.set_xlabel("Projected offset  inc - sys_inc (deg)")
    ax.set_ylabel("Probability density (per deg)")
    ax.set_title("Projected inclination scatter tracks sigma_i")
    ax.legend(loc="upper right", fontsize=8.5, title="sigma_i    sample std",
              title_fontsize=8.5, alignment="left")

    axb.hist(cos_i, bins=25, range=(0.0, 1.0), density=True,
             color=SERIES[4], edgecolor=SURFACE, linewidth=0.8)
    axb.axhline(1.0, color=INK_2, lw=1.2, dashes=(4, 2))
    axb.text(0.03, 1.13, "Uniform(0,1)", color=INK_2, fontsize=8.5)
    axb.set_ylim(0, 1.25)
    axb.set_xlabel("cos(inc)")
    axb.set_ylabel("Density")
    axb.set_title(f"isotropic=True\nKS vs uniform D = {d_uni:.4f}", fontsize=10)
    fig.savefig(OUT_DIR / "fig3_inclination_scatter.png", bbox_inches="tight")
    plt.close(fig)


def fig4_dichotomy(spm: dict, multi_frac: dict, stats: dict) -> None:
    stats["geo_spm"] = {(k if isinstance(k, str) else f"{k:g}"): float(v)
                        for k, v in spm.items()}
    stats["geo_multi_frac"] = {(k if isinstance(k, str) else f"{k:g}"):
                               float(v) for k, v in multi_frac.items()}

    xs = [k for k in spm if not isinstance(k, str)]
    fig, (ax, axf) = plt.subplots(
        2, 1, figsize=(7.2, 6.6), sharex=True,
        gridspec_kw={"height_ratios": [1.6, 1.0], "hspace": 0.14})

    ax.axhline(spm["isotropic"], color=MUTED, lw=1.4, dashes=(5, 3))
    ax.text(0.02, spm["isotropic"] * 1.02, f"isotropic  ({spm['isotropic']:.2f})",
            color=INK_2, fontsize=9, va="bottom")
    ax.plot(xs, [spm[k] for k in xs], color=SERIES[0], marker="o",
            markersize=6)
    ax.set_ylabel("Geometric singles-per-multi")
    ax.set_title("Transit-geometry dichotomy observable vs sigma_i\n"
                 f"({N_FIG4:,} systems per point; no detection model)")

    axf.axhline(multi_frac["isotropic"], color=MUTED, lw=1.4, dashes=(5, 3))
    axf.text(0.02, multi_frac["isotropic"] * 1.03,
             f"isotropic  ({multi_frac['isotropic']:.3f})",
             color=INK_2, fontsize=9, va="bottom")
    axf.plot(xs, [multi_frac[k] for k in xs], color=SERIES[1], marker="o",
             markersize=6)
    axf.set_ylabel("Multi fraction of transiting systems")
    axf.set_xlabel("Mutual-inclination Rayleigh scale sigma_i (deg)")

    for a in (ax, axf):
        a.set_xscale("symlog", linthresh=1.0)
        a.set_xticks([0, 1, 2, 3, 5, 10, 30])
        a.set_xticklabels(["0", "1", "2", "3", "5", "10", "30"])
        a.set_xlim(-0.15, 40)
    fig.savefig(OUT_DIR / "fig4_dichotomy_spm.png", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stats: dict = {"base_seed": BASE_SEED,
                   "n_systems": {"fig1_fig2_per_cell": N_FIG12,
                                 "fig3_per_cell": N_FIG3,
                                 "fig4_per_cell": N_FIG4}}
    t0 = time.time()

    print("Figures 1+2: radius cells", flush=True)
    radii, dlogr = gen_radius_cells()
    fig1_marginal(radii, stats)
    fig2_uniformity(dlogr, stats)

    print("Figure 3: inclination cells", flush=True)
    offsets, cos_i = gen_inclination_cells()
    fig3_inclinations(offsets, cos_i, stats)

    print("Figure 4: dichotomy cells", flush=True)
    spm, multi_frac = gen_dichotomy_cells()
    fig4_dichotomy(spm, multi_frac, stats)

    (OUT_DIR / "stats.json").write_text(json.dumps(stats, indent=2) + "\n")
    print(f"Done in {time.time() - t0:.0f} s -> {OUT_DIR}", flush=True)


if __name__ == "__main__":
    main()
