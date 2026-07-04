"""Publication-quality figures for the Kepler-dichotomy decomposition study.

Self-contained: reads only the cached JSONs (results/grid_inference.json,
results/mixture_inference.json), the raw per-universe grid cells
(results/grid/, results/grid_mixture/), and the real DR25 observables
(recomputed from data/ via exoverse). It re-derives the smoothed likelihood
surface with the *same* cubic-in-(ln sr, ln si) least-squares fit used by
analysis/grid_inference.py, and STOPS if the recomputed smooth best fit
does not match results/grid_inference.json "smooth_best" to 2 decimals.

Renders four figures into results/figures/paper/ as PNG (dpi=200) + PDF:
  fig1_surface   — smoothed Delta(2 lnL) surface, 68/95% contours, best-fit
                   star, simulated grid cells.
  fig2_marginals — sigma_i posteriors: sigma_R free vs fixed-uncorrelated,
                   with 68/95% upper-limit annotations.
  fig3_fit       — left: detected multiplicity N_k (real vs single vs
                   mixture); right: |dlogR| CDFs (real vs best vs uncorrelated).
  fig4_mixture   — Delta lnL and multiplicity deviance G2 vs f_hot.

Run: .venv/bin/python analysis/paper_figures.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as patheffects  # noqa: E402
from matplotlib.ticker import (  # noqa: E402
    FixedFormatter, FixedLocator, NullFormatter,
)

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from exoverse.kepler_data import load_koi_systems, load_stellar_targets  # noqa: E402
from exoverse.dichotomy import observables_from_koi  # noqa: E402

GRID = ROOT / "results" / "grid"
MIX = ROOT / "results" / "grid_mixture"
OUT = ROOT / "results" / "figures" / "paper"
KCATS = ["1", "2", "3", "4", "5", "6+"]
KLABELS = ["1", "2", "3", "4", "5", "6+"]
LVL_1SIG, LVL_2SIG = 2.30, 6.18   # Wilks Delta(2 lnL), 2 dof
VALLEY = 200.0

# --- style (mirrors analysis/phase1_validation_plots.py: dataviz palette) ---
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
    "font.family": "sans-serif",
    "font.size": 10,
    "text.color": INK,
    "axes.edgecolor": BASELINE,
    "axes.labelcolor": INK_2,
    "axes.titlecolor": INK,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "legend.frameon": False,
    "legend.fontsize": 9,
    "lines.linewidth": 2.0,
})


# ---------------------------------------------------------------------------
# Likelihood machinery (identical math to analysis/grid_inference.py)
# ---------------------------------------------------------------------------
def silverman_bw(x: np.ndarray) -> float:
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    return 0.9 * min(x.std(ddof=1), iqr / 1.34) * len(x) ** -0.2


def kde_loglike(real: np.ndarray, syn: np.ndarray, bw: float) -> float:
    z = (real[:, None] - syn[None, :]) / bw
    zr = (real[:, None] + syn[None, :]) / bw
    dens = (np.exp(-0.5 * z * z) + np.exp(-0.5 * zr * zr)).sum(axis=1)
    dens /= len(syn) * bw * math.sqrt(2 * math.pi)
    return float(np.log(np.maximum(dens, 1e-300)).sum())


def mult_loglike(real_nk: dict, counts: np.ndarray) -> float:
    p = (counts + 0.5) / (counts.sum() + 0.5 * len(KCATS))
    n = np.array([real_nk.get("6+" if k == "6+" else int(k), 0)
                  for k in KCATS], dtype=float)
    return float((n * np.log(p)).sum())


def pool_cells(cells: list):
    counts = np.zeros(len(KCATS))
    dlogr = []
    for d in cells:
        for i, k in enumerate(KCATS):
            counts[i] += d["n_k"].get(k, 0)
        dlogr.extend(d["dlogr"])
    return counts, np.array(dlogr)


def poly_design(x, y):
    return np.column_stack([np.ones_like(x), x, y, x * x, x * y, y * y,
                            x ** 3, x * x * y, x * y * y, y ** 3])


def load_cells(folder: Path, pattern: str) -> list:
    files = sorted(folder.glob(f"{pattern}_m*.json"),
                   key=lambda p: int(p.stem.rsplit("_m", 1)[1]))
    return [json.loads(p.read_text()) for p in files]


# ---------------------------------------------------------------------------
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    gi = json.loads((ROOT / "results" / "grid_inference.json").read_text())
    mi = json.loads((ROOT / "results" / "mixture_inference.json").read_text())

    axes = json.loads((GRID / "axes.json").read_text())
    sr_axis = [None if s == "None" else float(s) for s in axes["sigma_r"]]
    si_axis = [float(s) for s in axes["sigma_i"]]
    sr_free = [s for s in sr_axis if s is not None]

    real = observables_from_koi(load_koi_systems("data"),
                                n_targets=len(load_stellar_targets("data")))
    real_dlogr = np.array(real.dlogr)
    bw = silverman_bw(real_dlogr)
    real_total = int(real.n_systems)
    n_obs = np.array([real.n_k.get("6+" if k == "6+" else int(k), 0)
                      for k in KCATS], dtype=float)

    def cell_total(sr, si):
        tag = "None" if sr is None else sr
        counts, dl = pool_cells(load_cells(GRID, f"cell_sr{tag}_si{si}"))
        return mult_loglike(real.n_k, counts) + kde_loglike(real_dlogr, dl, bw)

    # --- recompute smoothed surface (grid_inference.py, exactly) -----------
    raw = np.array([[cell_total(sr, si) for si in si_axis]
                    for sr in sr_free])
    llmax_raw = raw.max()
    valley_keys = [(i, j) for i in range(len(sr_free))
                   for j in range(len(si_axis))
                   if 2 * (llmax_raw - raw[i, j]) < VALLEY]
    vx = np.log([sr_free[i] for i, _ in valley_keys])
    vy = np.log([si_axis[j] for _, j in valley_keys])
    vz = np.array([raw[i, j] for i, j in valley_keys])
    beta, *_ = np.linalg.lstsq(poly_design(vx, vy), vz, rcond=None)
    gx = np.linspace(vx.min(), vx.max(), 220)
    gy = np.linspace(vy.min(), vy.max(), 220)
    GX, GY = np.meshgrid(gx, gy, indexing="ij")
    smooth = (poly_design(GX.ravel(), GY.ravel()) @ beta).reshape(GX.shape)
    i_b, j_b = np.unravel_index(np.argmax(smooth), smooth.shape)
    best_sr, best_si = float(np.exp(gx[i_b])), float(np.exp(gy[j_b]))
    d2s = 2 * (smooth.max() - smooth)

    jb = gi["smooth_best"]
    if (round(best_sr, 2) != round(jb["sigma_r"], 2)
            or round(best_si, 2) != round(jb["sigma_i"], 2)):
        raise SystemExit(
            "STOP: recomputed smooth best "
            f"(sr={best_sr:.4f}, si={best_si:.4f}) does not match "
            f"grid_inference.json smooth_best "
            f"(sr={jb['sigma_r']:.4f}, si={jb['sigma_i']:.4f}). "
            "Not shipping fig1.")
    print(f"[fig1] smooth best sr={best_sr:.4f} si={best_si:.4f} "
          f"matches JSON ({jb['sigma_r']:.4f}, {jb['sigma_i']:.4f}); "
          f"{len(valley_keys)} valley cells")

    def savefig(fig, name):
        fig.savefig(OUT / f"{name}.png", dpi=200)
        fig.savefig(OUT / f"{name}.pdf")
        plt.close(fig)

    # =====================================================================
    # fig1: smoothed Delta(2 lnL) surface
    # =====================================================================
    si_mesh_grid, sr_mesh_grid = np.exp(gy), np.exp(gx)
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    pm = ax.pcolormesh(si_mesh_grid, sr_mesh_grid, np.minimum(d2s, 60),
                       cmap="viridis_r", shading="auto")
    cs = ax.contour(si_mesh_grid, sr_mesh_grid, d2s,
                    levels=[LVL_1SIG, LVL_2SIG],
                    colors="white", linewidths=[2.0, 1.6],
                    linestyles=["-", "--"])
    ax.clabel(cs, fmt={LVL_1SIG: "68%", LVL_2SIG: "95%"}, fontsize=9,
              inline=True)
    for i, j in valley_keys:
        ax.plot(si_axis[j], sr_free[i], ".", color="0.15", ms=3, alpha=0.45)
    ax.plot(best_si, best_sr, "*", ms=18, mfc="white", mec="k", mew=1.2,
            zorder=5)
    ax.annotate("best fit", xy=(best_si, best_sr),
                xytext=(10, 9), textcoords="offset points",
                color="k", fontsize=10, fontweight="bold", zorder=6,
                path_effects=[patheffects.withStroke(linewidth=2.5,
                                                     foreground="white")])
    ax.set_xscale("log")
    ax.set_yscale("log")
    xticks = [t for t in (0.2, 0.3, 0.5, 1.0, 2.0, 5.0)
              if si_mesh_grid.min() <= t <= si_mesh_grid.max()]
    yticks = [t for t in (1.0, 1.5, 2.0, 3.0)
              if sr_mesh_grid.min() <= t <= sr_mesh_grid.max()]
    ax.xaxis.set_major_locator(FixedLocator(xticks))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in xticks]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.yaxis.set_major_formatter(FixedFormatter([f"{t:g}" for t in yticks]))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$\sigma_i$  (deg, Rayleigh mutual inclination)")
    ax.set_ylabel(r"$\sigma_R$  (radius-correlation scale; smaller = stronger)")
    ax.set_title(r"DR25 joint fit: smoothed $\Delta(2\ln L)$ over the valley")
    cb = fig.colorbar(pm, ax=ax, label=r"$\Delta(2\ln L)$  (capped at 60)")
    cb.outline.set_edgecolor(BASELINE)
    fig.tight_layout()
    savefig(fig, "fig1_surface")

    # =====================================================================
    # fig2: sigma_i marginals (from the posterior block of the JSON)
    # =====================================================================
    post = gi["posterior"]
    si_mesh = np.array(post["si_mesh"])
    post_si_free = np.array(post["post_si_free"])
    si_mesh_fixed = np.array(post["si_mesh_fixed"])
    post_si_fixed = np.array(post["post_si_fixed"])
    mfree = gi["marginals"]["si_free"]
    mfix = gi["marginals"]["si_fixed_uncorrelated"]

    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    ax.plot(si_mesh, post_si_free / post_si_free.max(), "-",
            color=SERIES[0], label=r"$\sigma_R$ free (marginalized)")
    ax.plot(si_mesh_fixed, post_si_fixed / post_si_fixed.max(), "--",
            color=SERIES[5], label=r"$\sigma_R$ fixed: uncorrelated radii")
    for q, style, col, lab in [
            (mfree["q68"], (0, (5, 2)), SERIES[0], "68%"),
            (mfree["q95"], (0, (1, 1.5)), SERIES[0], "95%"),
            (mfix["q68"], (0, (5, 2)), SERIES[5], "68%"),
            (mfix["q95"], (0, (1, 1.5)), SERIES[5], "95%")]:
        ax.axvline(q, color=col, ls=style, lw=1.3, alpha=0.8, ymax=0.88)
    ax.set_ylim(0, 1.05)
    # annotate the upper limits in the empty centre-right region
    ax.text(0.62, 0.62, f"$\\sigma_R$ free:\n"
            f"UL$_{{68}}$ = {mfree['q68']:.2f} deg\n"
            f"UL$_{{95}}$ = {mfree['q95']:.2f} deg",
            transform=ax.transData, color=SERIES[0], fontsize=9,
            va="top", ha="left")
    ax.text(0.62, 0.36, f"$\\sigma_R$ fixed (uncorr.):\n"
            f"UL$_{{68}}$ = {mfix['q68']:.2f} deg\n"
            f"UL$_{{95}}$ = {mfix['q95']:.2f} deg",
            transform=ax.transData, color=SERIES[5], fontsize=9,
            va="top", ha="left")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sigma_i$  (deg)")
    ax.set_ylabel("relative posterior  (log-uniform prior)")
    ax.set_title(r"$\sigma_i$ marginal posterior: joint vs "
                 r"radius-independent model")
    ax.legend(loc="upper right")
    fig.tight_layout()
    savefig(fig, "fig2_marginals")

    # =====================================================================
    # fig3: multiplicity fit (left) + |dlogR| CDFs (right)
    # =====================================================================
    single_shape = np.array(gi["multiplicity_gof"]["model_n_k_shape_matched"])
    # mixture shape: pool n_k over the 8 mixture cells, normalize to real total
    mix_cells = load_cells(MIX, "cell_si0.2_f0.1_sih10.0")
    mix_counts, _ = pool_cells(mix_cells)
    mix_shape = mix_counts / mix_counts.sum() * real_total

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.2))
    x = np.arange(len(KCATS))
    w = 0.27
    a1.bar(x - w, n_obs, w, color=INK_2, label=f"DR25 KOI ({real_total})")
    a1.bar(x, single_shape, w, color=SERIES[0],
           label="single-pop best fit")
    a1.bar(x + w, mix_shape, w, color=SERIES[2],
           label=r"mixture ($f_{\rm hot}=0.1$, $\sigma_{i,{\rm hot}}=10^\circ$)")
    a1.set_yscale("log")
    a1.set_xticks(x)
    a1.set_xticklabels(KLABELS)
    a1.set_xlabel("detected multiplicity  $k$")
    a1.set_ylabel(r"systems  (shape-normalized to %d)" % real_total)
    a1.set_title("Detected multiplicity")
    a1.legend()
    a1.grid(axis="y", color=GRIDLINE, lw=0.8)
    a1.set_axisbelow(True)

    # right: |dlogR| CDFs
    _, best_dl = pool_cells(load_cells(GRID, "cell_sr1.5244_si0.4183"))
    unc_cells = load_cells(GRID, "cell_srNone_si0.2")
    _, unc_dl = pool_cells(unc_cells)
    print(f"[fig3] best-fit pairs={len(best_dl)} (from "
          f"{len(load_cells(GRID, 'cell_sr1.5244_si0.4183'))} seeds); "
          f"uncorrelated pairs={len(unc_dl)} from {len(unc_cells)} seeds")

    def ecdf(v):
        v = np.sort(v)
        return v, np.arange(1, len(v) + 1) / len(v)

    for v, col, lab in [
            (real_dlogr, INK_2, f"DR25 KOI ({len(real_dlogr)} pairs)"),
            (best_dl, SERIES[0], "best fit (correlated, pooled)"),
            (unc_dl, SERIES[5], r"uncorrelated best $\sigma_i$")]:
        xs, ys = ecdf(v)
        a2.plot(xs, ys, color=col, label=lab)
    a2.set_xlim(0, 1.2)
    a2.set_ylim(0, 1.0)
    a2.set_xlabel(r"$|\Delta \log_{10} R|$  (adjacent small-planet pairs)")
    a2.set_ylabel("cumulative fraction")
    a2.set_title("Size uniformity")
    a2.legend(loc="lower right")
    a2.grid(color=GRIDLINE, lw=0.8)
    a2.set_axisbelow(True)
    fig.tight_layout()
    savefig(fig, "fig3_fit")

    # =====================================================================
    # fig4: mixture Delta lnL and G2 vs f_hot
    # =====================================================================
    single0 = mi["single_table"]["0.2"]      # single-pop best (cold si=0.2)
    ll0 = single0["lnL"]
    g20 = single0["g2"]
    f_axis = mi["axes"]["f_hot"]
    fs = [0.0] + list(f_axis)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.2))
    for sih, col in [(10.0, SERIES[0]), (30.0, SERIES[2])]:
        ll = [ll0] + [mi["table"][f"0.2|{f}|{sih}"]["lnL"] for f in f_axis]
        g2 = [g20] + [mi["table"][f"0.2|{f}|{sih}"]["g2"] for f in f_axis]
        a1.plot(fs, np.array(ll) - ll0, "o-", color=col,
                label=rf"$\sigma_{{i,{{\rm hot}}}}={sih:g}^\circ$")
        a2.plot(fs, g2, "o-", color=col,
                label=rf"$\sigma_{{i,{{\rm hot}}}}={sih:g}^\circ$")
    a1.axhline(0, color=BASELINE, lw=0.8)
    a1.set_xlabel(r"$f_{\rm hot}$")
    a1.set_ylabel(r"$\Delta \ln L$  vs single population")
    a1.set_title("Mixture likelihood gain")
    a1.legend(loc="upper right")
    a1.grid(color=GRIDLINE, lw=0.8)
    a1.set_axisbelow(True)
    a2.axhline(7.81, color=SERIES[5], lw=1.2, ls="--",
               label=r"$\chi^2$ $p=0.05$ (3 dof)")
    a2.set_xlabel(r"$f_{\rm hot}$")
    a2.set_ylabel(r"multiplicity deviance  $G^2$")
    a2.set_title("Multiplicity goodness of fit")
    a2.legend(loc="upper right")
    a2.grid(color=GRIDLINE, lw=0.8)
    a2.set_axisbelow(True)
    fig.suptitle(r"Two-population inclination mixture at "
                 r"$\sigma_R=1.5244$, cold $\sigma_i=0.2^\circ$",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    savefig(fig, "fig4_mixture")

    print(f"\nWrote 8 files to {OUT}")
    for f in sorted(OUT.glob("fig*")):
        print(f"  {f.name}  {f.stat().st_size // 1024} KB")


if __name__ == "__main__":
    main()
