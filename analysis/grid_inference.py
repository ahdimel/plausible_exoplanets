"""Inference on the fine (sigma_r, sigma_i) grid: likelihood surface,
credible contours, and the sigma_i marginal with sigma_r free vs fixed.

Why a likelihood and not the raw distance: the two pre-registered
distance components live on incommensurable scales (at the optimum the
size term is ~15x the multiplicity term), so contours on their unweighted
sum are arbitrary. Both components map onto likelihoods:

  - multiplicity: the multinomial-KL of the N_k shape equals (up to a
    theta-independent constant) the exact multinomial log-likelihood of
    the real N_k under pooled synthetic frequencies (Laplace 0.5), so
    logL_mult = sum_k n_k_real * ln p_k(theta) is the same statistic on
    a likelihood scale, not a new metric.
  - size: the KS statistic has no likelihood interpretation; instead the
    real |dlogR| sample is scored under a Gaussian KDE (reflected at 0)
    of the pooled synthetic |dlogR|. The bandwidth is FIXED across cells
    (Silverman's rule on the real sample) so logL_size is comparable
    across theta.

Monte-Carlo noise: the per-cell seed noise on logL exceeds the Wilks
1-sigma level (Delta 2lnL = 2.30), so contours are NOT read off raw
cells. Instead a cubic polynomial in (ln sigma_r, ln sigma_i) is fit to
the valley cells (Delta 2lnL < 200) by least squares; the fit residual
RMS is checked against the seed-noise estimate (from half-vs-half seed
splits), and contours/marginals come from the smoothed surface. The
sigma_r=None row (fixed-uncorrelated comparison model) gets the same
treatment in 1D. Marginals use a log-uniform prior.

Run after grid_sweep.py --fine (and optionally grid_topup.py):
  .venv/bin/python analysis/grid_inference.py
Outputs: results/grid_inference.json, results/figures/dichotomy_*.png
"""
import json
import math
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from exoverse.kepler_data import load_koi_systems, load_stellar_targets
from exoverse.dichotomy import observables_from_koi

ROOT = Path(__file__).resolve().parent.parent
GRID = ROOT / "results" / "grid"
FIGS = ROOT / "results" / "figures"
KCATS = ["1", "2", "3", "4", "5", "6+"]
LVL_1SIG, LVL_2SIG = 2.30, 6.18   # Wilks Delta(2 lnL), 2 dof
VALLEY = 200.0                    # Delta(2 lnL) cut for smoothing region


def silverman_bw(x: np.ndarray) -> float:
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    return 0.9 * min(x.std(ddof=1), iqr / 1.34) * len(x) ** -0.2


def kde_loglike(real: np.ndarray, syn: np.ndarray, bw: float) -> float:
    """Sum of log reflected-Gaussian-KDE densities of `syn` at `real`."""
    z = (real[:, None] - syn[None, :]) / bw
    zr = (real[:, None] + syn[None, :]) / bw          # reflection at 0
    dens = (np.exp(-0.5 * z * z) + np.exp(-0.5 * zr * zr)).sum(axis=1)
    dens /= len(syn) * bw * math.sqrt(2 * math.pi)
    return float(np.log(np.maximum(dens, 1e-300)).sum())


def mult_loglike(real_nk: dict, syn_counts: np.ndarray) -> float:
    p = (syn_counts + 0.5) / (syn_counts.sum() + 0.5 * len(KCATS))
    n = np.array([real_nk.get("6+" if k == "6+" else int(k), 0)
                  for k in KCATS], dtype=float)
    return float((n * np.log(p)).sum())


def pooled_loglike(cells: list, real_nk: dict, real_dlogr: np.ndarray,
                   bw: float) -> dict:
    counts = np.zeros(len(KCATS))
    dlogr = []
    for d in cells:
        for i, k in enumerate(KCATS):
            counts[i] += d["n_k"].get(k, 0)
        dlogr.extend(d["dlogr"])
    lm = mult_loglike(real_nk, counts)
    ls = kde_loglike(real_dlogr, np.array(dlogr), bw)
    return {"mult": lm, "size": ls, "total": lm + ls,
            "counts": counts.tolist()}


def poly_design(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Cubic 2D polynomial design matrix in (x, y)."""
    return np.column_stack([np.ones_like(x), x, y, x * x, x * y, y * y,
                            x ** 3, x * x * y, x * y * y, y ** 3])


def upper_quantiles(x_axis, post, qs=(0.05, 0.32, 0.68, 0.95)):
    cdf = np.cumsum(post) / np.sum(post)
    lx = np.log(x_axis)
    return {f"q{int(q*100)}": float(np.exp(np.interp(q, cdf, lx)))
            for q in qs}


def main() -> None:
    axes = json.loads((GRID / "axes.json").read_text())
    sr_axis = [None if s == "None" else float(s) for s in axes["sigma_r"]]
    si_axis = [float(s) for s in axes["sigma_i"]]

    targets = load_stellar_targets("data")
    real = observables_from_koi(load_koi_systems("data"),
                                n_targets=len(targets))
    real_dlogr = np.array(real.dlogr)
    bw = silverman_bw(real_dlogr)
    print(f"real: {real.n_systems} systems, {len(real_dlogr)} adjacent "
          f"pairs, KDE bandwidth {bw:.4f}")

    def load_cells(sr, si):
        tag_sr = "None" if sr is None else sr
        files = sorted(GRID.glob(f"cell_sr{tag_sr}_si{si}_m*.json"),
                       key=lambda p: int(p.stem.rsplit("_m", 1)[1]))
        return [json.loads(p.read_text()) for p in files]

    result = {}
    for sr in sr_axis:
        for si in si_axis:
            cells = load_cells(sr, si)
            r = pooled_loglike(cells, real.n_k, real_dlogr, bw)
            half = len(cells) // 2
            for h, cs in enumerate((cells[:half], cells[half:])):
                r[f"total_h{h}"] = pooled_loglike(
                    cs, real.n_k, real_dlogr, bw)["total"]
            r["n_seeds"] = len(cells)
            multis = sum(sum(d["n_k"].get(k, 0) for k in KCATS[1:])
                         for d in cells)
            singles = sum(d["n_k"].get("1", 0) for d in cells)
            r["spm"] = singles / multis if multis else float("inf")
            result[f"{'None' if sr is None else sr}|{si}"] = r

    def cell(sr, si):
        return result[f"{'None' if sr is None else sr}|{si}"]

    sr_free = [s for s in sr_axis if s is not None]
    raw = np.array([[cell(sr, si)["total"] for si in si_axis]
                    for sr in sr_free])
    llmax_raw = raw.max()

    # --- seed-noise estimate over the valley ------------------------------
    valley_keys, diffs = [], []
    for i, sr in enumerate(sr_free):
        for j, si in enumerate(si_axis):
            if 2 * (llmax_raw - raw[i, j]) < VALLEY:
                valley_keys.append((i, j))
                r = cell(sr, si)
                diffs.append(r["total_h0"] - r["total_h1"])
    # std(h0-h1) = 2 * sigma_pooled (each half holds M/2 seeds)
    noise = float(np.std(diffs) / 2)
    print(f"valley: {len(valley_keys)} cells, seed noise on pooled logL "
          f"~{noise:.2f} (Wilks 1-sigma needs Delta lnL = 1.15)")

    # --- smoothed surface: cubic in (ln sr, ln si) over the valley --------
    vx = np.log([sr_free[i] for i, _ in valley_keys])
    vy = np.log([si_axis[j] for _, j in valley_keys])
    vz = np.array([raw[i, j] for i, j in valley_keys])
    beta, *_ = np.linalg.lstsq(poly_design(vx, vy), vz, rcond=None)
    resid = vz - poly_design(vx, vy) @ beta
    resid_rms = float(np.std(resid))
    print(f"cubic fit residual RMS {resid_rms:.2f} vs seed noise "
          f"{noise:.2f} -> smoothing is "
          f"{'consistent with pure seed noise' if resid_rms < 1.5 * noise else 'ABSORBING REAL STRUCTURE (check!)'}")

    # dense mesh over the valley bounding box
    gx = np.linspace(vx.min(), vx.max(), 220)
    gy = np.linspace(vy.min(), vy.max(), 220)
    GX, GY = np.meshgrid(gx, gy, indexing="ij")
    smooth = (poly_design(GX.ravel(), GY.ravel()) @ beta).reshape(GX.shape)
    i_b, j_b = np.unravel_index(np.argmax(smooth), smooth.shape)
    llmax = float(smooth[i_b, j_b])
    best_sr, best_si = float(np.exp(gx[i_b])), float(np.exp(gy[j_b]))
    d2s = 2 * (llmax - smooth)
    edge_max = i_b in (0, len(gx) - 1) or j_b in (0, len(gy) - 1)

    # degeneracy: local quadratic curvature at the smooth maximum
    dx, dy = gx[1] - gx[0], gy[1] - gy[0]
    ic = min(max(i_b, 2), len(gx) - 3)
    jc = min(max(j_b, 2), len(gy) - 3)
    hxx = (smooth[ic + 1, jc] - 2 * smooth[ic, jc] + smooth[ic - 1, jc]) / dx**2
    hyy = (smooth[ic, jc + 1] - 2 * smooth[ic, jc] + smooth[ic, jc - 1]) / dy**2
    hxy = (smooth[ic + 1, jc + 1] - smooth[ic + 1, jc - 1]
           - smooth[ic - 1, jc + 1] + smooth[ic - 1, jc - 1]) / (4 * dx * dy)
    corr = float(hxy / math.sqrt(hxx * hyy)) if hxx * hyy > 0 else float("nan")

    # degeneracy, robust to an edge maximum: posterior correlation of
    # (ln sigma_r, ln sigma_i) over the smoothed surface
    wn = np.exp(smooth - llmax)
    wn /= wn.sum()
    ex, ey = (wn * GX).sum(), (wn * GY).sum()
    cxy = (wn * (GX - ex) * (GY - ey)).sum()
    post_corr = float(cxy / math.sqrt(
        (wn * (GX - ex) ** 2).sum() * (wn * (GY - ey) ** 2).sum()))

    # --- marginals (log-uniform prior; mesh is uniform in log) ------------
    w = np.exp(smooth - llmax)
    post_si = w.sum(axis=0)
    post_sr = w.sum(axis=1)
    si_mesh, sr_mesh = np.exp(gy), np.exp(gx)

    # fixed-uncorrelated model: 1D cubic in ln si over its own valley
    ll_none = np.array([cell(None, si)["total"] for si in si_axis])
    m_none = 2 * (ll_none.max() - ll_none) < VALLEY
    xn = np.log(np.array(si_axis)[m_none])
    Xn = np.column_stack([np.ones_like(xn), xn, xn**2, xn**3])
    bn, *_ = np.linalg.lstsq(Xn, ll_none[m_none], rcond=None)
    gyn = np.linspace(xn.min(), xn.max(), 220)
    ll_none_s = np.column_stack(
        [np.ones_like(gyn), gyn, gyn**2, gyn**3]) @ bn
    post_si_fixed = np.exp(ll_none_s - ll_none_s.max())

    marg = {
        "si_free": upper_quantiles(si_mesh, post_si),
        "si_fixed_uncorrelated": upper_quantiles(np.exp(gyn),
                                                 post_si_fixed),
        "sr": upper_quantiles(sr_mesh, post_sr),
    }
    dll_corr = float(llmax - ll_none_s.max())

    # --- component decomposition along the two axes ------------------------
    j_col = int(np.argmin(np.abs(np.log(si_axis) - math.log(best_si))))
    i_row = int(np.argmin(np.abs(np.log(sr_free) - math.log(best_sr))))
    decomp = {
        "along_sigma_r_at_si": si_axis[j_col],
        "mult_vs_sr": [round(cell(sr, si_axis[j_col])["mult"], 1)
                       for sr in sr_free],
        "size_vs_sr": [round(cell(sr, si_axis[j_col])["size"], 1)
                       for sr in sr_free],
        "along_sigma_i_at_sr": sr_free[i_row],
        "mult_vs_si": [round(cell(sr_free[i_row], si)["mult"], 1)
                       for si in si_axis],
        "size_vs_si": [round(cell(sr_free[i_row], si)["size"], 1)
                       for si in si_axis],
    }

    # --- multiplicity goodness of fit at the best raw cell ----------------
    i_r, j_r = np.unravel_index(np.argmax(raw), raw.shape)
    best_cell = cell(sr_free[i_r], si_axis[j_r])
    counts = np.array(best_cell["counts"])
    p = (counts + 0.5) / (counts.sum() + 0.5 * len(KCATS))
    n_obs = np.array([real.n_k.get("6+" if k == "6+" else int(k), 0)
                      for k in KCATS], dtype=float)
    exp_nk = p * n_obs.sum()
    G2 = float(2 * np.sum(np.where(n_obs > 0,
                                   n_obs * np.log(n_obs / exp_nk), 0.0)))
    gof = {"raw_best": [sr_free[i_r], si_axis[j_r]],
           "real_n_k": n_obs.tolist(),
           "model_n_k_shape_matched": [round(float(v), 1) for v in exp_nk],
           "deviance_G2": G2, "dof": len(KCATS) - 1 - 2,
           "spm_model": best_cell["spm"], "spm_real": 3.399}

    out = {
        "axes": {"sigma_r": axes["sigma_r"], "sigma_i": si_axis},
        "kde_bandwidth": bw, "valley_delta2lnl": VALLEY,
        "seed_noise_logl": noise, "fit_residual_rms": resid_rms,
        "logL": result,
        "smooth_best": {"sigma_r": best_sr, "sigma_i": best_si,
                        "logL": llmax, "at_mesh_edge": edge_max},
        "delta_lnL_corr_vs_uncorr": dll_corr,
        "log_curvature_correlation": corr,
        "posterior_correlation": post_corr,
        "marginals": marg,
        "decomposition": decomp,
        "multiplicity_gof": gof,
        "convergence": {
            "raw_best": [sr_free[i_r], si_axis[j_r]],
            "smooth_best": [best_sr, best_si],
            "n_valley_cells": len(valley_keys)},
        "posterior": {
            "si_mesh": si_mesh.tolist(), "post_si_free": post_si.tolist(),
            "sr_mesh": sr_mesh.tolist(), "post_sr": post_sr.tolist(),
            "si_mesh_fixed": np.exp(gyn).tolist(),
            "post_si_fixed": post_si_fixed.tolist()},
        "wilks_levels": {"1sig": LVL_1SIG, "2sig": LVL_2SIG},
    }
    (ROOT / "results" / "grid_inference.json").write_text(
        json.dumps(out, indent=1))

    # --- report ------------------------------------------------------------
    print(f"\nsmooth best: sigma_r={best_sr:.2f}, sigma_i={best_si:.2f} deg"
          f"{' (AT MESH EDGE)' if edge_max else ''}; raw best "
          f"({sr_free[i_r]}, {si_axis[j_r]}), spm={best_cell['spm']:.2f}")
    print(f"Delta lnL (correlated - uncorrelated best) = {dll_corr:.1f}")
    print(f"sigma_i: free sigma_r 68/95% UL "
          f"{marg['si_free']['q68']:.2f}/{marg['si_free']['q95']:.2f} deg | "
          f"fixed uncorrelated "
          f"{marg['si_fixed_uncorrelated']['q68']:.2f}/"
          f"{marg['si_fixed_uncorrelated']['q95']:.2f} deg")
    print(f"sigma_r: 68% [{marg['sr']['q32']:.2f}, {marg['sr']['q68']:.2f}]"
          f", 95% [{marg['sr']['q5']:.2f}, {marg['sr']['q95']:.2f}]")
    print(f"degeneracy: curvature correlation at max {corr:.2f}; "
          f"posterior correlation (ln sr, ln si) {post_corr:.2f}")
    print(f"multiplicity GoF at raw best: G2={G2:.1f} on {gof['dof']} dof; "
          f"model N_k {[float(v) for v in gof['model_n_k_shape_matched']]} vs real "
          f"{n_obs.astype(int).tolist()}")
    print("\ndecomposition (logL relative to column max):")
    mv = np.array(decomp["mult_vs_sr"]); sv = np.array(decomp["size_vs_sr"])
    print(f"  vs sigma_r at si={decomp['along_sigma_r_at_si']}:")
    print("    sr:   " + " ".join(f"{s:7.2f}" for s in sr_free))
    print("    mult: " + " ".join(f"{v:7.1f}" for v in mv - mv.max()))
    print("    size: " + " ".join(f"{v:7.1f}" for v in sv - sv.max()))
    mv = np.array(decomp["mult_vs_si"]); sv = np.array(decomp["size_vs_si"])
    print(f"  vs sigma_i at sr={decomp['along_sigma_i_at_sr']}:")
    print("    si:   " + " ".join(f"{s:7.2f}" for s in si_axis))
    print("    mult: " + " ".join(f"{v:7.1f}" for v in mv - mv.max()))
    print("    size: " + " ".join(f"{v:7.1f}" for v in sv - sv.max()))

    # --- figures -------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGS.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    pm = ax.pcolormesh(np.exp(gy), np.exp(gx), np.minimum(d2s, 60),
                       cmap="viridis_r", shading="auto")
    cs = ax.contour(np.exp(gy), np.exp(gx), d2s,
                    levels=[LVL_1SIG, LVL_2SIG],
                    colors=["w", "w"], linestyles=["-", "--"])
    ax.clabel(cs, fmt={LVL_1SIG: "68%", LVL_2SIG: "95%"}, fontsize=9)
    ax.plot(best_si, best_sr, "w*", ms=14, mec="k")
    for i, j in valley_keys:
        ax.plot(si_axis[j], sr_free[i], "k.", ms=2, alpha=0.4)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\sigma_i$ (deg, Rayleigh mutual inclination)")
    ax.set_ylabel(r"$\sigma_R$ (radius correlation; low = strong)")
    ax.set_title(r"DR25 joint fit: smoothed $\Delta(2\ln L)$, valley region")
    fig.colorbar(pm, label=r"$\Delta(2\ln L)$ (capped at 60)")
    fig.tight_layout()
    fig.savefig(FIGS / "dichotomy_surface.png", dpi=150)

    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    ax.plot(si_mesh, post_si / post_si.max(), "-",
            label=r"$\sigma_R$ free (marginalized)")
    ax.plot(np.exp(gyn), post_si_fixed / post_si_fixed.max(), "--",
            label=r"$\sigma_R$ fixed: uncorrelated radii")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\sigma_i$ (deg)")
    ax.set_ylabel("relative posterior (log-uniform prior)")
    ax.legend()
    ax.set_title(r"$\sigma_i$ marginal: joint vs radius-independent model")
    fig.tight_layout()
    fig.savefig(FIGS / "dichotomy_sigma_i_marginal.png", dpi=150)

    best_cells = load_cells(sr_free[i_r], si_axis[j_r])
    syn_shape = exp_nk
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.2))
    x = np.arange(len(KCATS))
    a1.bar(x - 0.2, n_obs, 0.4, label="DR25 KOI (fiducial cuts)")
    a1.bar(x + 0.2, syn_shape, 0.4, label="best fit (shape-matched)")
    a1.set_xticks(x, KCATS); a1.set_yscale("log")
    a1.set_xlabel("detected multiplicity k"); a1.set_ylabel("systems")
    a1.legend(fontsize=8)
    syn_dlogr = np.sort(np.concatenate([d["dlogr"] for d in best_cells]))
    rd = np.sort(real_dlogr)
    a2.plot(rd, np.linspace(0, 1, len(rd)), label="DR25 KOI")
    a2.plot(syn_dlogr, np.linspace(0, 1, len(syn_dlogr)),
            label="best fit (pooled seeds)")
    a2.set_xlabel(r"$|\Delta \log R|$ adjacent pairs"); a2.set_ylabel("CDF")
    a2.set_xlim(0, 1.2); a2.legend(fontsize=8)
    fig.suptitle("Fit quality at best cell")
    fig.tight_layout()
    fig.savefig(FIGS / "dichotomy_fit_quality.png", dpi=150)
    print(f"\nfigures written to {FIGS}")


if __name__ == "__main__":
    main()
