"""Model comparison for the two-population inclination mixture at fixed
sigma_r = 1.5244 (joint best fit): single Rayleigh sigma_i vs
(1 - f_hot) Rayleigh(sigma_i) + f_hot Rayleigh(sigma_i_hot).

Likelihood machinery is imported from grid_inference (same real
observables, same fixed KDE bandwidth). Both models are pooled over the
SAME number of seeds (M=8) so KDE log-likelihoods are comparable: the
single-population baselines use fine-grid cells m=0..7 at sr=1.5244.

Comparison via AIC (mixture adds 2 parameters: f_hot, sigma_i_hot);
the LRT null is on the f_hot=0 boundary, so Wilks is not exact — AIC
plus the multiplicity goodness-of-fit deviance G2 carry the conclusion.
Per-cell seed noise is quoted from a half-split at the best mixture cell.

Run: .venv/bin/python analysis/mixture_inference.py
Outputs: results/mixture_inference.json, results/figures/dichotomy_mixture.png
"""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from grid_inference import (  # noqa: E402
    KCATS, pooled_loglike, silverman_bw,
)
from exoverse.kepler_data import load_koi_systems, load_stellar_targets
from exoverse.dichotomy import observables_from_koi

ROOT = Path(__file__).resolve().parent.parent
MIX = ROOT / "results" / "grid_mixture"
GRID = ROOT / "results" / "grid"
SIGMA_R = 1.5244
M = 8


def load_cells(pattern: str, folder: Path, ms=range(M)) -> list:
    out = []
    for m in ms:
        p = folder / f"cell_{pattern}_m{m}.json"
        out.append(json.loads(p.read_text()))
    return out


def g2_deviance(counts: np.ndarray, n_obs: np.ndarray) -> float:
    p = (counts + 0.5) / (counts.sum() + 0.5 * len(KCATS))
    exp_nk = p * n_obs.sum()
    return float(2 * np.sum(np.where(n_obs > 0,
                                     n_obs * np.log(n_obs / exp_nk), 0.0)))


def main() -> None:
    axes = json.loads((MIX / "axes.json").read_text())
    targets = load_stellar_targets("data")
    real = observables_from_koi(load_koi_systems("data"),
                                n_targets=len(targets))
    real_dlogr = np.array(real.dlogr)
    bw = silverman_bw(real_dlogr)
    n_obs = np.array([real.n_k.get("6+" if k == "6+" else int(k), 0)
                      for k in KCATS], dtype=float)

    def score(cells):
        r = pooled_loglike(cells, real.n_k, real_dlogr, bw)
        counts = np.array(r["counts"])
        multis = counts[1:].sum()
        r["spm"] = float(counts[0] / multis) if multis else float("inf")
        r["g2"] = g2_deviance(counts, n_obs)
        return r

    # single-population baselines at sr=1.5244, m=0..7 (seed-matched)
    single = {}
    for si in axes["sigma_i"]:
        single[si] = score(load_cells(f"sr{SIGMA_R}_si{si}", GRID))
    si_best_single = max(single, key=lambda s: single[s]["total"])
    s0 = single[si_best_single]

    # mixture grid
    mix = {}
    for si in axes["sigma_i"]:
        for f in axes["f_hot"]:
            for sih in axes["sigma_i_hot"]:
                mix[(si, f, sih)] = score(
                    load_cells(f"si{si}_f{f}_sih{sih}", MIX))
    (si_m, f_m, sih_m) = max(mix, key=lambda k: mix[k]["total"])
    m0 = mix[(si_m, f_m, sih_m)]

    # seed noise at the best mixture cell (half split, M/2 each)
    best_cells = load_cells(f"si{si_m}_f{f_m}_sih{sih_m}", MIX)
    h0 = pooled_loglike(best_cells[:M // 2], real.n_k, real_dlogr, bw)
    h1 = pooled_loglike(best_cells[M // 2:], real.n_k, real_dlogr, bw)
    noise = abs(h0["total"] - h1["total"]) / 2

    dlnl = m0["total"] - s0["total"]
    daic = 2 * 2 - 2 * dlnl   # mixture has +2 params; negative favors mixture

    out = {
        "sigma_r_fixed": SIGMA_R, "axes": axes,
        "single_best": {"sigma_i": si_best_single,
                        "lnL": s0["total"], "mult": s0["mult"],
                        "size": s0["size"], "g2": s0["g2"],
                        "spm": s0["spm"]},
        "mixture_best": {"sigma_i": si_m, "f_hot": f_m,
                         "sigma_i_hot": sih_m,
                         "lnL": m0["total"], "mult": m0["mult"],
                         "size": m0["size"], "g2": m0["g2"],
                         "spm": m0["spm"]},
        "delta_lnL": dlnl, "delta_AIC": daic,
        "seed_noise_lnL": noise,
        "real_spm": 3.399,
        "table": {f"{si}|{f}|{sih}": {
            "lnL": v["total"], "g2": v["g2"], "spm": v["spm"]}
            for (si, f, sih), v in mix.items()},
        "single_table": {str(si): {
            "lnL": v["total"], "g2": v["g2"], "spm": v["spm"]}
            for si, v in single.items()},
    }
    (ROOT / "results" / "mixture_inference.json").write_text(
        json.dumps(out, indent=1))

    print(f"single best:  sigma_i={si_best_single}  lnL={s0['total']:.1f}"
          f"  G2={s0['g2']:.1f}  spm={s0['spm']:.2f}")
    print(f"mixture best: sigma_i={si_m} f_hot={f_m} sigma_i_hot={sih_m}"
          f"  lnL={m0['total']:.1f}  G2={m0['g2']:.1f}  spm={m0['spm']:.2f}")
    print(f"Delta lnL = {dlnl:.1f}  Delta AIC = {daic:.1f} "
          f"(negative favors mixture)  seed noise ~{noise:.1f}")
    print(f"real spm = 3.399\n")
    for sih in axes["sigma_i_hot"]:
        print(f"lnL by (rows sigma_i_cold, cols f_hot) at "
              f"sigma_i_hot={sih}; rel. to mixture best")
        print("  si\\f  " + " ".join(f"{f:>7}" for f in axes["f_hot"]))
        for si in axes["sigma_i"]:
            row = [mix[(si, f, sih)]["total"] - m0["total"]
                   for f in axes["f_hot"]]
            print(f"  {si:5} " + " ".join(f"{v:7.1f}" for v in row))
        print("G2 (same layout):")
        for si in axes["sigma_i"]:
            row = [mix[(si, f, sih)]["g2"] for f in axes["f_hot"]]
            print(f"  {si:5} " + " ".join(f"{v:7.1f}" for v in row))

    # figure: lnL and G2 vs f_hot at the best cold sigma_i
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 4.2))
    fs = [0.0] + list(axes["f_hot"])
    for sih in axes["sigma_i_hot"]:
        ll = [s0["total"]] + [mix[(si_m, f, sih)]["total"]
                              for f in axes["f_hot"]]
        g2 = [single[si_m]["g2"]] + [mix[(si_m, f, sih)]["g2"]
                                     for f in axes["f_hot"]]
        a1.plot(fs, np.array(ll) - s0["total"], "o-",
                label=rf"$\sigma_{{i,hot}}={sih}^\circ$")
        a2.plot(fs, g2, "o-", label=rf"$\sigma_{{i,hot}}={sih}^\circ$")
    a1.axhline(0, color="k", lw=0.6)
    a1.set_xlabel(r"$f_{hot}$")
    a1.set_ylabel(r"$\Delta \ln L$ vs single population")
    a2.axhline(7.8, color="r", lw=0.8, ls="--", label="G2 p=0.05 (3 dof)")
    a2.set_xlabel(r"$f_{hot}$")
    a2.set_ylabel(r"multiplicity deviance $G^2$")
    a1.legend(fontsize=8); a2.legend(fontsize=8)
    fig.suptitle(rf"Mixture test at $\sigma_R={SIGMA_R}$, "
                 rf"$\sigma_i={si_m}^\circ$ (cold)")
    fig.tight_layout()
    fig.savefig(ROOT / "results" / "figures" / "dichotomy_mixture.png",
                dpi=150)
    print("\nfigure: results/figures/dichotomy_mixture.png")


if __name__ == "__main__":
    main()
