#!/usr/bin/env python
"""Three panels: what the experiment shows, what the score shows, and why.

A  The published relationship on this compound set. Functional efficacy tracks
   log residence time and not affinity. This is the claim the review makes.
B  The docking score against the same three experimental axes.
C  Redocking RMSD against sampling budget, which separates a search failure
   from a scoring failure.
"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dataset import PUBLISHED

WORK = "/home/kjb9412/a2a_work"
OUT = f"{WORK}/out"
FIGDIR = "/home/kjb9412/docking-score-audit/figures/a2a"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.linewidth": 0.8, "savefig.dpi": 600, "savefig.bbox": "tight",
    "legend.fontsize": 7, "legend.frameon": False,
})
BLU, RED, GRY, GRN = "#0077B6", "#E63946", "#6C757D", "#2D6A4F"


def fit_line(ax, x, y, color):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return None
    sl, ic, r, p, _ = stats.linregress(x, y)
    xs = np.linspace(x.min(), x.max(), 50)
    ax.plot(xs, sl * xs + ic, "-", color=color, lw=1.1, zorder=1)
    return r ** 2, p, len(x)


def scatter(ax, x, y, labels, xlabel, ylabel, title, color, annotate=()):
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[ok], y[ok], s=26, c=color, edgecolors="black",
               linewidths=0.4, zorder=3)
    st = fit_line(ax, x, y, color)
    for nm in annotate:
        if nm in labels:
            i = labels.index(nm)
            if ok[i]:
                ax.annotate(nm, (x[i], y[i]), textcoords="offset points",
                            xytext=(5, 4), fontsize=6, color="black")
    ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
    if st:
        r2, p, n = st
        ax.set_title(f"{title}\n$r^2$={r2:.2f}, p={p:.3f}, n={n}", pad=4)
    else:
        ax.set_title(title, pad=4)
    ax.spines[["top", "right"]].set_visible(False)


def main():
    os.makedirs(FIGDIR, exist_ok=True)
    A = json.load(open(f"{OUT}/analysis.json"))
    rows = A["rows"]
    name = [r["name"] for r in rows]
    logRT = [r["logRT"] for r in rows]
    logKi = [np.log10(r["Ki_nM"]) for r in rows]
    pKi = [r["pKi"] for r in rows]
    eff = [r["eff_impedance"] for r in rows]
    cnn = [r["cnn_best"] for r in rows]
    MARK = ("UK432097", "LUF5835", "CGS21680")

    # ---------------- Panel A: the experiment ----------------
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 2.9))
    scatter(axes[0], logRT, eff, name, "log residence time (min)",
            "efficacy (% of CGS21680)", "Efficacy tracks residence time",
            BLU, MARK)
    scatter(axes[1], logKi, eff, name, "log $K_i$ (nM)",
            "efficacy (% of CGS21680)", "Efficacy does not track affinity",
            GRY, MARK)
    for ax, k in zip(axes, ("eff_impedance_vs_logRT", "eff_impedance_vs_logKi")):
        ax.text(0.97, 0.05, f"Guo 2012: $r^2$={PUBLISHED[k]['r2']}",
                transform=ax.transAxes, fontsize=6.5, va="bottom", ha="right",
                color=GRY)
    fig.suptitle("A   Published experiment, A$_{2A}$ agonists (Guo et al. 2012)",
                 x=0.01, ha="left", fontsize=9, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{FIGDIR}/FigA_experiment.pdf"); plt.close(fig)

    # ---------------- Panel B: the score ----------------
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 2.9))
    scatter(axes[0], pKi, cnn, name, "experimental p$K_i$",
            "GNINA CNNaffinity", "Score vs affinity", RED, MARK)
    scatter(axes[1], logRT, cnn, name, "log residence time (min)",
            "GNINA CNNaffinity", "Score vs kinetics", RED, MARK)
    scatter(axes[2], eff, cnn, name, "efficacy (% of CGS21680)",
            "GNINA CNNaffinity", "Score vs function", RED, MARK)
    fig.suptitle("B   Docking score against the same three axes",
                 x=0.01, ha="left", fontsize=9, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{FIGDIR}/FigB_score.pdf"); plt.close(fig)

    # ---------------- Panel C: search vs scoring ----------------
    sweep = A.get("sampling_sweep") or []
    if sweep:
        fig, ax = plt.subplots(figsize=(5.2, 2.9))
        for pdb, col in zip(sorted({s["pdb"] for s in sweep}), (BLU, GRN, RED)):
            pts = sorted([s for s in sweep if s["pdb"] == pdb and not s.get("error")],
                         key=lambda s: s["exhaustiveness"])
            if not pts:
                continue
            ax.plot([p["exhaustiveness"] for p in pts],
                    [p["rmsd_best"] for p in pts], "o-", color=col, lw=1.1,
                    ms=4, label=f"{pdb} (best pose)")
            ax.plot([p["exhaustiveness"] for p in pts],
                    [p["rmsd_rank1"] for p in pts], "s--", color=col, lw=0.9,
                    ms=3.5, alpha=0.6, label=f"{pdb} (top-scored)")
        ax.axhline(2.0, color="black", lw=0.7, ls=":")
        ax.text(16.4, 2.2, "2 A success threshold", fontsize=6, color="black")
        ax.set_xscale("log", base=2)
        ax.set_xticks([16, 32, 64]); ax.set_xticklabels(["16", "32", "64"])
        ax.set_xlabel("sampling budget (exhaustiveness)")
        ax.set_ylabel("redocking RMSD (A)")
        ax.set_ylim(0, 8.6)
        ax.set_title("C   More search does not help", loc="left",
                     fontsize=9, fontweight="bold", pad=6)
        ax.legend(ncol=1, loc="center left", bbox_to_anchor=(1.02, 0.5))
        ax.spines[["top", "right"]].set_visible(False)
        fig.tight_layout()
        fig.savefig(f"{FIGDIR}/FigC_sampling.pdf"); plt.close(fig)

    print("wrote panels to", FIGDIR)
    for f in sorted(os.listdir(FIGDIR)):
        print("  ", f)


if __name__ == "__main__":
    main()
