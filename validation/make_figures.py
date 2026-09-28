"""
Figures and tables for the DockCert v2 manuscript, built only from validation/results/.

    python validation/make_figures.py

Style: Elsevier double-column width 190 mm, 8 pt sans-serif text, one fixed colour and marker per
RMSD method in every figure (validated categorical palette; identity never relies on colour alone).
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import beta  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
FIG = os.path.join(HERE, "figures")
TAB = os.path.join(HERE, "tables")
MM = 1 / 25.4
DOUBLE = 190 * MM

METHODS = {   # key -> (column, colour, marker, label)
    "rdkit": ("rmsd_rdkit", "#1baf7a", "D", "RDKit CalcRMS (reference)"),
    "dockcert": ("rmsd_dockcert", "#2a78d6", "o", "DockCert 1.1 (graph isomorphism)"),
    "naive": ("rmsd_naive", "#8a8984", "x", "no symmetry correction"),
    "hung_lb": ("rmsd_hung_lb", "#eda100", "v", "per-element Hungarian (lower bound)"),
    "hung_v100": ("rmsd_hung_v100", "#e87ba4", "^", "DockCert 1.0.0 Hungarian, as called"),
}
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"


def setup():
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "legend.fontsize": 6.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.color": INK2, "ytick.color": INK2, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.5, "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 1.2,
        "lines.markersize": 4, "savefig.dpi": 600, "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def panel(ax, letter):
    ax.text(-0.13, 1.03, f"({letter})", transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left", color=INK)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def clopper(k, n):
    lo = beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


def astex_success(poses):
    rows = []
    n = poses.id.nunique()
    top = poses[poses["rank"] == 1].set_index("id")
    for key, (col, *_rest, label) in METHODS.items():
        best = poses.groupby("id")[col].min()
        for crit, series in (("top-1", top[col]), ("best of 9", best)):
            k = int((series <= 2.0).sum())
            lo, hi = clopper(k, n)
            rows.append({"method": key, "label": label, "criterion": crit, "n": n, "success": k / n,
                         "ci_low": lo, "ci_high": hi, "k": k})
    return pd.DataFrame(rows)


def fig_astex(poses, succ):
    fig, (a, b) = plt.subplots(1, 2, figsize=(DOUBLE, 72 * MM), gridspec_kw={"width_ratios": [1, 1.15]})
    lim = 16
    a.plot([0, lim], [0, lim], color=INK2, lw=0.7, ls="--", zorder=1)
    a.axhline(2, color=INK2, lw=0.5, ls=":")
    a.axvline(2, color=INK2, lw=0.5, ls=":")
    for key in ("naive", "hung_v100", "hung_lb", "dockcert"):
        col, c, m, lab = METHODS[key]
        a.scatter(poses.rmsd_rdkit, poses[col], s=7 if key != "dockcert" else 9, color=c, marker=m,
                  lw=0.6 if m in "x+" else 0.2, edgecolor="white" if m not in "x+" else None, label=lab,
                  zorder=3 if key == "dockcert" else 2, alpha=0.9)
    a.set(xlim=(0, lim), ylim=(0, lim), xlabel="RDKit symmetry-corrected RMSD (Å)", ylabel="RMSD by method (Å)")
    a.set_aspect("equal")
    panel(a, "a")

    order = ["rdkit", "dockcert", "naive", "hung_lb", "hung_v100"]
    x = np.arange(2)
    w = 0.16
    for i, key in enumerate(order):
        col, c, m, lab = METHODS[key]
        d = succ[succ.method == key].set_index("criterion").loc[["top-1", "best of 9"]]
        xs = x + (i - 2) * w
        b.bar(xs, d.success, w * 0.9, color=c, label=lab, zorder=2)
        b.errorbar(xs, d.success, yerr=[d.success - d.ci_low, d.ci_high - d.success], fmt="none",
                   ecolor=INK2, elinewidth=0.6, capsize=1.5, zorder=3)
    b.set_xticks(x, ["top-ranked pose", "best of 9 poses"])
    b.set(ylabel="success rate (RMSD ≤ 2 Å)", ylim=(0, 1.0))
    panel(b, "b")
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=METHODS[k][1], marker=METHODS[k][2], ls="none", markersize=5,
                      markeredgewidth=0.8, label=METHODS[k][3]) for k in order]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0), columnspacing=1.5)
    fig.tight_layout(w_pad=3.0, rect=(0, 0.08, 1, 1))
    save(fig, "fig2_astex")


def table_astex(succ, cases):
    lines = [r"\begin{tabular}{lrr}", r"\toprule",
             r"RMSD definition & top-ranked pose & best of 9 poses \\", r"\midrule"]
    for key in ("rdkit", "dockcert", "naive", "hung_lb", "hung_v100"):
        d = succ[succ.method == key].set_index("criterion")
        cells = [f"{d.loc[c,'success']:.2f} [{d.loc[c,'ci_low']:.2f}, {d.loc[c,'ci_high']:.2f}]" for c in ("top-1", "best of 9")]
        lines.append(f"{METHODS[key][3]} & {cells[0]} & {cells[1]} \\\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    open(os.path.join(TAB, "table_astex.tex"), "w").write("\n".join(lines) + "\n")


def main():
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(TAB, exist_ok=True)
    setup()
    poses = pd.read_csv(os.path.join(RES, "astex_poses.csv"))
    cases = pd.read_csv(os.path.join(RES, "astex_cases.csv"))
    ok = set(cases[cases.ok].id)
    poses = poses[poses.id.isin(ok)]
    succ = astex_success(poses)
    succ.to_csv(os.path.join(RES, "astex_success.csv"), index=False)
    fig_astex(poses, succ)
    table_astex(succ, cases)
    print(succ[["method", "criterion", "k", "n", "success", "ci_low", "ci_high"]].to_string(index=False))


if __name__ == "__main__":
    main()
