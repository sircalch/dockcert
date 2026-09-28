"""
Graphical abstract for the DockCert v2 manuscript (Elsevier: 531 x 1328 px minimum, readable at
5 x 13 cm). Built only from validation/results/.
"""
import os

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
INK, INK2 = "#0b0b0b", "#52514e"


def main():
    succ = pd.read_csv(os.path.join(RES, "astex_success.csv"))
    cov = pd.read_csv(os.path.join(RES, "bootstrap_coverage.csv"))
    dl = pd.read_csv(os.path.join(RES, "auc_ci_coverage.csv"))
    matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                                "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                                "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
    fig = plt.figure(figsize=(13 / 2.54, 5 / 2.54))
    a = fig.add_axes([0.08, 0.2, 0.42, 0.55])
    b = fig.add_axes([0.63, 0.2, 0.35, 0.55])

    def get(m, c):
        return float(succ[(succ.method == m) & (succ.criterion == c)].success.iloc[0])

    bars = [("no\nsymmetry", get("naive", "top-1"), "#8a8984"),
            ("graph RMSD\n(DockCert)", get("dockcert", "top-1"), "#2a78d6"),
            ("Hungarian\n(low. bound)", get("hung_lb", "top-1"), "#eda100"),
            ("best of\n9 poses", get("dockcert", "best of 9"), "#1baf7a")]
    x = np.arange(len(bars))
    a.bar(x, [v for _, v, _ in bars], 0.7, color=[c for *_, c in bars])
    for xi, (_, v, _) in zip(x, bars):
        a.text(xi, v + 0.02, f"{v:.0%}", ha="center", va="bottom", fontsize=8, color=INK)
    a.set_xticks(x, [l for l, _, _ in bars], fontsize=6)
    a.set_ylim(0, 1.0)
    a.set_yticks([0, 0.5, 1.0], ["0", "50%", "100%"], fontsize=7)
    a.set_ylabel("RMSD ≤ 2 Å", fontsize=7)
    a.set_title("Same Vina poses, four ways to count success", fontsize=8, color=INK)

    for key, c, m, lab, src in (("auc", "#8a8984", "x", "AUC, bootstrap", cov),
                                ("auc_delong_logit", "#2a78d6", "o", "AUC, DeLong-logit", dl),
                                ("ef5", "#eda100", "v", "EF 5%, bootstrap", cov)):
        e = src[(src.mu == 2.0) & (src.metric == key)].sort_values("actives")
        b.plot(e.actives, e.coverage, color=c, marker=m, label=lab, ms=4)
    b.axhline(0.95, color=INK2, ls="--", lw=0.7)
    b.set_xscale("log")
    b.set_xticks([20, 50, 100], ["20", "50", "100"], fontsize=7)
    b.minorticks_off()
    b.set_ylim(0.80, 1.0)
    b.set_yticks([0.8, 0.85, 0.9, 0.95, 1.0], ["0.80", "0.85", "0.90", "0.95", "1.00"], fontsize=7)
    b.set_xlabel("actives in the library", fontsize=7)
    b.set_title("Coverage of 95% intervals", fontsize=8, color=INK)
    b.legend(fontsize=5.5, frameon=False, loc="lower right", borderaxespad=0.2)
    fig.text(0.5, 0.975, "DockCert: symmetry, ranking and uncertainty in docking benchmarks",
             ha="center", va="top", fontsize=9, fontweight="bold", color=INK)
    out = os.path.join(HERE, "figures", "graphical_abstract")
    fig.savefig(out + ".png", dpi=400)
    fig.savefig(out + ".pdf")
    from PIL import Image
    w, h = Image.open(out + ".png").size
    print("graphical abstract", w, "x", h, "px")


if __name__ == "__main__":
    main()
