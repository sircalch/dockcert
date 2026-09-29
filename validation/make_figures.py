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


def table_astex_per_complex(poses, run_status):
    """Supplementary table: one row per complex."""
    top = poses[poses["rank"] == 1].set_index("id")
    best = poses.groupby("id").rmsd_rdkit.min()
    npose = poses.groupby("id").size()
    head = (r"Complex & RDKit & DockCert & naive & Hung.\ LB & best & poses & score & heterogens \\")
    lines = [r"\begin{longtable}{lrrrrrrrp{2.6cm}}",
             r"\caption{Per-complex redocking results on the Astex Diverse Set (AutoDock Vina 1.2.7). "
             r"RMSD (\AA) of the top-ranked pose by four definitions, best RMSD over all poses, number of "
             r"poses, Vina score of the top-ranked pose (kcal/mol) and heterogens kept in the receptor.}"
             r"\label{tab:s_astex}\\",
             r"\toprule", head, r"\midrule", r"\endfirsthead",
             r"\toprule", head, r"\midrule", r"\endhead"]
    for cid in sorted(top.index):
        r = top.loc[cid]
        het = ", ".join(run_status.get(cid, {}).get("kept_heterogens", [])) or "--"
        cid_tex = cid.replace("_", r"\_")
        lines.append(f"{cid_tex} & {r.rmsd_rdkit:.2f} & {r.rmsd_dockcert:.2f} & {r.rmsd_naive:.2f} & "
                     f"{r.rmsd_hung_lb:.2f} & {best[cid]:.2f} & {npose[cid]} & {r.vina_kcal:.2f} & {het} " + r"\\")
    lines += [r"\bottomrule", r"\end{longtable}"]
    with open(os.path.join(TAB, "table_s_astex.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


COV = {   # metric -> (colour, marker, label)
    "auc": ("#8a8984", "x", "ROC-AUC, percentile bootstrap"),
    "auc_delong_logit": ("#2a78d6", "o", "ROC-AUC, DeLong (logit)"),
    "bedroc20": ("#1baf7a", "D", r"BEDROC ($\alpha$ = 20), bootstrap"),
    "ef1": ("#e87ba4", "^", r"EF$_{1\%}$, bootstrap"),
    "ef5": ("#eda100", "v", r"EF$_{5\%}$, bootstrap"),
}


def fig_coverage(boot, delong):
    from matplotlib.lines import Line2D
    from scipy.stats import norm
    d = pd.concat([boot[["mu", "actives", "metric", "coverage", "reps"]],
                   delong[["mu", "actives", "metric", "coverage", "reps"]]], ignore_index=True)
    fig, axes = plt.subplots(1, 2, figsize=(DOUBLE, 64 * MM), sharey=True)
    reps = int(d.reps.iloc[0])
    band = 1.96 * np.sqrt(0.95 * 0.05 / reps)
    for ax, mu, letter in zip(axes, (1.0, 2.0), "ab"):
        ax.axhspan(0.95 - band, 0.95 + band, color=INK2, alpha=0.08, lw=0)
        ax.axhline(0.95, color=INK2, lw=0.7, ls="--")
        for key, (c, m, lab) in COV.items():
            e = d[(d.mu == mu) & (d.metric == key)].sort_values("actives")
            ax.plot(e.actives, e.coverage, color=c, marker=m, label=lab, ms=4.5,
                    mew=0.9 if m in "x+" else 0.5, mec=c if m in "x+" else "white")
        ax.set(xscale="log", xlabel="actives in the library (50 decoys per active)", ylim=(0.82, 1.0),
               title=rf"separation $\mu$ = {mu:g} (true ROC-AUC = {norm.cdf(mu / np.sqrt(2)):.2f})")
        ax.set_xticks([20, 50, 100], ["20", "50", "100"])
        ax.minorticks_off()
        panel(ax, letter)
    axes[0].set_ylabel("coverage (nominal 0.95)")
    handles = [Line2D([], [], color=c, marker=m, label=l, ms=5, mew=0.9 if m in "x+" else 0.5,
                      mec=c if m in "x+" else "white") for c, m, l in COV.values()]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0), columnspacing=1.5)
    fig.tight_layout(w_pad=2.0, rect=(0, 0.12, 1, 1))
    save(fig, "fig3_coverage")


def table_coverage(boot, delong):
    d = pd.concat([boot, delong], ignore_index=True)
    names = {"auc": "ROC-AUC (bootstrap)", "auc_delong_logit": "ROC-AUC (DeLong, logit)",
             "bedroc20": r"BEDROC ($\alpha=20$)", "ef1": r"EF$_{1\%}$", "ef5": r"EF$_{5\%}$"}
    lines = [r"\begin{tabular}{llrrrrr}", r"\toprule",
             r"$\mu$ & metric & true value & 20 actives & 50 actives & 100 actives & width (50 actives) \\",
             r"\midrule"]
    for mu in (1.0, 2.0):
        for key, name in names.items():
            e = d[(d.mu == mu) & (d.metric == key)].set_index("actives")
            cov = " & ".join(f"{e.loc[a, 'coverage']:.3f}" for a in (20, 50, 100))
            lines.append(f"{mu:g} & {name} & {e['true'].iloc[0]:.3f} & {cov} & {e.loc[50, 'median_width']:.3f} " + r"\\")
        if mu == 1.0:
            lines.append(r"\midrule")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_coverage.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def table_dude(m):
    """DUD-E enrichment metrics with 95% intervals (DeLong for ROC-AUC, bootstrap otherwise)."""
    names = {"ampc": "AmpC", "inha": "InhA"}
    lines = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
             r"Target & actives / decoys & ROC-AUC & BEDROC ($\alpha=20$) & BEDROC ($\alpha=80.5$) & EF$_{1\%}$ (max) & EF$_{5\%}$ \\",
             r"\midrule"]
    for _, r in m.iterrows():
        lines.append(
            f"{names.get(r.target, r.target)} & {int(r.actives)} / {int(r.decoys)} & "
            f"{r.auc:.2f} [{r.auc_lo:.2f}, {r.auc_hi:.2f}] & "
            f"{r.bedroc20:.2f} [{r.bedroc20_lo:.2f}, {r.bedroc20_hi:.2f}] & "
            f"{r.bedroc80:.2f} [{r.bedroc80_lo:.2f}, {r.bedroc80_hi:.2f}] & "
            f"{r.ef1:.1f} [{r.ef1_lo:.1f}, {r.ef1_hi:.1f}] ({r.ef1_max:.0f}) & "
            f"{r.ef5:.1f} [{r.ef5_lo:.1f}, {r.ef5_hi:.1f}] " + r"\\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    with open(os.path.join(TAB, "table_dude.tex"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


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
    import json, zipfile
    status = {}
    zp = os.path.join(RES, "astex_vina_outputs.zip")
    if os.path.exists(zp):
        with zipfile.ZipFile(zp) as z:
            for n in z.namelist():
                if n.endswith("status.json"):
                    st = json.loads(z.read(n))
                    status[st["id"]] = st
    table_astex_per_complex(poses, status)
    bp, dp = os.path.join(RES, "bootstrap_coverage.csv"), os.path.join(RES, "auc_ci_coverage.csv")
    if os.path.exists(bp) and os.path.exists(dp):
        boot, delong = pd.read_csv(bp), pd.read_csv(dp)
        fig_coverage(boot, delong)
        table_coverage(boot, delong)
    mp = os.path.join(RES, "dude_metrics.csv")
    if os.path.exists(mp):
        table_dude(pd.read_csv(mp))
    print(succ[["method", "criterion", "k", "n", "success", "ci_low", "ci_high"]].to_string(index=False))


if __name__ == "__main__":
    main()
