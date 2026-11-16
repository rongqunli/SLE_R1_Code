import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
from pathlib import Path
import hashlib, json
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "work"
WORK.mkdir(parents=True, exist_ok=True)
TFROOT = ROOT.parent / "tf_recheck_20261115"
GROOT = ROOT.parent / "targeted_recheck_20261115" / "association_results"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 8.5,
        "axes.titlesize": 9,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    }
)
SOURCES = []


def read(p):
    SOURCES.append({"path": str(p), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()})
    return pd.read_csv(p, sep="\t")


def save(fig, name):
    for ext in ["png", "pdf", "svg"]:
        fig.savefig(WORK / f"{name}.{ext}", dpi=450)
    plt.close(fig)


tf = read(TFROOT / "four_tf_primary_summary.tsv")
fig, axes = plt.subplots(2, 1, figsize=(5.12, 5.0))
fig.subplots_adjust(left=0.16, right=0.75, top=0.91, bottom=0.12, hspace=0.63)
for ax, model, tag, title in zip(
    axes,
    ["all_adjusted", "sle_adjusted"],
    ["A", "B"],
    ["All pediatric donors (n = 44)", "SLE complete cases (n = 31)"],
):
    rows = tf[tf.model.eq(model)].set_index("TF").loc[["IRF7", "STAT1", "STAT2", "IRF9"]]
    for j, (name, r) in enumerate(rows.iterrows()):
        y = 3 - j
        ax.errorbar(
            r.rho,
            y,
            xerr=[[r.rho - r.bootstrap_low], [r.bootstrap_high - r.rho]],
            fmt="o",
            color="#346A8C",
            ms=4.6,
            lw=1.25,
            capsize=2,
        )
        ax.text(
            1.035,
            y,
            f"{r.rho:.3f}   {r.permutation_fdr:.3g}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=8,
        )
    ax.text(1.035, 1.08, "ρ       q", transform=ax.transAxes, fontsize=8, weight="bold")
    ax.set(
        yticks=[3, 2, 1, 0],
        yticklabels=rows.index,
        xlim=(-0.5, 1),
        ylim=(-0.5, 3.55),
        xticks=[-0.5, 0, 0.5, 1],
        xlabel="Partial rank correlation with IFN score",
    )
    ax.axvline(0, color="#888888", ls="--", lw=0.6)
    ax.grid(axis="x", color="#E8E8E8", lw=0.5)
    ax.tick_params(axis="y", length=0)
    ax.set_title(title, loc="left", pad=12, weight="bold")
    ax.text(-0.24, 1.12, tag, transform=ax.transAxes, weight="bold", fontsize=11)
save(fig, "Figure6")
tf.to_csv(ROOT / "figures" / "Figure6_plotted_values.tsv", sep="\t", index=False)
tabs = [
    read(GROOT / n)
    for n in ["gse121239_patient_means_corrected.tsv", "gse121239_longitudinal_corrected.tsv"]
]
fig, axes = plt.subplots(2, 1, figsize=(5.12, 5.0))
fig.subplots_adjust(left=0.32, right=0.95, top=0.91, bottom=0.19, hspace=0.72)
labels = ["GZMK state", "B cell / BCR", "Immature neutrophil"]
for ax, tab, tag, title in zip(
    axes,
    tabs,
    ["A", "B"],
    ["Patient means: 65 patients", "Within patients: 65 patients, 292 visits"],
):
    for j, relation in enumerate(["GZMK", "B_cell", "Neutrophil"]):
        y = 2 - j
        for variant, offset, color, marker in [
            ("noCCL5", 0.12, "#2166AC", "o"),
            ("full", -0.12, "#949494", "s"),
        ]:
            r = tab[(tab.relation == relation) & (tab.score_variant == variant)].iloc[0]
            ax.errorbar(
                r.rho,
                y + offset,
                xerr=[[r.rho - r.ci_lo], [r.ci_hi - r.rho]],
                fmt=marker,
                ms=4,
                color=color,
                lw=1.2,
                capsize=2,
                zorder=3,
            )
    ax.axvline(0, color="#888888", lw=0.65, ls="--")
    ax.set(
        yticks=[2, 1, 0],
        yticklabels=labels,
        xlim=(-0.65, 0.92),
        ylim=(-0.45, 2.5),
        xticks=[-0.5, 0, 0.5],
        xlabel="IFN-adjusted rank association with cytotoxicity",
    )
    ax.set_title(title, loc="left", pad=11, weight="bold", fontsize=8.6)
    ax.text(-0.47, 1.14, tag, transform=ax.transAxes, weight="bold", fontsize=11)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#E8E8E8", lw=0.5)
legend = [
    Line2D(
        [0],
        [0],
        marker="o",
        color="#2166AC",
        lw=1.2,
        label="CCL5 removed from both GZMK and cytotoxicity",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        color="#949494",
        lw=1.2,
        label="Original complete scores (sensitivity)",
    ),
]
fig.legend(
    handles=legend,
    loc="lower center",
    bbox_to_anchor=(0.5, 0.027),
    ncol=1,
    frameon=False,
    fontsize=7.7,
)
save(fig, "Figure9")
pd.concat(tabs).to_csv(ROOT / "figures" / "Figure9_plotted_values.tsv", sep="\t", index=False)
(ROOT / "figures" / "corrected_estimate_sources.json").write_text(
    json.dumps(SOURCES, indent=2), encoding="utf-8"
)
print("Rendered Figure6 and Figure9 from frozen result tables.")
