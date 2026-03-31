"""
Supplementary Figure 21a — compact marker heatmap used in the manuscript.

Source logic:
- Suppl21/cardio_marker_heatmaps.ipynb
- Manuscript assembly in NCOMMS-25-09704A_NTVE-2.pdf page 48

This panel is not the large iPSC->CM overview heatmap from the notebook.
It is the compact NTVE-only day 0-9 marker panel shown in the manuscript.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


ROOT = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent / "suppl21a_plots"
OUT_DIR.mkdir(exist_ok=True)
CARDIO_PKL = ROOT / "Suppl21" / "cardio_data.pkl"

MESODERM_MARKERS = ["TBXT", "PDGFRA", "MESP1", "HAND1", "CDH2", "HAND2"]
CARDIAC_MARKERS = ["NKX2-5", "ACTN2", "TNNT2"]
DAYS = list(range(10))
PSEUDOCOUNT = 0.1


def save(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, dpi=300, bbox_inches="tight")
        print(f"Saved: {out}")


def load_marker_matrix() -> pd.DataFrame:
    print("Loading cardio_data.pkl ...")
    with open(CARDIO_PKL, "rb") as handle:
        data = pickle.load(handle)

    tpm = data["tpm_matrix_pc"].copy()
    sample_meta = data["sample_meta_df"].copy()
    gene_id_to_name = data["gene_id_to_gene_name"]
    gene_name_to_id = {v: k for k, v in gene_id_to_name.items()}

    rows = []
    labels = []
    for marker in MESODERM_MARKERS + CARDIAC_MARKERS:
        gene_id = gene_name_to_id.get(marker)
        if gene_id is None or gene_id not in tpm.index:
            continue
        values = []
        for day in DAYS:
            samples = sample_meta[
                (sample_meta["cell_type"] == "cm")
                & (sample_meta["sample_type"] == "sn")
                & (sample_meta["day"] == float(day))
            ].index.tolist()
            values.append(float(tpm.loc[gene_id, samples].mean()))
        rows.append(values)
        labels.append(marker)

    marker_df = pd.DataFrame(rows, index=labels, columns=DAYS)
    print("Markers included:", list(marker_df.index))
    log_df = np.log2(marker_df + PSEUDOCOUNT)
    return log_df.sub(log_df.mean(axis=1), axis=0)


def draw_heatmap(ax: plt.Axes, data: pd.DataFrame, title: str, show_xlabels: bool) -> plt.AxesImage:
    vmin, vmax = -1.5, 1.5
    im = ax.imshow(data.values, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)
    ax.set_yticks(range(len(data.index)))
    ax.set_yticklabels(data.index, fontsize=8, fontstyle="italic")
    ax.set_xticks(range(len(DAYS)))
    ax.set_xticklabels([f"d{d}" for d in DAYS] if show_xlabels else [], fontsize=8)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=10, weight="bold", pad=8)
    for x in range(len(DAYS) + 1):
        ax.axvline(x - 0.5, color="#999999", lw=0.6, ls=(0, (2, 2)), zorder=0)
    for spine in ax.spines.values():
        spine.set_edgecolor("#cccccc")
    return im


def main() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none"})
    marker_df = load_marker_matrix()
    meso = marker_df.loc[[g for g in MESODERM_MARKERS if g in marker_df.index]]
    cardio = marker_df.loc[[g for g in CARDIAC_MARKERS if g in marker_df.index]]

    fig = plt.figure(figsize=(8.2, 3.8))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[0.08, 1.0], height_ratios=[0.9, 0.55, 0.35], wspace=0.18, hspace=0.35)
    ax_scale = fig.add_subplot(gs[0:2, 0])
    ax_top = fig.add_subplot(gs[0, 1])
    ax_bottom = fig.add_subplot(gs[1, 1], sharex=ax_top)
    ax_timeline = fig.add_subplot(gs[2, 1], sharex=ax_top)

    im = draw_heatmap(ax_top, meso, "Mesodermal marker genes quantified by NTVE", show_xlabels=False)
    draw_heatmap(ax_bottom, cardio, "Cardiac marker genes quantified by NTVE", show_xlabels=False)

    cbar = plt.colorbar(im, cax=ax_scale)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label("log2(TPM+0.1) - mean", fontsize=8)

    ax_timeline.set_xlim(-0.5, len(DAYS) - 0.5)
    ax_timeline.set_ylim(0, 1)
    ax_timeline.spines["left"].set_visible(False)
    ax_timeline.spines["right"].set_visible(False)
    ax_timeline.spines["bottom"].set_visible(False)
    ax_timeline.set_yticks([])
    ax_timeline.set_xticks(range(len(DAYS)))
    ax_timeline.set_xticklabels([f"d{d}" for d in DAYS], fontsize=8)
    ax_timeline.tick_params(length=0)
    for x in range(len(DAYS) + 1):
        ax_timeline.axvline(x - 0.5, color="#999999", lw=0.6, ls=(0, (2, 2)), zorder=0)
    ax_timeline.set_title("Timepoint (days)", fontsize=10, pad=2)
    ax_timeline.annotate(
        "mesoderm\ninduction",
        xy=(0.0, 0.04),
        xytext=(0.55, 0.88),
        ha="center",
        va="bottom",
        fontsize=8,
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black"),
    )
    ax_timeline.annotate(
        "cardiac\ninduction",
        xy=(2.0, 0.04),
        xytext=(2.45, 0.88),
        ha="center",
        va="bottom",
        fontsize=8,
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black"),
    )
    ax_timeline.annotate(
        "first\ncontraction",
        xy=(6.0, 0.04),
        xytext=(6.2, 0.88),
        ha="center",
        va="bottom",
        fontsize=8,
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black"),
    )

    fig.tight_layout()
    save(fig, "suppl21a_marker_heatmap_compact")
    plt.close(fig)
    print("Done.")


if __name__ == "__main__":
    main()
