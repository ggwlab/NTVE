"""
Suppl Fig 19a-c — compact lineage marker heatmaps

Manuscript-faithful refactor of the NTVE marker part of Supplementary Figure 19:
- a: ectoderm markers (SOX2, PAX6)
- b: mesoderm markers (CDH5, PDGFRB)
- c: endoderm markers (CXCR4, SOX17)

Important scope note:
- The manuscript panel also contains endpoint flow-cytometry plots.
- Those are intentionally not reproduced here because the corresponding raw FACS
  source data are not available in this repo.
- This script only regenerates the NTVE heatmap subpanels.
- It prefers the refactored raw-data rebuild with iPSC baseline:
    Figure6/figure6_three_lineage_data.pkl
- and falls back to the original notebook pickle if that file is absent.
"""

from pathlib import Path
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def find_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in (here.parent, here.parent.parent):
        if (candidate / "resources").exists():
            return candidate
    return here.parent.parent


def first_existing_path(*candidates: Path) -> Path:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


ROOT = find_root()
OUT_DIR = Path(__file__).parent / "suppl19_plots"
OUT_DIR.mkdir(exist_ok=True)

REFAC_PICKLE = ROOT / "Figure6" / "figure6_three_lineage_data.pkl"
LEGACY_PICKLE = first_existing_path(
    ROOT / "resources" / "figure6_three_lineage_data.pkl",
)

# The manuscript top subpanels show these two markers per lineage.
PANEL_MARKERS = {
    "a": {
        "title": "NTVE ectoderm",
        "genes": ["SOX2", "PAX6"],
    },
    "b": {
        "title": "NTVE mesoderm",
        "genes": ["CDH5", "PDGFRB"],
    },
    "c": {
        "title": "NTVE endoderm",
        "genes": ["CXCR4", "SOX17"],
    },
}

LINEAGE_ORDER = ["ectoderm", "endoderm", "mesoderm"]
LINEAGE_LABELS = {"ectoderm": "Ecto", "endoderm": "Endo", "mesoderm": "Meso"}
# The original notebook uses all available days in the metadata. With the
# refactored loader that means iPSC/day0 plus 1, 2, 5, 6, 7.
DAYS_TO_PLOT = [0, 1, 2, 5, 6, 7]
DAY_LABELS = {0: "iPSC", 1: "1", 2: "2", 5: "5", 6: "6", 7: "7"}


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    pickle_path = REFAC_PICKLE if REFAC_PICKLE.exists() else LEGACY_PICKLE
    with pickle_path.open("rb") as f:
        data = pickle.load(f)

    tpm_df = data["tpm_matrix_pc"].copy()
    sample_meta = data["sample_meta_df"].copy()
    gene_id_to_name = data.get("gene_id_to_gene_name", {})

    return tpm_df, sample_meta, gene_id_to_name


def build_gene_name_to_id(gene_id_to_name: dict[str, str]) -> dict[str, str]:
    gene_name_to_id: dict[str, str] = {}
    for gene_id, gene_name in gene_id_to_name.items():
        if gene_name and gene_name not in gene_name_to_id:
            gene_name_to_id[gene_name] = gene_id
    return gene_name_to_id


def average_by_lineage_and_day(
    tpm_df: pd.DataFrame,
    sample_meta: pd.DataFrame,
    gene_id: str,
) -> pd.DataFrame:
    values = []

    for lineage in LINEAGE_ORDER:
        row = []
        for day in DAYS_TO_PLOT:
            if day == 0:
                samples = sample_meta[sample_meta["day"] == 0].index.tolist()
            else:
                samples = sample_meta[
                    (sample_meta["lineage"] == lineage)
                    & (sample_meta["day"] == day)
                ].index.tolist()

            if samples:
                expr = pd.to_numeric(tpm_df.loc[gene_id, samples], errors="coerce")
                row.append(float(np.nanmean(expr)))
            else:
                row.append(np.nan)

        values.append(row)

    return pd.DataFrame(values, index=LINEAGE_ORDER, columns=DAYS_TO_PLOT)


def log_gene_matrix(df: pd.DataFrame) -> pd.DataFrame:
    return np.log10(df.astype(float) + 1.0)


def plot_panel(
    panel_key: str,
    title: str,
    gene_mats: list[tuple[str, pd.DataFrame]],
) -> None:
    vmax = max(np.nanmax(mat.to_numpy()) for _, mat in gene_mats)

    fig, axes = plt.subplots(
        nrows=len(gene_mats),
        ncols=1,
        figsize=(3.0, 2.5),
        sharex=True,
    )
    if len(gene_mats) == 1:
        axes = [axes]

    fig.suptitle(title, fontsize=9, fontweight="bold", y=0.98)

    for ax, (gene, mat) in zip(axes, gene_mats):
        im = ax.imshow(
            mat.to_numpy(),
            aspect="auto",
            cmap="Reds",
            vmin=0,
            vmax=vmax,
            interpolation="nearest",
        )
        ax.set_yticks(range(len(LINEAGE_ORDER)))
        ax.set_yticklabels([LINEAGE_LABELS[x] for x in LINEAGE_ORDER], fontsize=7)
        ax.set_xticks(range(len(DAYS_TO_PLOT)))
        ax.set_xticklabels([DAY_LABELS[d] for d in DAYS_TO_PLOT], fontsize=7)
        ax.set_ylabel("Lineage", fontsize=7)
        ax.set_title(gene, fontsize=8, pad=2)

        for x in range(1, len(DAYS_TO_PLOT)):
            ax.axvline(x - 0.5, color="white", lw=0.8, alpha=0.8)
        for y in range(1, len(LINEAGE_ORDER)):
            ax.axhline(y - 0.5, color="white", lw=0.8, alpha=0.8)

    axes[-1].set_xlabel("day", fontsize=7)
    cbar = fig.colorbar(im, ax=axes, fraction=0.04, pad=0.04)
    cbar.set_label("log10(TPM+1)", fontsize=7)
    cbar.ax.tick_params(labelsize=7)

    fig.tight_layout(rect=[0, 0, 0.96, 0.95])

    stem = f"suppl19{panel_key}_marker_heatmap"
    for ext in ("png", "svg"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)


def plot_combined(panels: list[tuple[str, str, list[tuple[str, pd.DataFrame]]]]) -> None:
    all_mats = [mat for _, _, gene_mats in panels for _, mat in gene_mats]
    vmax = max(np.nanmax(mat.to_numpy()) for mat in all_mats)

    fig, axes = plt.subplots(
        nrows=2,
        ncols=3,
        figsize=(7.2, 3.6),
        sharex=True,
        sharey=True,
    )

    for col_idx, (panel_key, title, gene_mats) in enumerate(panels):
        for row_idx, (gene, mat) in enumerate(gene_mats):
            ax = axes[row_idx, col_idx]
            im = ax.imshow(
                mat.to_numpy(),
                aspect="auto",
                cmap="Reds",
                vmin=0,
                vmax=vmax,
                interpolation="nearest",
            )
            if row_idx == 0:
                ax.set_title(f"{panel_key}    {title}", fontsize=9, loc="left", pad=8)
            ax.set_yticks(range(len(LINEAGE_ORDER)))
            ax.set_yticklabels([LINEAGE_LABELS[x] for x in LINEAGE_ORDER], fontsize=7)
            ax.set_xticks(range(len(DAYS_TO_PLOT)))
            ax.set_xticklabels([DAY_LABELS[d] for d in DAYS_TO_PLOT], fontsize=7)
            if col_idx == 0:
                ax.set_ylabel("Lineage", fontsize=7)
            ax.set_xlabel("day", fontsize=7)
            ax.text(
                0.5,
                1.02,
                gene,
                transform=ax.transAxes,
                ha="center",
                va="bottom",
                fontsize=8,
                fontstyle="italic",
            )
            for x in range(1, len(DAYS_TO_PLOT)):
                ax.axvline(x - 0.5, color="white", lw=0.8, alpha=0.8)
            for y in range(1, len(LINEAGE_ORDER)):
                ax.axhline(y - 0.5, color="white", lw=0.8, alpha=0.8)

    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.03)
    cbar.set_label("log10(TPM+1)", fontsize=7)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout(rect=[0, 0, 0.96, 1.0])

    for ext in ("png", "svg"):
        out = OUT_DIR / f"suppl19abc_marker_heatmaps.{ext}"
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)


def main() -> None:
    print("Loading three-lineage data ...")
    tpm_df, sample_meta, gene_id_to_name = load_data()
    gene_name_to_id = build_gene_name_to_id(gene_id_to_name)

    print(f"TPM matrix: {tpm_df.shape}")
    print(f"Sample metadata: {sample_meta.shape}")
    print(f"Available days: {sorted(sample_meta['day'].unique().tolist())}")
    print(f"Available lineages: {sorted(sample_meta['lineage'].unique().tolist())}")

    panels = []

    for panel_key, panel_cfg in PANEL_MARKERS.items():
        gene_mats: list[tuple[str, pd.DataFrame]] = []
        print(f"\nBuilding panel {panel_key}: {panel_cfg['title']}")

        for gene in panel_cfg["genes"]:
            gene_id = gene_name_to_id.get(gene)
            if gene_id is None:
                print(f"  Missing gene in mapping: {gene}")
                continue

            avg_df = average_by_lineage_and_day(tpm_df, sample_meta, gene_id)
            logged = log_gene_matrix(avg_df)
            gene_mats.append((gene, logged))
            print(f"  {gene}: using gene_id {gene_id}")

        if not gene_mats:
            print(f"  No genes available for panel {panel_key}; skipping")
            continue

        plot_panel(panel_key, panel_cfg["title"], gene_mats)
        panels.append((panel_key, panel_cfg["title"], gene_mats))

    if panels:
        plot_combined(panels)

    print("\nDone.")


if __name__ == "__main__":
    plt.rcParams.update({"font.size": 8, "svg.fonttype": "none"})
    main()
