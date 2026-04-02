"""
Figure 6d/e/f — static reproduction from the original standalone scatter export.

This script uses the compressed data payload from the legacy interactive Figure 6
standalone export, but renders simple static manuscript-style panels and exports
tidy numerical tables for the showcased genes and scatter coordinates.

The showcase genes were inferred from the manuscript panel labels.
"""

from __future__ import annotations

import base64
import gzip
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd


def find_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in (here.parent, here.parent.parent):
        if (candidate / "resources").exists():
            return candidate
    return here.parent.parent


ROOT = find_root()
OUT_DIR = Path(__file__).parent / "6def_plots"
CSV_DIR = Path(__file__).parent / "6def_csv"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR.mkdir(exist_ok=True)

DATA_JS = ROOT / "resources" / "fig6_6def" / "data_compressed.js"

LINEAGE_ORDER = ["ectoderm", "mesoderm", "endoderm"]
LINEAGE_COLORS = {
    "ectoderm": "#e74c3c",
    "mesoderm": "#3498db",
    "endoderm": "#2ecc71",
}
LINEAGE_PANEL_LETTERS = {"ectoderm": "d", "mesoderm": "e", "endoderm": "f"}
LINEAGE_TITLES = {"ectoderm": "ectoderm", "mesoderm": "mesoderm", "endoderm": "endoderm"}
SHOWCASE_GENES = {
    "ectoderm": ["ZMIZ1", "PAX6", "CLDN4"],
    "mesoderm": ["CARM1", "MSX2", "MESP1"],
    "endoderm": ["GPX4", "CKB", "SUCLG1"],
}

LRT_X_THRESHOLD = 3.841
LRT_Y_THRESHOLD = 7.815


def decode_payload() -> tuple[dict, dict]:
    text = DATA_JS.read_text(encoding="utf-8")
    scatter_match = re.search(r'window\.SCATTER_B64="([^"]+)"', text)
    ts_match = re.search(r'window\.TS_B64="([^"]+)"', text)
    if not scatter_match or not ts_match:
        raise RuntimeError(f"Could not find compressed payloads in {DATA_JS}")

    scatter_payload = json.loads(gzip.decompress(base64.b64decode(scatter_match.group(1))))
    ts_payload = json.loads(gzip.decompress(base64.b64decode(ts_match.group(1))))
    return scatter_payload, ts_payload


def load_scatter_tables(scatter_payload: dict) -> dict[str, pd.DataFrame]:
    tables: dict[str, pd.DataFrame] = {}
    for lineage, payload in scatter_payload.items():
        df = pd.DataFrame(
            payload["customdata"],
            columns=[
                "Gene",
                "GeneName",
                "LRT_stat",
                "stat_full",
                "LRT_day2",
                "LRT_day7",
                "log2FoldChange",
            ],
        )
        df["x_lrt"] = pd.to_numeric(payload["x_lrt"])
        df["y_stat_full"] = pd.to_numeric(payload["y"])
        for col in ["LRT_stat", "stat_full", "LRT_day2", "LRT_day7", "log2FoldChange"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        df["panel_lineage"] = lineage
        tables[lineage] = df
    return tables


def log10_mean_sem(mean: float, sem: float) -> tuple[float, float]:
    y = float(np.log10(mean + 1.0))
    yerr = float(np.log10(mean + sem + 1.0) - y) if sem is not None else np.nan
    return y, max(yerr, 0.0)


def extract_showcase_timeseries(
    lineage: str,
    scatter_df: pd.DataFrame,
    ts_payload: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    ts_lineage = ts_payload[lineage]
    genes = ts_lineage["genes"]
    gene_to_idx = {gene_id: idx for idx, gene_id in enumerate(genes)}
    rows = []
    selected_rows = []

    showcase_df = scatter_df[scatter_df["GeneName"].isin(SHOWCASE_GENES[lineage])].copy()
    showcase_df["showcase_rank"] = showcase_df["GeneName"].map(
        {gene: rank for rank, gene in enumerate(SHOWCASE_GENES[lineage])}
    )
    showcase_df = showcase_df.sort_values("showcase_rank")

    missing = [g for g in SHOWCASE_GENES[lineage] if g not in showcase_df["GeneName"].tolist()]
    if missing:
        raise FileNotFoundError(
            f"Missing showcase genes for {lineage}: {', '.join(missing)}. "
            f"Check resources/fig6_6def/data_compressed.js."
        )

    for _, row in showcase_df.iterrows():
        gene_id = row["Gene"]
        gene_name = row["GeneName"]
        gene_idx = gene_to_idx.get(gene_id)
        if gene_idx is None:
            raise KeyError(f"{gene_name} ({gene_id}) not found in timeseries payload for {lineage}")

        selected_rows.append(
            {
                "panel_lineage": lineage,
                "Gene": gene_id,
                "GeneName": gene_name,
                "LRT_stat": row["LRT_stat"],
                "stat_full": row["stat_full"],
                "LRT_day2": row["LRT_day2"],
                "LRT_day7": row["LRT_day7"],
                "log2FoldChange": row["log2FoldChange"],
            }
        )

        for plotted_lineage in LINEAGE_ORDER:
            for day in ts_lineage["days"]:
                mean_list = ts_lineage["ts"][plotted_lineage][str(day)]
                sem_list = ts_lineage["ts"][plotted_lineage][f"{day}_sem"]
                mean_cpm = mean_list[gene_idx]
                sem_cpm = sem_list[gene_idx]
                if mean_cpm is None:
                    continue
                log10_mean, log10_sem = log10_mean_sem(float(mean_cpm), float(sem_cpm or 0.0))
                rows.append(
                    {
                        "panel_lineage": lineage,
                        "showcase_gene_name": gene_name,
                        "showcase_gene_id": gene_id,
                        "plotted_lineage": plotted_lineage,
                        "day": int(day),
                        "mean_cpm": float(mean_cpm),
                        "sem_cpm": float(sem_cpm or 0.0),
                        "log10_cpm_plus_1": log10_mean,
                        "log10_sem_approx": log10_sem,
                    }
                )

    return pd.DataFrame(rows), pd.DataFrame(selected_rows)


def get_log2fc_norm(scatter_tables: dict[str, pd.DataFrame]) -> colors.Normalize:
    all_values = pd.concat(
        [df["log2FoldChange"].abs() for df in scatter_tables.values()],
        ignore_index=True,
    ).dropna()
    vmax = max(float(all_values.quantile(0.995)), 1.0)
    return colors.TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)


def draw_lineage_panel(
    fig: plt.Figure,
    outer_spec,
    lineage: str,
    scatter_df: pd.DataFrame,
    showcase_ts_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    norm: colors.Normalize,
) -> None:
    subgrid = GridSpecFromSubplotSpec(
        2,
        3,
        subplot_spec=outer_spec,
        height_ratios=[1.0, 1.7],
        hspace=0.35,
        wspace=0.35,
    )

    gene_order = SHOWCASE_GENES[lineage]
    for idx, gene_name in enumerate(gene_order):
        ax = fig.add_subplot(subgrid[0, idx])
        subset = showcase_ts_df[showcase_ts_df["showcase_gene_name"] == gene_name].copy()
        for plotted_lineage in LINEAGE_ORDER:
            line_df = subset[subset["plotted_lineage"] == plotted_lineage].sort_values("day")
            ax.errorbar(
                line_df["day"],
                line_df["log10_cpm_plus_1"],
                yerr=line_df["log10_sem_approx"],
                color=LINEAGE_COLORS[plotted_lineage],
                marker="o",
                lw=1.6 if plotted_lineage == lineage else 1.2,
                ms=4,
                capsize=2,
                alpha=1.0 if plotted_lineage == lineage else 0.7,
            )
        ax.set_title(gene_name, fontsize=8, fontstyle="italic", pad=2)
        ax.set_xticks([2, 4, 6])
        ax.set_xlim(1, 7)
        if idx == 0:
            ax.set_ylabel("log$_{10}$ CPM(+1)", fontsize=8)
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if idx == 1:
            ax.set_xlabel("days", fontsize=8)

    ax = fig.add_subplot(subgrid[1, :])
    scatter = ax.scatter(
        scatter_df["LRT_stat"],
        scatter_df["stat_full"],
        c=scatter_df["log2FoldChange"],
        cmap="PiYG",
        norm=norm,
        s=3.5,
        alpha=0.55,
        linewidths=0,
        rasterized=True,
    )
    selected_points = scatter_df[scatter_df["GeneName"].isin(gene_order)].copy()
    ax.scatter(
        selected_points["LRT_stat"],
        selected_points["stat_full"],
        s=20,
        facecolors="none",
        edgecolors="black",
        linewidths=0.8,
        zorder=4,
    )
    for _, row in selected_points.iterrows():
        ax.annotate(
            row["GeneName"],
            (row["LRT_stat"], row["stat_full"]),
            xytext=(3, 3),
            textcoords="offset points",
            fontsize=7,
            fontstyle="italic",
        )
    ax.axvline(LRT_X_THRESHOLD, color="#aa2222", lw=0.8, ls="--")
    ax.axhline(LRT_Y_THRESHOLD, color="#aa2222", lw=0.8, ls="--")
    ax.set_xlabel("LRT of expression over all timepoints", fontsize=8)
    ax.set_ylabel("LRT of tissue-specific expression trajectory", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title(
        f"{LINEAGE_PANEL_LETTERS[lineage]}    {LINEAGE_TITLES[lineage]}",
        loc="left",
        fontsize=11,
        fontweight="bold",
        color=LINEAGE_COLORS[lineage],
        pad=4,
    )


def save_panel_outputs(
    lineage: str,
    scatter_df: pd.DataFrame,
    showcase_ts_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    norm: colors.Normalize,
) -> None:
    fig = plt.figure(figsize=(5.0, 4.8), constrained_layout=False)
    draw_lineage_panel(fig, GridSpec(1, 1, figure=fig)[0], lineage, scatter_df, showcase_ts_df, selected_df, norm)
    cax = fig.add_axes([0.92, 0.18, 0.018, 0.26])
    sm = ScalarMappable(norm=norm, cmap="PiYG")
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("avg log2FC", fontsize=7)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout(rect=[0.0, 0.0, 0.9, 1.0])

    stem = f"6{LINEAGE_PANEL_LETTERS[lineage]}_{lineage}_static"
    for ext in ("png", "svg"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, dpi=250, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)


def save_combined_output(
    scatter_tables: dict[str, pd.DataFrame],
    showcase_ts_tables: dict[str, pd.DataFrame],
    selected_tables: dict[str, pd.DataFrame],
    norm: colors.Normalize,
) -> None:
    fig = plt.figure(figsize=(12.8, 4.8), constrained_layout=False)
    outer = GridSpec(1, 3, figure=fig, wspace=0.35)
    for idx, lineage in enumerate(LINEAGE_ORDER):
        draw_lineage_panel(
            fig,
            outer[idx],
            lineage,
            scatter_tables[lineage],
            showcase_ts_tables[lineage],
            selected_tables[lineage],
            norm,
        )
    cax = fig.add_axes([0.985, 0.18, 0.01, 0.26])
    sm = ScalarMappable(norm=norm, cmap="PiYG")
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label("avg log2FC", fontsize=7)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout(rect=[0.0, 0.0, 0.98, 1.0])

    for ext in ("png", "svg"):
        out = OUT_DIR / f"6def_combined_static.{ext}"
        fig.savefig(out, dpi=250, bbox_inches="tight")
        print(f"Saved: {out}")
    plt.close(fig)


def main() -> None:
    if not DATA_JS.exists():
        raise FileNotFoundError(
            f"Missing Figure 6 standalone payload at {DATA_JS}. "
            "Add it under resources/fig6_6def/data_compressed.js."
        )

    scatter_payload, ts_payload = decode_payload()
    scatter_tables = load_scatter_tables(scatter_payload)
    norm = get_log2fc_norm(scatter_tables)

    showcase_ts_tables: dict[str, pd.DataFrame] = {}
    selected_tables: dict[str, pd.DataFrame] = {}
    all_showcase_ts = []
    all_selected = []
    all_scatter = []

    for lineage in LINEAGE_ORDER:
        scatter_df = scatter_tables[lineage].copy()
        showcase_ts_df, selected_df = extract_showcase_timeseries(lineage, scatter_df, ts_payload)
        showcase_ts_tables[lineage] = showcase_ts_df
        selected_tables[lineage] = selected_df
        all_showcase_ts.append(showcase_ts_df)
        all_selected.append(selected_df)
        all_scatter.append(scatter_df)
        save_panel_outputs(lineage, scatter_df, showcase_ts_df, selected_df, norm)

    save_combined_output(scatter_tables, showcase_ts_tables, selected_tables, norm)

    scatter_out = pd.concat(all_scatter, ignore_index=True)
    scatter_out.to_csv(CSV_DIR / "6def_scatter_points.csv", index=False)
    print(f"Saved: {CSV_DIR / '6def_scatter_points.csv'}")

    selected_out = pd.concat(all_selected, ignore_index=True)
    selected_out.to_csv(CSV_DIR / "6def_showcase_gene_stats.csv", index=False)
    print(f"Saved: {CSV_DIR / '6def_showcase_gene_stats.csv'}")

    ts_out = pd.concat(all_showcase_ts, ignore_index=True)
    ts_out.to_csv(CSV_DIR / "6def_showcase_timeseries.csv", index=False)
    print(f"Saved: {CSV_DIR / '6def_showcase_timeseries.csv'}")

    print("Done.")


if __name__ == "__main__":
    main()
