"""
Supplementary Figure 13b — TES-aligned coverage profiles

Notebook-faithful reproduction of the SPAM outputs from
Suppl12_13/coverage_absolute_tes_refactored_cleaned.ipynb.

Exports:
  - SPAM lysate line plot
  - SPAM supernatant line plot
  - SPAM SN/L ratio plot
  - combined 3-panel summary sheet
"""

from collections import defaultdict
from pathlib import Path
import re
import sqlite3
import warnings

import lz4.frame
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
    return here.parent


ROOT = find_root()
OUT_DIR = Path(__file__).parent / "suppl13b_plots"
OUT_DIR.mkdir(exist_ok=True)
CSV_DIR = Path(__file__).parent / "suppl13_csv"
CSV_DIR.mkdir(exist_ok=True)

warnings.filterwarnings("ignore")
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 8

DB_PATH = ROOT / "resources" / "transcript_coverage.db"
LIMS_CSV_PATH = ROOT / "resources" / "Project_1716_lims_simplified_Spam2_deleted.csv"
SALMON_QUANT_FILE = ROOT / "resources" / "quantfiles_filtered_pipeline" / "homo_sapiens" / "25L008000.sf.gz"
PLOT_WIDTH = 3.5
PLOT_HEIGHT = 3.5

NAME_MAPPING = {
    "Benchmarking-S1": "NTVE",
    "Benchmarking-S2": "HIV-gag WT",
    "Benchmarking-S3": "TRACE-seq",
    "Benchmarking-S4": "MMLV-gag WT",
    "Benchmarking-S5": "HIV-gag:MCP",
    "SPAM-S1": "NTVE",
    "SPAM-S2": "HIV-gag:sPAM2 (p6)",
    "SPAM-S3": "HIV-gag:sPAM2 (NC)",
    "SPAM-S4": "HIV-gag WT",
    "SPAM-S5": "mini-gag:PABP",
    "SPAM-S6": "NLuc",
    "SPAM-S7": "NTVE (G2A)",
}
COLOR_MAPPING = {
    "NTVE": "#BBE5F0FF",
    "HIV-gag WT": "#9FAED5FF",
    "TRACE-seq": "#E28E7AFF",
    "MMLV-gag WT": "#93AE8FFF",
    "HIV-gag:MCP": "#C4D9C1F8",
    "HIV-gag:sPAM2 (p6)": "#FF6B6BFF",
    "HIV-gag:sPAM2 (NC)": "#4ECDC4FF",
    "mini-gag:PABP": "#F2D4E2",
    "NLuc": "#FFA500FF",
    "NTVE (G2A)": "#A8E6CFFF",
}


def get_db_connection(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def decompress_lz4(compressed_data: bytes) -> np.ndarray:
    return np.frombuffer(lz4.frame.decompress(compressed_data), dtype=np.uint16)


def batch_load_transcripts(transcript_ids: list[str]) -> dict:
    conn = get_db_connection(DB_PATH)
    cursor = conn.cursor()
    results = {tid: {} for tid in transcript_ids}
    placeholders = ",".join("?" * len(transcript_ids))
    cursor.execute(
        f"""
        SELECT transcript_id, sample_id, compressed_coverage
        FROM transcript_coverage
        WHERE transcript_id IN ({placeholders})
        """,
        transcript_ids,
    )
    for row in cursor.fetchall():
        results[row[0]][row[1]] = decompress_lz4(row[2])
    conn.close()
    return results


def parse_spam_sample_name(name: str):
    match = re.match(r"(SN|L)_(\d+)_(\d+)", str(name))
    if match:
        return {"Type": match.group(1), "Sample": int(match.group(2)), "Replicate": int(match.group(3)), "Experiment": "SPAM"}
    return None


def parse_benchmarking_sample_name(name: str):
    match = re.match(r"(SN|Lysate)_S(\d+)_(\d+)", str(name))
    if match:
        return {"Type": match.group(1), "Sample": int(match.group(2)), "Replicate": int(match.group(3)), "Experiment": "Benchmarking"}
    return None


def select_top_transcripts(quant_file: Path, top_n: int = 6000) -> dict[str, float]:
    df = pd.read_csv(quant_file, sep="\t")
    df = df[df["Name"] != "mGreenLantern"]
    df = df.sort_values("NumReads", ascending=False).head(top_n)
    return {row["Name"]: float(row["NumReads"]) for _, row in df.iterrows()}


def get_coverage_positions_relative_to_tes(exonic_length: int) -> np.ndarray:
    tes_index = exonic_length - 1
    return np.arange(exonic_length) - tes_index


def bin_coverage_by_absolute_position(coverage: np.ndarray, positions_relative_to_tes: np.ndarray, bin_edges: np.ndarray) -> np.ndarray:
    n_bins = len(bin_edges) - 1
    bin_indices = np.searchsorted(bin_edges, positions_relative_to_tes, side="right") - 1
    valid_mask = (bin_indices >= 0) & (bin_indices < n_bins)
    valid_indices = bin_indices[valid_mask]
    valid_coverage = coverage[valid_mask]
    binned_sums = np.bincount(valid_indices, weights=valid_coverage, minlength=n_bins)
    bin_counts = np.bincount(valid_indices, minlength=n_bins)
    result = np.zeros(n_bins)
    mask = bin_counts > 0
    result[mask] = binned_sums[mask] / bin_counts[mask]
    return result


def normalize_binned_coverage(binned: np.ndarray) -> np.ndarray:
    mean_val = np.mean(binned)
    return binned / mean_val if mean_val > 0 else binned


def apply_name_mapping(experiment: str, sample_num: int, sample_type: str) -> str:
    mapped_name = NAME_MAPPING.get(f"{experiment}-S{sample_num}", f"S{sample_num}")
    return f"{mapped_name}-{sample_type}"


def get_color_for_sample(experiment: str, sample_num: int) -> str:
    return COLOR_MAPPING.get(NAME_MAPPING.get(f"{experiment}-S{sample_num}", ""), "#808080FF")


def save_plot(fig: plt.Figure, stem: str) -> None:
    for ext in ("svg", "png"):
        out = OUT_DIR / f"{stem}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
        print(f"Saved: {out}")


def export_curve_tables(experiment: str, coverage_stats_by_type: dict[str, dict[int, dict]], ratio_stats: dict[int, dict]) -> None:
    rows = []
    for sample_type, stats in coverage_stats_by_type.items():
        for sample_num in sorted(stats):
            data = stats[sample_num]
            label = apply_name_mapping(experiment, sample_num, sample_type)
            for pos, mean, lo, hi in zip(bin_centers, data["means"], data["ci_lower"], data["ci_upper"]):
                rows.append(
                    {
                        "experiment": experiment,
                        "panel": "coverage",
                        "sample_type": sample_type,
                        "sample_num": sample_num,
                        "label": label,
                        "position_relative_to_tes_bp": pos,
                        "mean": mean,
                        "ci_lower": lo,
                        "ci_upper": hi,
                    }
                )

    for sample_num in sorted(ratio_stats):
        data = ratio_stats[sample_num]
        label = apply_name_mapping(experiment, sample_num, "SN")
        for pos, mean, lo, hi in zip(bin_centers, data["means"], data["ci_lower"], data["ci_upper"]):
            rows.append(
                {
                    "experiment": experiment,
                    "panel": "sn_lysate_ratio",
                    "sample_type": "ratio",
                    "sample_num": sample_num,
                    "label": label,
                    "position_relative_to_tes_bp": pos,
                    "mean": mean,
                    "ci_lower": lo,
                    "ci_upper": hi,
                }
            )

    pd.DataFrame(rows).to_csv(CSV_DIR / "suppl13b_spam_tes_coverage_curves.csv", index=False)


def compute_ratio_stats(experiment: str) -> dict[int, dict]:
    lysate_type = "Lysate" if experiment == "Benchmarking" else "L"
    ratios: dict[int, dict] = {}

    for sample_num in PANEL_B_SAMPLES:
        sn_key = (experiment, "SN", sample_num)
        lysate_key = (experiment, lysate_type, sample_num)
        if sn_key not in sample_groups or lysate_key not in sample_groups:
            continue

        sn_replicate_ids = sample_groups[sn_key]
        lysate_replicate_ids = sample_groups[lysate_key]
        sn_binned_list = []
        lysate_binned_list = []

        for transcript_id in preprocessed_data.keys():
            binned_dict = preprocessed_data[transcript_id]
            sn_rep_binned = [binned_dict[sid] for sid in sn_replicate_ids if sid in binned_dict]
            lys_rep_binned = [binned_dict[sid] for sid in lysate_replicate_ids if sid in binned_dict]
            if sn_rep_binned and lys_rep_binned:
                sn_binned_list.extend(sn_rep_binned)
                lysate_binned_list.extend(lys_rep_binned)

        if not sn_binned_list or not lysate_binned_list:
            continue

        min_len = min(
            min(len(arr) for arr in sn_binned_list),
            min(len(arr) for arr in lysate_binned_list),
        )
        sn_binned = np.array([arr[:min_len] for arr in sn_binned_list])
        lysate_binned = np.array([arr[:min_len] for arr in lysate_binned_list])

        sn_means = np.mean(sn_binned, axis=0)
        lysate_means = np.mean(lysate_binned, axis=0)
        ratio_means = np.divide(
            sn_means,
            lysate_means,
            where=lysate_means > 0,
            out=np.zeros_like(sn_means),
        )

        sn_std = np.std(sn_binned, axis=0)
        lysate_std = np.std(lysate_binned, axis=0)
        sn_sem = sn_std / np.sqrt(sn_binned.shape[0])
        lysate_sem = lysate_std / np.sqrt(lysate_binned.shape[0])
        ratio_sem = np.abs(ratio_means) * np.sqrt(
            np.divide(sn_sem**2, sn_means**2, where=sn_means > 0, out=np.zeros_like(sn_sem))
            + np.divide(
                lysate_sem**2,
                lysate_means**2,
                where=lysate_means > 0,
                out=np.zeros_like(lysate_sem),
            )
        )
        ci = 1.96 * ratio_sem

        ratios[sample_num] = {
            "means": ratio_means,
            "ci_lower": ratio_means - ci,
            "ci_upper": ratio_means + ci,
        }

    return ratios


print("Loading metadata and coverage database...")
lims_df = pd.read_csv(LIMS_CSV_PATH)
sample_metadata = {}
lims_df["parsed"] = lims_df["Sample_NameLIMS"].apply(parse_spam_sample_name)
for _, row in lims_df[lims_df["parsed"].notna()].iterrows():
    sample_metadata[row["Sample_Name"]] = row["parsed"]
conn = get_db_connection(DB_PATH)
cursor = conn.cursor()
cursor.execute("SELECT DISTINCT sample_id FROM transcript_coverage")
all_samples = [row[0] for row in cursor.fetchall()]
conn.close()
for sample_id in all_samples:
    if sample_id not in sample_metadata:
        parsed = parse_benchmarking_sample_name(sample_id)
        if parsed:
            sample_metadata[sample_id] = parsed

top_candidate_transcripts = select_top_transcripts(SALMON_QUANT_FILE, top_n=6000)
transcript_ids = list(top_candidate_transcripts)
all_coverage_data = {}
for i in range(0, len(transcript_ids), 500):
    batch_data = batch_load_transcripts(transcript_ids[i:i + 500])
    valid_batch_data = {tid: data for tid, data in batch_data.items() if data}
    all_coverage_data.update(valid_batch_data)
    if len(all_coverage_data) >= 5000:
        all_coverage_data = dict(list(all_coverage_data.items())[:5000])
        break

bin_size = 10
upstream_bp = 3000
bin_edges = np.arange(-upstream_bp, 1, bin_size)
bin_centers = bin_edges[:-1] + bin_size / 2
preprocessed_data = {}
for transcript_id, transcript_coverage_dict in all_coverage_data.items():
    exonic_length = len(next(iter(transcript_coverage_dict.values())))
    positions_relative_to_tes = get_coverage_positions_relative_to_tes(exonic_length)
    preprocessed_data[transcript_id] = {}
    for sample_id, coverage in transcript_coverage_dict.items():
        binned = bin_coverage_by_absolute_position(coverage, positions_relative_to_tes, bin_edges)
        preprocessed_data[transcript_id][sample_id] = normalize_binned_coverage(binned)

all_sample_ids = set()
for transcript_coverage in all_coverage_data.values():
    all_sample_ids.update(transcript_coverage.keys())
sample_groups = defaultdict(list)
for sample_id in all_sample_ids:
    if sample_id in sample_metadata:
        meta = sample_metadata[sample_id]
        sample_groups[(meta["Experiment"], meta["Type"], meta["Sample"])] .append(sample_id)

PANEL_B_SAMPLES = [1, 2, 4, 5]  # NTVE, HIV:sPAM p6, HIV-gag WT, mini-gag:PABP

coverage_stats = {}
for experiment in ["Benchmarking", "SPAM"]:
    coverage_stats[experiment] = {}
    sample_types = ["SN", "Lysate"] if experiment == "Benchmarking" else ["SN", "L"]
    for sample_type in sample_types:
        coverage_stats[experiment][sample_type] = {}
        for sample_num in PANEL_B_SAMPLES:
            group_key = (experiment, sample_type, sample_num)
            if group_key not in sample_groups:
                continue
            replicate_ids = sample_groups[group_key]
            all_binned = []
            for transcript_id in preprocessed_data:
                rep_binned = [preprocessed_data[transcript_id][sid] for sid in replicate_ids if sid in preprocessed_data[transcript_id]]
                if rep_binned:
                    all_binned.extend(rep_binned)
            if not all_binned:
                continue
            min_len = min(len(arr) for arr in all_binned)
            all_binned = np.array([arr[:min_len] for arr in all_binned])
            means = np.mean(all_binned, axis=0)
            sem = np.std(all_binned, axis=0) / np.sqrt(all_binned.shape[0])
            ci = 1.96 * sem
            coverage_stats[experiment][sample_type][sample_num] = {
                "means": means,
                "ci_lower": np.maximum(means - ci, 0),
                "ci_upper": means + ci,
            }

spam_ratio_stats = compute_ratio_stats("SPAM")


def plot_coverage_lineplot(experiment: str, sample_type: str, stem: str) -> None:
    fig, ax = plt.subplots(figsize=(PLOT_WIDTH, PLOT_HEIGHT))
    for sample_num in sorted(coverage_stats[experiment][sample_type]):
        data = coverage_stats[experiment][sample_type][sample_num]
        label = apply_name_mapping(experiment, sample_num, sample_type)
        color = get_color_for_sample(experiment, sample_num)
        ax.plot(bin_centers, data["means"], color=color, linewidth=2.5, label=label)
        ax.fill_between(bin_centers, data["ci_lower"], data["ci_upper"], color=color, alpha=0.2)
    ax.set_xlabel("Transcript Distance from TES (bp)")
    ax.set_ylabel("Mean Normalized Coverage")
    ax.set_title(f"{experiment}: {sample_type} coverage (TES-aligned)")
    ax.axvline(x=0, color="red", linestyle="--", linewidth=1, alpha=0.5, label="TES")
    ax.set_xlim(-upstream_bp, 0)
    ax.set_ylim(0, None)
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.legend(fontsize=7, loc="best", framealpha=0.95)
    fig.tight_layout()
    save_plot(fig, stem)
    plt.close(fig)


def plot_sn_lysate_ratio(experiment: str, ratio_stats: dict[int, dict], stem: str) -> None:
    fig, ax = plt.subplots(figsize=(PLOT_WIDTH, PLOT_HEIGHT))
    for sample_num in sorted(ratio_stats):
        data = ratio_stats[sample_num]
        label = apply_name_mapping(experiment, sample_num, "SN")
        color = get_color_for_sample(experiment, sample_num)
        ax.plot(bin_centers, data["means"], color=color, linewidth=2.5, label=label)
        ax.fill_between(bin_centers, data["ci_lower"], data["ci_upper"], color=color, alpha=0.2)
    ax.set_xlabel("Transcript Distance from TES (bp)")
    ax.set_ylabel("Coverage Ratio (SN / Lysate)")
    ax.set_title(f"{experiment}: SN/Lysate coverage ratio (TES-aligned)")
    ax.axhline(y=1, color="gray", linestyle=":", linewidth=1, alpha=0.5, label="Ratio = 1")
    ax.axvline(x=0, color="red", linestyle="--", linewidth=1, alpha=0.5, label="TES")
    ax.set_xlim(-upstream_bp, 0)
    ax.set_ylim(0, None)
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.legend(fontsize=7, loc="best", framealpha=0.95)
    fig.tight_layout()
    save_plot(fig, stem)
    plt.close(fig)


def _draw_coverage_panel(ax: plt.Axes, stats: dict[int, dict], experiment: str, sample_type: str, title: str) -> None:
    for sample_num in sorted(stats):
        data = stats[sample_num]
        label = apply_name_mapping(experiment, sample_num, sample_type)
        color = get_color_for_sample(experiment, sample_num)
        ax.plot(bin_centers, data["means"], color=color, linewidth=2.2, label=label)
        ax.fill_between(bin_centers, data["ci_lower"], data["ci_upper"], color=color, alpha=0.2)
    ax.set_xlabel("Transcript Distance from TES (bp)")
    ax.set_ylabel("Mean Normalized Coverage")
    ax.set_title(title, fontsize=8)
    ax.axvline(x=0, color="red", linestyle="--", linewidth=1, alpha=0.5)
    ax.set_xlim(-upstream_bp, 0)
    ax.set_ylim(0, None)
    ax.grid(True, alpha=0.3, linestyle="--")


def _draw_ratio_panel(ax: plt.Axes, stats: dict[int, dict], experiment: str, title: str) -> None:
    for sample_num in sorted(stats):
        data = stats[sample_num]
        label = apply_name_mapping(experiment, sample_num, "SN")
        color = get_color_for_sample(experiment, sample_num)
        ax.plot(bin_centers, data["means"], color=color, linewidth=2.2, label=label)
        ax.fill_between(bin_centers, data["ci_lower"], data["ci_upper"], color=color, alpha=0.2)
    ax.set_xlabel("Transcript Distance from TES (bp)")
    ax.set_ylabel("Coverage Ratio (SN / Lysate)")
    ax.set_title(title, fontsize=8)
    ax.axhline(y=1, color="gray", linestyle=":", linewidth=1, alpha=0.5)
    ax.axvline(x=0, color="red", linestyle="--", linewidth=1, alpha=0.5)
    ax.set_xlim(-upstream_bp, 0)
    ax.set_ylim(0, None)
    ax.grid(True, alpha=0.3, linestyle="--")


def plot_spam_combined(ratio_stats: dict[int, dict]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), sharex=True)
    _draw_coverage_panel(axes[0], coverage_stats["SPAM"]["L"], "SPAM", "L", "Lysate")
    _draw_coverage_panel(axes[1], coverage_stats["SPAM"]["SN"], "SPAM", "SN", "Supernatant")
    _draw_ratio_panel(axes[2], ratio_stats, "SPAM", "ratio")
    handles, labels = axes[2].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.99, 0.5), fontsize=7, framealpha=0.95)
    fig.tight_layout(rect=[0, 0, 0.92, 1])
    save_plot(fig, "suppl13b_spam_combined_absolute_tes")
    plt.close(fig)


plot_coverage_lineplot("SPAM", "L", "suppl13b_spam_l_lineplot_absolute_tes")
plot_coverage_lineplot("SPAM", "SN", "suppl13b_spam_sn_lineplot_absolute_tes")
plot_sn_lysate_ratio("SPAM", spam_ratio_stats, "suppl13b_spam_sn_lysate_ratio_absolute_tes")
plot_spam_combined(spam_ratio_stats)
export_curve_tables("SPAM", coverage_stats["SPAM"], spam_ratio_stats)
print("Done.")
