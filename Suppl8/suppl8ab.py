"""
Supplementary Figure 8a, 8b — Dox vs No-Dox DESeq2 volcano plots.

Source notebook: Suppl8/DESeq2_Four_Lines_Dox_vs_NoDox-all.ipynb

Pipeline:
  Step 1 — Load Salmon quantification files → gene-level read-count matrix
            (protein-coding genes only, ntvetools GTF annotation)
  Step 2 — Run DESeq2 + apeglm shrinkage for each of the 4 lines via Docker
            (image: deseq2-with-shrinkage:latest)
  Step 3 — Load results CSVs → volcano plots (shared axes)

Default manuscript-faithful outputs go to suppl8ab_plots/:
  deseq2_results/<strain>/deseq2_<strain>_Dox_vs_NoDox.csv
  deseq2_results/<strain>/deseq2_<strain>_Dox_vs_NoDox_significant.csv
  deseq2_results/<strain>/normalized_counts_<strain>.csv
  deseq2_results/<strain>/size_factors_<strain>.csv
  deseq2_results/summary_statistics_all_lines.csv
  volcano_ISBM_64.svg / .png
  volcano_DE_WT.svg / .png
  

The notebook can emit additional helper outputs for the other two lines and a
summary barplot, but those are disabled by default because they are not shown
in the final paper.
"""

import matplotlib
matplotlib.use("Agg")

import gzip
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import seaborn as sns

try:
    import adjustText
except ImportError:
    print("Warning: adjustText not available, gene labels may overlap")
    adjustText = None

# ── Paths ──────────────────────────────────────────────────────────────────────
def find_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in (here.parent, here.parent.parent):
        if (candidate / "resources").exists():
            return candidate
    return here.parent.parent


REPO_ROOT = find_root()
RESOURCES = REPO_ROOT / "resources"
QUANTIF_DIR = RESOURCES / "quantfiles_filtered_pipeline" / "homo_sapiens"
GTF_CSV = RESOURCES / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"

OUT_DIR = Path(__file__).resolve().parent / "suppl8ab_plots"
DESEQ2_OUT = OUT_DIR / "deseq2_results"
OUT_DIR.mkdir(exist_ok=True)
DESEQ2_OUT.mkdir(exist_ok=True)

# ── Style ──────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 1,
})
sns.set_style("whitegrid")

# ── Sample groups ──────────────────────────────────────────────────────────────
SAMPLE_GROUPS = {
    "ISBM-64": {
        "no_dox": ["25L000325", "25L000326", "25L000327"],
        "dox":    ["25L000328", "25L000329", "25L000330"],
    },
    "ISBM-83": {
        "no_dox": ["25L000331", "25L000332", "25L000333"],
        "dox":    ["25L000334", "25L000335", "25L000336"],
    },
    "DE-WT": {
        "no_dox": ["25L000337", "25L000338"],
        "dox":    ["25L000340", "25L000341", "25L000342"],
    },
    "ISBM-17": {
        "no_dox": ["25L000343", "25L000344", "25L000345"],
        "dox":    ["25L000346", "25L000347", "25L000348"],
    },
}

ALL_SAMPLES = [s for g in SAMPLE_GROUPS.values() for grp in g.values() for s in grp]
MANUSCRIPT_STRAINS = ["ISBM-64", "DE-WT"]


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1 — Build gene-level read-count matrix
# ══════════════════════════════════════════════════════════════════════════════

def load_gtf_mappings():
    from ntvetools import load_gtf_df
    gtf = load_gtf_df(str(GTF_CSV))
    gtf_df = gtf["gtf_df"]

    gene_biotype = gtf_df.query("feature=='gene'").set_index("gene_id")["gene_biotype"].to_dict()
    gene_biotype["Gag_common_core"] = "protein_coding"

    t2g = defaultdict(str, gtf_df.set_index("transcript_id")["gene_id"].to_dict())
    g2name = defaultdict(str, gtf["gene_id_to_name"])

    # synthetic genes
    t2g["Gag_common_core"] = "Gag_common_core"
    g2name["Gag_common_core"] = "Gag_common_core"

    biotype_dd = defaultdict(lambda: "NA", gene_biotype)
    return t2g, g2name, biotype_dd


def load_sample(sf_path: Path, t2g: dict) -> tuple[dict, dict]:
    with gzip.open(sf_path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t")
    df["GeneID"] = [t2g[n.split(".")[0]] for n in df["Name"]]
    agg = df.groupby("GeneID")[["TPM", "NumReads"]].sum()
    return agg["TPM"].to_dict(), agg["NumReads"].to_dict()


def build_count_matrix(t2g, g2name, biotype_dd) -> tuple[pd.DataFrame, dict]:
    numreads_data = {}
    for sample_id in ALL_SAMPLES:
        sf = QUANTIF_DIR / f"{sample_id}.sf.gz"
        if not sf.exists():
            raise FileNotFoundError(f"Quantification file missing: {sf}")
        _, nr = load_sample(sf, t2g)
        numreads_data[sample_id] = nr
        print(f"  loaded {sample_id}")

    all_genes = sorted(g for d in numreads_data.values() for g in d if g)
    all_genes = sorted(set(all_genes))

    matrix = pd.DataFrame(
        {sid: [numreads_data[sid].get(g, 0) for g in all_genes] for sid in numreads_data},
        index=all_genes,
    )

    is_pc = [biotype_dd[g] == "protein_coding" for g in matrix.index]
    count_pc = matrix[is_pc].astype(int)
    print(f"  protein-coding genes: {count_pc.shape[0]}, samples: {count_pc.shape[1]}")

    # build gene_id → name mapping for DESeq2 R script
    g2name_subset = {g: g2name[g] if g2name[g] else g for g in count_pc.index}
    return count_pc, g2name_subset


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2 — Run DESeq2 via Docker for each strain
# ══════════════════════════════════════════════════════════════════════════════

DESEQ2_R_TEMPLATE = r"""#!/usr/bin/env Rscript
library(DESeq2)

cat("\n=== DESeq2: {strain} Dox vs No Dox ===\n")

setwd("/workspace/input")
output_dir <- "/workspace/output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

count_matrix    <- read.csv("count_matrix.csv",    row.names = 1, check.names = FALSE)
sample_metadata <- read.csv("sample_metadata.csv", stringsAsFactors = FALSE)
gene_id_to_name <- read.csv("gene_id_to_name.csv", row.names = 1)
gene_id_to_name_map <- setNames(gene_id_to_name$GeneName, rownames(gene_id_to_name))

sample_metadata <- sample_metadata[order(match(sample_metadata$Sample, colnames(count_matrix))), ]
count_matrix    <- count_matrix[, sample_metadata$Sample]

sample_metadata$Treatment <- factor(sample_metadata$Treatment, levels = c("No_Dox", "Dox"))
rownames(sample_metadata) <- sample_metadata$Sample

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(count_matrix),
  colData   = sample_metadata,
  design    = ~ Treatment
)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

res <- results(dds, contrast = c("Treatment", "Dox", "No_Dox"), alpha = 0.05)
cat(paste("  Significant genes (padj < 0.05):", sum(res$padj < 0.05, na.rm=TRUE), "\n"))

coef_name <- resultsNames(dds)[2]
tryCatch({
  res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")
  cat("  apeglm shrinkage applied\n")
}, error = function(e) {
  cat("  apeglm failed, using unshrunk results\n")
  res_shrunk <<- res
})

results_df <- as.data.frame(res_shrunk)
results_df$Gene     <- rownames(results_df)
results_df$GeneName <- gene_id_to_name_map[results_df$Gene]
results_df <- results_df[, c("Gene", "GeneName", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")]
results_df <- results_df[!is.na(results_df$padj), ]

strain_safe <- gsub("-", "_", "{strain}")
write.csv(results_df,
          file.path(output_dir, paste0("deseq2_", strain_safe, "_Dox_vs_NoDox.csv")),
          row.names = FALSE)
write.csv(results_df[results_df$padj < 0.05, ],
          file.path(output_dir, paste0("deseq2_", strain_safe, "_Dox_vs_NoDox_significant.csv")),
          row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
norm_df     <- data.frame(Gene = rownames(norm_counts), norm_counts)
write.csv(norm_df,
          file.path(output_dir, paste0("normalized_counts_", strain_safe, ".csv")),
          row.names = FALSE)

sf_df <- data.frame(Sample = names(sizeFactors(dds)), SizeFactor = as.numeric(sizeFactors(dds)))
write.csv(sf_df,
          file.path(output_dir, paste0("size_factors_", strain_safe, ".csv")),
          row.names = FALSE)

cat("=== Done ===\n")
"""


def run_deseq2(strain_name: str, count_matrix: pd.DataFrame, g2name: dict) -> bool:
    groups = SAMPLE_GROUPS[strain_name]
    samples = groups["no_dox"] + groups["dox"]

    meta_rows = []
    for s in groups["no_dox"]:
        meta_rows.append({"Sample": s, "Strain": strain_name, "Treatment": "No_Dox",
                          "Replicate": groups["no_dox"].index(s) + 1})
    for s in groups["dox"]:
        meta_rows.append({"Sample": s, "Strain": strain_name, "Treatment": "Dox",
                          "Replicate": groups["dox"].index(s) + 1})
    meta_df = pd.DataFrame(meta_rows)

    strain_count = count_matrix[samples]
    strain_safe = strain_name.replace("-", "_")
    strain_out = DESEQ2_OUT / strain_safe
    strain_out.mkdir(exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        strain_count.to_csv(tmp / "count_matrix.csv")
        meta_df.to_csv(tmp / "sample_metadata.csv", index=False)
        pd.DataFrame(
            {"GeneName": [g2name.get(g, g) for g in strain_count.index]},
            index=strain_count.index,
        ).to_csv(tmp / "gene_id_to_name.csv")

        r_script = DESEQ2_R_TEMPLATE.replace("{strain}", strain_name)
        (tmp / "run_deseq2.R").write_text(r_script)

        cmd = [
            "docker", "run", "--rm",
            "-v", f"{tmp.absolute()}:/workspace/input",
            "-v", f"{strain_out.absolute()}:/workspace/output",
            "deseq2-with-shrinkage:latest",
            "Rscript", "/workspace/input/run_deseq2.R",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        print(result.stdout)
        if result.returncode != 0:
            print(f"ERROR (Docker, {strain_name}):\n{result.stderr}")
            return False
        return True


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3 — Load results and plot
# ══════════════════════════════════════════════════════════════════════════════

def load_all_results() -> dict[str, pd.DataFrame | None]:
    all_results = {}
    for strain_name in SAMPLE_GROUPS:
        strain_safe = strain_name.replace("-", "_")
        csv = DESEQ2_OUT / strain_safe / f"deseq2_{strain_safe}_Dox_vs_NoDox.csv"
        if csv.exists():
            all_results[strain_name] = pd.read_csv(csv)
            print(f"  loaded {strain_name}: {len(all_results[strain_name])} genes")
        else:
            print(f"  MISSING: {csv}")
            all_results[strain_name] = None
    return all_results


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT_DIR / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(OUT_DIR / f"{stem}.png", dpi=1200, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {stem}.svg / .png")


def create_volcano_plot(
    results_df: pd.DataFrame, strain_name: str,
    xlim=None, ylim=None,
    padj_threshold=0.05, lfc_threshold=1.0, n_labels=20,
) -> plt.Figure:
    df = results_df.dropna(subset=["padj", "log2FoldChange"]).copy()
    df["-log10_padj"] = -np.log10(df["padj"] + 1e-300)

    df["significant"] = "Not Sig"
    df.loc[(df["padj"] < padj_threshold) & (df["log2FoldChange"] >  lfc_threshold), "significant"] = "Up"
    df.loc[(df["padj"] < padj_threshold) & (df["log2FoldChange"] < -lfc_threshold), "significant"] = "Down"

    colors = {"Up": "#d62728", "Down": "#1f77b4", "Not Sig": "#7f7f7f"}
    fig, ax = plt.subplots(figsize=(2.5, 2.5))

    for cat in ["Not Sig", "Down", "Up"]:
        sub = df[df["significant"] == cat]
        ax.scatter(sub["log2FoldChange"], sub["-log10_padj"],
                   c=colors[cat], alpha=0.6, s=10, label=cat,
                   edgecolors="none", rasterized=True)

    ax.axhline(y=-np.log10(padj_threshold), color="black", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(x= lfc_threshold,  color="black", linestyle="--", linewidth=0.8, alpha=0.7)
    ax.axvline(x=-lfc_threshold,  color="black", linestyle="--", linewidth=0.8, alpha=0.7)

    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)

    for spine in ax.spines.values():
        spine.set_linewidth(1)

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title(f"{strain_name}: Dox vs No Dox", pad=10)
    ax.tick_params(axis="both", which="major", labelsize=8)

    n_up   = (df["significant"] == "Up").sum()
    n_down = (df["significant"] == "Down").sum()
    ax.text(0.02, 0.98,
            f"Total DEGs: {n_up + n_down}\nUp: {n_up}\nDown: {n_down}",
            transform=ax.transAxes, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="gray"))

    if n_labels > 0 and adjustText is not None:
        df_sig = df[df["significant"].isin(["Up", "Down"])].copy()
        df_sig["score"] = df_sig["-log10_padj"] * df_sig["log2FoldChange"].abs()
        df_sig = df_sig.nlargest(n_labels, "score")
        texts = []
        for _, row in df_sig.iterrows():
            gene_name = row.get("GeneName", row["Gene"])
            if pd.isna(gene_name) or gene_name == "":
                gene_name = row["Gene"]
            gene_name = str(gene_name)[:30]
            texts.append(ax.text(row["log2FoldChange"], row["-log10_padj"],
                                 gene_name, fontweight="bold"))
        try:
            adjustText.adjust_text(
                texts,
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
                ax=ax,
                expand=(1.3, 1.3),
            )
        except Exception:
            pass

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3], loc="upper right", frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.2, linestyle=":", linewidth=0.5)
    fig.tight_layout()
    return fig


def create_summary_barplot(summary_df: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(summary_df))
    w = 0.35
    b1 = ax.bar(x - w / 2, summary_df["Upregulated"],   w, label="Upregulated",   color="#d62728", alpha=0.8)
    b2 = ax.bar(x + w / 2, summary_df["Downregulated"], w, label="Downregulated", color="#1f77b4", alpha=0.8)

    ax.set_xlabel("Strain", fontsize=10, fontweight="bold")
    ax.set_ylabel("Number of DEGs (padj < 0.05)", fontsize=10, fontweight="bold")
    ax.set_title("Differential Expression: Dox vs No Dox — All Lines",
                 fontsize=11, fontweight="bold", pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(summary_df["Strain"], rotation=0)
    ax.legend(loc="upper right", frameon=True, fancybox=True, shadow=True)
    ax.grid(True, alpha=0.2, linestyle=":", linewidth=0.5, axis="y")

    for rects in (b1, b2):
        for rect in rects:
            h = rect.get_height()
            if h > 0:
                ax.annotate(f"{int(h)}",
                            xy=(rect.get_x() + rect.get_width() / 2, h),
                            xytext=(0, 3), textcoords="offset points",
                            ha="center", va="bottom", fontsize=9, fontweight="bold")

    fig.tight_layout()
    return fig


# ── Main ───────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    # ── Step 1: build count matrix ─────────────────────────────────────────────
    print("\n=== Step 1: Loading quantification data ===")
    t2g, g2name, biotype_dd = load_gtf_mappings()
    count_matrix, g2name_subset = build_count_matrix(t2g, g2name, biotype_dd)

    # ── Step 2: run DESeq2 via Docker ──────────────────────────────────────────
    print("\n=== Step 2: Running DESeq2 (Docker) ===")
    for strain in SAMPLE_GROUPS:
        print(f"\n── {strain} ──")
        ok = run_deseq2(strain, count_matrix, g2name_subset)
        if not ok:
            print(f"  WARNING: DESeq2 failed for {strain}")

    # ── Step 3: load results + plot ────────────────────────────────────────────
    print("\n=== Step 3: Plotting ===")
    all_results = load_all_results()

    # compute global shared axis limits
    all_lfc, all_logp = [], []
    for df in all_results.values():
        if df is not None:
            df2 = df.dropna(subset=["padj", "log2FoldChange"]).copy()
            df2["-log10_padj"] = -np.log10(df2["padj"] + 1e-300)
            all_lfc.extend(df2["log2FoldChange"].tolist())
            all_logp.extend(df2["-log10_padj"].tolist())

    lfc_min, lfc_max = np.min(all_lfc), np.max(all_lfc)
    lfc_pad = (lfc_max - lfc_min) * 0.1
    xlim = (lfc_min - lfc_pad, lfc_max + lfc_pad)

    logp_max = np.max(all_logp)
    ylim = (0, logp_max * 1.1)

    print(f"  shared xlim={xlim}, ylim={ylim}")

    # manuscript volcano plots only
    for strain, df in all_results.items():
        if df is not None:
            if strain in MANUSCRIPT_STRAINS:
                print(f"  volcano {strain} …")
                fig = create_volcano_plot(df, strain, xlim=xlim, ylim=ylim)
                safe = strain.replace("-", "_")
                save(fig, f"volcano_{safe}")

    # summary table is still useful provenance for all analyzed lines
    summary_rows = []
    for strain, df in all_results.items():
        if df is not None:
            n_up   = ((df["padj"] < 0.05) & (df["log2FoldChange"] >  0)).sum()
            n_down = ((df["padj"] < 0.05) & (df["log2FoldChange"] < 0)).sum()
            summary_rows.append({"Strain": strain, "Upregulated": n_up, "Downregulated": n_down})

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(DESEQ2_OUT / "summary_statistics_all_lines.csv", index=False)

    # Additional notebook helper outputs are intentionally disabled by default.
    #
    # for strain, df in all_results.items():
    #     if df is not None and strain not in MANUSCRIPT_STRAINS:
    #         fig = create_volcano_plot(df, strain, xlim=xlim, ylim=ylim)
    #         safe = strain.replace("-", "_")
    #         save(fig, f"volcano_{safe}")
    #
    # save(create_summary_barplot(summary_df), "comparison_summary_barplot")

    print(f"\nDone. All outputs in {OUT_DIR}")
