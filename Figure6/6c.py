"""
Figure 6c — Trilineage UMAP (SNs + projected Lysates)
Standalone reproduction of Figure6/UMAP.ipynb
Outputs to: refactoring_roadmap/6c_plots/
"""

import joblib
import anndata as ad
import pandas as pd
import numpy as np
import re
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from pathlib import Path
import sys
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Setup paths
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = SCRIPT_DIR / "6c_plots"
OUT_DIR.mkdir(exist_ok=True)

# Data paths
DATA_DICT_PATH = ROOT / "resources" / "data_dict.joblib"
BARCODE_PATH = ROOT / "resources" / "Sample_barcode.csv"
GTF_PATH = ROOT / "resources" / "merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz"

def main():
    print("Loading data...")
    # 1. Load data
    data_dict = joblib.load(DATA_DICT_PATH)
    sample_df = pd.read_csv(BARCODE_PATH, sep=";")
    
    # 2. Map barcodes to sample names
    barcode_to_name = sample_df.set_index("Barcode")["Sample_name"].to_dict()
    for key in ["numreads_df", "tpm_df_protein_coding_genes"]:
        if key in data_dict:
            data_dict[key] = data_dict[key].rename(columns=barcode_to_name)
    
    # 3. Filter high-read protein-coding columns
    total_reads = data_dict["numreads_df"].sum(axis=0)
    high_read_cols = total_reads[total_reads > 0].index
    tpm_df = data_dict["tpm_df_protein_coding_genes"][high_read_cols]
    
    # 4. Create AnnData
    print("Preparing AnnData...")
    adata = ad.AnnData(X=tpm_df.values.T)
    adata.var_names = tpm_df.index
    adata.obs_names = tpm_df.columns
    
    # 5. Filter for Trilineage + iPSC patterns
    pattern = r"^(SN_(mes|en|ec)_d\d+_\d+|L_(mes|en|ec|IPSC)_\d+)$"
    adata = adata[adata.obs_names.str.match(pattern)].copy()
    
    # 6. Parse metadata
    lineage = []
    timepoint = []
    replicate = []
    lysate = []
    
    for sample_name in adata.obs_names:
        if sample_name.startswith("SN_"):
            parts = sample_name.split("_")
            lineage.append(parts[1].title())
            timepoint.append(int(parts[2].replace("d", "")))
            replicate.append(int(parts[3]))
            lysate.append(False)
        elif sample_name.startswith("L_"):
            parts = sample_name.split("_")
            lin = parts[1]
            lineage.append("Ipsc" if lin == "IPSC" else lin.title())
            timepoint.append(0)
            replicate.append(int(parts[2]))
            lysate.append(True)
        else:
            lineage.append(None); timepoint.append(None); replicate.append(None); lysate.append(None)
            
    adata.obs["lineage"] = lineage
    adata.obs["timepoint"] = timepoint
    adata.obs["replicate"] = replicate
    adata.obs["lysate"] = lysate
    
    # 7. Gene symbol conversion (matches notebook)
    print("Converting to gene symbols...")
    from ntvetools import load_gtf_df
    gtf_dict = load_gtf_df()
    ensg2symbol = gtf_dict["gene_id_to_name"]
    adata.var["gene_symbols"] = adata.var_names.map(ensg2symbol)
    adata.var_names = adata.var["gene_symbols"].fillna(
        pd.Series(adata.var_names, index=adata.var.index)
    )
    adata.var_names_make_unique()

    # 8. UMAP Analysis — fit on SNs only, no log1p (matches notebook)
    print("Running UMAP analysis...")
    adata_sn = adata[~adata.obs['lysate']].copy()

    sc.tl.pca(adata_sn)
    sc.pp.neighbors(adata_sn)
    sc.tl.umap(adata_sn)

    # 9. Custom Plotting
    print("Generating custom UMAP plots...")
    plot_custom_umap(adata, adata_sn)

def plot_custom_umap(adata, adata_sn):
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'sans-serif'

    cm_to_inch = 1 / 2.54
    subplot_width_inch = 4 * cm_to_inch
    subplot_height_inch = 5 * cm_to_inch
    left_margin, right_margin, bottom_margin, top_margin = 0.65, 0.6, 0.5, 0.3
    fwidth = left_margin + subplot_width_inch + right_margin
    fheight = bottom_margin + subplot_height_inch + top_margin

    markers = {6: '^', 9: 's', 16: 'o'}
    lineage_base_colors = {
        'Ec': np.array([0.8, 0.2, 0.2]),
        'En': np.array([0.2, 0.6, 0.2]),
        'Mes': np.array([0.2, 0.4, 0.8]),
        'Ipsc': np.array([0.6, 0.6, 0.6]),
    }

    unique_tps = sorted(adata_sn.obs['timepoint'].unique())
    tp_min, tp_max = min(unique_tps), max(unique_tps)

    def get_color(base_color, tp):
        norm_tp = 0.5 if tp_max == tp_min else (tp - tp_min) / (tp_max - tp_min)
        sat = 0.3 + norm_tp * 0.7
        white = np.array([1.0, 1.0, 1.0])
        return np.clip(white + sat * (base_color - white), 0, 1)

    def get_size(tp):
        norm_tp = 0.5 if tp_max == tp_min else (tp - tp_min) / (tp_max - tp_min)
        return 5 + norm_tp * 65

    fig = plt.figure(figsize=(fwidth, fheight), facecolor='white')
    axw = subplot_width_inch / fwidth
    axh = subplot_height_inch / fheight
    ax = fig.add_axes([left_margin / fwidth, bottom_margin / fheight, axw, axh])

    coords_sn = adata_sn.obsm['X_umap']
    for idx, row in adata_sn.obs.iterrows():
        x, y = coords_sn[adata_sn.obs.index.get_loc(idx)]
        lin_key = row['lineage'][:3] if row['lineage'] != 'Ipsc' else 'Ipsc'
        c = get_color(lineage_base_colors[lin_key], row['timepoint'])
        ax.scatter(x, y, marker=markers[row['replicate']], s=get_size(row['timepoint']),
                   c=[c], alpha=0.6, edgecolors='black', linewidths=0.5)

    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.set_title('Differentiation of 3 hiPSC clones', fontsize=8)

    # Timepoint colorbar
    import matplotlib.cm as mcm
    sm = plt.cm.ScalarMappable(cmap='Greys', norm=plt.Normalize(vmin=tp_min, vmax=tp_max))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.4, aspect=10, pad=0.02)
    cbar.set_label('Timepoint', fontsize=7)
    cbar.set_ticks([tp_min, tp_max])

    # Legend: one entry per lineage × clone combination
    legend_elements = []
    for lin, label in [('Ec', 'ectoderm'), ('Mes', 'mesoderm'), ('En', 'endoderm')]:
        base = lineage_base_colors[lin]
        for clone, marker in markers.items():
            c = get_color(base, tp_max)
            legend_elements.append(
                Line2D([0], [0], marker=marker, color='w', markerfacecolor=c,
                       markersize=6, markeredgecolor='black', markeredgewidth=0.5,
                       label=f'{label} clone {clone}')
            )
    ax.legend(handles=legend_elements, loc='best', fontsize=5, ncol=3)

    fig.savefig(OUT_DIR / "trilineage_umap.png", dpi=300, bbox_inches='tight')
    fig.savefig(OUT_DIR / "trilineage_umap.svg", format='svg', bbox_inches='tight')
    plt.close(fig)
    print(f"Saved UMAP plots to {OUT_DIR}")

if __name__ == "__main__":
    main()
