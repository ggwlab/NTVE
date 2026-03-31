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
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

OUT_DIR = Path(__file__).parent / "6c_plots"
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
    
    # 7. UMAP Analysis (Fit on SNs, Project Lysates)
    print("Running UMAP analysis and projection...")
    adata_sn = adata[~adata.obs['lysate']].copy()
    adata_lysate = adata[adata.obs['lysate']].copy()
    
    sc.pp.log1p(adata_sn)
    sc.tl.pca(adata_sn)
    sc.pp.neighbors(adata_sn)
    sc.tl.umap(adata_sn)
    
    if len(adata_lysate) > 0:
        sc.pp.log1p(adata_lysate)
        sc.tl.ingest(adata_lysate, adata_sn, obs='lineage')
        
        # Combine
        lysate_mask = adata.obs['lysate'].values
        X_umap_combined = np.zeros((len(adata), 2))
        X_umap_combined[~lysate_mask] = adata_sn.obsm['X_umap']
        X_umap_combined[lysate_mask] = adata_lysate.obsm['X_umap']
        adata.obsm['X_umap'] = X_umap_combined
    else:
        adata.obsm['X_umap'] = adata_sn.obsm['X_umap']

    # 8. Custom Plotting
    print("Generating custom UMAP plots...")
    plot_custom_umap(adata, adata_sn)

def plot_custom_umap(adata, adata_sn):
    plt.rcParams['font.size'] = 8
    plt.rcParams['font.family'] = 'sans-serif'
    
    cm_to_inch = 1 / 2.54
    subplot_width_inch = 4 * cm_to_inch
    subplot_height_inch = 5 * cm_to_inch
    left_margin, right_margin, bottom_margin, top_margin, spacing = 0.65, 0.6, 0.5, 0.3, 0.5
    fwidth = left_margin + 2 * subplot_width_inch + spacing + right_margin
    fheight = bottom_margin + subplot_height_inch + top_margin
    
    markers = {6: '^', 9: 's', 16: 'o'}
    lineage_base_colors = {
        'Ec': np.array([0.8, 0.2, 0.2]),    # Red
        'En': np.array([0.2, 0.6, 0.2]),    # Green
        'Mes': np.array([0.2, 0.4, 0.8]),   # Blue
        'Ipsc': np.array([0.6, 0.6, 0.6])   # Gray
    }
    
    unique_tps = sorted(adata.obs['timepoint'].unique())
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
    axw, axh = subplot_width_inch / fwidth, subplot_height_inch / fheight
    ax1 = fig.add_axes([left_margin / fwidth, bottom_margin / fheight, axw, axh])
    ax2 = fig.add_axes([(left_margin + subplot_width_inch + spacing) / fwidth, bottom_margin / fheight, axw, axh])
    
    # Left Panel: SNs only
    coords_sn = adata_sn.obsm['X_umap']
    for idx, row in adata_sn.obs.iterrows():
        x, y = coords_sn[adata_sn.obs.index.get_loc(idx)]
        c = get_color(lineage_base_colors[row['lineage'][:3] if row['lineage'] != 'Ipsc' else 'Ipsc'], row['timepoint'])
        ax1.scatter(x, y, marker=markers[row['replicate']], s=get_size(row['timepoint']), 
                   c=[c], alpha=0.6, edgecolors='black', linewidths=0.5)
    
    # Right Panel: Combined
    coords_all = adata.obsm['X_umap']
    for idx, row in adata.obs.iterrows():
        x, y = coords_all[adata.obs.index.get_loc(idx)]
        lin_key = row['lineage'][:3] if row['lineage'] != 'Ipsc' else 'Ipsc'
        c = get_color(lineage_base_colors[lin_key], row['timepoint'])
        s = get_size(row['timepoint'])
        if row['lysate']:
            ax2.scatter(x, y, marker=markers[row['replicate']], s=s, facecolors='none', 
                       edgecolors=c, alpha=0.6, linewidths=1.5)
        else:
            ax2.scatter(x, y, marker=markers[row['replicate']], s=s, c=[c], 
                       alpha=0.6, edgecolors='black', linewidths=0.5)
    
    for ax, title in zip([ax1, ax2], ['SNs only', 'SNs + projected lysates']):
        ax.set_xlabel('UMAP1'); ax.set_ylabel('UMAP2'); ax.set_title(title)
        ax.set_xlim(coords_sn[:, 0].min()-1, coords_sn[:, 0].max()+1)
        ax.set_ylim(coords_sn[:, 1].min()-1, coords_sn[:, 1].max()+1)

    # Legends
    legend_elements_replicates = [
        Line2D([0], [0], marker=m, color='w', markerfacecolor='gray', markersize=6, 
               markeredgecolor='black', markeredgewidth=0.5, label=f'Clone {c}')
        for c, m in markers.items()
    ]
    ax1.legend(handles=legend_elements_replicates, loc='best', fontsize=6, title='Replicate', title_fontsize=6)
    
    legend_elements_lineage = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=lineage_base_colors['Ec'], markersize=6, label='Ectoderm'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=lineage_base_colors['En'], markersize=6, label='Endoderm'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=lineage_base_colors['Mes'], markersize=6, label='Mesoderm'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=lineage_base_colors['Ipsc'], markersize=6, label='iPSC'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=6, label='SN'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='none', markeredgecolor='gray', markersize=6, markeredgewidth=1.5, label='Lysate'),
    ]
    ax2.legend(handles=legend_elements_lineage, loc='best', fontsize=6, title='Lineage & Type', title_fontsize=6)

    plt.savefig(OUT_DIR / "trilineage_umap_two_panels.png", dpi=300, bbox_inches='tight')
    plt.savefig(OUT_DIR / "trilineage_umap_two_panels.svg", format='svg', bbox_inches='tight')
    print(f"Saved UMAP plots to {OUT_DIR}")

if __name__ == "__main__":
    main()
