"""
Figure 4c, 4d, 4e — Neuron TP2 vs TP3 DESeq2 volcano plots
Standalone reproduction of Figure4/neuron_TP2_vs_TP3_DESeq2-Wald.ipynb
Outputs to: refactoring_roadmap/4cde_plots/

Manuscript mapping used here:
  Fig 4c = sample2
  Fig 4d = sample3
  Fig 4e = sample4

sample1 remains available in the upstream DESeq2 output directory, but is not
emitted by default because it does not appear in the final paper.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Setup paths
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT))

OUT_DIR = Path(__file__).parent / "4cde_plots"
OUT_DIR.mkdir(exist_ok=True)

# Configuration from notebook
PADJ_THRESHOLD = 0.05
LFC_THRESHOLD = 1.0
SPECIAL_GENES = {'Creb1', 'Mapk8', 'Dusp1', 'Cflar'}

def cm_to_inch(cm):
    return cm / 2.54

def create_volcano_plot(results_df, sample_num, output_dir, 
                       padj_threshold=PADJ_THRESHOLD, lfc_threshold=LFC_THRESHOLD,
                       xlim=None, ylim=None, n_top_labels=10,
                       special_genes=SPECIAL_GENES):
    """
    Create publication-quality volcano plot.
    Matches the logic in Figure4/neuron_TP2_vs_TP3_DESeq2-Wald.ipynb
    """
    df = results_df.copy()
    df = df.dropna(subset=['padj', 'log2FoldChange'])
    df['-log10_padj'] = -np.log10(np.maximum(df['padj'], 1e-300))
    
    df['category'] = 'Not significant'
    df.loc[(df['padj'] < padj_threshold) & (df['log2FoldChange'] > lfc_threshold), 'category'] = 'Upregulated'
    df.loc[(df['padj'] < padj_threshold) & (df['log2FoldChange'] < -lfc_threshold), 'category'] = 'Downregulated'
    
    # Select top genes to label: significant, |log2FC| > 1, and top by p-value
    significant_genes = df[df['padj'] < padj_threshold].copy()
    significant_genes['abs_log2FC'] = np.abs(significant_genes['log2FoldChange'])
    significant_genes = significant_genes[significant_genes['abs_log2FC'] > 1]
    significant_genes = significant_genes.sort_values(['padj', 'abs_log2FC'], ascending=[True, False])
    top_genes_to_label = significant_genes.head(n_top_labels)
    
    # Find special genes in the dataset
    special_genes_to_label = df[df['GeneName'].isin(special_genes)]
    
    # Use the same figure size as the notebook
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
    plt.rcParams['font.size'] = 10
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    fig, ax = plt.subplots(figsize=(cm_to_inch(8), cm_to_inch(10)), dpi=300)
    
    colors = {
        'Upregulated': '#E64B35',
        'Downregulated': '#4DBBD5',
        'Not significant': '#D3D3D3'
    }
    
    for category in ['Not significant', 'Downregulated', 'Upregulated']:
        subset = df[df['category'] == category]
        ax.scatter(subset['log2FoldChange'], subset['-log10_padj'],
                  c=colors[category], alpha=0.6 if category == 'Not significant' else 0.8,
                  s=8 if category == 'Not significant' else 12,
                  label=category, edgecolors='none', rasterized=False)
    
    # Add labels for top genes (regular color)
    for _, row in top_genes_to_label.iterrows():
        gene_name = row['GeneName'] if pd.notna(row['GeneName']) else row['Gene']
        if gene_name in special_genes:
            continue
        ax.annotate(gene_name,
                   xy=(row['log2FoldChange'], row['-log10_padj']),
                   xytext=(5, 5),
                   textcoords='offset points',
                   fontsize=7,
                   alpha=0.9,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='none'),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=0.5, alpha=0.5))
    
    # Add labels for special genes in a different color (orange)
    special_color = '#FF8C00'
    for _, row in special_genes_to_label.iterrows():
        gene_name = row['GeneName'] if pd.notna(row['GeneName']) else row['Gene']
        ax.annotate(gene_name,
                   xy=(row['log2FoldChange'], row['-log10_padj']),
                   xytext=(5, 5),
                   textcoords='offset points',
                   fontsize=7,
                   fontweight='bold',
                   alpha=0.95,
                   color=special_color,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8, edgecolor=special_color, linewidth=1.5),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0', lw=1.0, alpha=0.7, color=special_color))
    
    ax.axhline(y=-np.log10(padj_threshold), color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=lfc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(x=-lfc_threshold, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('log2 fold change (shrunken)', fontsize=10, fontweight='bold')
    ax.set_ylabel('-log10 (adjusted p-value)', fontsize=10, fontweight='bold')
    ax.set_title(f'Sample {sample_num}: TP2 vs TP3', fontsize=11, fontweight='bold', pad=10)
    
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    
    n_up = (df['category'] == 'Upregulated').sum()
    n_down = (df['category'] == 'Downregulated').sum()
    ax.text(0.02, 0.98, f'Up: {n_up}\nDown: {n_down}',
           transform=ax.transAxes, fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    ax.legend(fontsize=8, loc='upper right', frameon=True, edgecolor='black')
    ax.grid(True, alpha=0.15, linestyle=':', linewidth=0.5)
    
    plt.tight_layout(pad=0.3)
    
    for ext in ['svg', 'png', 'pdf']:
        save_path = output_dir / f'volcano_sample{sample_num}.{ext}'
        plt.savefig(save_path, format=ext, bbox_inches='tight')
    plt.close()
    print(f"Saved volcano plots for Sample {sample_num} to {output_dir}")

def main():
    results_dir = Path(__file__).parent / "4cde_output"
    all_results = {}
    
    for sample_num in [2, 3, 4]:
        csv_path = results_dir / f"sample{sample_num}" / "deseq2_results.csv"
        if csv_path.exists():
            all_results[sample_num] = pd.read_csv(csv_path)
        else:
            print(f"Warning: Results for Sample {sample_num} not found at {csv_path}")

    if not all_results:
        print("No results found. Exiting.")
        return

    # Calculate consistent axis limits as in notebook
    all_log2fc = []
    for df in all_results.values():
        df_clean = df.dropna(subset=['log2FoldChange', 'padj'])
        all_log2fc.extend(df_clean['log2FoldChange'].values)
    
    log2fc_99_9 = np.percentile(all_log2fc, 99.9)
    log2fc_min = -log2fc_99_9 * 1.5
    log2fc_max = log2fc_99_9 * 1.5
    neg_log10_padj_max = 8 # Matches notebook's manual setting
    
    volcano_xlim = (log2fc_min, log2fc_max)
    volcano_ylim = (0, neg_log10_padj_max)

    for sample_num, df in sorted(all_results.items()):
        create_volcano_plot(df, sample_num, OUT_DIR, xlim=volcano_xlim, ylim=volcano_ylim)

    # Notebook helper output intentionally disabled by default:
    # sample1 can still be rendered from refactoring_roadmap/4cde_output/sample1
    # if needed for debugging or comparison.

if __name__ == "__main__":
    main()
