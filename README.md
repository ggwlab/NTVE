# NTVE — Analysis Code and Data for the NTVE Manuscript

**Non-destructive transcriptomics by vesicular export (NTVE), a platform for multi-time-point monitoring of RNA expression dynamics in living cells.**
<img width="2048" height="1136" alt="Graphical_Abstract" src="https://github.com/user-attachments/assets/44a5a4f7-7f13-48b1-adf5-9f478b16c912" />

_NTVE-enabled cells package a part of their mRNA into VLPs of ~100nm diameter and secrete it into the surrounding compartment. This enables us to measure the transcriptome without having to sacrifice the cell._

---

This repository contains notebooks, helper tools, and bundled data snapshots
to reproduce all analyses and figures for the NTVE manuscript.

## 0. Full Text and Complete Repository
The bioRxiv preprint is at:

[https://www.biorxiv.org/content/10.1101/2024.11.11.622832v2](https://www.biorxiv.org/content/10.1101/2024.11.11.622832v2)


doi: [https://doi.org/10.1101/2024.11.11.622832](https://doi.org/10.1101/2024.11.11.622832) 

**IMPORTANT:** Due to file size restrictions, the indices for minimap and the resource folders are not in this repo. 
There is a complete mirror repo including those files at Zenodo:

[https://zenodo.org/records/18634314](https://zenodo.org/records/18634314)



Deeper step-by-step reproduction instructions live in
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

---

## 1. System requirements

| Requirement | Details |
|---|---|
| **Operating system** | Linux (tested on Ubuntu 22.04 and 24.04). macOS should work for the Python notebooks; Docker is required for all R-based workflows regardless of OS. |
| **Python** | 3.14 (via Conda / Mamba) |
| **Docker** | Required for DESeq2 and fgsea R workflows (R 4.4.0 inside containers) |
| **Disk space** | ~15 GB for the repository including bundled resources |
| **Hardware** | No non-standard hardware required |

### Key Python packages

pandas, numpy, scipy, matplotlib, seaborn, scanpy, scikit-learn, statsmodels,
adjustText, biopython, joblib, tqdm, natsort, umap-learn, anndata, leidenalg

All pinned versions are recorded in `resources/ntve_environment.yml`.

### Key R packages (inside Docker)

DESeq2, fgsea, apeglm, ashr, msigdbr, pheatmap, ggplot2, dplyr

### Optional

A Rust toolchain is only needed if you want to recompile
`resources/bam-filtering-rust` from source.

---

## 2. Installation guide

### 2a. Create the Python / Jupyter environment

```bash
conda env create -f resources/ntve_environment.yml   # ~5-15 min
conda activate ntve
```

### 2b. Install the local helper library

```bash
pip install -e resources/ntvetools
```

### 2c. Build the Docker images for R workflows

```bash
docker build -t deseq2-cubicsplines:latest   resources/docker_cubic_splines
docker build -t deseq2-with-shrinkage:latest resources/docker_DeSeq_with_Shrinkage
docker build -t fgsea:latest                 resources/docker_fgsea
```

Each image takes ~15 min to build (~45 min total).

**Typical total install time on a normal desktop: ~60-75 min.**

---

## 3. Demo

A good lightweight starting point is the enrichment-factor distribution
notebook, which loads only bundled resources and produces a single figure:

```bash
conda activate ntve
python Figure1/1g.py
```

- **Expected output**: SVG plot written to `Figure1/1g_plots/`.
- **Expected run time**: < 2 min on a normal desktop.


---

## 4. Instructions for use

### Running individual figure notebooks

There is a python file for each subfigure (or in rare cases a python file that creates two panels). Run them from the repo root with the `ntve` conda environment, for example `python Figure1/1e.py` or `python Suppl16/suppl16ab.py`.


### Docker-based R analyses

Several notebooks call Docker containers to run DESeq2 or fgsea.
Make sure the Docker images are built (step 2c) and that the Docker daemon is
running before executing those notebooks.

### Full reproduction workflow

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for the complete
ordered list of steps, expected inputs/outputs, and path-portability guidance.


---

## Repository layout

| Directory | Contents |
|---|---|
| `Figure1` – `Figure7` | Main figure analyses and outputs |
| `Suppl5` – `Suppl21` | Supplementary figure analyses and outputs |
| `resources/` | Shared data, Docker build contexts, `ntvetools`, and utility code |
| `docs/` | Detailed reproduction and portability documentation |

---

## Core resources bundled in this repository

### Gene annotations
- `resources/Homo_sapiens.GRCh38.108.gtf.gz` — Ensembl v108 (CC-BY-4.0)
- `0_RNA_Seq_Pipelines/indices_for_minimap/homo_sapiens_cDNA_ncRNA_additional_sequences.fa.gz` — Ensembl v108 (CC-BY-4.0)
- `0_RNA_Seq_Pipelines/indices_for_minimap/homo_sapiens_mus_musculus_cDNA_ncRNA_additional_sequences.gtf.gz` — Ensembl v108/v109 (CC-BY-4.0)
- `0_RNA_Seq_Pipelines/indices_for_minimap/mus_musculus_cDNA_ncRNA_additional_sequences.gtf.gz` — Ensembl v109 (CC-BY-4.0)
- `resources/merged_gtf_homosapiens_v108_musmusculus_v109.csv.gz` — derived from Ensembl v108/v109 (CC-BY-4.0)

### Pathway gene sets (MSigDB, CC-BY-4.0)
- `resources/h.all.v2024.1.Hs.json` — Hallmark gene sets
- `resources/c5.go.v2025.1.Hs.json` — GO gene sets
- `resources/c2.all.v2025.1.Hs.json` — Curated gene sets (includes KEGG and BioCarta subsets; see License section)

### Project metadata, counts, and quantifications
- `resources/Project_1716_lims_simplified.csv`
- `resources/Project_1716_lims_simplified_Spam2_without_missing.csv`
- `resources/harmonized_harmonized_gene_counts_rv_stranded.txt`
- `resources/comparison_harmonized_gene_counts.txt.gz`
- `resources/transcript_coverage.db`
- `resources/dominant_isoform_genes.csv.gz`
- `resources/data_dict.joblib`
- `resources/salmon_harmonized/*`

### Cell-cycle gene list
- `Figure6/regev_lab_cell_cycle_genes.txt` — S-phase and G2/M marker genes (Tirosh et al. 2016, via [scanpy_usage](https://github.com/scverse/scanpy_usage), BSD-3-Clause)

### Trilineage analysis pipeline
- `Figure6/pipeline/` — 5 notebooks that load quantification data and run DESeq2-based analyses (LRT, cubic splines, per-day comparisons) feeding into the interactive scatter plots. See `Figure6/pipeline/README.md` for execution order and Docker requirements.

### Interactive plots
- [`https://magro.codeberg.page/Trilinplots/`](https://magro.codeberg.page/Trilinplots/) — standalone HTML scatter plots bundling Plotly.js v3.3.1 (MIT, Plotly Inc.)

### Helper tools
- `resources/ntvetools/` — Python utility library (MIT)
- `resources/bam-filtering-rust/` — Rust BAM-filtering tool (MIT)


## License and attribution

### Code

All original code in this repository (notebooks, `ntvetools`, `bam-filtering-rust`,
and scripts) is released under the **MIT License**. See [`LICENSE`](LICENSE).

### Bundled data

| Resource | Source | License |
|---|---|---|
| Ensembl GTF / derived CSVs | Ensembl (EMBL-EBI) | CC-BY-4.0 |
| MSigDB Hallmark (`h.all`) | Liberzon et al., MSigDB | CC-BY-4.0 |
| MSigDB GO (`c5.go`) | MSigDB / Gene Ontology Consortium | CC-BY-4.0 |
| Cell-cycle gene list (`regev_lab_cell_cycle_genes.txt`) | Tirosh et al. 2016 / scanpy_usage (Theis Lab) | BSD-3-Clause |
| Plotly.js v3.3.1 (`plotly.min.js`) | Plotly, Inc. | MIT |
| Project data (counts, metadata, coverage, quantifications) | This study | MIT |

**Note on `c2.all.v2025.1.Hs.json`:** This file contains pathways not compatible with CC-BY-4.0. Please download it from MSigDB

### Attribution

If you use Ensembl annotations, please cite:

> Cunningham F, et al. "Ensembl 2022." *Nucleic Acids Research* (2022).

If you use MSigDB gene sets, please cite:

> Liberzon A, et al. "The Molecular Signatures Database Hallmark Gene Set
> Collection." *Cell Systems* (2015).
>
> Subramanian A, et al. "Gene set enrichment analysis: A knowledge-based
> approach." *PNAS* (2005).
