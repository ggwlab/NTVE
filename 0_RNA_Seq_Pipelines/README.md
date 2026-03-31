# RNA-seq Harmonized Submission Pipelines 

## 1. What these pipelines do

There are two independent workflows:

1. `Minimap2 -> name-sort -> filter_bam_pairs -> Salmon` (SLURM jobs)
- Locates and pairs FASTQ files under a user-provided raw-data root.
- For each sample, processes **three transcriptome references**:
  - `homo_sapiens`
  - `mus_musculus`
  - `HS_MM_combined`
- Deletes intermediate BAM files from scratch storage.
- Keeps only final Salmon outputs.

2. `STAR -> featureCounts` on harmonized references
- Locates and pairs FASTQ files under a user-provided raw-data root.
- Aligns each sample with STAR to a harmonized STAR genome index.
- Writes BAM outputs under `star_alignments/Project_<PROJECT>/star_<SPECIES>/`.
- Runs `featureCounts` twice on those BAMs:
  - unstranded (`-s 0`)
  - reversely stranded (`-s 2`)

## 2. Runtime and software requirements

- OS: Linux (HPC environment with SLURM required for the minimap pipeline as originally run).
- Python: 3.10.x (notebook metadata used 3.10.11 and 3.10.14).
- Required Python packages:
  - none (only Python standard library is used in the standalone script)
- Command-line tools:
  - `STAR`
  - `minimap2` (notebook path used `minimap2-2.30_x64-linux`)
  - `samtools` (notebook path used `samtools-1.16.1`)
  - `salmon`
  - `featureCounts` (Subread, notebook path used `subread-2.0.3`)
  - `filter_bam_pairs` (Rust binary; set `FILTER_PATH` to your local binary location if it is not in `PATH`)
- Scheduler tools:
  - Optional for STAR/minimap submissions: `sbatch`, `squeue`

## 3. Required input files and directories

## 3.1 Raw FASTQs

The minimap workflow recursively scans a base directory for:
- `*.fastq.gz`
- `*.fq.gz`

R1/R2 pairing is derived from filename suffixes. Supported patterns:
- `_R1` / `_R2`
- `_1` / `_2`
- `.R1` / `.R2`

Sample names are inferred from the shared prefix (everything before the R1 token).

## 3.2 Transcriptome reference files (FASTA)

The minimap workflow expects exactly these files in `INDEX_DIR`:
- `homo_sapiens_cDNA_ncRNA_additional_sequences.fa.gz`
- `mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz`
- `homo_sapiens_mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz`

These three files are expected in:
- `indices_minimap/`

## 3.3 Inputs for direct STAR -> featureCounts workflow

Required:
- Raw FASTQ root containing `*.fastq.gz` / `*.fq.gz` files
- STAR genome index directory (for harmonized reference)
- Annotation GTF for featureCounts (typically merged harmonized GTF)

Default BAM location produced by STAR stage:
- `./star_alignments/Project_<PROJECT>/star_<SPECIES>/*.bam`

## 3.4 Local index assets included in this repository

STAR-related source files:
- `indices_for_star/gag_common_core_and_mGL.fa`
- `indices_for_star/synthetic_ISBM.gtf`
- `indices_for_star/build_instructions.txt`

Minimap-related FASTA references:
- `indices_minimap/harmonized_additional.fa.gz` (included in this repository)
- `indices_minimap/homo_sapiens_cDNA_ncRNA_additional_sequences.fa.gz` (included in this repository)
- `indices_minimap/mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz` (included in this repository)
- `indices_minimap/homo_sapiens_mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz` (included in this repository)

## 4. Default parameters reproduced from notebooks

## 4.1 Minimap/Salmon workflow defaults

- Minimap preset: `-ax sr`
- Threads per tool stage: `8`
- Keep properly paired reads during SAM->BAM: `samtools view -f 0x2`
- Filtering:
  - complexity cutoff: `0.8`
  - minimum mapped bases: `80`
- Salmon quant:
  - `-l A`
  - `--numBootstraps 100`
  - `--seqBias --gcBias`
  - `-p 8`
- SLURM options:
  - `--time=20:00:00 --mem=16G -n 8 --job-name harmonized_rna`
- Scratch directory for temporary BAMs:
  - `/scratch/<USER>/rna_harmonized`

## 4.2 Direct STAR -> featureCounts workflow defaults

STAR stage:
- `--readFilesCommand zcat`
- Threads: `8`
- Output BAM type: `--outSAMtype BAM SortedByCoordinate`
- Quant mode: `--quantMode GeneCounts TranscriptomeSAM`
- Notebook SLURM options (when using scheduler): `--time=03:00:00 --mem=110G -n 8 --job-name STAR_RNA`

- Grouping: `-g gene_id`
- Feature type: `-t exon`
- Paired-end mode: `-p`
- Threads: `-T 32`
- Strandedness runs:
  - unstranded: `-s 0`
  - reversely stranded: `-s 2`

## 5. Building the STAR and minimap reference indices

This section documents how to rebuild the references from source inputs and the local synthetic additions.

## 5.1 Rebuild STAR harmonized genome index (`indices_for_star`)

The original instructions are in:
- `indices_for_star/build_instructions.txt`

Required external source files (not bundled here):
- Human genome FASTA: `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
- Mouse genome FASTA: `Mus_musculus.GRCm39.dna.primary_assembly.fa.gz`
- Human GTF: `Homo_sapiens.GRCh38.108.gtf.gz`
- Mouse GTF: `Mus_musculus.GRCm39.109.gtf.gz`

Run:

```bash
cd indices_for_star

# Example variable setup (edit to your local file locations)
HS_DNA=/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
MM_DNA=/path/to/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
HS_GTF=/path/to/Homo_sapiens.GRCh38.108.gtf.gz
MM_GTF=/path/to/Mus_musculus.GRCm39.109.gtf.gz
STAR_BIN=/path/to/STAR

# Prefix chromosome/contig headers to keep species separable
gzip -dc "$HS_DNA" | sed 's/^>/>hs_/' > human_modified.fa
gzip -dc "$MM_DNA" | sed 's/^>/>mm_/' > mouse_modified.fa

# Prefix first-column seqnames in GTF (same rule as original notes)
gzip -dc "$HS_GTF" | sed 's/^[0-9XYM]/hs_&/' > human_modified.gtf
gzip -dc "$MM_GTF" | sed 's/^[0-9XYM]/mm_&/' > mouse_modified.gtf

# Merge with repository-provided synthetic references
cat human_modified.fa mouse_modified.fa gag_common_core_and_mGL.fa > merged.fa
cat human_modified.gtf mouse_modified.gtf synthetic_ISBM.gtf > merged.gtf

# Build STAR genome index
mkdir -p h38m39108112
"$STAR_BIN" \
  --runMode genomeGenerate \
  --genomeFastaFiles merged.fa \
  --genomeDir h38m39108112 \
  --sjdbGTFfile merged.gtf \
  --runThreadN 32
```

Outputs used downstream:
- `merged.fa`
- `merged.gtf`
- STAR index directory (example name from original run): `h38m39108112/`

## 5.2 Rebuild minimap transcriptome references (`indices_minimap`)

`indices_minimap/` stores the final gzipped transcriptome FASTAs used by the minimap workflow. The notebook aligns directly to these `.fa.gz` files (no prebuilt `.mmi` required).

Required external source files (not bundled here):
- Human cDNA FASTA (Ensembl release-matched)
- Human ncRNA FASTA (Ensembl release-matched)
- Mouse cDNA FASTA (Ensembl release-matched)
- Mouse ncRNA FASTA (Ensembl release-matched)

Repository-provided additional transcripts:
- `indices_minimap/harmonized_additional.fa.gz`

Run:

```bash
cd indices_minimap

# Example variable setup (edit to your local file locations)
HS_CDNA=/path/to/Homo_sapiens.GRCh38.cdna.all.fa.gz
HS_NCRNA=/path/to/Homo_sapiens.GRCh38.ncrna.fa.gz
MM_CDNA=/path/to/Mus_musculus.GRCm39.cdna.all.fa.gz
MM_NCRNA=/path/to/Mus_musculus.GRCm39.ncrna.fa.gz

# Build gzipped FASTA references via gzip-stream concatenation
cat "$HS_CDNA" "$HS_NCRNA" harmonized_additional.fa.gz \
  > homo_sapiens_cDNA_ncRNA_additional_sequences.fa.gz

cat "$MM_CDNA" "$MM_NCRNA" harmonized_additional.fa.gz \
  > mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz

cat "$HS_CDNA" "$HS_NCRNA" "$MM_CDNA" "$MM_NCRNA" harmonized_additional.fa.gz \
  > homo_sapiens_mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz
```

Optional pre-indexing for faster startup (not required by notebook logic):

```bash
minimap2 -d homo_sapiens_cDNA_ncRNA_additional_sequences.mmi \
  homo_sapiens_cDNA_ncRNA_additional_sequences.fa.gz
minimap2 -d mus_musculus_cDNA_ncRNA_additional_sequences.mmi \
  mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz
minimap2 -d homo_sapiens_mus_musculus_cDNA_ncRNA_additional_sequences.mmi \
  homo_sapiens_mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz
```

Validation:

```bash
for f in *_additional_sequences.fa.gz; do
  echo "===== $f"
  gzip -dc "$f" | sed -n '1,4p'
done
```

## 6. Standalone execution: Minimap -> filter -> Salmon

Run from repository root (or adapt paths explicitly):

```bash
python3 - <<'PY'
import os
import re
import shutil
from pathlib import Path
import subprocess
import tempfile
import json

# ------------------------------
# Configuration (EDIT AS NEEDED)
# ------------------------------
BASE_RAW_DATA = Path("/path/to/raw_fastq_root")
BASE_OUTPUT = Path("./outputs/harmonized_submission")
INDEX_DIR = Path("./indices_minimap")

MINIMAP_PATH = "minimap2"
SAMTOOLS_PATH = "samtools"
SALMON_PATH = "salmon"
FILTER_PATH = "filter_bam_pairs"

COMPLEXITY_CUTOFF = 0.8
MIN_MAPPED_BASES = 80
SBATCH_OPTIONS = "--time=20:00:00 --mem=16G -n 8 --job-name harmonized_rna"
SCRATCH_DIR = "/scratch/tmp/rna_harmonized"

PAIR_PATTERNS = [
    re.compile(r"(.+)_R1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
    re.compile(r"(.+)_1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
    re.compile(r"(.+)\.R1(\.fastq\.gz|\.fq\.gz)$", re.IGNORECASE),
]


def get_all_transcriptome_indices(index_dir):
    return {
        'homo_sapiens': index_dir / "homo_sapiens_cDNA_ncRNA_additional_sequences.fa.gz",
        'mus_musculus': index_dir / "mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz",
        'HS_MM_combined': index_dir / "homo_sapiens_mus_musculus_cDNA_ncRNA_additional_sequences.fa.gz"
    }


def construct_pipeline_command(sample_name, r1_fastq, r2_fastq, output_base,
                               minimap_path, samtools_path, salmon_path, filter_path,
                               transcriptome_indices, complexity_cutoff=0.8, min_mapped=0,
                               scratch_dir="/scratch/tmp/rna_harmonized"):
    pipeline_cmds = []
    for index_name, transcriptome_path in transcriptome_indices.items():
        scratch_sample_dir = Path(scratch_dir) / index_name / sample_name
        salmon_outdir = output_base / "salmon" / index_name / sample_name / "quant"
        salmon_outdir.parent.mkdir(parents=True, exist_ok=True)

        raw_bam = scratch_sample_dir / f"{sample_name}_raw.bam"
        namesorted_bam = scratch_sample_dir / f"{sample_name}_namesorted.bam"
        filtered_bam = scratch_sample_dir / f"{sample_name}_filtered.bam"

        cmd = (
            f"echo 'Processing {index_name}...' && "
            f"mkdir -p {scratch_sample_dir} && "
            f"{minimap_path} -ax sr -t 8 {transcriptome_path} {r1_fastq} {r2_fastq} | "
            f"{samtools_path} view -bS -f 0x2 - > {raw_bam} && "
            f"{samtools_path} sort -n -@ 8 -o {namesorted_bam} {raw_bam} && "
            f"rm {raw_bam} && "
            f"{filter_path} -i {namesorted_bam} -o {filtered_bam} -c {complexity_cutoff} -m {min_mapped} && "
            f"rm {namesorted_bam} && "
            f"{salmon_path} quant -t {transcriptome_path} -l A -a {filtered_bam} "
            f"-o {salmon_outdir} -p 8 --numBootstraps 100 --seqBias --gcBias && "
            f"gzip -f {salmon_outdir}/quant.sf && "
            f"rm {filtered_bam} && "
            f"rm -rf {scratch_sample_dir}"
        )
        pipeline_cmds.append(cmd)

    return " && ".join(pipeline_cmds) + " && echo 'All indices processed successfully!'"


def submit_sbatch_job(pipeline_cmd, sample_name, sbatch_options):
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".sh") as tmp:
        tmp.write("#!/bin/bash\n")
        tmp.write(f"#SBATCH {sbatch_options}\n")
        tmp.write(f"# Sample: {sample_name}\n")
        tmp.write("set -e\n")
        tmp.write(f"{pipeline_cmd}\n")
        tmp_name = tmp.name

    try:
        result = subprocess.run(['sbatch', tmp_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=10)
        if result.returncode == 0:
            return result.stdout.split()[-1]
        return None
    finally:
        os.remove(tmp_name)


def derive_r2_filename(r1_name):
    if "_R1" in r1_name:
        return r1_name.replace("_R1", "_R2", 1)
    if "_1" in r1_name:
        return r1_name.replace("_1", "_2", 1)
    if ".R1" in r1_name:
        return r1_name.replace(".R1", ".R2", 1)
    return None


def sample_name_from_r1(r1_name):
    for rx in PAIR_PATTERNS:
        m = rx.match(r1_name)
        if m:
            return m.group(1)
    return None


def discover_fastq_pairs(base_raw_data):
    fastq_files = sorted(base_raw_data.glob("**/*.fastq.gz")) + sorted(base_raw_data.glob("**/*.fq.gz"))
    by_name = {p.name: p for p in fastq_files}

    pairs = {}
    unmatched_r1 = []
    for p in fastq_files:
        sample = sample_name_from_r1(p.name)
        if sample is None:
            continue
        r2_name = derive_r2_filename(p.name)
        if r2_name and r2_name in by_name:
            pairs[sample] = (p, by_name[r2_name])
        else:
            unmatched_r1.append(p.name)
    return pairs, unmatched_r1


def binary_exists(cmd):
    return Path(cmd).exists() or shutil.which(cmd) is not None


# ----------
# Pre-flight
# ----------
required_paths = [BASE_RAW_DATA, INDEX_DIR]
missing_paths = [str(p) for p in required_paths if not Path(p).exists()]
if missing_paths:
    raise FileNotFoundError("Missing required paths:\n" + "\n".join(missing_paths))

required_bins = [MINIMAP_PATH, SAMTOOLS_PATH, SALMON_PATH, FILTER_PATH]
missing_bins = [b for b in required_bins if not binary_exists(b)]
if missing_bins:
    raise FileNotFoundError("Missing required tools:\n" + "\n".join(missing_bins))

# ----------
# Discover FASTQ pairs from input directory
# ----------
sample_to_fastq, unmatched_r1 = discover_fastq_pairs(BASE_RAW_DATA)
print(f"Matched FASTQ pairs: {len(sample_to_fastq)}")
print(f"Unmatched R1 files: {len(unmatched_r1)}")

# ----------
# Build and submit jobs
# ----------
transcriptome_indices = get_all_transcriptome_indices(INDEX_DIR)
for name, path in transcriptome_indices.items():
    if not path.exists():
        raise FileNotFoundError(f"Missing transcriptome index: {name} -> {path}")

job_commands = []
for sample_name, (r1, r2) in sample_to_fastq.items():
    cmd = construct_pipeline_command(
        sample_name, r1, r2, BASE_OUTPUT,
        MINIMAP_PATH, SAMTOOLS_PATH, SALMON_PATH, FILTER_PATH,
        transcriptome_indices, COMPLEXITY_CUTOFF, MIN_MAPPED_BASES,
        scratch_dir=SCRATCH_DIR,
    )
    job_commands.append((sample_name, cmd))

submitted_jobs = []
failed_jobs = []
for sample_name, cmd in job_commands:
    job_id = submit_sbatch_job(cmd, sample_name, SBATCH_OPTIONS)
    if job_id:
        submitted_jobs.append((sample_name, job_id))
        print(f"SUBMITTED {sample_name} -> {job_id}")
    else:
        failed_jobs.append(sample_name)
        print(f"FAILED {sample_name}")

BASE_OUTPUT.mkdir(parents=True, exist_ok=True)
tracking = {
    'submitted_jobs': submitted_jobs,
    'failed_jobs': failed_jobs,
    'unmatched_r1_files': unmatched_r1,
    'total_submitted': len(submitted_jobs),
    'total_failed': len(failed_jobs),
    'output_directory': str(BASE_OUTPUT),
    'complexity_cutoff': COMPLEXITY_CUTOFF,
    'min_mapped_bases': MIN_MAPPED_BASES,
}
with open(BASE_OUTPUT / 'job_tracking.json', 'w') as f:
    json.dump(tracking, f, indent=2)

print('\nDone.')
print(f"Job tracking file: {BASE_OUTPUT / 'job_tracking.json'}")
PY
```

## 7. Standalone execution: Direct STAR -> featureCounts workflow

- `run_harmonized_submission_direct.py`

It supports three stages:
- `--mode star`: STAR only
- `--mode featurecounts`: featureCounts only
- `--mode all`: STAR then featureCounts

### 7.1 Run complete workflow locally (portable default)

```bash
python3 run_harmonized_submission_direct.py \
  --mode all \
  --workflow-root . \
  --project NTVE \
  --species harmonized \
  --raw-data-root /path/to/raw_fastq_root \
  --star-index /path/to/harmonized_star_index \
  --annotation-gtf /path/to/merged.gtf \
  --star-bin STAR \
  --featurecounts-bin featureCounts
```

### 7.2 Run STAR with SLURM, then featureCounts

```bash
python3 run_harmonized_submission_direct.py \
  --mode all \
  --executor slurm \
  --wait-for-star \
  --workflow-root . \
  --project NTVE \
  --species harmonized \
  --raw-data-root /path/to/raw_fastq_root \
  --star-index /path/to/harmonized_star_index \
  --annotation-gtf /path/to/merged.gtf \
  --slurm-options "--time=03:00:00 --mem=110G -n 8 --job-name STAR_RNA"
```

### 7.3 Run featureCounts only on existing BAMs

Using default STAR output location:

```bash
python3 run_harmonized_submission_direct.py \
  --mode featurecounts \
  --workflow-root . \
  --project NTVE \
  --species harmonized \
  --annotation-gtf /path/to/merged.gtf \
  --featurecounts-bin featureCounts
```

Using explicit BAM glob:

```bash
python3 run_harmonized_submission_direct.py \
  --mode featurecounts \
  --annotation-gtf /path/to/merged.gtf \
  --bam-glob "/path/to/bams/*.bam"
```

Use `--dry-run` in any mode to print discovered inputs/commands without execution.

## 8. Expected output structure

## 8.1 Minimap/Salmon workflow

Final outputs under:
- `${BASE_OUTPUT}/salmon/homo_sapiens/<SAMPLE>/quant/`
- `${BASE_OUTPUT}/salmon/mus_musculus/<SAMPLE>/quant/`
- `${BASE_OUTPUT}/salmon/HS_MM_combined/<SAMPLE>/quant/`

Expected key files per sample/index:
- `quant.sf.gz` (compressed by pipeline)
- Salmon auxiliary files/logs

Tracking file:
- `${BASE_OUTPUT}/job_tracking.json`

Scratch path `${SCRATCH_DIR}` is temporary and should be cleaned per sample when jobs finish successfully.

## 8.2 Direct STAR -> featureCounts workflow

STAR outputs:
- `star_alignments/Project_<PROJECT>/star_<SPECIES>/*_Aligned.sortedByCoord.out.bam`
- STAR log files next to each sample output prefix

featureCounts outputs:
- `featurecounts_results/harmonized_harmonized_gene_counts_unstranded.txt`
- `featurecounts_results/harmonized_harmonized_gene_counts_rv_stranded.txt`

Tracking file:
- `direct_workflow_tracking.json`

## 9. Validation checks

After minimap workflow submission:

```bash
# Submitted/failed summary
jq '.total_submitted, .total_failed' /path/to/output/job_tracking.json

# Check queue state (first <=20 jobs example)
squeue -j $(jq -r '.submitted_jobs[:20][][1]' /path/to/output/job_tracking.json | paste -sd, -)
```

After minimap jobs finish:

```bash
# Count quant files by index
find /path/to/output/salmon/homo_sapiens -name quant.sf.gz | wc -l
find /path/to/output/salmon/mus_musculus -name quant.sf.gz | wc -l
find /path/to/output/salmon/HS_MM_combined -name quant.sf.gz | wc -l
```

After direct workflow:

```bash
find star_alignments/Project_<PROJECT>/star_<SPECIES> -name '*.bam' | wc -l
ls -lh featurecounts_results/*.txt
head -n 2 featurecounts_results/harmonized_harmonized_gene_counts_unstranded.txt
head -n 2 featurecounts_results/harmonized_harmonized_gene_counts_rv_stranded.txt
jq '.mode, .star_pairs_discovered, .featurecounts_bam_count' direct_workflow_tracking.json
```

## 10. Portability checklist (required when moving to a different system)

Update these path groups before running:
- Raw data root: `BASE_RAW_DATA`
- Output root: `BASE_OUTPUT`
- Index root: `INDEX_DIR`
- STAR index root: `--star-index`
- Local index asset dirs (if rebuilding): `indices_for_star/`, `indices_minimap/`
- Tool binaries: `STAR`, `featureCounts`, `MINIMAP_PATH`, `SAMTOOLS_PATH`, `SALMON_PATH`, `FILTER_PATH`
- Workflow labels: `--project`, `--species`
- Scratch root: `SCRATCH_DIR`
- SLURM partition and resources: minimap `SBATCH_OPTIONS`, STAR `--slurm-options`

Also verify:
- `sbatch` available on execution host.
- Scratch path exists and is writable from worker nodes.
- Output path is writable.
- FASTQ filenames follow one of the supported R1 naming patterns.

