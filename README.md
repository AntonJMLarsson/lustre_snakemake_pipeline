# Configuration

Describe how to configure the workflow (using config.yaml and maybe additional files).
All of them need to be present with example entries inside of the config folder.

---

# System Requirements

## Software Dependencies

| Software | Version tested | Notes |
|---|---|---|
| [Snakemake](https://snakemake.readthedocs.io) | ≥ 6.0 (tested with 7.32) | Workflow manager |
| [Python](https://www.python.org/) | ≥ 3.8 | Required for all scripts |
| [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) | 2.2.1 | Read alignment; install via bioconda |
| [fastp](https://github.com/OpenGene/fastp) | ≥ 0.23 | Read trimming and QC |
| [samtools](https://www.htslib.org/) | ≥ 1.15 | BAM processing |
| [bedtools](https://bedtools.readthedocs.io) | ≥ 2.30 | Interval operations |
| [split_barcoded_bam](https://github.com/AntonJMLarsson/split_barcoded_bam) | latest | Cell barcode demultiplexing (Rust binary, must be compiled) |

### Python packages (main pipeline)

| Package | pip / conda name | Used by |
|---|---|---|
| `pysam` | `pysam` | BAM I/O throughout |
| `pandas` | `pandas` | Statistics and tables |
| `numpy` | `numpy` | Numerical operations |
| `joblib` | `joblib` | Parallel processing |
| `pyfaidx` | `pyfaidx` | FASTA indexing |
| `matplotlib` | `matplotlib` | QC plots |
| `seaborn` | `seaborn` | QC plots |
| `mgzip` | `mgzip` | Multi-threaded gzip I/O |

### Python packages (auxiliary scripts only)

These are only needed if using scripts outside the main Snakemake workflow (e.g. `demultiplex_file.py`, `rupture_and_polyA.py`):

| Package | pip install name |
|---|---|
| `Levenshtein` | `python-levenshtein` |
| `fastqandfurious` | `fastq-and-furious` |
| `mgzip` | `mgzip` |
| `networkx` | `networkx` |
| `ruptures` | `ruptures` |

## Operating Systems

Tested on:
- Linux (Ubuntu 20.04, 22.04)
- macOS (tested for development only; production use on Linux recommended)

## Required Hardware

- **CPU:** Multi-core strongly recommended. The pipeline is configured via `threads` in `config.yaml`; 48 threads were used in production.
- **RAM:** ≥ 64 GB recommended for whole-genome datasets.
- **Storage:** ≥ 500 GB free disk space for intermediate BAM files and results (varies with dataset size).
- **No non-standard hardware required.** AVX2 support is recommended for bwa-mem2 performance but not mandatory (alternative binaries available).

---

# Installation Guide

## Instructions

1. **Clone the repository:**
   ```bash
   git clone https://github.com/AntonJMLarsson/lustre_snakemake_pipeline.git
   cd lustre_snakemake_pipeline
   ```

2. **Create the conda environment** with all required tools and Python packages:
   ```bash
   conda create -n lustre -c bioconda -c conda-forge \
     snakemake samtools bedtools fastp bwa-mem2 \
     "python=3.10" pysam pandas numpy joblib \
     pyfaidx matplotlib seaborn -y
   conda activate lustre
   ```

3. **Install remaining Python packages via pip:**
   ```bash
   pip install mgzip pyfaidx
   ```
   For auxiliary scripts (demultiplexing, transduction analysis):
   ```bash
   pip install fastq-and-furious python-levenshtein mgzip networkx ruptures
   ```

4. **Install split_barcoded_bam** (Rust tool — requires [Rust/Cargo](https://rustup.rs/)):
   ```bash
   git clone https://github.com/AntonJMLarsson/split_barcoded_bam.git
   cd split_barcoded_bam && cargo build --release
   cd ..
   ```
   Note the path to the compiled binary (`split_barcoded_bam/target/release/split_barcoded_bam`) — it is required in `config/config.yaml`.

5. **Download reference files** (hg38 genome, L1/KR/KNR BED files) and update paths in `config/config.yaml`.

6. **Create `config/cells_config.yaml`** — this file is required by the Snakefile:
   ```bash
   # For a new project, populate with your barcode→donor mapping:
   #   cbc_to_donor:
   #     BARCODE1: DONOR_A
   #     BARCODE2: DONOR_B
   # For the demo dataset:
   cp demo/cells_config.yaml config/cells_config.yaml
   ```
   See `config/cells_config.yaml` in this repository for a documented example.

## Typical Install Time

On a standard desktop computer with a fast internet connection: **30–60 minutes**, including dependency downloads and compilation of the Rust binary.

---

# Demo

## Running on Simulated Data

A self-contained simulation script is provided at [demo/simulate_lustre_demo.py](../demo/simulate_lustre_demo.py). It fetches real hg38 flanking sequences from the UCSC REST API and generates synthetic paired-end FASTQ files — no local reference genome is required for this step.

The simulation creates 100 cells across three insertion types:

| Insertion | Site | Cells | Expected pipeline output |
|---|---|---|---|
| KR (reference L1HS) | chr1:30,568,600 | all 100 | `results/KR_bam/` (concordant) |
| KNR (non-reference L1) | chr1:11,239,472 | all 100 | `results/KNR_bam/` (discordant) |
| Somatic (novel) | chr2:100,000,000 | 50 (DONOR_A) | `results/UNK_bam/` → insertion table |

**Step 1 – Generate demo FASTQ and config files (~2 min, internet required):**
```bash
cd <repo_root>
python demo/simulate_lustre_demo.py
```

**Step 2 – Set genome and tool paths in the generated config:**
```bash
# Fill in fasta_file and split_barcoded_bam in demo/demo_config.yaml
```

**Step 3 – Copy config files and run:**
```bash
cp demo/demo_config.yaml    config/config.yaml
cp demo/cells_config.yaml   config/cells_config.yaml
cp demo/DONOR_A_cells.txt   config/DONOR_A_cells.txt
cp demo/DONOR_B_cells.txt   config/DONOR_B_cells.txt

snakemake --snakefile workflow/Snakefile --cores 4
```

## Expected Output

The pipeline produces the following key output files under `results/`:

| Path | Description |
|---|---|
| `results/fastq_files_trimmed/` | Trimmed FASTQ files |
| `results/mapped/` | BWA-MEM2 aligned BAM |
| `results/KNR_bam/`, `results/KR_bam/`, `results/UNK_bam/` | Per-cell BAMs sorted into KNR / KR / UNK categories |
| `results/KNR_insertions_final.bed`, `results/UNK_insertions_final.bed` | Called insertion sites |
| `results/regular_stats/` | Per-donor read count matrices and split-read statistics |
| `results/insertion_tables/` | Insertion tables (CSV) with per-insertion metadata |
| `results/detailed_stats/` | Detailed per-insertion statistics including polyA signal |

## Expected Run Time (Demo)

| Stage | Time |
|---|---|
| `simulate_lustre_demo.py` (UCSC fetch + FASTQ generation) | ~2 min |
| Full pipeline on simulated dataset (~7,400 read pairs, 4 cores) | 15–45 min |
| Full production dataset (hundreds of millions of reads, 48 cores) | 12–48 hours |

---

# Instructions for Use

## Running the Pipeline on Your Data

### 1. Prepare inputs

- **FASTQ files:** Concatenated R1/R2 FASTQ files from a LUSTRE sequencing run. Use `demultiplex_fastq` (separate repository) to embed cell barcodes into the query names prior to running this pipeline.
- **Reference genome:** hg38 FASTA (indexed with bwa-mem2).
- **BED files:** KNR, KR, and nonref L1 BED files (examples provided in `resources/`).
- **Samplesheet:** CSV linking barcodes to samples/donors.
- **Expected barcodes:** Text file with one expected cell barcode per line.

### 2. Configure `config/config.yaml`

Update all paths and parameters:

```yaml
project: <project_name>           # Unique identifier for your run
samplesheet: config/<project>_samplesheet.csv
expected_barcodes: config/<project>_expected_barcodes.txt
fasta_file: /path/to/hg38.fna
KNR_file: /path/to/nonref.collection.hg38.chr.L1.bed
old_KNR_file: /path/to/KNR_chr.bed
KR_file: /path/to/KR_chr.bed
KNR_remap: /path/to/KNR_full.bed
KR_remap: /path/to/KR_full.bed
r1: /path/to/reads_1.fq.gz
r2: /path/to/reads_2.fq.gz
threads: 48                        # Adjust to available cores

# Tool paths
samtools: samtools
bedtools: bedtools
bwa: /path/to/bwa-mem2.avx
fastp: fastp
split_barcoded_bam: /path/to/split_barcoded_bam

# Per-donor custom remapping configuration
custom_references:
  DONOR_A:
    file_list: ['results/regular_stats/<project>_UNK_DONOR_A_stats.csv']
    cell_file: ['config/DONOR_A_cells.txt']
```

### 3. Configure `config/cells_config.yaml`

This file is **required** — the Snakefile loads it at startup (`configfile: "config/cells_config.yaml"`). Provide a `cbc_to_donor` mapping from cell barcodes to donor identifiers. This is used during custom remapping to assign each cell to its donor-specific reference.

```yaml
cbc_to_donor:
  ACGTACGTACGT: DONOR_A
  TTGGCCAATTGG: DONOR_A
  GCTAGCTAGCTA: DONOR_B
```

The donor IDs must match the keys under `custom_references` in `config/config.yaml`. A blank `cbc_to_donor: {}` template is provided in `config/cells_config.yaml`.

### 4. Run the workflow

```bash
# Dry run to verify the DAG:
snakemake --snakefile workflow/Snakefile --cores <N> -n

# Full run:
snakemake --snakefile workflow/Snakefile --cores <N>
```

### 5. Pipeline stages

The workflow runs the following stages in order:

1. **Trimming & mapping** (`trim_and_map.smk`) – fastp trimming followed by bwa-mem2 alignment.
2. **Preprocessing** (`preprocess.smk`) – barcode extraction, per-cell demultiplexing, and sorting into KR / KNR / UNK read categories.
3. **Calling** (`calling.smk`) – peak grouping and artefact filtering to produce `*_insertions_final.bed` files.
4. **Reporting (regular)** (`reporting_regular.smk`) – per-insertion split-read analysis and per-cell read count matrices.
5. **Custom remapping** (`custom_remapping.smk`) – extraction of polyA-rich partial reads, construction of donor-specific references, and re-alignment for increased sensitivity.
6. **Reporting (remapped)** (`reporting_remapped.smk`) – quality statistics after custom mapping.
7. **Table generation** (`make_tables.smk`) – final insertion tables and detailed per-insertion statistics.
