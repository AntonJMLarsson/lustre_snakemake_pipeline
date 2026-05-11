#!/usr/bin/env python3
"""
LUSTRE paired-end read simulator for pipeline demo.

Generates synthetic paired-end FASTQ reads for three L1 insertion types
across 100 simulated single cells.  Flanking genomic sequences are fetched
directly from the UCSC REST API (hg38) — no local reference genome is needed
to run this script.

Insertion types
---------------
  KR       – chr1:30568722 (primer site within reference L1HS at 30567813-30568750, + strand)
             Reference sequence here starts with the 31-bp primer + L1 poly-A tail,
             so Read 2 passes filter_bam.py AND maps concordantly (no soft-clip)
  KNR      – chr1:11239472 (Ewing-2010b-S4, + strand, from KNR_chr.bed)
             Read 2 = PRIMER + polyA + flank → large 5′ soft-clip → discordant
  Somatic  – chr2:100000000 (novel, absent from KR/KNR files)
             Same structure as KNR; present in 50 out of 100 cells

Read structure
--------------
  All insertions (KR, KNR, Somatic):
    Read 2:  TGCACATGTACCCTAAAACTTAGAGTATAAT  (31 bp primer) + [continuation]
    Read 1:  reverse-complement of genomic sequence ~insert_size downstream

  KR    – [continuation] = reference sequence after primer site in L1 (poly-A + flank);
           primer maps concordantly to the reference L1 → concordant pair
  KNR / Somatic – [continuation] = A×n (20-55 nt synthetic poly-A) + flanking genome;
           primer + poly-A soft-clip on alignment → discordant pair

Cells
-----
  100 cells; 12-bp barcodes embedded in query name as  sim########-BC:<barcode>
  move_barcode.py will later extract this to a BAM BC tag.
  Donors: first 50 cells → DONOR_A, last 50 → DONOR_B.
  Somatic insertion: present only in the first 50 cells (DONOR_A).

Output (written to demo/)
--------------------------
  demo_R1.fq.gz       paired FASTQ (Read 1)
  demo_R2.fq.gz       paired FASTQ (Read 2)
  samplesheet.csv     barcode → Donor table (index = barcode)
  expected_barcodes.txt
  cells_config.yaml   cbc_to_donor mapping (copy to config/ before running pipeline)
  DONOR_A_cells.txt   barcodes for DONOR_A cells (copy to config/)
  DONOR_B_cells.txt   barcodes for DONOR_B cells (copy to config/)
  demo_config.yaml    ready-to-use pipeline config (fill in fasta_file and tool paths)

Usage
-----
  cd <repo_root>
  python demo/simulate_lustre_demo.py

  Then follow the printed instructions to wire the outputs into the pipeline.

Requirements
------------
  Python ≥ 3.8, numpy
"""

import gzip
import json
import random
import sys
import urllib.request
from pathlib import Path
from typing import List, Optional, Set, Tuple

import numpy as np

# ─── Tuneable parameters ──────────────────────────────────────────────────────
PRIMER         = "TGCACATGTACCCTAAAACTTAGAGTATAAT"   # 31 bp
READ_LEN       = 150
N_CELLS        = 100
N_TAG_SITES    = 4     # normal tagmentation sites per cell per insertion
N_SPLIT_SITES  = 1     # short-insert split-read sites per cell (KNR/Somatic only)
READS_PER_TAG  = 50    # Poisson mean PCR copies per tagmentation site
SPLIT_INSERT_MAX = 100 # max insert size (bp of genomic flank) for split-read sites
FETCH_SIZE     = 2000  # bp fetched around each insertion site
SEED           = 42

# KR: primer site within L1HS chr1:30567813-30568750 (+).
# The 31-bp primer TGCACATGTACCCTAAAACTTAGAGTATAAT is at chr1:30568722 in hg38,
# followed by the L1 poly-A tail.  Read 2 uses this reference sequence directly,
# so the primer maps concordantly (no soft-clip → concordant pair).
KR_CHROM, KR_POS    = "chr1", 30568722

# KNR: chr1:11239472-11239522 from KNR_chr.bed (+).
KNR_CHROM, KNR_POS  = "chr1", 11239472

# Somatic: chr2:100000000 — confirmed >10 kb from any KR/KNR entry.
SOM_CHROM, SOM_POS  = "chr2", 100_000_000

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ─── Sequence utilities ───────────────────────────────────────────────────────
def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def fake_qual(length: int, phred: int = 40) -> str:
    return chr(phred + 33) * length


# ─── UCSC REST API ────────────────────────────────────────────────────────────
def fetch_ucsc(chrom: str, start: int, end: int, genome: str = "hg38") -> str:
    """Return upper-case sequence for chrom:start-end (0-based, half-open)."""
    url = (
        f"https://api.genome.ucsc.edu/getData/sequence"
        f"?genome={genome}&chrom={chrom}&start={start}&end={end}"
    )
    print(f"  Fetching {chrom}:{start:,}-{end:,} from UCSC …", flush=True)
    try:
        with urllib.request.urlopen(url, timeout=60) as resp:
            payload = json.loads(resp.read())
        return payload["dna"].upper()
    except Exception as exc:
        sys.exit(f"UCSC fetch failed for {chrom}:{start}-{end}: {exc}")


# ─── Cell barcode generation ──────────────────────────────────────────────────
def make_barcodes(n: int, length: int = 12) -> List[str]:
    rnd = random.Random(SEED)
    pool: Set[str] = set()
    while len(pool) < n:
        pool.add("".join(rnd.choices("ACGT", k=length)))
    return sorted(pool)


# ─── Read-pair generators ─────────────────────────────────────────────────────
def kr_pair(
    seq: str, insert: int, jitter: int
) -> Tuple[str, str]:
    """
    KR concordant pair.
    Read 2 = actual genomic sequence from within the reference L1 body.
    BWA will align this without a large soft-clip → concordant.
    """
    r2 = seq[jitter : jitter + READ_LEN]
    r1 = revcomp(seq[jitter + insert - READ_LEN : jitter + insert])
    return r1, r2


def discordant_pair(
    seq: str, insert: int, polya_len: int, jitter: int,
    l1_polya_seq: str = "",
) -> Tuple[str, str]:
    """
    KNR / Somatic discordant pair.

    Read 2 = PRIMER + poly-A + genomic flank (large primer+polyA soft-clip on
    alignment marks the pair as discordant).

    Read 1 depends on insert size relative to READ_LEN:

    Normal insert (insert >= READ_LEN):
        Read 1 = RC(flank[insert-READ_LEN : insert])  — entirely in 3' flank.

    Short insert (insert < READ_LEN):
        The fragment's right end is inside the poly-A tail, so Read 1 includes
        the L1 poly-A sequence placed at the 5' end of the reverse read
        (= rightmost genomic position).  This produces:

            cigar = [(READ_LEN-insert) S][insert M]   (soft-clip at cigar[0])

        analyze_split_reads.py detects this (cigartuples[0][0] == 4 for reverse
        reads).  get_softclip_sequence reverses the cigar[0] soft-clip, exposing
        the poly-A content which drives seq_pos_A > 7.

        l1_polya_seq must supply ≥ (READ_LEN-insert) bp of L1 poly-A sequence
        (starts just after the primer site in the reference L1).  Using the real
        L1 sequence (e.g. AAAAAAAAAT...) rather than pure poly-A prevents bwa
        from mapping the soft-clip region to a genomic poly-A repeat (e.g.
        chr12), ensuring the primary alignment stays at the insertion site.
    """
    geno_len = READ_LEN - len(PRIMER) - polya_len
    r2 = PRIMER + "A" * polya_len + seq[jitter : jitter + geno_len]

    if insert >= READ_LEN:
        r1 = revcomp(seq[jitter + insert - READ_LEN : jitter + insert])
    else:
        softclip_len = READ_LEN - insert
        # Short-insert split read.  Orientation matters:
        #
        # The original REVERSE molecule (5'→3' on − strand) =
        #   [flank RC (X bp, rightmost)]  +  [L1 poly-A RC (n bp, leftmost)]
        #
        # BAM stores the RC of this = [L1 poly-A (n bp, leftmost)] + [flank (X bp, rightmost)]
        # Cigar (left→right) = [n S][X M]  → cigar[0] = soft-clip ← read_is_split=True
        # get_softclip_sequence reverses cigar[0] bases → L1 poly-A → seq_pos_A high.
        l1_src = l1_polya_seq[:softclip_len] if l1_polya_seq else "A" * softclip_len
        r1 = revcomp(seq[jitter : jitter + insert]) + revcomp(l1_src)

    return r1, r2


# ─── Per-insertion read simulation ───────────────────────────────────────────
def simulate_insertion(
    seq: str,
    insertion_type: str,
    barcodes: List[str],
    active_barcodes: Optional[Set[str]],
    rng: np.random.Generator,
    counter_start: int = 0,
    l1_polya_seq: str = "",
) -> List[Tuple[str, str, str]]:
    """
    Returns list of (read_name, r1_seq, r2_seq).

    Biological model
    ----------------
    Each cell × insertion generates N_TAG_SITES unique tagmentation events.
    A tagmentation event = a specific genomic cut site downstream of the
    insertion, characterised by its insert size (distance from primer to cut).
    PCR amplifies each event ~READS_PER_TAG times.

    Read 2 always starts at the primer position (position 0 of seq) — no
    jitter — because the primer anneals to a fixed site in the L1.
    Read 1 varies only via the insert size (tagmentation site location).

    Parameters
    ----------
    active_barcodes
        Cells that carry this insertion.  ``None`` means every cell.
    """
    pairs: List[Tuple[str, str, str]] = []
    idx = counter_start
    for bc in barcodes:
        if active_barcodes is not None and bc not in active_barcodes:
            continue

        # Draw N_TAG_SITES normal-insert tagmentation sites for this cell.
        tag_inserts = [
            max(READ_LEN + 5, min(int(rng.lognormal(np.log(390), 0.23)), FETCH_SIZE - 55))
            for _ in range(N_TAG_SITES)
        ]

        # For KNR/Somatic, add N_SPLIT_SITES short-insert sites to produce split
        # reads at the poly-A/flank boundary (needed by prepare_insertion_table.py).
        if insertion_type != "KR":
            split_inserts = [
                int(rng.integers(30, SPLIT_INSERT_MAX + 1))
                for _ in range(N_SPLIT_SITES)
            ]
            tag_inserts = tag_inserts + split_inserts

        # For KNR/Somatic, pick one poly-A length per cell (L1 3′-end property).
        polya_len = int(rng.integers(20, 56))

        for insert in tag_inserts:
            n_copies = int(rng.poisson(READS_PER_TAG))
            for _ in range(n_copies):
                if insertion_type == "KR":
                    r1, r2 = kr_pair(seq, insert, 0)
                else:
                    r1, r2 = discordant_pair(seq, insert, polya_len, 0, l1_polya_seq)
                pairs.append((f"sim{idx:08d}-BC:{bc}", r1, r2))
                idx += 1
    return pairs


# ─── Writers ─────────────────────────────────────────────────────────────────
def write_fastqs(
    pairs: List[Tuple[str, str, str]], out_dir: Path
) -> Tuple[Path, Path]:
    r1_path = out_dir / "demo_R1.fq.gz"
    r2_path = out_dir / "demo_R2.fq.gz"
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for name, r1, r2 in pairs:
            f1.write(f"@{name}\n{r1}\n+\n{fake_qual(len(r1))}\n")
            f2.write(f"@{name}\n{r2}\n+\n{fake_qual(len(r2))}\n")
    return r1_path, r2_path


def write_samplesheet(barcodes: List[str], out_dir: Path) -> Path:
    path = out_dir / "samplesheet.csv"
    with open(path, "w") as f:
        f.write("barcode,Donor\n")
        for i, bc in enumerate(barcodes):
            donor = "DONOR_A" if i < N_CELLS // 2 else "DONOR_B"
            f.write(f"{bc},{donor}\n")
    return path


def write_expected_barcodes(barcodes: List[str], out_dir: Path) -> Path:
    path = out_dir / "expected_barcodes.txt"
    path.write_text("\n".join(barcodes) + "\n")
    return path


def write_cells_config(barcodes: List[str], out_dir: Path) -> Path:
    path = out_dir / "cells_config.yaml"
    with open(path, "w") as f:
        f.write("cbc_to_donor:\n")
        for i, bc in enumerate(barcodes):
            donor = "DONOR_A" if i < N_CELLS // 2 else "DONOR_B"
            f.write(f"  {bc}: {donor}\n")
    return path


def write_cell_lists(barcodes: List[str], out_dir: Path) -> Tuple[Path, Path]:
    donor_a = [bc for i, bc in enumerate(barcodes) if i < N_CELLS // 2]
    donor_b = [bc for i, bc in enumerate(barcodes) if i >= N_CELLS // 2]
    a_path = out_dir / "DONOR_A_cells.txt"
    b_path = out_dir / "DONOR_B_cells.txt"
    a_path.write_text("\n".join(donor_a) + "\n")
    b_path.write_text("\n".join(donor_b) + "\n")
    return a_path, b_path


def write_demo_config(
    out_dir: Path,
    r1_path: Path,
    r2_path: Path,
    ss_path: Path,
    bc_path: Path,
    repo_root: Path,
    somatic_site: Tuple[str, int],
) -> Path:
    resources = (repo_root / "resources").resolve()
    cfg_path = out_dir / "demo_config.yaml"
    som_chrom, som_pos = somatic_site
    with open(cfg_path, "w") as f:
        f.write(f"""\
# Demo config for the LUSTRE snakemake pipeline.
# Fill in fasta_file and tool paths before running.

project: demo

samplesheet: {ss_path.resolve()}
expected_barcodes: {bc_path.resolve()}

# Path to hg38 reference genome FASTA (must be BWA-indexed)
fasta_file: /path/to/hg38.fna   # UPDATE THIS

KNR_file: {resources / "nonref.collection.hg38.chr.L1.bed"}
old_KNR_file: {resources / "KNR_chr.bed"}
KR_file: {resources / "KR_chr.bed"}
KNR_remap: {resources / "nonref.collection.hg38.chr.L1.bed"}
KR_remap: {resources / "KR_chr.bed"}

r1: {r1_path.resolve()}
r2: {r2_path.resolve()}

threads: 4   # Adjust to available cores

samtools: samtools
bedtools: bedtools
bwa: bwa-mem2            # or full path to bwa-mem2 binary
fastp: fastp
demultiplex_fastq:       # not needed for this demo
split_barcoded_bam: /path/to/split_barcoded_bam  # UPDATE THIS

# Simulated somatic insertion: {som_chrom}:{som_pos:,} (present in DONOR_A cells)
custom_references:
  DONOR_A:
    file_list: ['results/regular_stats/demo_UNK_DONOR_A_stats.csv']
    cell_file: ['config/DONOR_A_cells.txt']
  DONOR_B:
    file_list: ['results/regular_stats/demo_UNK_DONOR_B_stats.csv']
    cell_file: ['config/DONOR_B_cells.txt']
""")
    return cfg_path


# ─── Main ─────────────────────────────────────────────────────────────────────
def main() -> None:
    repo_root = Path(__file__).resolve().parent.parent
    out_dir = repo_root / "demo"
    out_dir.mkdir(exist_ok=True)

    rng = np.random.default_rng(SEED)

    print("Fetching flanking sequences from UCSC hg38 REST API …")
    kr_seq  = fetch_ucsc(KR_CHROM,  KR_POS,  KR_POS  + FETCH_SIZE)
    knr_seq = fetch_ucsc(KNR_CHROM, KNR_POS, KNR_POS + FETCH_SIZE)
    som_seq = fetch_ucsc(SOM_CHROM, SOM_POS, SOM_POS + FETCH_SIZE)

    print("\nGenerating cell barcodes …")
    barcodes = make_barcodes(N_CELLS)
    somatic_cells = set(barcodes[: N_CELLS // 2])   # DONOR_A = first 50 cells

    # L1 poly-A sequence used as the soft-clip in short-insert split reads.
    # Taken from the reference L1 just after the primer site (kr_seq[31:]).
    # Using the actual L1 sequence (rather than pure poly-A) prevents bwa from
    # mapping the soft-clip to a genomic poly-A repeat (e.g. chr12).
    l1_polya_seq = kr_seq[len(PRIMER):]   # AAAAAAAAATAATAAATCAA...

    print("\nSimulating reads …")
    kr_pairs  = simulate_insertion(kr_seq,  "KR",      barcodes, None,          rng, 0)
    knr_pairs = simulate_insertion(knr_seq, "KNR",     barcodes, None,          rng, len(kr_pairs),
                                   l1_polya_seq=l1_polya_seq)
    som_pairs = simulate_insertion(som_seq, "Somatic", barcodes, somatic_cells, rng, len(kr_pairs) + len(knr_pairs),
                                   l1_polya_seq=l1_polya_seq)

    all_pairs = kr_pairs + knr_pairs + som_pairs

    # Shuffle so reads from different insertions are interleaved (realistic)
    perm = rng.permutation(len(all_pairs)).tolist()
    all_pairs = [all_pairs[i] for i in perm]

    print(f"\n  KR      reads: {len(kr_pairs):>7,}  (100 cells × {N_TAG_SITES} tag sites × ~{READS_PER_TAG} reads @ {KR_CHROM}:{KR_POS:,})")
    print(f"  KNR     reads: {len(knr_pairs):>7,}  (100 cells × {N_TAG_SITES} tag sites × ~{READS_PER_TAG} reads @ {KNR_CHROM}:{KNR_POS:,})")
    print(f"  Somatic reads: {len(som_pairs):>7,}  ( 50 cells × {N_TAG_SITES} tag sites × ~{READS_PER_TAG} reads @ {SOM_CHROM}:{SOM_POS:,})")
    print(f"  Total:         {len(all_pairs):>7,}")

    print("\nWriting output files …")
    r1_path, r2_path = write_fastqs(all_pairs, out_dir)
    ss_path  = write_samplesheet(barcodes, out_dir)
    bc_path  = write_expected_barcodes(barcodes, out_dir)
    write_cells_config(barcodes, out_dir)
    write_cell_lists(barcodes, out_dir)
    cfg_path = write_demo_config(
        out_dir, r1_path, r2_path, ss_path, bc_path, repo_root, (SOM_CHROM, SOM_POS)
    )

    print(
        "\n"
        + "─" * 60
        + "\nGenerated files:\n"
        + "\n".join(f"  {p.relative_to(repo_root)}" for p in sorted(out_dir.iterdir()))
    )
    print(
        "\n"
        "Next steps\n"
        "──────────\n"
        f"1. Edit  {cfg_path.relative_to(repo_root)}\n"
        "     • Set  fasta_file  to your hg38 FASTA (must be BWA-indexed)\n"
        "     • Set  split_barcoded_bam  to the compiled Rust binary path\n"
        "     • Adjust tool paths / threads as needed\n"
        "\n"
        "2. Copy config files into the pipeline config/ directory:\n"
        "     cp demo/demo_config.yaml    config/config.yaml\n"
        "     cp demo/cells_config.yaml   config/cells_config.yaml\n"
        "     cp demo/DONOR_A_cells.txt   config/DONOR_A_cells.txt\n"
        "     cp demo/DONOR_B_cells.txt   config/DONOR_B_cells.txt\n"
        "\n"
        "3. Dry-run:\n"
        "     snakemake --snakefile workflow/Snakefile --cores 4 -n\n"
        "\n"
        "4. Full run:\n"
        "     snakemake --snakefile workflow/Snakefile --cores 4\n"
        "\n"
        "Expected signals after the pipeline completes\n"
        "──────────────────────────────────────────────\n"
        f"  KR    insertion  → results/KR_bam/         (concordant reads, all 100 cells)\n"
        f"  KNR   insertion  → results/KNR_bam/        (discordant, all 100 cells)\n"
        f"  Somatic          → results/UNK_bam/        (discordant, 50 DONOR_A cells)\n"
        f"  Somatic site in  → results/insertion_tables/demo_UNK_DONOR_A_insertion_table.csv\n"
        "─" * 60
    )


if __name__ == "__main__":
    main()
