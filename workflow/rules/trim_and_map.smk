configfile: "config.yaml"

(SAMPLES,) = glob_wildcards("fastq_lustre_cat_files/{sample}_read_1.fq.gz")
REF = config["fasta_file"]

rule all:
    input:
        mapped = expand("mapped/{sample}.BC.bam", sample=SAMPLES),
        demx = expand("{sample}_demx.done", sample=SAMPLES)
rule fastp:
    input:
        r1 = "fastq_lustre_cat_files/{sample}_read_1.fq.gz",
        r2 = "fastq_lustre_cat_files/{sample}_read_2.fq.gz"
    output:
        out_r1 = "fastq_files_trimmed/{sample}_read_1.fq.gz",
        out_r2 = "fastq_files_trimmed/{sample}_read_2.fq.gz"
    threads: config["threads"]
    shell:
        "{config[fastp]} --in1 {input.r1} --in2 {input.r2} --out1 {output.out_r1} --out2 {output.out_r2} -g -Q -c --dont_eval_duplication -w {threads}"

rule bwa_mem:
    input:
        reads=["fastq_files_trimmed/{sample}_read_1.fq.gz", "fastq_files_trimmed/{sample}_read_2.fq.gz"]
    output:
        temp("mapped/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index=REF,
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
        sort_extra=""            # Extra args for samtools/picard.
    threads: config["threads"]
    wrapper:
        "0.49.0/bio/bwa/mem"

rule move_barcode:
    input: "mapped/{sample}.bam"
    output: "mapped/{sample}.BC.bam"
    shell: "python3 ../scripts/preprocessing/move_barcode.py {input} {output}"

rule demx:
    input: "mapped/{sample}.BC.bam"
    output: touch("{sample}_demx.done")
    shell:"""
    mkdir mapped_qname
    {config[split_barcoded_bam]} {input} mapped_qname/ BC {config[expected_barcodes]}
    """

