
ALL_SAMPLES, = glob_wildcards("polyA_rich_fastqs/{sample}.fastq.gz")
SPECIFIC_SAMPLES = set([line.rstrip() for line in open(config["cell_file"])])

SAMPLES = list(SPECIFIC_SAMPLES.intersection(ALL_SAMPLES))

rule all:
    input: expand("results/polyA_rich_mapped_custom/{sample}.bam", sample=SAMPLES)
rule bwa_mem:
    input:
        reads="polyA_rich_fastqs/{sample}.fastq.gz"
    output: "results/polyA_rich_mapped_custom/{sample}.bam"
    log:
        "logs/bwa_mem_extra/{sample}.no_alt.txt"
    params:
        index=config["custom_ref"],
        extra=r"-R '@RG\tBC:{sample}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="TMP_DIR=tmp/"            # Extra args for samtools/picard.
    threads: 1
    wrapper:
        "0.49.0/bio/bwa/mem"
