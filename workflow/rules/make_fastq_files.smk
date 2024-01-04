SAMPLES, = glob_wildcards("mapped_qname_r1/{sample}.bam")

rule all:
    input:
        expand("polyA_rich_fastqs/{sample}.fastq.gz", sample=SAMPLES)

rule make_fastqs:
    input:
        "mapped_qname_r1/{sample}.bam"
    output:
        "polyA_rich_fastqs/{sample}.fastq.gz"
    shell:
        "python3 ../scripts/refine/extract_polyA_rich_reads.py -i {input} -o {output}"
