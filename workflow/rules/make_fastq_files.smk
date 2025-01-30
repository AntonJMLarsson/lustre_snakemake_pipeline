

rule all:
    input:
        expand("polyA_rich_fastqs/{cell}.fastq.gz", cell=CELLS)

rule make_fastqs:
    input: "results/mapped_qname_r1/{cell}.bam"
    output: "results/polyA_rich_fastqs/{cell}.fastq.gz"
    shell: "python3 workflow/scripts/refine/extract_polyA_rich_reads.py -i {input} -o {output}"
