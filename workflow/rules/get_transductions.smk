configfile: "config.yaml"
SAMPLES, = glob_wildcards("mapped_qname_r1/{sample}.bam")
REF = "/home/antonl/meta/reference_genomes/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
rule all:
    input:
        expand("Transduction_reads.sorted.bam.bai")

rule make_fastqs:
    input:
        "mapped_qname_r1/{sample}.bam"
    output:
        temp("possible_transduction_fastqs/{sample}.fastq.gz")
    shell:
        "python3 ../scripts/refine/extract_possible_transduction_reads.py -i {input} -o {output}"
rule bwa_mem:
    input:
        reads="possible_transduction_fastqs/{sample}.fastq.gz"
    output:
        temp("possible_transduction_mapped/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.no_alt.txt"
    params:
        index=REF,
        extra=r"-R '@RG\tID:{sample}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="TMP_DIR=tmp/"            # Extra args for samtools/picard.
    threads: 1
    wrapper:
        "0.49.0/bio/bwa/mem"
rule cat:
    input: expand("possible_transduction_mapped/{sample}.bam", sample=SAMPLES)
    output: temp("Transduction_reads.bam")
    shell: "samtools cat -o {output} possible_transduction_mapped/*.bam"
rule sort:
    input: "Transduction_reads.bam"
    output:  "Transduction_reads.sorted.bam"
    threads: 30
    shell: "samtools sort -t {threads} -m 2G -o {output} {input}"
rule index:
    input: "Transduction_reads.sorted.bam"
    output: "Transduction_reads.sorted.bam.bai"
    threads: 30
    shell: "samtools index -@ {threads} {input}"

rule find:
    input: bam = "Transduction_reads.sorted.bam", bai = "Transduction_reads.sorted.bam.bai", UNK = "{project}_UNK.csv".format(project=config["project"])
    output: "{project}.transductions.csv".format(project=config["project"])
    shell: "python3 ../scripts/extra/find_three_prime_transductions.py -i {input.bam} -o {output} --KR {config[KR_file]} --KNR {config[KNR_file]} --UNK {input.UNK}" 