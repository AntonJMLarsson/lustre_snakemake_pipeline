
rule all:
    input: "results/{project}.transductions.csv".format(project=config["project"])

rule make_fastqs:
    input:
        "results/mapped_qname_r1/{cell}.bam"
    output:
        temp("results/possible_transduction_fastqs/{cell}.fastq.gz")
    shell:
        "python3 workflow/scripts/refine/extract_possible_transduction_reads.py -i {input} -o {output}"
rule bwa_mem:
    input:
        reads="results/possible_transduction_fastqs/{cell}.fastq.gz"
    output:
        temp("results/possible_transduction_mapped/{cell}.bam")
    log:
        "results/logs/bwa_mem/{cell}.no_alt.txt"
    params:
        index=REF,
        extra=r"-R '@RG\tID:{cell}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="TMP_DIR=tmp/"            # Extra args for samtools/picard.
    threads: 1
    wrapper:
        "0.49.0/bio/bwa/mem"
rule cat:
    input: expand("results/possible_transduction_mapped/{cell}.bam", cell=cellS)
    output: temp("results/Transduction_reads.bam")
    shell: "samtools cat -o {output} results/possible_transduction_mapped/*.bam"
rule sort:
    input: "results/Transduction_reads.bam"
    output:  "results/Transduction_reads.sorted.bam"
    threads: 30
    shell: "samtools sort -t {threads} -m 1G -o {output} {input}"
rule index:
    input: "results/Transduction_reads.sorted.bam"
    output: "results/Transduction_reads.sorted.bam.bai"
    threads: 30
    shell: "samtools index -@ {threads} {input}"

rule find:
    input: bam = "results/Transduction_reads.sorted.bam", bai = "results/Transduction_reads.sorted.bam.bai", UNK = "results/{project}_UNK.csv".format(project=config["project"])
    output: "results/{project}.transductions.csv".format(project=config["project"])
    shell: "python3 workflow/scripts/extra/find_three_prime_transductions.py -i {input.bam} -o {output} --KR {config[KR_file]} --KNR {config[KNR_file]} --UNK {input.UNK}" 