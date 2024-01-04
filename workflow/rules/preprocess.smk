configfile: "config.yaml"
SAMPLES, = glob_wildcards("mapped_qname/{sample}.bam")

rule all:
    input:
        UNK = expand("UNK_bam/{sample}.bam", sample=SAMPLES),
        KNR = expand("KNR_bam/{sample}.bam", sample=SAMPLES),
        KR = expand("KR_bam/{sample}.bam", sample=SAMPLES)


rule filter_bam:
    input:
        "mapped_qname/{sample}.bam"
    output:
        bam = temp("mapped_qname_filtered/{sample}.bam"),
        report = "filter_stats/{sample}.csv"
    params:
        sample = lambda w: w.sample
    shell:
        "python3 ../scripts/preprocessing/filter_bam.py -i {input} -o {output.bam} -r {output.report} -s {params.sample}"

rule split_files:
    input:
        "mapped_qname_filtered/{sample}.bam"
    output:
        r1 = "mapped_qname_r1/{sample}.bam",
        r2 = "mapped_qname_r2/{sample}.bam"
    shell:
        "python3 ../scripts/preprocessing/split_into_read1_and_read2.py -i {input} -r1 {output.r1} -r2 {output.r2}"

rule get_KR:
    input:
        "mapped_qname_r1/{sample}.bam"
    output:
        temp("KR_bam_unsorted/{sample}.bam")
    shell:
        "{config[bedtools]} window -u -w 10000 -b {config[KR_file]} -abam {input} > {output}"

rule sort_KR:
    input:
        "KR_bam_unsorted/{sample}.bam"
    output:
        "KR_bam/{sample}.bam"
    shell:
        "{config[samtools]} sort -o {output} {input}"

rule get_KNR:
    input:
        "mapped_qname_r1/{sample}.bam"
    output:
        temp("KNR_bam_unsorted/{sample}.bam")
    shell:
        "{config[bedtools]} window -u -w 10000 -b {config[KNR_file]} -abam {input} > {output}"

rule remove_concordant_KNR:
    input:
        "KNR_bam_unsorted/{sample}.bam"
    output:
        discon = temp("KNR_bam_filtered/{sample}.bam"),
        concord = "KNR_bam_concord/{sample}.bam"
    shell:
        "python3 ../scripts/preprocessing/remove_concordant.py -i {input} -o1 {output.discon} -o2 {output.concord}"

rule sort_KNR:
    input:
        "KNR_bam_filtered/{sample}.bam"
    output:
        "KNR_bam/{sample}.bam"
    shell:
        "{config[samtools]} sort -o {output} {input}"

rule filter_KR:
    input:
        "mapped_qname_r1/{sample}.bam"
    output:
        temp("KR_filter/{sample}.bam")
    shell:
        "{config[bedtools]} window -v -w 10000 -b {config[KR_file]} -abam {input} > {output}"

rule filter_KNR:
    input:
        "KR_filter/{sample}.bam"
    output:
        temp("KNR_filter/{sample}.bam")
    shell:
        "{config[bedtools]} window -v -w 10000 -b {config[KNR_file]} -abam {input} > {output}"
        
rule remove_concordant_UNK:
    input:
        "KNR_filter/{sample}.bam"
    output:
        discon = temp("UNK_bam_unsorted/{sample}.bam"),
        concord = "UNK_bam_concord/{sample}.bam"
    shell:
        "python3 ../scripts/preprocessing/remove_concordant.py -i {input} -o1 {output.discon} -o2 {output.concord}"
        
rule sort_UNK:
    input:
        "UNK_bam_unsorted/{sample}.bam"
    output:
        "UNK_bam/{sample}.bam"
    shell: 
        "{config[samtools]} sort -o {output} {input}"
