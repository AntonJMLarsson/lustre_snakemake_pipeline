

rule filter_bam:
    input: bam = "results/mapped_qname/{cell,[A-Z]+}.bam"
    output:
        bam = temp("results/mapped_qname_filtered/{cell,[A-Z]+}.bam"),
        report = "results/filter_stats/{cell,[A-Z]+}.csv"
    params:
        cell = lambda w: w.cell
    shell:
        "python3 workflow/scripts/preprocessing/filter_bam.py -i {input} -o {output.bam} -r {output.report} -s {params.cell}"

rule split_files:
    input: "results/mapped_qname_filtered/{cell,[A-Z]+}.bam"
    output:
        r1 = "results/mapped_qname_r1/{cell,[A-Z]+}.bam",
        r2 = "results/mapped_qname_r2/{cell,[A-Z]+}.bam"
    shell:
        "python3 workflow/scripts/preprocessing/split_into_read1_and_read2.py -i {input} -r1 {output.r1} -r2 {output.r2}"

rule get_KR:
    input:
        "results/mapped_qname_r1/{cell,[A-Z]+}.bam"
    output:
        temp("results/KR_bam_unsorted/{cell,[A-Z]+}.bam")
    shell:
        "{config[bedtools]} window -u -w 10000 -b {config[KR_file]} -abam {input} > {output}"

rule sort_KR:
    input:
        "results/KR_bam_unsorted/{cell,[A-Z]+}.bam"
    output:
        "results/KR_bam/{cell,[A-Z]+}.bam"
    shell:
        "{config[samtools]} sort -o {output} {input}"

rule get_KNR:
    input:
        "results/mapped_qname_r1/{cell,[A-Z]+}.bam"
    output:
        temp("results/KNR_bam_unsorted/{cell,[A-Z]+}.bam")
    shell:
        "{config[bedtools]} window -u -w 10000 -b {config[KNR_file]} -abam {input} > {output}"

rule remove_concordant_KNR:
    input: "results/KNR_bam_unsorted/{cell,[A-Z]+}.bam"
    output:
        discon = temp("results/KNR_bam_filtered/{cell,[A-Z]+}.bam"),
        concord = "results/KNR_bam_concord/{cell,[A-Z]+}.bam"
    shell:
        "python3 workflow/scripts/preprocessing/remove_concordant.py -i {input} -o1 {output.discon} -o2 {output.concord}"

rule sort_KNR:
    input: "results/KNR_bam_filtered/{cell,[A-Z]+}.bam"
    output: "results/KNR_bam/{cell,[A-Z]+}.bam"
    shell: "{config[samtools]} sort -o {output} {input}"

rule filter_KR:
    input: "results/mapped_qname_r1/{cell,[A-Z]+}.bam"
    output: temp("results/KR_filter/{cell,[A-Z]+}.bam")
    shell: "{config[bedtools]} window -v -w 10000 -b {config[KR_file]} -abam {input} > {output}"

rule filter_KNR:
    input: "results/KR_filter/{cell,[A-Z]+}.bam"
    output: temp("results/KNR_filter/{cell,[A-Z]+}.bam")
    shell: "{config[bedtools]} window -v -w 10000 -b {config[KNR_file]} -abam {input} > {output}"
        
rule remove_concordant_UNK:
    input: "results/KNR_filter/{cell,[A-Z]+}.bam"
    output: discon = temp("results/UNK_bam_unsorted/{cell,[A-Z]+}.bam"), concord = "results/UNK_bam_concord/{cell,[A-Z]+}.bam"
    shell: "python3 workflow/scripts/preprocessing/remove_concordant.py -i {input} -o1 {output.discon} -o2 {output.concord}"
        
rule sort_UNK:
    input: "results/UNK_bam_unsorted/{cell,[A-Z]+}.bam"
    output: "results/UNK_bam/{cell,[A-Z]+}.bam"
    shell:  "{config[samtools]} sort -o {output} {input}"
