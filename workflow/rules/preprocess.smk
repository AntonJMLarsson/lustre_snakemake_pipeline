

rule filter_bam:
    input: bam = "results/mapped_qname/{cell}"
    output:
        bam = temp("results/mapped_qname_filtered/{cell}.bam"),
        report = "results/filter_stats/{cell}.csv"
    params:
        cell = lambda w: w.cell
    shell:
        "python3 workflow/scripts/preprocessing/filter_bam.py -i {input} -o {output.bam} -r {output.report} -s {params.cell}"

rule split_files:
    input: "results/mapped_qname_filtered/{cell}.bam"
    output:
        r1 = "results/mapped_qname_r1/{cell}.bam",
        r2 = "results/mapped_qname_r2/{cell}.bam"
    shell:
        "python3 workflow/scripts/preprocessing/split_into_read1_and_read2.py -i {input} -r1 {output.r1} -r2 {output.r2}"

rule get_KR:
    input: "results/mapped_qname_r1/{cell}.bam"
    output:
        temp("results/KR_bam_unsorted/{cell}.bam")
    shell:
        "{config[bedtools]} window -u -w 10000 -b {config[KR_file]} -abam {input} > {output}"

rule sort_KR:
    input: "results/KR_bam_unsorted/{cell}.bam"
    output: "results/KR_bam/{cell}"
    shell: "{config[samtools]} sort -o {output} {input}"

rule get_KNR:
    input: "results/mapped_qname_r1/{cell}.bam"
    output: temp("results/KNR_bam_unsorted/{cell}.bam")
    shell:
        "{config[bedtools]} window -u -w 10000 -b {config[KNR_file]} -abam {input} > {output}"

rule remove_concordant_KNR:
    input: "results/KNR_bam_unsorted/{cell}.bam"
    output:
        discon = temp("results/KNR_bam_filtered/{cell}.bam"),
        concord = "results/KNR_bam_concord/{cell}.bam"
    shell:
        "python3 workflow/scripts/preprocessing/remove_concordant.py -i {input} -o1 {output.discon} -o2 {output.concord}"

rule sort_KNR:
    input: "results/KNR_bam_filtered/{cell}"
    output: "results/KNR_bam/{cell}"
    shell: "{config[samtools]} sort -o {output} {input}"

rule filter_KR:
    input: "results/mapped_qname_r1/{cell}"
    output: temp("results/KR_filter/{cell}")
    shell: "{config[bedtools]} window -v -w 10000 -b {config[KR_file]} -abam {input} > {output}"

rule filter_KNR:
    input: "results/KR_filter/{cell}"
    output: temp("results/KNR_filter/{cell}")
    shell: "{config[bedtools]} window -v -w 10000 -b {config[KNR_file]} -abam {input} > {output}"
        
rule remove_concordant_UNK:
    input: "results/KNR_filter/{cell}"
    output: discon = temp("results/UNK_bam_unsorted/{cell}"), concord = "results/UNK_bam_concord/{cell}"
    shell: "python3 workflow/scripts/preprocessing/remove_concordant.py -i {input} -o1 {output.discon} -o2 {output.concord}"
        
rule sort_UNK:
    input: "results/UNK_bam_unsorted/{cell}"
    output: "results/UNK_bam/{cell}"
    shell:  "{config[samtools]} sort -o {output} {input}"
