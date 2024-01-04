configfile: "config.yaml"

rule all:
    input: KNR = "KNR_regular.done", UNK = "UNK_regular.done"

rule analyze_split_reads_UNK:
    input: bam = "UNK_discond_merged.sorted.bam", bed = "UNK_insertions_final.bed"
    output: "UNK_split_report.csv"
    threads: config["threads"]
    shell: "python3 ../scripts/calling/analyze_split_reads.py -bam {input.bam} -bed {input.bed} -o {output} -t {threads}"

rule analyze_split_reads_KNR:
    input: bam = "KNR_discond_merged.sorted.bam", bed = "KNR_insertions_final.bed"
    output: "KNR_split_report.csv"
    threads: config["threads"]
    shell: "python3 ../scripts/calling/analyze_split_reads.py -bam {input.bam} -bed {input.bed} -o {output} -t {threads}"

rule count_cell_detection_UNK:
    input: bam =  "UNK_discond_merged.sorted.bam", stats = "UNK_split_report.csv"
    output: touch("UNK_regular.done") 
    shell: "python3 ../scripts/count_cell_detection.py -bam {input.bam} -stats {input.stats} -sample {config[samplesheet]} -p regular_stats/{config[project]}_UNK -ct BC -t 50"

rule count_cell_detection_KNR:
    input: bam =  "KNR_discond_merged.sorted.bam", stats = "KNR_split_report.csv"
    output: touch("KNR_regular.done") 
    shell: "python3 ../scripts/count_cell_detection.py -bam {input.bam} -stats {input.stats} -sample {config[samplesheet]} -p regular_stats/{config[project]}_KNR -ct BC -t 50"


