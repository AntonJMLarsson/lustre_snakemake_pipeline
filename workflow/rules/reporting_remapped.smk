configfile: "config.yaml"

rule all:
    input: KNR = "KNR_remapped.done", UNK = "UNK_remapped.done"

rule count_cell_detection_UNK:
    input: bam =  "polyA_rich_mapped_custom_concat.transformed.sorted.bam", stats = "UNK_split_report.csv"
    output: touch("results/UNK_remapped.done") 
    shell: "python3 workflow/scripts/count_cell_detection.py -bam {input.bam} -stats {input.stats} -sample {config[samplesheet]} -p remapped_stats/{config[project]}_UNK -ct BC -t 50"

rule count_cell_detection_KNR:
    input: bam =  "polyA_rich_mapped_custom_concat.transformed.sorted.bam", stats = "KNR_split_report.csv"
    output: touch("results/KNR_remapped.done") 
    shell: "python3 workflow/scripts/count_cell_detection.py -bam {input.bam} -stats {input.stats} -sample {config[samplesheet]} -p remapped_stats/{config[project]}_KNR -ct BC -t 50"


