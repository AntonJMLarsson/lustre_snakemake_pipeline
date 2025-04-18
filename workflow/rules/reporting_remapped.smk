rule count_cell_detection_UNK_remapped:
    input: bam =  "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam", bai = "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam.bai", stats = "results/UNK_split_report.csv"
    output: touch("results/UNK_remapped.done") 
    shell: "python3 workflow/scripts/count_cell_detection.py -bam {input.bam} -stats {input.stats} -sample {config[samplesheet]} -p results/remapped_stats/{config[project]}_UNK -ct BC -t 50"

rule count_cell_detection_KNR_remapped:
    input: bam =  "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam", bai = "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam.bai", stats = "results/KNR_split_report.csv"
    output: touch("results/KNR_remapped.done") 
    shell: "python3 workflow/scripts/count_cell_detection.py -bam {input.bam} -stats {input.stats} -sample {config[samplesheet]} -p results/remapped_stats/{config[project]}_KNR -ct BC -t 50"


