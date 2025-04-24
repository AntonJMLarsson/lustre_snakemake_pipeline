rule analyze_split_reads_UNK:
    input: bam = "results/UNK_discond_merged.sorted.bam", bed = "results/UNK_insertions_final.bed"
    output: "results/UNK_split_report.csv"
    threads: config["threads"]
    shell: "python3 workflow/scripts/calling/analyze_split_reads.py -bam {input.bam} -bed {input.bed} -o {output} -t {threads}"

rule analyze_split_reads_KNR:
    input: bam = "results/KNR_discond_merged.sorted.bam", bed = "results/KNR_insertions_final.bed"
    output: "results/KNR_split_report.csv"
    threads: config["threads"]
    shell: "python3 workflow/scripts/calling/analyze_split_reads.py -bam {input.bam} -bed {input.bed} -o {output} -t {threads}"

def flatten_file_list(config):
    output = []
    for c in config["custom_references"].values():
        output.extend(c["file_list"])
    return output

rule count_cell_detection_UNK:
    input: bam =  "results/UNK_discond_merged.sorted.bam", stats = "results/UNK_split_report.csv"
    output: flatten_file_list(config), touch("results/UNK_regular.done")
    shell: "python3 workflow/scripts/count_cell_detection.py -bam {input.bam} -stats {input.stats} -sample {config[samplesheet]} -p results/regular_stats/{config[project]}_UNK -ct BC -t 50"

rule count_cell_detection_KNR:
    input: bam =  "results/KNR_discond_merged.sorted.bam", stats = "results/KNR_split_report.csv"
    output: touch("results/KNR_regular.done") 
    shell: "python3 workflow/scripts/count_cell_detection.py -bam {input.bam} -stats {input.stats} -sample {config[samplesheet]} -p results/regular_stats/{config[project]}_KNR -ct BC -t 50"


