

rule viz_split_reads_analysis_UNK:
    input: "UNK_split_report.csv"
    output: touch("viz_split_reads_UNK.done")
    shell: "mkdir -p results/plots && python3 workflow/scripts/visualization/viz_split_reads_analysis.py -i {input} -p results/plots/{config[project]}_UNK"

rule viz_split_reads_analysis_KNR:
    input: "KNR_split_report.csv"
    output: touch("viz_split_reads_KNR.done")
    shell: "mkdir -p results/plots && python3 workflow/scripts/visualization/viz_split_reads_analysis.py -i {input} -p results/plots/{config[project]}_KNR"

rule summarize_filtering_stats:
    input: "filter_stats/"
    output: "{project}_stats_filter.csv".format(project=config["project"])
    shell: "mkdir -p results/plots && python3 workflow/scripts/visualization/summarize_filtering_stats.py -f {input} -s {config[samplesheet]} -o {output}"

rule viz_reads_filter:
    input: "{project}_stats_filter.csv".format(project=config["project"])
    output: touch("viz_filtering.done")
    shell: "mkdir -p results/plots && python3 workflow/scripts/visualization/viz_reads_filter.py -i {input} -g Sample -p results/plots/{config[project]}_filtering"