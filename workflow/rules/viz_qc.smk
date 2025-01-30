configfile: "config.yaml"

rule all:
    input: split_reads_UNK = "viz_split_reads_UNK.done", split_reads_KNR = "viz_split_reads_KNR.done", filtering = "viz_filtering.done"

rule viz_split_reads_analysis_UNK:
    input: "UNK_split_report.csv"
    output: touch("viz_split_reads_UNK.done")
    shell: "python3 workflow/scripts/visualization/viz_split_reads_analysis.py -i {input} -p {config[project]}_UNK"

rule viz_split_reads_analysis_KNR:
    input: "KNR_split_report.csv"
    output: touch("viz_split_reads_KNR.done")
    shell: "python3 workflow/scripts/visualization/viz_split_reads_analysis.py -i {input} -p {config[project]}_KNR"

rule summarize_filtering_stats:
    input: "filter_stats/"
    output: "{project}_stats_filter.csv".format(project=config["project"])
    shell: "python3 workflow/scripts/summarize_filtering_stats.py -f {input} -s {config[samplesheet]} -o {output}"

rule viz_reads_filter:
    input: "{project}_stats_filter.csv".format(project=config["project"])
    output: touch("viz_filtering.done")
    shell: "python3 workflow/scripts/visualization/viz_reads_filter.py -i {input} -g Sample"