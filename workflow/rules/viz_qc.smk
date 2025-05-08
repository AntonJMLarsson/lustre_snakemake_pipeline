rule viz_split_reads_analysis_UNK:
    input: "results/UNK_split_report.csv"
    output: touch("results/viz_split_reads_UNK.done")
    shell: "mkdir -p results/plots && python3 workflow/scripts/visualization/viz_split_reads_analysis.py -i {input} -p results/plots/{config[project]}_UNK"

rule viz_split_reads_analysis_KNR:
    input: "results/KNR_split_report.csv"
    output: touch("results/viz_split_reads_KNR.done")
    shell: "mkdir -p results/plots && python3 workflow/scripts/visualization/viz_split_reads_analysis.py -i {input} -p results/plots/{config[project]}_KNR"

rule summarize_filtering_stats:
    input: lambda wildcards: expand("results/filter_stats/{cell}.csv", cell=get_cells(wildcards))
    output: "results/{project}_stats_filter.csv".format(project=config["project"])
    params: folder = "results/filter_stats/"
    shell: "mkdir -p results/plots && python3 workflow/scripts/visualization/summarize_filtering_stats.py -f {params.folder} -s {config[samplesheet]} -o {output}"

rule viz_reads_filter:
    input: "results/{project}_stats_filter.csv".format(project=config["project"])
    output: touch("results/viz_filtering.done")
    shell: "mkdir -p results/plots && python3 workflow/scripts/visualization/viz_reads_filter.py -i {input} -g Sample -p results/plots/{config[project]}_filtering"