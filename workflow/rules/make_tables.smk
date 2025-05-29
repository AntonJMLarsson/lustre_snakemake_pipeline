
DONORS = list(config["custom_references"].keys())
PROJECTS = [config["project"]]

for donor in config["custom_references"].keys():
    rule:
        name: "prepare_insertion_table_KNR_{}".format(donor)
        input: "results/KNR_remapped.done", "results/KNR_regular.done"
        output: "results/insertion_tables/{project}_KNR_{donor}_insertion_table.csv".format(project=config["project"], donor=donor)
        params: donor = donor, regular_stats = "results/regular_stats/{project}_KNR_{donor}_stats.csv".format(project=config["project"], donor=donor), regular_read_mat = "results/regular_stats/{project}_KNR_{donor}_read_depth.csv".format(project=config["project"], donor=donor), remapped_stats = "results/regular_stats/{project}_KNR_{donor}_stats.csv".format(project=config["project"], donor=donor), remapped_read_mat = "results/regular_stats/{project}_KNR_{donor}_read_depth.csv".format(project=config["project"], donor=donor)
        shell: "mkdir -p results/insertion_tables && python3 workflow/scripts/refine/prepare_insertion_table.py --prefix results/insertion_tables/{config[project]}_KNR_{params.donor} --name {config[project]} --regular {params.regular_stats} --remapped {params.remapped_stats} --read-mat {params.regular_read_mat} --read-mat-remapped {params.remapped_read_mat}"

    rule:
        name: "prepare_insertion_table_UNK_{}".format(donor)
        input: "results/UNK_remapped.done", "results/UNK_regular.done"
        output: "results/insertion_tables/{project}_UNK_{donor}_insertion_table.csv".format(project=config["project"], donor=donor)
        params: donor = donor, regular_stats = "results/regular_stats/{project}_UNK_{donor}_stats.csv".format(project=config["project"], donor=donor), regular_read_mat = "results/regular_stats/{project}_UNK_{donor}_read_depth.csv".format(project=config["project"], donor=donor), remapped_stats = "results/regular_stats/{project}_UNK_{donor}_stats.csv".format(project=config["project"], donor=donor), remapped_read_mat = "results/regular_stats/{project}_UNK_{donor}_read_depth.csv".format(project=config["project"], donor=donor)
        shell: "mkdir -p results/insertion_tables && python3 workflow/scripts/refine/prepare_insertion_table.py --prefix results/insertion_tables/{config[project]}_UNK_{params.donor} --name {config[project]} --regular {params.regular_stats} --remapped {params.remapped_stats} --read-mat {params.regular_read_mat} --read-mat-remapped {params.remapped_read_mat}"
    
    rule:
        name: "make_detailed_insertion_stats_UNK_{}".format(donor)
        input: insertion_table = "results/insertion_tables/{project}_UNK_{donor}_insertion_table.csv".format(project=config["project"], donor=donor), bam = "results/UNK_discond_merged.sorted.bam", bam_polyA =  "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam", bai_polyA = "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam.bai"
        output: "results/detailed_stats/detailed_insertion_stats_UNK_{}.csv".format(donor)
        params: donor = donor
        threads: config["threads"]
        shell: "python3 workflow/scripts/extra/detailed_insertion_stats.py -bam {input.bam} -bam-polyA {input.bam_polyA} --stats {input.insertion_table} --samplesheet {config[samplesheet]} --donor {params.donor} -p {config[project]} --tag BC --threads {threads} --output {output}"