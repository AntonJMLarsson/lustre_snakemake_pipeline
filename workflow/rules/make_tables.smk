
DONORS = list(config["custom_references"].keys())
PROJECTS = [config["project"]]

for donor in config["custom_references"].keys():
    rule:
        name: "prepare_insertion_table_KNR_{}".format(donor)
        input: "results/KNR_remapped.done", "results/KNR_regular.done"
        output: "results/insertion_tables/{project}_KNR_{donor}_insertion_table.csv".format(project=config["project"], donor=donor)
        params: regular_stats = "results/regular_stats/{project}_KNR_{donor}_stats.csv".format(project=config["project"], donor=donor), regular_read_mat = "results/regular_stats/{project}_KNR_{donor}_read_depth.csv".format(project=config["project"], donor=donor), remapped_stats = "results/regular_stats/{project}_KNR_{donor}_stats.csv".format(project=config["project"], donor=donor), remapped_read_mat = "results/regular_stats/{project}_KNR_{donor}_read_depth.csv".format(project=config["project"], donor=donor)
        shell: "python3 workflow/scripts/refine/prepare_insertion_table.py --prefix insertion_tables/{config[project]}_KNR_{wildcards.donor} --name {config[project]} --regular {input.regular_stats} --remapped {input.remapped_stats} --read-mat {input.regular_read_mat} --read-mat-remapped {input.remapped_read_mat}"

    rule:
        name: "prepare_insertion_table_UNK_{}".format(donor)
        input: "results/UNK_remapped.done", "results/UNK_regular.done"
        output: "results/insertion_tables/{project}_UNK_{donor}_insertion_table.csv".format(project=config["project"], donor=donor)
        params: regular_stats = "results/regular_stats/{project}_UNK_{donor}_stats.csv".format(project=config["project"], donor=donor), regular_read_mat = "results/regular_stats/{project}_UNK_{donor}_read_depth.csv".format(project=config["project"], donor=donor), remapped_stats = "results/regular_stats/{project}_UNK_{donor}_stats.csv".format(project=config["project"], donor=donor), remapped_read_mat = "results/regular_stats/{project}_UNK_{donor}_read_depth.csv".format(project=config["project"], donor=donor)
        shell: "python3 workflow/scripts/refine/prepare_insertion_table.py --prefix insertion_tables/{config[project]}_UNK_{wildcards.donor} --name {config[project]} --regular {input.regular_stats} --remapped {input.remapped_stats} --read-mat {input.regular_read_mat} --read-mat-remapped {input.remapped_read_mat}"