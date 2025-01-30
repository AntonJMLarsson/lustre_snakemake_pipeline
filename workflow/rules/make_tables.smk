configfile: "config.yaml"
DONORS, = glob_wildcards("custom_mapping_{donor}.smk")
PROJECTS, = glob_wildcards("{project}_demx.done")

rule all:
    input: KNR = expand("insertion_tables/{project}_KNR_{donor}_insertion_table.csv", project=PROJECTS ,donor=DONORS), UNK = expand("insertion_tables/{project}_UNK_{donor}_insertion_table.csv", project=PROJECTS, donor=DONORS)

rule prepare_insertion_table_KNR:
    input: regular_stats = "regular_stats/{project}_KNR_{donor}_stats.csv", regular_read_mat = "regular_stats/{project}_KNR_{donor}_read_depth.csv",remapped_stats = "remapped_stats/{project}_KNR_{donor}_stats.csv", remapped_read_mat = "remapped_stats/{project}_KNR_{donor}_read_depth.csv" 
    output: "insertion_tables/{project}_KNR_{donor}_insertion_table.csv"
    shell: "python3 workflow/scripts/refine/prepare_insertion_table.py --prefix insertion_tables/{config[project]}_KNR_{wildcards.donor} --name {config[project]} --regular {input.regular_stats} --remapped {input.remapped_stats} --read-mat {input.regular_read_mat} --read-mat-remapped {input.remapped_read_mat}"

rule prepare_insertion_table_UNK:
    input: regular_stats = "regular_stats/{project}_UNK_{donor}_stats.csv" , regular_read_mat = "regular_stats/{project}_UNK_{donor}_read_depth.csv" ,remapped_stats = "remapped_stats/{project}_UNK_{donor}_stats.csv" , remapped_read_mat = "remapped_stats/{project}_UNK_{donor}_read_depth.csv"
    output: "insertion_tables/{project}_UNK_{donor}_insertion_table.csv" 
    shell: "python3 workflow/scripts/refine/prepare_insertion_table.py --prefix insertion_tables/{config[project]}_UNK_{wildcards.donor} --name {config[project]} --regular {input.regular_stats} --remapped {input.remapped_stats} --read-mat {input.regular_read_mat} --read-mat-remapped {input.remapped_read_mat}"