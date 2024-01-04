configfile: "config.yaml"
DONORS, = glob_wildcards("custom_mapping_{donor}.smk")

rule all:
    input: expand("custom_{donor}.done", donor=DONORS)

rule custom_mapping:
    input: smk = "custom_mapping_{donor}.smk", txt = "{donor}_cells.txt"
    output: touch("custom_{donor}.done")
    threads: config["threads"]
    shell: "snakemake -s {input.smk} -j {threads}"
