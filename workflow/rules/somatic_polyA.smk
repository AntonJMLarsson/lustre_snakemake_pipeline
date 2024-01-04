configfile: "config.yaml"
LOCI, = glob_wildcards("text_files/{loci}.txt")

rule all:
    input: expand("polyA_pkls/{loci}_polyA.pkl", loci=LOCI)

rule make_pkl:
    input: "text_files/{loci}.txt"
    output: "polyA_pkls/{loci}_polyA.pkl"
    shell: "python3 ../scripts/rupture/rupture_and_polyA.py {wildcards.loci} {config[bamfile]} {config[cellfile]}"
