configfile: "config.yaml"
LOCI, = glob_wildcards("text_files/{loci}.txt")

rule all:
    input:
        expand("polyA_pkls/{loci}_polyA.pkl", loci=LOCI),
        "polyA_pkls/polyA_pickles.csv"

rule make_pkl:
    input: "text_files/{loci}.txt"
    output: "polyA_pkls/{loci}_polyA.pkl"
    shell: "python3 {workflow.basedir}/../../workflow/scripts/rupture/rupture_and_polyA.py {wildcards.loci} {config[bamfile]} {config[cellfile]}"


rule polyA_pickles:
    input:
        expand("polyA_pkls/{loci}_polyA.pkl", loci=LOCI)
    output:
        "polyA_pkls/polyA_pickles.csv"
    run:
        insertion_csv = config.get("somatic_insertion_csv") or config.get("insertion_csv")
        samplesheet_csv = (
            config.get("somatic_samplesheet_csv")
            or config.get("samplesheet_csv")
            or config.get("samplesheet")
        )

        if not insertion_csv:
            raise ValueError(
                "Missing insertion CSV path. Set somatic_insertion_csv (or insertion_csv) in the config file."
            )
        if not samplesheet_csv:
            raise ValueError(
                "Missing samplesheet CSV path. Set somatic_samplesheet_csv (or samplesheet_csv/samplesheet) in the config file."
            )

        shell(
            "python3 {workflow.basedir}/../../workflow/scripts/rupture/polyA_pickles.py "
            "{insertion_csv:q} {samplesheet_csv:q} {output:q}"
        )
