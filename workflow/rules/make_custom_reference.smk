configfile: "custom_reference.yaml"

rule all:
    input: "{prefix}.fa.bwt".format(prefix = config["prefix"])

rule make_bed_file:
    input: config["file_list"]
    output: "{prefix}.UNK_unmerged.bed".format(prefix = config["prefix"])
    shell: "python3 ../scripts/refine/make_bed_file.py -i {input} -o {output}"
rule sort_bed_file:
    input: "{prefix}.UNK_unmerged.bed".format(prefix = config["prefix"])
    output: "{prefix}.UNK_unmerged.sorted.bed".format(prefix = config["prefix"])
    shell: "sort -k1,1 -k2,2n {input} > {output}"
rule merge_bed_file:
    input: "{prefix}.UNK_unmerged.sorted.bed".format(prefix = config["prefix"])
    output: "{prefix}.UNK_merged.bed".format(prefix = config["prefix"])
    shell: "{config[bedtools]} merge -i {input} -d 2000 > {output}"
rule add_germline_intervals:
    input: "{prefix}.UNK_merged.bed".format(prefix = config["prefix"])
    output: "{prefix}.all_merged.bed".format(prefix = config["prefix"])
    shell: "cat {input} {config[KNR_bed]} {config[KR_bed]} > {output}"    
rule make_fasta:
    input: "{prefix}.all_merged.bed".format(prefix = config["prefix"])
    output: "{prefix}.fa".format(prefix = config["prefix"])
    shell: "python3 ../scripts/refine/make_fasta_from_bed_file.py -i {input} -f {config[fasta_file]} -o {output}"
rule index_fasta:
    input: "{prefix}.fa".format(prefix = config["prefix"])
    output: "{prefix}.fa.bwt".format(prefix = config["prefix"])
    shell: "bwa index {input}"