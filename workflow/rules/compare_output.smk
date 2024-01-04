SAMPLES, = glob_wildcards("insertion_tables/{sample}_insertion_table.csv")

rule all:
    input: "samples_multiinter.bed"

rule make_bed_files:
    input: "insertion_tables/{sample}_insertion_table.csv"
    output: "UNK_bed_files/{sample}.bed"
    shell: "python3 ../scripts/refine/make_bed_file.py --input {input} --output {output}"

rule sort_bed_files:
    input: "UNK_bed_files/{sample}.bed"
    output: "UNK_bed_files/{sample}.sorted.bed"
    shell: "sort -k 1,1 -k2,2n {input} > {output}"

rule multiinter:
    input: expand("UNK_bed_files/{sample}.sorted.bed", sample=SAMPLES)
    output: "samples_multiinter.bed"
    shell: "bedtools multiinter -header -i {input} > {output}"