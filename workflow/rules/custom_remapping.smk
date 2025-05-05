def get_cells_donor(wildcards, l, donor):
        checkpoint_output = checkpoints.demx.get(**wildcards).output[0]
        print(checkpoint_output)
        res = [f.replace(".bam", "") for f in os.listdir(checkpoint_output) if f.endswith(".bam") and f.split('.')[0] in l]
        print(donor, len(res))
        return res

rule:
    name: "make_fastqs"
    input: "results/mapped_qname_r1/{sample}.bam"
    output: "results/polyA_rich_fastqs/{sample}.fastq.gz"
    shell: "python3 workflow/scripts/refine/extract_polyA_rich_reads.py -i {input} -o {output}"

for donor, donor_config in config['custom_references'].items():
    TAG = "BC"
    rule:
        name: "make_bed_file_{}".format(donor)
        input: donor_config["file_list"]
        output: "results/custom_references/{prefix}.UNK_unmerged.bed".format(prefix = donor)
        shell: "python3 workflow/scripts/refine/make_bed_file.py -i {input} -o {output}"
    rule:
        name: "sort_bed_file_{}".format(donor)
        input: "results/custom_references/{prefix}.UNK_unmerged.bed".format(prefix = donor)
        output: "results/custom_references/{prefix}.UNK_unmerged.sorted.bed".format(prefix = donor)
        shell: "sort -k1,1 -k2,2n {input} > {output}"
    rule:
        name: "merge_bed_file_{}".format(donor)
        input: "results/custom_references/{prefix}.UNK_unmerged.sorted.bed".format(prefix = donor)
        output: "results/custom_references/{prefix}.UNK_merged.bed".format(prefix = donor)
        shell: "{config[bedtools]} merge -i {input} -d 2000 > {output}"
    rule:
        name: "add_germline_intervals_{}".format(donor)
        input: "results/custom_references/{prefix}.UNK_merged.bed".format(prefix = donor)
        output: "results/custom_references/{prefix}.all_merged.bed".format(prefix = donor)
        shell: "cat {input} {config[KNR_remap]} {config[KR_remap]} > {output}"    
    rule:
        name: "make_fasta_{}".format(donor)
        input: "results/custom_references/{prefix}.all_merged.bed".format(prefix = donor)
        output: "results/custom_references/{prefix}.fa".format(prefix = donor)
        shell: "python3 workflow/scripts/refine/make_fasta_from_bed_file.py -i {input} -f {config[fasta_file]} -o {output}"
    rule:
        name: "index_fasta_{}".format(donor)
        input: "results/custom_references/{prefix}.fa".format(prefix = donor)
        output: "results/custom_references/{prefix}.fa.ann".format(prefix = donor)
        shell: "{config[bwa]} index {input}"
rule:
    name: "bwa_mem"
    input:
        reads="results/polyA_rich_fastqs/{sample}.fastq.gz", fa = lambda wildcards: "results/custom_references/{prefix}.fa".format(prefix = config['cbc_to_donor'][wildcards.sample]), faidx = lambda wildcards: "results/custom_references/{prefix}.fa.ann".format(prefix = config['cbc_to_donor'][wildcards.sample])
    output: "results/polyA_rich_mapped_custom/{sample}.bam" 
    log: "logs/bwa_mem_extra/{sample}.no_alt.txt"
    params:
        index= lambda wildcards: "results/custom_references/{prefix}.fa".format(prefix = config['cbc_to_donor'][wildcards.sample]),
        extra=r"-R '@RG\tBC:{sample,[A-Z]+}'",
        sort="none",             # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="TMP_DIR=tmp/"            # Extra args for samtools/picard.
    threads: 1
    shell: "{config[bwa]} mem {params.index} {input.reads} | {config[samtools]} view -bS > {output}"

rule:
    name: "transform_bam"
    input:
        polyA_bam = "results/polyA_rich_mapped_custom/{sample}.bam",
        header_bam = "results/UNK_discond_merged.sorted.bam"
    output: "results/polyA_rich_mapped_custom_tagged_transformed/{sample}.bam"
    shell: "python3 workflow/scripts/refine/convert_coordinates_custom_mapping.py -i {input.polyA_bam} -head {input.header_bam} -o {output}"

rule:
    name: "concat_bams"
    input: lambda wildcards: expand("results/polyA_rich_mapped_custom_tagged_transformed/{sample}.bam", sample=get_cells(wildcards))
    output: "results/polyA_rich_mapped_custom_concat.transformed.bam"
    shell: "samtools cat -o {output} {input}"

rule sort_bam_remap:
    input: "results/polyA_rich_mapped_custom_concat.transformed.bam"
    output: "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam"
    threads: config["threads"]
    shell: "samtools sort -@ {threads} -o {output} {input}"

rule index_bam_remap:
    input: "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam"
    output: "results/polyA_rich_mapped_custom_concat.transformed.sorted.bam.bai"
    threads: config["threads"]
    shell: "samtools index -@ {threads} {input}"




    
