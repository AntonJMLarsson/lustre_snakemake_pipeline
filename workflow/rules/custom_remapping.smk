rule make_fastqs:
    input: "results/mapped_qname_r1/{cell}.bam"
    output: "results/polyA_rich_fastqs/{cell}.fastq.gz"
    shell: "python3 workflow/scripts/refine/extract_polyA_rich_reads.py -i {input} -o {output}"

for donor, donor_config in config['custom_references'].items():
    SPECIFIC_SAMPLES = set([line.rstrip() for line in open(donor_config["cell_file"][0])])
    SAMPLES = lambda wildcards: set(SPECIFIC_SAMPLES.intersection(set(get_cells(wildcards))))
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
        output: "results/custom_references/{prefix}.fa.bwt".format(prefix = donor)
        shell: "bwa index {input}"
    rule:
        name: "bwa_mem_{}".format(donor)
        input:
            reads="polyA_rich_fastqs/{sample}.fastq.gz", fa = "results/custom_references/{prefix}.fa".format(prefix = donor), faidx = "results/custom_references/{prefix}.fa.bwt".format(prefix = donor)
        output: "results/polyA_rich_mapped_custom/{sample}.bam"
        log: "logs/bwa_mem_extra/{sample}.no_alt.txt"
        params:
            index="results/custom_references/{prefix}.fa".format(prefix = donor),
            extra=r"-R '@RG\tBC:{sample}'",
            sort="none",             # Can be 'none', 'samtools' or 'picard'.
            sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
            sort_extra="TMP_DIR=tmp/"            # Extra args for samtools/picard.
        threads: 1
        wrapper:
            "0.49.0/bio/bwa/mem"

    rule:
        name: "move_custom_{}".format(donor)
        input: "results/polyA_rich_mapped_custom/{sample}.bam"
        output: "results/polyA_rich_mapped_custom_tagged/{sample}.bam"
        params: tag = TAG
        shell: "python3 workflow/scripts/refine/move_BC_from_header_to_tag.py -i {input} -o {output} --tag {params.tag}"

    rule:
        name: "transform_bam_{}".format(donor)
        input:
            polyA_bam = "results/polyA_rich_mapped_custom_tagged/{sample}.bam",
            header_bam = "results/UNK_discond_merged.sorted.bam"
        output: "polyA_rich_mapped_custom_tagged_transformed/{sample}.bam"
        shell: "python3 workflow/scripts/refine/convert_coordinates_custom_mapping.py -i {input.polyA_bam} -head {input.header_bam} -o {output}"

    rule:
        name: "concat_bams_{}".format(donor)
        input: expand("results/polyA_rich_mapped_custom/{sample}.bam", sample=SAMPLES)
        output: "results/polyA_rich_mapped_custom_concat.transformed.{donor}.bam".format(donor=donor)
        shell: "samtools cat -o {output} {input}"

rule concat_donor_bams:
    input: expand("results/polyA_rich_mapped_custom_concat.transformed.{donor}.bam", donor=config['custom_references'].keys())
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




    