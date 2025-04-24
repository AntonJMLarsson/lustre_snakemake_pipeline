rule make_fastqs:
    input: "results/mapped_qname_r1/{cell}.bam"
    output: "results/polyA_rich_fastqs/{cell}.fastq.gz"
    shell: "python3 workflow/scripts/refine/extract_polyA_rich_reads.py -i {input} -o {output}"

checkpoint fastq:
    input: lambda wildcards: expand("results/polyA_rich_fastqs/{fq}.fastq.gz", fq=get_cells(wildcards))
    output: directory("results/polyA_rich_fastqs/")

def check_fastqs(wildcards):
    fq_output = checkpoints.fastq.get(**wildcards).output[0]
    FQ, = glob_wildcards(os.path.join(fq_output, "{fq}.fastq.gz"))
    return expand("{fq}", fq=FQ)

for donor, donor_config in config['custom_references'].items():
    SPECIFIC_SAMPLES = set([line.rstrip() for line in open(donor_config["cell_file"][0])])
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
        name: "bwa_mem_{}".format(donor)
        input:
            reads="results/polyA_rich_fastqs/{cell}.fastq.gz", fa = "results/custom_references/{prefix}.fa".format(prefix = donor), faidx = "results/custom_references/{prefix}.fa.ann".format(prefix = donor)
        output: "results/polyA_rich_mapped_custom/{cell}_" + "{prefix}.bam".format(prefix = donor)
        log: "logs/bwa_mem_extra/{cell}.no_alt.txt"
        params:
            index="results/custom_references/{prefix}.fa".format(prefix = donor),
            extra=r"-R '@RG\tBC:{sample}'",
            sort="none",             # Can be 'none', 'samtools' or 'picard'.
            sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
            sort_extra="TMP_DIR=tmp/"            # Extra args for samtools/picard.
        threads: 1
        shell: "{config[bwa]} mem {params.index} {input.reads} | {config[samtools]} view -bS > {output}"

    rule:
        name: "transform_bam_{}".format(donor)
        input:
            polyA_bam = "results/polyA_rich_mapped_custom/{cell}_" + "{prefix}.bam".format(prefix = donor),
            header_bam = "results/UNK_discond_merged.sorted.bam"
        output: "results/polyA_rich_mapped_custom_tagged_transformed/{cell}.bam"
        shell: "python3 workflow/scripts/refine/convert_coordinates_custom_mapping.py -i {input.polyA_bam} -head {input.header_bam} -o {output}"

    rule:
        name: "concat_bams_{}".format(donor)
        input: lambda wildcards: expand("results/polyA_rich_mapped_custom_tagged_transformed/{cell}_" + "{prefix}.bam".format(prefix = donor), cell=set(SPECIFIC_SAMPLES.intersection(set(check_fastqs(wildcards)))))
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




    
