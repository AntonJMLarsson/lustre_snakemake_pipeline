SAMPLES, = glob_wildcards("polyA_rich_mapped_custom_tagged/{sample}.bam")

rule all:
    input: "polyA_rich_mapped_custom_concat.transformed.sorted.bam"

rule transform_bam:
    input:
        polyA_bam = "polyA_rich_mapped_custom_tagged/{sample}.bam",
        header_bam = "UNK_discond_merged.sorted.bam"
    output: "polyA_rich_mapped_custom_tagged_transformed/{sample}.bam"
    shell: "python3 workflow/scripts/refine/convert_coordinates_custom_mapping.py -i {input.polyA_bam} -head {input.header_bam} -o {output}"

rule concat_bams:
    input: expand("polyA_rich_mapped_custom_tagged_transformed/{sample}.bam", sample=SAMPLES)
    output: "polyA_rich_mapped_custom_concat.transformed.bam"
    shell: "samtools cat -o {output} polyA_rich_mapped_custom_tagged_transformed/*.bam"

rule sort_bam:
    input: "polyA_rich_mapped_custom_concat.transformed.bam"
    output: "polyA_rich_mapped_custom_concat.transformed.sorted.bam"
    shell:"""
    samtools sort -o {output} {input}
    samtools index {output}
    """
