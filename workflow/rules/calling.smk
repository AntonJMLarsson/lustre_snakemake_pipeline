configfile: "config.yaml"
SAMPLES, = glob_wildcards("KNR_bam_concord/{sample}.bam")

rule all:
    input:
        UNK = "UNK_insertions_final.bed",
        KNR = "KNR_insertions_final.bed"

  
rule get_concord_intervals_KNR:
    input:
        "KNR_bam_concord/{sample}.bam"
    output:
        temp("KNR_concord_bed/{sample}.bed")
    shell:
        "{config[bedtools]} bamtobed -i {input} | sort -k1,1 -k2,2n | {config[bedtools]} merge -c 5 -o mean -i - > {output}"
rule get_concord_intervals_UNK:
    input:
        "UNK_bam_concord/{sample}.bam"
    output:
        temp("UNK_concord_bed/{sample}.bed")
    shell:
        "{config[bedtools]} bamtobed -i {input} | sort -k1,1 -k2,2n | {config[bedtools]} merge -c 5 -o mean -i - > {output}"

rule merge_concord_intervals_KNR:
    input:
        expand("KNR_concord_bed/{sample}.bed", sample=SAMPLES)
    output:
        "KNR_concord_merged.bed"
    shell:
        "cat KNR_concord_bed/*.bed | sort -k1,1 -k2,2n | {config[bedtools]} merge -c 4 -o mean -i - > {output}"
rule concat_concord_bam_KNR:
    input:
        expand( "KNR_bam_concord/{sample}.bam", sample=SAMPLES)
    output:
        "KNR_concord_merged.bam"
    shell:
        "{config[samtools]} cat -o {output} KNR_bam_concord/*.bam"
rule sort_concord_bam_KNR:
    input:
        "KNR_concord_merged.bam"
    output:
        "KNR_concord_merged.sorted.bam"
    threads:
        30
    shell:
        """
        {config[samtools]} sort -@ {threads} -m 2G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """  
rule concat_concord_bam_UNK:
    input:
        expand( "UNK_bam_concord/{sample}.bam", sample=SAMPLES)
    output:
        "UNK_concord_merged.bam"
    shell:
        "{config[samtools]} cat -o {output} UNK_bam_concord/*.bam"
rule sort_concord_bam_UNK:
    input:
        "UNK_concord_merged.bam"
    output:
        "UNK_concord_merged.sorted.bam"
    threads:
        30
    shell:
        """
        {config[samtools]} sort -@ {threads} -m 2G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """
rule concat_discond_bam_KNR:
    input:
        expand( "KNR_bam/{sample}.bam", sample=SAMPLES)
    output:
        "KNR_discond_merged.bam"
    shell:
        "{config[samtools]} cat -o {output} KNR_bam/*.bam"
rule sort_discond_bam_KNR:
    input:
        "KNR_discond_merged.bam"
    output:
        "KNR_discond_merged.sorted.bam"
    threads:
        30
    shell:
        """
        {config[samtools]} sort -@ {threads} -m 2G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """
rule concat_discond_bam_UNK:
    input:
        expand( "UNK_bam/{sample}.bam", sample=SAMPLES)
    output:
        "UNK_discond_merged.bam"
    shell:
        "{config[samtools]} cat -o {output} UNK_bam/*.bam"
rule sort_discond_bam_UNK:
    input:
        "UNK_discond_merged.bam"
    output:
        "UNK_discond_merged.sorted.bam"
    threads:
        30
    shell:
        """
        {config[samtools]} sort -@ {threads} -m 2G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """
rule merge_concord_intervals_UNK:
    input:
        expand("UNK_concord_bed/{sample}.bed", sample=SAMPLES)
    output:
        "UNK_concord_merged.bed"
    shell:
        "cat UNK_concord_bed/*.bed | sort -k1,1 -k2,2n | {config[bedtools]} merge -c 4 -o mean -i - > {output}"
rule filter_concord_intervals_KNR:
    input:
        bed = "KNR_concord_merged.bed",
        ic = "KNR_concord_merged.sorted.bam",
        idi = "KNR_discond_merged.sorted.bam"
    output:
        "KNR_concord_merged.filtered.bed"
    shell:
        "python3 ../scripts/calling/filter_concord_bed_file.py -i {input.bed} -ic {input.ic} -id {input.idi} -o {output}"
rule filter_concord_intervals_UNK:
    input:
        bed = "UNK_concord_merged.bed",
        ic = "UNK_concord_merged.sorted.bam",
        idi = "UNK_discond_merged.sorted.bam"
    output:
        "UNK_concord_merged.filtered.bed"
    shell:
        "python3 ../scripts/calling/filter_concord_bed_file.py -i {input.bed} -ic {input.ic} -id {input.idi} -o {output}"

rule concord_filter_KNR:
     input:
          bam = "KNR_bam/{sample}.bam",
          bed = "KNR_concord_merged.filtered.bed"
     output:
          "KNR_bam_concord_filtered/{sample}.bam"
     shell:
          "{config[bedtools]} intersect -v -abam {input.bam} -b {input.bed} > {output}"

rule concord_filter_UNK:
     input:
          bam = "UNK_bam/{sample}.bam",
          bed = "UNK_concord_merged.filtered.bed"
     output:
          "UNK_bam_concord_filtered/{sample}.bam"
     shell:
          "{config[bedtools]} intersect -v -abam {input.bam} -b {input.bed} > {output}"

rule remove_supplemental_KNR:
    input:
        "KNR_bam_concord_filtered/{sample}.bam"
    output:
        "KNR_bam_concord_filtered_no_supp/{sample}.bam"
    shell:
        "{config[samtools]} view -o {output} -F 0x800 {input}"

rule remove_supplemental_UNK:
    input:
        "UNK_bam_concord_filtered/{sample}.bam"
    output:
        temp("UNK_bam_concord_filtered_no_supp/{sample}.bam")
    shell:
        "{config[samtools]} view -o {output} -F 0x800 {input}"

rule get_KNR_intervals:
    input:
        "KNR_bam_concord_filtered_no_supp/{sample}.bam"
    output:
        "KNR_intervals_bed/{sample}.bed"
    shell:
        "{config[bedtools]} bamtobed -i {input} | sort -k1,1 -k2,2n | {config[bedtools]} merge -d 1000  -c 5 -o mean -i - > {output}"


rule get_UNK_intervals:
    input:
        "UNK_bam_concord_filtered_no_supp/{sample}.bam"
    output:
        "UNK_intervals_bed/{sample}.bed"
    shell:
        "{config[bedtools]} bamtobed -i {input} | sort -k1,1 -k2,2n | {config[bedtools]} merge -d 1000 -c 5 -o mean -i - > {output}"

rule merge_KNR_intervals:
    input:
        expand("KNR_intervals_bed/{sample}.bed", sample=SAMPLES)
    output:
        "KNR_intervals_merged.bed"
    shell:
        "cat KNR_intervals_bed/*.bed | sort -k1,1 -k2,2n | {config[bedtools]} merge -d 1000 -c 1,4 -o count,mean -i - > {output}"
        
rule merge_UNK_intervals:
    input:
        expand("UNK_intervals_bed/{sample}.bed", sample=SAMPLES)
    output:
        "UNK_intervals_merged.bed"
    shell:
        "cat UNK_intervals_bed/*.bed | sort -k1,1 -k2,2n | {config[bedtools]} merge -d 1000 -c 1,4 -o count,mean -i - > {output}"


rule KR_filter_KNR:
    input:
        "KNR_intervals_merged.bed"
    output:
        "KNR_intervals_merged.KR_filtered.bed"
    shell:
        "{config[bedtools]} window -w 10000 -a KNR_intervals_merged.bed -b {config[KR_file]} -v > KNR_intervals_merged.KR_filtered.bed"

rule move_KNR_from_UNK:
    input:
        UNK = "UNK_intervals_merged.bed",
        KNR = "KNR_intervals_merged.KR_filtered.bed"
    output:
        touch("moved_KNR.done")
    shell:
        "{config[bedtools]} intersect -a {input.UNK} -b {config[old_KNR_file]} -wa | {config[bedtools]} intersect -v -a - -b {config[KR_file]} >> {input.KNR}"

rule remove_KNR_from_UNK:
    input:
        UNK = "UNK_intervals_merged.bed",
        KNR = "KNR_intervals_merged.KR_filtered.bed",
        moved_KNR = "moved_KNR.done"
    output:
        "UNK_intervals_merged.KNR_filtered.bed"
    shell:
        "{config[bedtools]} intersect -a {input.UNK} -b {input.KNR} -v > {output}"
rule filter_UNK:
    input:
        bed = "UNK_intervals_merged.KNR_filtered.bed",
        bam = "UNK_discond_merged.sorted.bam"
    output:
        "UNK_insertions_unmerged.bed"
    shell:
        "python3 ../scripts/calling/tag_site_filtering.py -i {input.bed} -b {input.bam} -o {output}"
rule filter_KNR:
    input:
        bed = "KNR_intervals_merged.KR_filtered.bed",
        bam = "KNR_discond_merged.sorted.bam"
    output:
        "KNR_insertions_unmerged.bed"
    shell:
        "python3 ../scripts/calling/tag_site_filtering.py -i {input.bed} -b {input.bam} -o {output}"
rule merge_UNK:
    input:
        "UNK_insertions_unmerged.bed"
    output:
        "UNK_insertions_final.bed"
    shell:
        "sort -k1,1 -k2,2n {input} | {config[bedtools]} merge -c 4,5,6 -o sum,max,sum -d 2000 -i - > {output}"
rule merge_KNR:
    input:
        "KNR_insertions_unmerged.bed"
    output:
        "KNR_insertions_final.bed"
    shell:
        "sort -k1,1 -k2,2n {input} | {config[bedtools]} merge -c 4,5,6 -o sum,max,sum -d 2000 -i - > {output}"