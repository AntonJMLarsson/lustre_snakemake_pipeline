rule get_concord_intervals_KNR:
    input: "results/KNR_bam_concord/{cell}.bam"
    output: temp("results/KNR_concord_bed/{cell}.bed")
    shell: "{config[bedtools]} bamtobed -i {input} | sort -k1,1 -k2,2n | {config[bedtools]} merge -c 5 -o mean -i - > {output}"

rule get_concord_intervals_UNK:
    input:"results/UNK_bam_concord/{cell}.bam"
    output: temp("results/UNK_concord_bed/{cell}.bed")
    shell: "{config[bedtools]} bamtobed -i {input} | sort -k1,1 -k2,2n | {config[bedtools]} merge -c 5 -o mean -i - > {output}"

rule merge_concord_intervals_KNR:
    input: lambda wildcards:expand("results/KNR_concord_bed/{cell}.bed", cell=get_cells(wildcards))
    output:"results/KNR_concord_merged.bed"
    shell: "cat results/KNR_concord_bed/*.bed | sort -k1,1 -k2,2n | {config[bedtools]} merge -c 4 -o mean -i - > {output}"

rule merge_concord_intervals_UNK:
    input: lambda wildcards:expand("results/UNK_concord_bed/{cell}.bed", cell=get_cells(wildcards))
    output: "results/UNK_concord_merged.bed"
    shell: "cat results/UNK_concord_bed/*.bed | sort -k1,1 -k2,2n | {config[bedtools]} merge -c 4 -o mean -i - > {output}"

rule concat_concord_bam_KNR:
    input: lambda wildcards:expand("results/KNR_bam_concord/{cell}.bam", cell=get_cells(wildcards))
    output: "results/KNR_concord_merged.bam"
    shell: "{config[samtools]} cat -o {output} results/KNR_bam_concord/*.bam"

rule sort_concord_bam_KNR:
    input: "results/KNR_concord_merged.bam"
    output: "results/KNR_concord_merged.sorted.bam"
    threads: config["threads"]
    shell:"""
        {config[samtools]} sort -@ {threads} -m 1G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """ 

rule concat_concord_bam_UNK:
    input: lambda wildcards:expand("results/UNK_bam_concord/{cell}.bam", cell=get_cells(wildcards))
    output: "results/UNK_concord_merged.bam"
    shell: "{config[samtools]} cat -o {output} results/UNK_bam_concord/*.bam"
rule sort_concord_bam_UNK:
    input: "results/UNK_concord_merged.bam"
    output: "results/UNK_concord_merged.sorted.bam"
    threads: config["threads"]
    shell:"""
        {config[samtools]} sort -@ {threads} -m 1G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """

rule concat_discond_bam_KNR:
    input: lambda wildcards:expand("results/KNR_bam/{cell}.bam", cell=get_cells(wildcards))
    output: "results/KNR_discond_merged.bam"
    shell: "{config[samtools]} cat -o {output} results/KNR_bam/*.bam"
rule sort_discond_bam_KNR:
    input: "results/KNR_discond_merged.bam"
    output: "results/KNR_discond_merged.sorted.bam"
    threads: config["threads"]
    shell:"""
        {config[samtools]} sort -@ {threads} -m 1G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """

rule concat_bam_KR:
    input: lambda wildcards: expand("results/KR_bam/{cell}.bam", cell=get_cells(wildcards))
    output: temp("results/KR_merged.bam")
    shell: "{config[samtools]} cat -o {output} results/KR_bam/*.bam"

rule sort_bam_KR:
    input: "results/KR_merged.bam"
    output: "results/KR_merged.sorted.bam"
    threads: config["threads"]
    shell:"""
        {config[samtools]} sort -@ {threads} -m 1G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """

rule concat_discond_bam_UNK:
    input: lambda wildcards:expand("results/UNK_bam/{cell}.bam", cell=get_cells(wildcards))
    output: "results/UNK_discond_merged.bam"
    shell: "{config[samtools]} cat -o {output} results/UNK_bam/*.bam"
rule sort_discond_bam_UNK:
    input: "results/UNK_discond_merged.bam"
    output: "results/UNK_discond_merged.sorted.bam"
    threads: config["threads"]
    shell:"""
        {config[samtools]} sort -@ {threads} -m 1G -o {output} {input}
        {config[samtools]} index -@ {threads} {output}
        """

rule filter_concord_intervals_KNR:
    input: bed = "results/KNR_concord_merged.bed", ic = "results/KNR_concord_merged.sorted.bam", idi = "results/KNR_discond_merged.sorted.bam"
    output: "results/KNR_concord_merged.filtered.bed"
    shell: "python3 workflow/scripts/calling/filter_concord_bed_file.py -i {input.bed} -ic {input.ic} -id {input.idi} -o {output}"

rule filter_concord_intervals_UNK:
    input: bed = "results/UNK_concord_merged.bed", ic = "results/UNK_concord_merged.sorted.bam", idi = "results/UNK_discond_merged.sorted.bam"
    output: "results/UNK_concord_merged.filtered.bed"
    shell: "python3 workflow/scripts/calling/filter_concord_bed_file.py -i {input.bed} -ic {input.ic} -id {input.idi} -o {output}"

rule concord_filter_KNR:
     input: bam = "results/KNR_bam/{cell}.bam", bed = "results/KNR_concord_merged.filtered.bed"
     output: "results/KNR_bam_concord_filtered/{cell}.bam"
     shell: "{config[bedtools]} intersect -v -abam {input.bam} -b {input.bed} > {output}"

rule concord_filter_UNK:
     input: bam = "results/UNK_bam/{cell}.bam", bed = "results/UNK_concord_merged.filtered.bed"
     output: "results/UNK_bam_concord_filtered/{cell}.bam"
     shell: "{config[bedtools]} intersect -v -abam {input.bam} -b {input.bed} > {output}"

rule remove_supplemental_KNR:
    input: "results/KNR_bam_concord_filtered/{cell}.bam"
    output: "results/KNR_bam_concord_filtered_no_supp/{cell}.bam"
    shell: "{config[samtools]} view -o {output} -F 0x800 {input}"

rule remove_supplemental_UNK:
    input: "results/UNK_bam_concord_filtered/{cell}.bam"
    output: temp("results/UNK_bam_concord_filtered_no_supp/{cell}.bam")
    shell: "{config[samtools]} view -o {output} -F 0x800 {input}"

rule get_KNR_intervals:
    input: "results/KNR_bam_concord_filtered_no_supp/{cell}.bam"
    output: "results/KNR_intervals_bed/{cell}".replace('.bam', '.bed')
    shell: "{config[bedtools]} bamtobed -i {input} | sort -k1,1 -k2,2n | {config[bedtools]} merge -d 1000  -c 5 -o mean -i - > {output}"

rule get_UNK_intervals:
    input: "results/UNK_bam_concord_filtered_no_supp/{cell}.bam"
    output: "results/UNK_intervals_bed/{cell}".replace('.bam', '.bed')
    shell: "{config[bedtools]} bamtobed -i {input} | sort -k1,1 -k2,2n | {config[bedtools]} merge -d 1000 -c 5 -o mean -i - > {output}"

rule merge_KNR_intervals:
    input: lambda wildcards:expand("results/KNR_intervals_bed/{cell}.bed", cell=get_cells(wildcards))
    output: "results/KNR_intervals_merged.bed"
    shell: "cat results/KNR_intervals_bed/*.bed | sort -k1,1 -k2,2n | {config[bedtools]} merge -d 1000 -c 1,4 -o count,mean -i - > {output}"
        
rule merge_UNK_intervals:
    input: lambda wildcards:expand("results/UNK_intervals_bed/{cell}.bed", cell=get_cells(wildcards))
    output: "results/UNK_intervals_merged.bed"
    shell: "cat results/UNK_intervals_bed/*.bed | sort -k1,1 -k2,2n | {config[bedtools]} merge -d 1000 -c 1,4 -o count,mean -i - > {output}"

rule KR_filter_KNR:
    input: "results/KNR_intervals_merged.bed"
    output: "results/KNR_intervals_merged.KR_filtered.bed"
    shell: "{config[bedtools]} window -w 10000 -a {input} -b {config[KR_file]} -v > {output}"

rule move_KNR_from_UNK:
    input: UNK = "results/UNK_intervals_merged.bed", KNR = "results/KNR_intervals_merged.KR_filtered.bed"
    output: touch("results/moved_KNR.done")
    shell: "{config[bedtools]} intersect -a {input.UNK} -b {config[old_KNR_file]} -wa | {config[bedtools]} intersect -v -a - -b {config[KR_file]} >> {input.KNR}"

rule remove_KNR_from_UNK:
    input: UNK = "results/UNK_intervals_merged.bed", KNR = "results/KNR_intervals_merged.KR_filtered.bed", moved_KNR = "results/moved_KNR.done"
    output: "results/UNK_intervals_merged.KNR_filtered.bed"
    shell: "{config[bedtools]} intersect -a {input.UNK} -b {input.KNR} -v > {output}"

rule filter_UNK:
    input: bed = "results/UNK_intervals_merged.KNR_filtered.bed", bam = "results/UNK_discond_merged.sorted.bam"
    output: "results/UNK_insertions_unmerged.bed"
    shell: "python3 workflow/scripts/calling/tag_site_filtering.py -i {input.bed} -b {input.bam} -o {output}"

rule filter_KNR_1:
    input: bed = "results/KNR_intervals_merged.KR_filtered.bed", bam = "results/KNR_discond_merged.sorted.bam"
    output: "results/KNR_insertions_unmerged.bed"
    shell: "python3 workflow/scripts/calling/tag_site_filtering.py -i {input.bed} -b {input.bam} -o {output}"

rule merge_UNK:
    input: "results/UNK_insertions_unmerged.bed"
    output: "results/UNK_insertions_final.bed"
    shell: "sort -k1,1 -k2,2n {input} | {config[bedtools]} merge -c 4,5,6 -o sum,max,sum -d 2000 -i - > {output}"

rule merge_KNR:
    input: "results/KNR_insertions_unmerged.bed"
    output: "results/KNR_insertions_final.bed"
    shell: "sort -k1,1 -k2,2n {input} | {config[bedtools]} merge -c 4,5,6 -o sum,max,sum -d 2000 -i - > {output}"