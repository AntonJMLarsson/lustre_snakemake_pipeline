configfile: "config.yaml"

rule all:
    input:
        UNK = "UNK_insertions_final.bed",
        KNR = "KNR_insertions_final.bed"


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
        "python3 {workflow.basedir}/tag_site_filtering.py -i {input.bed} -b {input.bam} -o {output}"
rule filter_KNR:
    input:
        bed = "KNR_intervals_merged.KR_filtered.bed",
        bam = "KNR_discond_merged.sorted.bam"
    output:
        "KNR_insertions_unmerged.bed"
    shell:
        "python3 {workflow.basedir}/tag_site_filtering.py -i {input.bed} -b {input.bam} -o {output}"
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