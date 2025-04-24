
rule fastp:
    input:
        r1 = config["r1"],
        r2 = config["r2"]
    output:
        out_r1 = "results/fastq_files_trimmed/{project}_read_1.fq.gz".format(project=config["project"]),
        out_r2 = "results/fastq_files_trimmed/{project}_read_2.fq.gz".format(project=config["project"])
    threads: config["threads"]
    shell:
        "{config[fastp]} --in1 {input.r1} --in2 {input.r2} --out1 {output.out_r1} --out2 {output.out_r2} -g -Q -c --dont_eval_duplication -w {threads}"


rule bwa_mem:
    input: r1 = "results/fastq_files_trimmed/{project}_read_1.fq.gz".format(project=config["project"]), r2 = "results/fastq_files_trimmed/{project}_read_2.fq.gz".format(project=config["project"])
    output: temp("results/mapped/{project}.bam".format(project=config["project"]))
    log: "logs/bwa_mem/{project}.log".format(project=config["project"])
    params: index=REF, extra="-R '@RG\tID:{project}\tSM:{project}'".format(project=config["project"]), bwa = int(3*(config["threads"]/4)), samtools = int(config["threads"]/4)
    threads: config["threads"]
    shell: "{config[bwa]} mem -t {params.bwa} {params.extra} {params.index} {input.r1} {input.r2} | samtools view -@ {params.samtools} -bS -  > {output}" 
    

rule move_barcode:
    input: "results/mapped/{project}.bam".format(project=config["project"])
    output: "results/mapped/{project}.BC.bam".format(project=config["project"])
    shell: "python3 workflow/scripts/preprocessing/move_barcode.py {input} {output}"

checkpoint demx:
    input: "results/mapped/{project}.BC.bam".format(project=config["project"])
    output: directory("results/mapped_qname")
    shell:"""
    mkdir results/mapped_qname
    {config[split_barcoded_bam]} {input} results/mapped_qname/ BC {config[expected_barcodes]}
    """

def get_cells(wildcards):
    checkpoint_output = checkpoints.demx.get(**wildcards).output[0]
    print(checkpoint_output)
    return [f.replace(".bam", "") for f in os.listdir(checkpoint_output) if f.endswith(".bam")]
