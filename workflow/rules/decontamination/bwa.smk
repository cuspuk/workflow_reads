rule bwa__build_index:
    input:
        "{reference_dir}/{fasta}.fa",
    output:
        idx=protected(multiext("{reference_dir}/bwa_index/{fasta}", ".amb", ".ann", ".bwt", ".pac", ".sa")),
    params:
        prefix=lambda wildcards, output: os.path.splitext(output.idx[0])[0],
        approach="bwtsw",
    log:
        "{reference_dir}/bwa_index/logs/{fasta}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/bwa/index"


rule bwa__filter_reads:
    input:
        unpack(infer_fastqs_for_decontamination),
        index=get_bwa_indexes(),
    output:
        r1=temp_decontamination("results/reads/decontamination/{sample}_R1.fastq.gz"),
        r2=temp_decontamination("results/reads/decontamination/{sample}_R2.fastq.gz"),
    params:
        indices=lambda w, input: [os.path.splitext(i)[0] for i in input.index],
        keep_param="-F 2" if config["reads__decontamination__bwa"]["filter_mode"] == "exclude" else "-f 2",
    threads: min(config["threads"]["reads__decontamination"], config["max_threads"])
    log:
        "logs/decontamination/bwa__filter_reads/{sample}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/bwa/filter"


# rule bwa__map_reads:
#     input:
#         reads=get_fastq_for_mapping,
#         index=infer_bwa_index_for_mapping,
#         read_group="results/reads/.read_groups/{sample}.txt",
#     output:
#         bam=temp("results/reads/decontamination/{sample}.original.bam"),
#     threads: min(config["threads"]["reads__decontamination"], config["max_threads"])
#     resources:
#         mem_mb=get_mem_mb_for_mapping, # TODO
#     log:
#         "logs/reads/decontamination/{sample}.log",
#     wrapper:
#         "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/bwa/map"
# rule samtools__get_decontamination_bam:
#     input:
#         "{sample}.sam",
#     output:
#         bam="{sample}.bam",
#         idx="{sample}.bai",
#     log:
#         "{sample}.log",
#     params:
#         extra="-f 4",
#         region="",
#     threads: min(config["threads"]["reads__decontamination"], config["max_threads"])
#     wrapper:
#         "v3.12.1/bio/samtools/view"
# rule samtools__get_decontamination_fastq:
#     input:
#         "mapped/{sample}.bam",
#     output:
#         r1=temp_decontamination("results/reads/decontamination/{sample}_R1.fastq.gz"),
#         r2=temp_decontamination("results/reads/decontamination/{sample}_R2.fastq.gz"),
#     log:
#         "{sample}.separate.log",
#     params:
#         collate="",
#         fastq="-n",
#     threads: min(config["threads"]["reads__decontamination"], config["max_threads"])
#     wrapper:
#         "v3.12.1/bio/samtools/fastq/separate"
# rule bamtobed:
#     input:
#         "{sample}.bam",
#     output:
#         r1=temp_decontamination("results/reads/decontamination/{sample}_R1.fastq.gz"),
#         r2=temp_decontamination("results/reads/decontamination/{sample}_R2.fastq.gz"),
#     log:
#         "logs/bamtobed/{sample}.log",
#     shell:
#         "bedtools bamtofastq -i {input} -fq {output.r1} -fq2 {output.r2} > {log} 2>&1"
