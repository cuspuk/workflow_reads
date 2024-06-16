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
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.14.0/wrappers/bwa/index"


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
        fastq_param="-t -n -i",
    threads: min(config["threads"]["reads__decontamination"], config["max_threads"])
    log:
        "logs/decontamination/bwa__filter_reads/{sample}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.14.0/wrappers/bwa/filter"
