rule cutadapt__trim_reads_pe:
    input:
        unpack(infer_fastqs_for_trimming),
    output:
        fastq1=temp_trimming("results/reads/trimming/{sample}_R1.fastq.gz"),
        fastq2=temp_trimming("results/reads/trimming/{sample}_R2.fastq.gz"),
        qc="results/reads/trimming/{sample}.qc.txt",
    params:
        extra=get_cutadapt_extra_pe(),
    resources:
        mem_mb=get_mem_mb_for_trimming,
    threads: min(config["threads"]["reads__trimming"], config["max_threads"])
    log:
        "logs/cutadapt/trim_reads_pe/{sample}.log",
    wrapper:
        "v3.13.8/bio/cutadapt/pe"
