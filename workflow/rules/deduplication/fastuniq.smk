rule fastuniq__deduplicate_reads:
    input:
        unpack(infer_fastqs_for_deduplication),
    output:
        r1=temp_deduplication("results/reads/deduplication/{sample}_R1.fastq.gz"),
        r2=temp_deduplication("results/reads/deduplication/{sample}_R2.fastq.gz"),
        unzipped_out_r1=temp("results/reads/deduplication/{sample}_R1.fastq"),
        unzipped_out_r2=temp("results/reads/deduplication/{sample}_R2.fastq"),
        pair_description=temp("results/reads/deduplication/{sample}.txt"),
    threads: min(config["threads"]["reads__deduplication"], config["max_threads"])
    log:
        "logs/fastuniq/{sample}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/fastuniq/paired"
