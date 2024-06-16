rule fastqc__quality_report:
    input:
        read=infer_fastq_path,
    output:
        html=report(
            "results/reads/{step}/fastqc/{sample}_{pair}.html",
            category="Reads Quality Control",
            labels={
                "Sample": "{sample}",
                "Type": "FastQC {pair} - {step}",
            },
        ),
        zip="results/reads/{step}/fastqc/{sample}_{pair}.zip",
        qc_data="results/reads/{step}/fastqc/{sample}_{pair}/fastqc_data.txt",
        summary_txt="results/reads/{step}/fastqc/{sample}_{pair}/summary.txt",
    threads: min(config["threads"]["reads__fastqc"], config["max_threads"])
    resources:
        mem_mb=get_mem_mb_for_fastqc,
    log:
        "logs/fastqc/{step}/{sample}_{pair}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.14.0/wrappers/fastqc/quality"
