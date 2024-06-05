rule curl__download_kraken_db:
    output:
        protected("{prefix_dir}/hash.k2d"),
    params:
        url=lambda wildcards, output: "https://genome-idx.s3.amazonaws.com/kraken/{tag}.tar.gz".format(
            tag=os.path.basename(os.path.dirname(output[0]))
        ),
        dirpath=lambda wildcards, output: os.path.dirname(output[0]),
    retries: 1
    log:
        "{prefix_dir}/logs/download.log",
    conda:
        "../../envs/curl.yaml"
    shell:
        "(mkdir -p {params.dirpath} && curl -SL {params.url} | tar zxvf - -C {params.dirpath}) > {log} 2>&1"


rule kraken__analysis:
    input:
        unpack(infer_fastqs_for_decontamination),
        kraken_tax=os.path.join(config["reads__decontamination__kraken"]["kraken_dir"], "hash.k2d"),
    output:
        kraken_output=temp("results/kraken/{sample}.kraken"),
        report="results/kraken/{sample}.kreport2",
    params:
        save_memory="--memory-mapping" if config["reads__decontamination__kraken"]["save_memory"] else "",
        db_dir=lambda wildcards, input: os.path.dirname(input.kraken_tax),
    threads: min(config["threads"]["reads__decontamination"], config["max_threads"])
    log:
        "logs/kraken/analysis/{sample}.log",
    conda:
        "../../envs/kraken2.yaml"
    shell:
        "(kraken2 --db {params.db_dir} --threads {threads} --paired --gzip-compressed"
        " {params.save_memory} --report {output.report} {input.r1} {input.r2} 1> {output.kraken_output}) 2> {log}"


rule kraken__decontaminate:
    input:
        unpack(infer_fastqs_for_decontamination),
        kraken_output="results/kraken/{sample}.kraken",
        kraken_report="results/kraken/{sample}.kreport2",
    output:
        r1=temp_decontamination("results/reads/decontamination/{sample}_R1.fastq.gz"),
        r2=temp_decontamination("results/reads/decontamination/{sample}_R2.fastq.gz"),
        temp_r1=temp("results/reads/decontamination/{sample}_R1.fastq"),
        temp_r2=temp("results/reads/decontamination/{sample}_R2.fastq"),
        std_out=temp("results/reads/decontamination/{sample}_decontamination.out"),
    params:
        taxid=" ".join(str(taxa_id) for taxa_id in config["reads__decontamination__kraken"]["taxa_ids"]),
        extra=get_kraken_decontamination_params(),
        exclude=config["reads__decontamination__kraken"]["filter_mode", True],  # true or false
    log:
        "logs/kraken/decontaminate/{sample}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.6/wrappers/kraken/decontaminate_pe"


rule krona__update_taxonomy:
    output:
        protected("{prefix_dir}/taxonomy.tab"),
    params:
        tax_dir=lambda wildcards, output: os.path.dirname(output[0]),
    retries: 1
    log:
        "{prefix_dir}/logs/update_taxonomy.log",
    conda:
        "../../envs/kraken2.yaml"
    localrule: True
    shell:
        "ktUpdateTaxonomy.sh {params.tax_dir} 1> {log} 2>&1"


rule kraken__krona_chart:
    input:
        kraken_output="results/kraken/{sample}.kreport2",
        tax_tab=infer_krona_tab,
    output:
        report(
            "results/kraken/kronas/{sample}.html",
            category="Reads Quality Control",
            labels={
                "Sample": "{sample}",
                "Type": "Krona chart",
            },
        ),
    params:
        extra="-m 3 -t 5",
        tax_dir=lambda wildcards, input: os.path.dirname(input.tax_tab),
    log:
        "logs/kraken/krona_chart/{sample}.log",
    conda:
        "../../envs/kraken2.yaml"
    shell:
        "ktImportTaxonomy {params.extra} -tax {params.tax_dir} -o {output} {input.kraken_output} 1> {log} 2>&1"
