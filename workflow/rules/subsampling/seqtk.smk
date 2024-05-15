rule seqtk__subsample_reads_pe:
    input:
        unpack(infer_fastqs_for_subsampling),
    output:
        r1=temp_subsampling("results/reads/subsampling/{sample}_R1.fastq.gz"),
        r2=temp_subsampling("results/reads/subsampling/{sample}_R2.fastq.gz"),
    params:
        seed=config["reads__subsampling__seqtk"]["seed"],
        n_reads=config["reads__subsampling__seqtk"]["n_reads"],
        reduce_memory=config["reads__subsampling__seqtk"]["reduce_memory"],
    log:
        "logs/seqtk/subsample_reads_pe/{sample}.log",
    wrapper:
        "https://github.com/cuspuk/workflow_wrappers/raw/v1.13.4/wrappers/seqtk/subsample_paired"
