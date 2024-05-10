rule multiqc__report:
    input:
        **get_multiqc_inputs(),
        config=f"{workflow.basedir}/resources/multiqc.yaml",
    output:
        report(
            "results/_aggregation/multiqc.html",
            category="summary",
            labels={
                "Type": "MultiQC",
            },
        ),
    params:
        use_input_files_only=True,
    log:
        "logs/multiqc/all.log",
    wrapper:
        "v3.10.2/bio/multiqc"
