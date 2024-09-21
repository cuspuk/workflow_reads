rule multiqc__report:
    input:
        **get_multiqc_inputs(),
        config=f"{workflow.basedir}/resources/multiqc.yaml",
    output:
        report(
            "results/_aggregation/multiqc.html",
            category="Summary",
            labels={
                "Type": "MultiQC",
            },
        ),
    params:
        use_input_files_only=True,
    log:
        "logs/multiqc/all.log",
    wrapper:
        "v4.5.0/bio/multiqc"
