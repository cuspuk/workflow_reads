from snakemake.utils import validate


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml", set_default=False)


pepfile: config.get("pepfile", "config/pep/config.yaml")


validate(pep.sample_table, "../schemas/samples.schema.yaml")


### Layer for adapting other workflows  ###############################################################################


### Data input handling independent of wildcards ######################################################################


def get_sample_names():
    return list(pep.sample_table["sample_name"].values)


def get_constraints():
    constraints = {
        "sample": "|".join(get_sample_names()),
    }
    return constraints


def get_one_fastq_file(sample: str, read_pair="fq1"):
    return pep.sample_table.loc[sample][[read_pair]]


def get_fastq_paths(sample: str):
    return pep.sample_table.loc[sample][["fq1", "fq2"]]


def get_read_processing_steps() -> list[str]:
    return [step for step in config["reads"] if config["reads"][step] and not step.startswith("_")]


def get_last_step() -> str | None:
    steps = get_read_processing_steps()
    return steps[-1] if steps else None


def get_previous_step_from_step(current_step: str):
    previous = "original"
    for step in get_read_processing_steps():
        if step == current_step:
            return previous
        previous = step


def get_reads_for_step(step: str, sample: str):
    if step != "original":
        return {
            "r1": f"results/reads/{step}/{sample}_R1.fastq.gz",
            "r2": f"results/reads/{step}/{sample}_R2.fastq.gz",
        }
    paths = get_fastq_paths(sample)
    return {
        "r1": paths["fq1"],
        "r2": paths["fq2"],
    }


### Contract for other workflows ######################################################################################


def get_final_fastq_for_sample(sample: str):
    if last_step := get_last_step():
        return expand(f"results/reads/{last_step}/{sample}_{{pair}}.fastq.gz", pair=["R1", "R2"])
    return get_fastq_paths(sample)


def get_multiqc_inputs():
    outs = {}
    if config["reads"]["trimming"] == "cutadapt":
        outs["cutadapt"] = expand("results/reads/trimming/{sample}.qc.txt", sample=get_sample_names())

    if fastqc_steps := config["reads"]["_generate_fastqc_for"]:
        outs["fastqc"] = expand(
            f"results/reads/{{step}}/fastqc/{{sample}}_{{pair}}/fastqc_data.txt",
            sample=get_sample_names(),
            pair=["R1", "R2"],
            step=fastqc_steps,
        )
    if config["reads"]["decontamination"] == "kraken":
        outs["kraken"] = expand("results/kraken/{sample}.kreport2", sample=get_sample_names())
    return outs


### Global rule-set stuff #############################################################################################


def get_outputs():
    outputs = {}
    sample_names = get_sample_names()
    if last_step := get_last_step():
        outputs["reads"] = expand(
            f"results/reads/{last_step}/{{sample}}_{{pair}}.fastq.gz", sample=sample_names, pair=["R1", "R2"]
        )

    if config["reads"]["_generate_fastqc_for"]:
        outputs["fastqc_report"] = expand(
            "results/reads/{step}/fastqc/{sample}_{pair}.html",
            step=config["reads"]["_generate_fastqc_for"],
            sample=sample_names,
            pair=["R1", "R2"],
        )

    if config["reads"]["decontamination"] == "kraken":
        if config["reads__decontamination__kraken"]["generate_krona"]:
            outputs["krona_reports"] = expand(
                "results/kraken/kronas/{sample}.html",
                sample=sample_names,
            )
    return outputs


def get_standalone_outputs():
    # outputs that will be produced if the module is run as a standalone workflow, not as a part of a larger workflow
    return {
        "multiqc_report": "results/_aggregation/multiqc.html",
    }


def temp_decontamination(output_file):
    if get_last_step() == "decontamination":
        return output_file
    return temp(output_file)


def temp_trimming(output_file):
    if get_last_step() == "trimming":
        return output_file
    return temp(output_file)


def temp_deduplication(output_file):
    if get_last_step() == "deduplication":
        return output_file
    return temp(output_file)


def temp_subsampling(output_file):
    if get_last_step() == "subsampling":
        return output_file
    return temp(output_file)


def infer_fastq_path(wildcards):
    if wildcards.step != "original":
        return "results/reads/{step}/{sample}_{pair}.fastq.gz"
    if "pair" not in wildcards or wildcards.pair == "R1":
        return get_one_fastq_file(wildcards.sample, read_pair="fq1").iloc[0]
    elif wildcards.pair == "R2":
        return get_one_fastq_file(wildcards.sample, read_pair="fq2").iloc[0]


def infer_fastqs_for_decontamination(wildcards):
    return get_reads_for_step(get_previous_step_from_step("decontamination"), wildcards.sample)


def infer_fastqs_for_deduplication(wildcards):
    return get_reads_for_step(get_previous_step_from_step("deduplication"), wildcards.sample)


def infer_fastqs_for_trimming(wildcards):
    return get_reads_for_step(get_previous_step_from_step("trimming"), wildcards.sample)


def infer_fastqs_for_subsampling(wildcards):
    return get_reads_for_step(get_previous_step_from_step("subsampling"), wildcards.sample)


### Rule-granularity stuff ############################################################################################


def parse_adapter_removal_params(adapter_config):
    args_lst = [
        f"--action {adapter_config['action']}",
        f"--overlap {adapter_config['overlap']}",
        f"--times {adapter_config['times']}",
        f"--error-rate {adapter_config['error_rate']}",
    ]
    if adapter_config["keep_trimmed_only"]:
        args_lst.append("--discard-untrimmed")

    if path := adapter_config["adapters_anywhere_file"]:
        if not os.path.exists(path):
            raise ValueError(
                f"Adapter removal is enabled, but the adapter file specified in the adapters_anywhere_file element = {path} does not exist."
            )
        args_lst.append(f"--anywhere file:{path} -B file:{path}")

    if end3_file := adapter_config["adapters_3_end_file"]:
        if not os.path.exists(end3_file):
            raise ValueError(
                f"Adapter removal is enabled, but the adapter file specified in the adapters_3_end_file element = {end3_file} does not exist."
            )
        args_lst.append(f"--adapter file:{path} -A file:{path}")

    if end5_file := adapter_config["adapters_5_end_file"]:
        if not os.path.exists(end5_file):
            raise ValueError(
                f"Adapter removal is enabled, but the adapter file specified in the adapters_5_end_file element = {end5_file} does not exist."
            )
        args_lst.append(f"--front file:{path} -G file:{path}")

    return args_lst


def get_cutadapt_extra(cutadapt_config) -> list[str]:
    args_lst = []

    if (value := cutadapt_config["shorten_to_length"]) is not None:
        args_lst.append(f"--length {value}")
    if (value := cutadapt_config["cut_from_start_r1"]) is not None:
        args_lst.append(f"--cut {value}")
    if (value := cutadapt_config["cut_from_start_r2"]) is not None:
        args_lst.append(f"-U {value}")
    if (value := cutadapt_config["cut_from_end_r1"]) is not None:
        args_lst.append(f"--cut -{value}")
    if (value := cutadapt_config["cut_from_end_r2"]) is not None:
        args_lst.append(f"-U -{value}")

    if (value := cutadapt_config["max_n_bases"]) is not None:
        args_lst.append(f"--max-n {value}")
    if (value := cutadapt_config["max_expected_errors"]) is not None:
        args_lst.append(f"--max-expected-errors {value}")
    if value := cutadapt_config["trim_N_bases_on_ends"]:
        args_lst.append(f"--trim-n")
    if cutadapt_config["nextseq_trimming_mode"]:
        value = config["reads__trimming"]["quality_cutoff_from_3_end_r1"]
        args_lst.append(f"--nextseq-trim={value}")

    if cutadapt_config["do_adapter_removal"]:
        args_lst += parse_adapter_removal_params(cutadapt_config["adapter_removal"])

    return args_lst


def parse_paired_cutadapt_param(pe_config, param1, param2, arg_name) -> str:
    if pe_config.get(param1, None) is not None:
        if pe_config.get(param2, None) is not None:
            return f"{arg_name} {pe_config[param1]}:{pe_config[param2]}"
        else:
            return f"{arg_name} {pe_config[param1]}:"
    elif pe_config.get(param2, None) is not None:
        return f"{arg_name} :{pe_config[param2]}"
    return ""


def parse_cutadapt_comma_param(cutadapt_config, param1, param2, arg_name) -> str:
    if cutadapt_config.get(param1) is not None:
        if cutadapt_config.get(param2) is not None:
            return f"{arg_name} {cutadapt_config[param2]},{cutadapt_config[param1]}"
        else:
            return f"{arg_name} {cutadapt_config[param1]}"
    elif cutadapt_config.get(param2) is not None:
        return f"{arg_name} {cutadapt_config[param2]},0"
    return ""


def get_cutadapt_extra_pe() -> str:
    cutadapt_config = config["reads__trimming__cutadapt"]

    args_lst = get_cutadapt_extra(cutadapt_config)

    if parsed_arg := parse_paired_cutadapt_param(cutadapt_config, "max_length_r1", "max_length_r2", "--maximum-length"):
        args_lst.append(parsed_arg)
    if parsed_arg := parse_paired_cutadapt_param(cutadapt_config, "min_length_r1", "min_length_r2", "--minimum-length"):
        args_lst.append(parsed_arg)
    if qual_cut_arg_r1 := parse_cutadapt_comma_param(
        cutadapt_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "--quality-cutoff"
    ):
        args_lst.append(qual_cut_arg_r1)
    if qual_cut_arg_r2 := parse_cutadapt_comma_param(
        cutadapt_config, "quality_cutoff_from_3_end_r1", "quality_cutoff_from_5_end_r2", "-Q"
    ):
        args_lst.append(qual_cut_arg_r2)
    return " ".join(args_lst)


def get_kraken_decontamination_params():
    extra = []
    if config["reads__decontamination__kraken"]["filter_children"]:
        extra.append("--include-children")
    if config["reads__decontamination__kraken"]["filter_ancestors"]:
        extra.append("--include-parents")
    return " ".join(extra)


def infer_krona_tab(wildcards):
    return os.path.join(config["reads__decontamination__kraken"]["krona_dir"], "taxonomy.tab")


def get_all_relevant_extra_params():
    extra = ""
    if config["reads"]["trimming"] == "cutadapt":
        extra += f"Cutadapt: {get_cutadapt_extra_pe()}\n"
    return extra


### Resource handling #################################################################################################


def get_mem_mb_for_trimming(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["reads__trimming_mem_mb"] * attempt)


def get_mem_mb_for_fastqc(wildcards, attempt):
    return min(config["max_mem_mb"], config["resources"]["reads__fastqc_mem_mb"] * attempt)
