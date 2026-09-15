# Rule 2 fastqc

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import subprocess

rule run_fastqc:
    wildcard_constraints:
        sample=preparation_pattern(2, rule_layout_paired)
    input:
        trimmed_fastq1=lambda wildcards: preparation_trimmed_input(wildcards.sample, 1),
        trimmed_fastq2=lambda wildcards: preparation_trimmed_input(wildcards.sample, 2)
    output:
        report1=f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.trimmed_fastqc.html" if rule_layout_paired else f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}.trimmed_fastqc.html",
        report2=f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.trimmed_fastqc.html" if rule_layout_paired else [],
        report3=f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.trimmed_fastqc.zip" if rule_layout_paired else f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}.trimmed_fastqc.zip",
        report4=f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.trimmed_fastqc.zip" if rule_layout_paired else []
    params:
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['qc_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}"
    threads:
        Threads_Per_Rule['2']
    resources:
        mem_mb = Memory_Per_Rule['2'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['2']
    run:
        def run_fastqc(threads, input, outputfolder):
            log_once(logfile, "step2.header", "Generating FastQC reports...", f"EXECUTING STEP {master_config['qc_rule_num']}")
            log_once(logfile, "step2.inputfolder", f"Input folder: {params.inputfolder}")
            log_once(logfile, "step2.outputfolder", f"Output folder: {params.outputfolder}")
            tracking = begin_step_sample(master_config['qc_rule_num'], wildcards.sample, "run_fastqc")

            fastqc_version = subprocess.check_output(["fastqc", "--version"])
            log_once(logfile, "step2.fastqc_version", "\n"+fastqc_version.decode("utf-8"), "FASTQC VERSION")

            sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['qc_rule_num']-1], "step2.sanity")

            if sample_is_paired(wildcards.sample):
                record_step_note(master_config['qc_rule_num'], wildcards.sample, "running_fastqc_paired_end")
                fastqc_command = f"""
                    fastqc -t {threads} -o {outputfolder} {input.trimmed_fastq1} {input.trimmed_fastq2}
                """
            else:
                record_step_note(master_config['qc_rule_num'], wildcards.sample, "running_fastqc_single_end")
                fastqc_command = f"""
                    fastqc -t {threads} -o {outputfolder} {input.trimmed_fastq1}
                """

            fastqc_command = " ".join(fastqc_command.split())
            record_step_command(master_config['qc_rule_num'], wildcards.sample, fastqc_command)

            # Run the FastQC command
            try:
                shell(fastqc_command)
                if library_runtime:
                    reads = [str(input.trimmed_fastq1)] + ([str(input.trimmed_fastq2)] if sample_is_paired(wildcards.sample) else [])
                    reports = [(str(output.report1), str(output.report3))] + ([(str(output.report2), str(output.report4))] if sample_is_paired(wildcards.sample) else [])
                    for fastq, destinations in zip(reads, reports):
                        stem = re.sub(r"\.(fastq|fq)(\.gz)?$", "", os.path.basename(fastq))
                        for extension, destination in zip(("html", "zip"), destinations):
                            produced = os.path.join(outputfolder, f"{stem}_fastqc.{extension}")
                            if produced != destination:
                                os.replace(produced, destination)
                finish_step_sample(master_config['qc_rule_num'], wildcards.sample, "run_fastqc", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['qc_rule_num'], wildcards.sample, "run_fastqc", tracking["start_time"], "FAIL")
                raise

        # Call the function with parameters
        run_fastqc(threads, input, params.outputfolder)
