# Rule 2 fastqc

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os

rule run_fastqc:
    input:
        trimmed_fastq1= f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2= f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R2.trimmed.fastq.gz" if config['PAIRED'] else []
    output:
        report1=f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.trimmed_fastqc.html" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}.trimmed_fastqc.html",
        report2=f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.trimmed_fastqc.html" if config["PAIRED"] else [],
        report3=f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.trimmed_fastqc.zip" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}.trimmed_fastqc.zip",
        report4=f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.trimmed_fastqc.zip" if config["PAIRED"] else []
    params:
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['qc_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}"
    threads:
        Threads_Per_Rule['2']
    resources:
        mem_mb = Memory_Per_Rule['2']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['qc_rule_num']-1]}/benchmarks/{{sample}}_fastqc_benchmark.tsv"
    run:
        def run_fastqc(threads, input, outputfolder):
            log_it(logfile, "Generating FastQC reports...", f"EXECUTING STEP {master_config['qc']}")
            log_it(logfile, f"Input folder: {params.inputfolder}")
            log_it(logfile, f"Output folder: {params.outputfolder}")

            fastqc_version = subprocess.check_output("module load fastqc && fastqc --version", shell=True, executable='/bin/bash')
            log_it(logfile, "\n"+fastqc_version.decode("utf-8"), "FASTQC VERSION")

            sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['qc_rule_num']-1])

            if config['PAIRED']:
                log_it(logfile, "Running FastQC in paired end mode...")
                fastqc_command = f"""
                    module load fastqc && \
                    fastqc -t {threads} -o {outputfolder} {input.trimmed_fastq1} {input.trimmed_fastq2}
                """
            else:
                log_it(logfile, "Running FastQC in single end mode...")
                fastqc_command = f"""
                    module load fastqc && \
                    fastqc -t {threads} -o {outputfolder} {input.trimmed_fastq1}
                """

            # Run the FastQC command
            shell(fastqc_command)

        # Call the function with parameters
        run_fastqc(threads, input, params.outputfolder)