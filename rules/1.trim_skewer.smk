# Rule 1 option skewer

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob
import subprocess

rule run_skewer:
    input:
        fastq1=f"{experiment_dir}/{master_config['input_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R1.fastq.gz",
        fastq2=f"{experiment_dir}/{master_config['input_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R2.fastq.gz" if config["PAIRED"] else []
    output:
        trimmed_fastq1=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R2.trimmed.fastq.gz" if config["PAIRED"] else []
    params:
        seq_type=config["THETYPE"],
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['trim_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}"
    threads: 
        Threads_Per_Rule['1']
    resources:
        mem_mb = Memory_Per_Rule['1']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/benchmarks/{{sample}}_skewer_benchmark.tsv"
    run:
        def run_skewer(logfile, trim_tool, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample):
            log_it(logfile, "Trimming reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_it(logfile, f"Input folder: {inputfolder}")
            log_it(logfile, f"Output folder: {outputfolder}")
            log_it(logfile, f"Trim Tool: {trim_tool}")

            skewer_version = subprocess.check_output(["skewer", "--version"])
            log_it(logfile, "\n"+skewer_version.decode("utf-8"), "SKEWER VERSION")

            sanity_check_dir(logfile, inputfolder,  master_config['input_file_types'][master_config['trim_rule_num']-1])

            if seq_type == "ATAC":
                adapter_option = "-x CTGTCTCTTATACACATCT -y AGATGTGTATAAGAGACAG" if config["PAIRED"] else "-x CTGTCTCTTATACACATCT"
            else:
                adapter_option = ""

            if fastq2:
                log_it(logfile, "Running skewer in Paired End mode.")
                skewer_command = f"""
                    skewer --quiet {adapter_option} -m pe -q 15 -Q 15 -z -t {threads} -o "{outputfolder}/{sample}" {fastq1} {fastq2}
                """
            else:
                log_it(logfile, "Running skewer in Single End mode.")
                skewer_command = f"""
                    skewer --quiet {adapter_option} -m any -q 15 -Q 15 -z -t {threads} -o "{outputfolder}/{sample}" {fastq1}
                """

            # Run the skewer command
            shell(skewer_command, bench_record=bench_record)

            # Rename R1 trimmed files
            log_it(logfile, "Renaming trimmed results...")
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*pair1.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('-trimmed-pair1.fastq.gz', '_R1.trimmed.fastq.gz')) 
                os.rename(file_path, new_name)

            # Rename R2 trimmed files
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*pair2.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('-trimmed-pair2.fastq.gz', '_R2.trimmed.fastq.gz'))
                os.rename(file_path, new_name)

            # Rename unpaired trimmed files
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*-trimmed.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('-trimmed.fastq.gz', '.trimmed.fastq.gz'))
                os.rename(file_path, new_name)

            
            os.remove((os.path.join(f"{outputfolder}", f"{sample}" + "-trimmed.log" )))

        # Call the function with parameters
        run_skewer(logfile, config["THETRIMTOOL"], params.seq_type, threads, input.fastq1, input.fastq2, params.inputfolder, params.outputfolder, wildcards.sample)