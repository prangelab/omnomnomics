# Rule 6 Index BAM

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import shlex
import subprocess
rule index_bam:
    input:
        filtered_BAM= f"{experiment_dir}/{master_config['input_folders'][master_config['index_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.filtered.bam",
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.bai" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.filtered.bam.bai"
    threads:
        Threads_Per_Rule['6']
    params:
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['index_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['index_rule_num']-1]}"     
    resources:
        mem_mb = Memory_Per_Rule['6'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['6']
    run:
        log_it(logfile, f"Index BAM files...",f"EXECUTING STEP {master_config['index_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")
        log_it(logfile, f"Indexing {wildcards.sample}")
        samtools_version = subprocess.check_output(["samtools", "--version"], stderr=subprocess.STDOUT).decode("utf-8").splitlines()[:2]
        log_it(logfile, "\n" + "\n".join(samtools_version) + "\n", "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['index_rule_num']-1])
        index_command = f"samtools index -@ {threads} {shlex.quote(input.filtered_BAM)}"
        log_it(logfile, index_command, "SAMTOOLS INDEX COMMAND")
        shell(index_command)
