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
        log_once(logfile, "step6.header", "Index BAM files...", f"EXECUTING STEP {master_config['index_rule_num']}")
        log_once(logfile, "step6.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step6.outputfolder", f"Output folder: {params.outputfolder}")
        samtools_version = subprocess.check_output(["samtools", "--version"], stderr=subprocess.STDOUT).decode("utf-8").splitlines()[:2]
        log_once(logfile, "step6.samtools_version", "\n" + "\n".join(samtools_version) + "\n", "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['index_rule_num']-1])
        index_command = f"samtools index -@ {threads} {shlex.quote(input.filtered_BAM)}"
        tracking = begin_step_sample(master_config['index_rule_num'], wildcards.sample, "index_bam")
        record_step_note(master_config['index_rule_num'], wildcards.sample, "indexing_bam")
        record_step_command(master_config['index_rule_num'], wildcards.sample, index_command)
        try:
            shell(index_command)
            finish_step_sample(master_config['index_rule_num'], wildcards.sample, "index_bam", tracking["start_time"], "OK")
        except Exception:
            finish_step_sample(master_config['index_rule_num'], wildcards.sample, "index_bam", tracking["start_time"], "FAIL")
            raise
