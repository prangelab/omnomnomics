# Rule 9 Create Wiggles

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import shlex
import shutil
import subprocess
import tempfile

rule create_wiggles:
    input:
        input_bam=f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.filtered.bam",
        input_bai=f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.bai" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.filtered.bam.bai"
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num']-1]}/{{sample}}.extra_8.tmp"
    params:
        thetype=config['THETYPE'],
        inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num']-1]}"
    threads:
        Threads_Per_Rule['8']
    resources:
        mem_mb=Memory_Per_Rule['8'],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule['8']
    run:
        log_once(logfile, "step8.header", "Creating BigWigs...", f"EXECUTING STEP {master_config['wig_rule_num']}")
        log_once(logfile, "step8.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step8.outputfolder", f"Output folder: {params.outputfolder}")

        bamcoverage_version = subprocess.check_output(["bamCoverage", "--version"], stderr=subprocess.STDOUT)
        log_once(logfile, "step8.bamcoverage_version", "\n" + bamcoverage_version.decode("utf-8"), "DEEPTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['wig_rule_num'] - 1], "step8.sanity")

        def create_wig(input_bam, input_bai, outputfolder, sample):
            log_it(logfile, f"Sample {sample}: creating BigWig...")
            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample}.", dir=tmpdir_root if tmpdir_root else None)
            log_it(logfile, f"Scratch directory: {local_workdir}")

            def quote(path):
                return shlex.quote(path)

            def stage_input(path):
                local_path = os.path.join(local_workdir, os.path.basename(path))
                command = f"cp {quote(path)} {quote(local_path)}"
                log_it(logfile, f"Staging {os.path.basename(path)} to scratch...")
                shell(command)
                return local_path

            try:
                local_bam = stage_input(input_bam)
                local_bai = stage_input(input_bai)
                local_bw = os.path.join(local_workdir, f"{sample}.bw")
                output_bw = os.path.join(outputfolder, f"{sample}.bw")

                bamcoverage_command = f"""
                    bamCoverage -b {quote(local_bam)} -o {quote(local_bw)} \
                    --numberOfProcessors {threads} --binSize 10 --normalizeUsing CPM
                """
                bamcoverage_command = " ".join(bamcoverage_command.split())
                log_it(logfile, bamcoverage_command, "BAMCOVERAGE COMMAND")
                shell(bamcoverage_command)

                copy_bw_command = f"cp {quote(local_bw)} {quote(output_bw)}"
                log_it(logfile, f"Copying BigWig for {sample} back to project space...")
                shell(copy_bw_command)
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

            shell(f"""echo "necessity file for create wiggles. can delete this." > {outputfolder}/{sample}.extra_8.tmp""")

        create_wig(input.input_bam, input.input_bai, params.outputfolder, wildcards.sample)
