# Rule 1 option fastp

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

rule run_fastp:
    input:
        fastq1=lambda wildcards: resolve_fastq_input(
            wildcards.sample,
            "R1",
            master_config['input_folders'][master_config['trim_rule_num'] - 1],
        ),
        fastq2=(
            lambda wildcards: resolve_fastq_input(
                wildcards.sample,
                "R2",
                master_config['input_folders'][master_config['trim_rule_num'] - 1],
            )
        ) if config["PAIRED"] else []
    output:
        trimmed_fastq1=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R2.trimmed.fastq.gz" if config["PAIRED"] else []
    params:
        seq_type=config["THETYPE"],
        inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['trim_rule_num']-1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}"
    threads:
        Threads_Per_Rule['1']
    resources:
        mem_mb=Memory_Per_Rule['1'],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule['1']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/benchmarks/{{sample}}_fastp_benchmark.tsv"
    run:
        def run_fastp(logfile, trim_tool, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample, trimmed_fastq1, trimmed_fastq2):
            log_once(logfile, "step1.header", "Trimming reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_once(logfile, "step1.inputfolder", f"Input folder: {inputfolder}")
            log_once(logfile, "step1.outputfolder", f"Output folder: {outputfolder}")
            log_once(logfile, "step1.trimtool", f"Trim Tool: {trim_tool}")
            log_it(logfile, f"Sample {sample}: trimming reads...")

            fastp_version = subprocess.check_output(["fastp", "--version"], stderr=subprocess.STDOUT)
            log_once(logfile, "step1.fastp_version", "\n" + fastp_version.decode("utf-8"), "FASTP VERSION")

            sanity_check_dir(logfile, inputfolder, master_config['input_file_types'][master_config['trim_rule_num'] - 1], "step1.sanity")

            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample}.", dir=tmpdir_root if tmpdir_root else None)
            log_it(logfile, f"Scratch directory: {local_workdir}")

            def quote(path):
                return shlex.quote(path)

            def stage_input(path):
                local_path = os.path.join(local_workdir, os.path.basename(path))
                stage_command = f"cp {quote(path)} {quote(local_path)}"
                log_it(logfile, f"Staging {os.path.basename(path)} to scratch...")
                shell(stage_command)
                log_it(logfile, f"Staged {path} to {local_path}")
                return local_path

            try:
                local_fastq1 = stage_input(fastq1)
                local_fastq2 = stage_input(fastq2) if fastq2 else ""

                local_trimmed_fastq1 = os.path.join(local_workdir, os.path.basename(trimmed_fastq1))
                local_trimmed_fastq2 = os.path.join(local_workdir, os.path.basename(trimmed_fastq2)) if fastq2 else ""
                html_report = os.path.join(local_workdir, f"{sample}.fastp.html")
                json_report = os.path.join(local_workdir, f"{sample}.fastp.json")

                if fastq2:
                    log_it(logfile, "Running fastp in Paired End mode.")
                    fastp_command = f"""
                        fastp --in1 {quote(local_fastq1)} --in2 {quote(local_fastq2)} \
                        --out1 {quote(local_trimmed_fastq1)} --out2 {quote(local_trimmed_fastq2)} \
                        --thread {threads} --html {quote(html_report)} --json {quote(json_report)} \
                        --detect_adapter_for_pe
                    """
                else:
                    log_it(logfile, "Running fastp in Single End mode.")
                    fastp_command = f"""
                        fastp --in1 {quote(local_fastq1)} --out1 {quote(local_trimmed_fastq1)} \
                        --thread {threads} --html {quote(html_report)} --json {quote(json_report)}
                    """

                fastp_command = " ".join(fastp_command.split())
                log_it(logfile, fastp_command, "FASTP COMMAND")
                shell(fastp_command, bench_record=bench_record)
                log_it(logfile, f"fastp completed for {sample}")

                copy_fastq1_command = f"cp {quote(local_trimmed_fastq1)} {quote(trimmed_fastq1)}"
                log_it(logfile, f"Copying trimmed R1 for {sample} back to project space...")
                shell(copy_fastq1_command)
                log_it(logfile, f"Copied {local_trimmed_fastq1} to {trimmed_fastq1}")

                if fastq2 and trimmed_fastq2:
                    copy_fastq2_command = f"cp {quote(local_trimmed_fastq2)} {quote(trimmed_fastq2)}"
                    log_it(logfile, f"Copying trimmed R2 for {sample} back to project space...")
                    shell(copy_fastq2_command)
                    log_it(logfile, f"Copied {local_trimmed_fastq2} to {trimmed_fastq2}")
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

        run_fastp(
            logfile,
            config["THETRIMTOOL"],
            params.seq_type,
            threads,
            input.fastq1,
            input.fastq2,
            params.inputfolder,
            params.outputfolder,
            wildcards.sample,
            output.trimmed_fastq1,
            output.trimmed_fastq2,
        )
