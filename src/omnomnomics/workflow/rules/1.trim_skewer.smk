# Rule 1 option skewer

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob
import shlex
import shutil
import subprocess
import tempfile

rule run_skewer:
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
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['trim_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}"
    threads: 
        Threads_Per_Rule['1']
    resources:
        mem_mb = Memory_Per_Rule['1'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['1']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/benchmarks/{{sample}}_skewer_benchmark.tsv"
    run:
        def run_skewer(logfile, trim_tool, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample, trimmed_fastq1, trimmed_fastq2):
            log_once(logfile, "step1.header", "Trimming reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_once(logfile, "step1.inputfolder", f"Input folder: {inputfolder}")
            log_once(logfile, "step1.outputfolder", f"Output folder: {outputfolder}")
            log_once(logfile, "step1.trimtool", f"Trim Tool: {trim_tool}")
            log_it(logfile, f"Sample {sample}: trimming reads...")

            skewer_version = subprocess.check_output(["skewer", "--version"])
            log_once(logfile, "step1.skewer_version", "\n"+skewer_version.decode("utf-8"), "SKEWER VERSION")
            pigz_version = subprocess.check_output(["pigz", "--version"], stderr=subprocess.STDOUT)
            log_once(logfile, "step1.pigz_version", "\n"+pigz_version.decode("utf-8"), "PIGZ VERSION")

            sanity_check_dir(logfile, inputfolder,  master_config['input_file_types'][master_config['trim_rule_num']-1], "step1.sanity")

            if seq_type == "ATAC":
                adapter_option = "-x CTGTCTCTTATACACATCT -y AGATGTGTATAAGAGACAG" if config["PAIRED"] else "-x CTGTCTCTTATACACATCT"
            else:
                adapter_option = ""

            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample}.", dir=tmpdir_root if tmpdir_root else None)
            log_it(logfile, f"Scratch directory: {local_workdir}")

            def quote(path):
                return shlex.quote(path)

            def local_fastq_path(path):
                local_gz_path = os.path.join(local_workdir, os.path.basename(path))
                local_fastq_path = local_gz_path[:-3] if local_gz_path.endswith(".gz") else local_gz_path
                return local_gz_path, local_fastq_path

            def compress_to_output(source_path, target_path):
                command = f'pigz -p {threads} -c {quote(source_path)} > {quote(target_path)}'
                log_it(logfile, command, "PIGZ COMPRESS COMMAND")
                shell(command)
                log_it(logfile, f"Compressed {source_path} to {target_path}")

            try:
                local_fastq1_gz, local_fastq1 = local_fastq_path(fastq1)
                stage_fastq1_command = f'cp {quote(fastq1)} {quote(local_fastq1_gz)}'
                log_it(logfile, stage_fastq1_command, "STAGE INPUT COMMAND")
                shell(stage_fastq1_command)
                log_it(logfile, f"Staged {fastq1} to {local_fastq1_gz}")

                decompress_fastq1_command = f'pigz -d -p {threads} -c {quote(local_fastq1_gz)} > {quote(local_fastq1)}'
                log_it(logfile, decompress_fastq1_command, "PIGZ DECOMPRESS COMMAND")
                shell(decompress_fastq1_command)
                log_it(logfile, f"Decompressed {local_fastq1_gz} to {local_fastq1}")

                local_fastq2 = ""
                if fastq2:
                    local_fastq2_gz, local_fastq2 = local_fastq_path(fastq2)
                    stage_fastq2_command = f'cp {quote(fastq2)} {quote(local_fastq2_gz)}'
                    log_it(logfile, stage_fastq2_command, "STAGE INPUT COMMAND")
                    shell(stage_fastq2_command)
                    log_it(logfile, f"Staged {fastq2} to {local_fastq2_gz}")

                    decompress_fastq2_command = f'pigz -d -p {threads} -c {quote(local_fastq2_gz)} > {quote(local_fastq2)}'
                    log_it(logfile, decompress_fastq2_command, "PIGZ DECOMPRESS COMMAND")
                    shell(decompress_fastq2_command)
                    log_it(logfile, f"Decompressed {local_fastq2_gz} to {local_fastq2}")

                sample_prefix = os.path.join(local_workdir, sample)

                if fastq2:
                    log_it(logfile, "Running skewer in Paired End mode.")
                    skewer_command = f"""
                        skewer --quiet {adapter_option} -m pe -q 15 -Q 15 -t {threads} -o "{sample_prefix}" {local_fastq1} {local_fastq2}
                    """
                else:
                    log_it(logfile, "Running skewer in Single End mode.")
                    skewer_command = f"""
                        skewer --quiet {adapter_option} -m any -q 15 -Q 15 -t {threads} -o "{sample_prefix}" {local_fastq1}
                    """

                skewer_command = " ".join(skewer_command.split())
                log_it(logfile, skewer_command, "SKEWER COMMAND")

                # Run the skewer command
                shell(skewer_command, bench_record=bench_record)
                log_it(logfile, f"Skewer completed for {sample}")

                log_it(logfile, "Compressing and staging trimmed results...")

                pair1_files = glob.glob(sample_prefix + '*pair1.fastq')
                if fastq2 and len(pair1_files) != 1:
                    raise FileNotFoundError(f"Expected one paired-end R1 output for {sample}, found {len(pair1_files)} in {local_workdir}")
                if pair1_files:
                    compress_to_output(pair1_files[0], trimmed_fastq1)

                pair2_files = glob.glob(sample_prefix + '*pair2.fastq')
                if fastq2 and len(pair2_files) != 1:
                    raise FileNotFoundError(f"Expected one paired-end R2 output for {sample}, found {len(pair2_files)} in {local_workdir}")
                if pair2_files and trimmed_fastq2:
                    compress_to_output(pair2_files[0], trimmed_fastq2)

                for file_path in glob.glob(sample_prefix + '*-trimmed.fastq'):
                    base_name = os.path.basename(file_path)
                    new_name = os.path.join(outputfolder, base_name.replace('-trimmed.fastq', '.trimmed.fastq.gz'))
                    compress_to_output(file_path, new_name)
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

        # Call the function with parameters
        run_skewer(
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
