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
    wildcard_constraints:
        sample=lane_sample_wildcard_pattern
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
        trimmed_fastq2=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R2.trimmed.fastq.gz" if config["PAIRED"] else [],
        trim_metrics=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}.trim_metrics.tsv"
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
    run:
        def run_skewer(logfile, trim_tool, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample, trimmed_fastq1, trimmed_fastq2, trim_metrics_output):
            log_once(logfile, "step1.header", "Trimming reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_once(logfile, "step1.inputfolder", f"Input folder: {inputfolder}")
            log_once(logfile, "step1.outputfolder", f"Output folder: {outputfolder}")
            log_once(logfile, "step1.trimtool", f"Trim Tool: {trim_tool}")
            tracking = begin_step_sample(master_config['trim_rule_num'], sample, "run_skewer")

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
            record_step_note(master_config['trim_rule_num'], sample, f"scratch_dir={local_workdir}")

            def quote(path):
                return shlex.quote(path)

            def local_fastq_path(path):
                local_gz_path = os.path.join(local_workdir, os.path.basename(path))
                local_fastq_path = local_gz_path[:-3] if local_gz_path.endswith(".gz") else local_gz_path
                return local_gz_path, local_fastq_path

            def compress_to_output(source_path, target_path, count_path=None):
                if count_path:
                    command = (
                        f"bash -lc \"cat {quote(source_path)} | tee >(awk 'END {{print NR/4}}' > {quote(count_path)}) "
                        f"| pigz -p {threads} -c > {quote(target_path)}\""
                    )
                else:
                    command = f'pigz -p {threads} -c {quote(source_path)} > {quote(target_path)}'
                record_step_note(master_config['trim_rule_num'], sample, f"compressing {os.path.basename(source_path)}")
                shell(command)

            try:
                local_fastq1_gz, local_fastq1 = local_fastq_path(fastq1)
                stage_fastq1_command = f'cp {quote(fastq1)} {quote(local_fastq1_gz)}'
                record_step_note(master_config['trim_rule_num'], sample, f"staging {os.path.basename(fastq1)}")
                shell(stage_fastq1_command)

                raw_count_r1 = os.path.join(local_workdir, f"{sample}.raw_r1.count.txt")
                decompress_fastq1_command = (
                    f"bash -lc \"pigz -d -p {threads} -c {quote(local_fastq1_gz)} "
                    f"| tee >(awk 'END {{print NR/4}}' > {quote(raw_count_r1)}) > {quote(local_fastq1)}\""
                )
                record_step_note(master_config['trim_rule_num'], sample, f"decompressing {os.path.basename(local_fastq1_gz)}")
                shell(decompress_fastq1_command)

                local_fastq2 = ""
                raw_count_r2 = None
                if fastq2:
                    local_fastq2_gz, local_fastq2 = local_fastq_path(fastq2)
                    stage_fastq2_command = f'cp {quote(fastq2)} {quote(local_fastq2_gz)}'
                    record_step_note(master_config['trim_rule_num'], sample, f"staging {os.path.basename(fastq2)}")
                    shell(stage_fastq2_command)

                    raw_count_r2 = os.path.join(local_workdir, f"{sample}.raw_r2.count.txt")
                    decompress_fastq2_command = (
                        f"bash -lc \"pigz -d -p {threads} -c {quote(local_fastq2_gz)} "
                        f"| tee >(awk 'END {{print NR/4}}' > {quote(raw_count_r2)}) > {quote(local_fastq2)}\""
                    )
                    record_step_note(master_config['trim_rule_num'], sample, f"decompressing {os.path.basename(local_fastq2_gz)}")
                    shell(decompress_fastq2_command)

                sample_prefix = os.path.join(local_workdir, sample)

                if fastq2:
                    record_step_note(master_config['trim_rule_num'], sample, "running_skewer_paired_end")
                    skewer_command = f"""
                        skewer --quiet {adapter_option} -m pe -q 15 -Q 15 -t {threads} -o "{sample_prefix}" {local_fastq1} {local_fastq2}
                    """
                else:
                    record_step_note(master_config['trim_rule_num'], sample, "running_skewer_single_end")
                    skewer_command = f"""
                        skewer --quiet {adapter_option} -m any -q 15 -Q 15 -t {threads} -o "{sample_prefix}" {local_fastq1}
                    """

                skewer_command = " ".join(skewer_command.split())
                record_step_command(master_config['trim_rule_num'], sample, skewer_command)

                # Run the skewer command
                shell(skewer_command)
                record_step_note(master_config['trim_rule_num'], sample, "compressing_and_staging_trimmed_results")

                pair1_files = glob.glob(sample_prefix + '*pair1.fastq')
                if fastq2 and len(pair1_files) != 1:
                    raise FileNotFoundError(f"Expected one paired-end R1 output for {sample}, found {len(pair1_files)} in {local_workdir}")
                trimmed_count_r1 = os.path.join(local_workdir, f"{sample}.trimmed_r1.count.txt")
                if pair1_files:
                    compress_to_output(pair1_files[0], trimmed_fastq1, trimmed_count_r1)

                pair2_files = glob.glob(sample_prefix + '*pair2.fastq')
                if fastq2 and len(pair2_files) != 1:
                    raise FileNotFoundError(f"Expected one paired-end R2 output for {sample}, found {len(pair2_files)} in {local_workdir}")
                trimmed_count_r2 = None
                if pair2_files and trimmed_fastq2:
                    trimmed_count_r2 = os.path.join(local_workdir, f"{sample}.trimmed_r2.count.txt")
                    compress_to_output(pair2_files[0], trimmed_fastq2, trimmed_count_r2)

                for file_path in glob.glob(sample_prefix + '*-trimmed.fastq'):
                    base_name = os.path.basename(file_path)
                    new_name = os.path.join(outputfolder, base_name.replace('-trimmed.fastq', '.trimmed.fastq.gz'))
                    compress_to_output(file_path, new_name)

                def read_count(count_path):
                    if count_path and os.path.exists(count_path):
                        with open(count_path, "r") as handle:
                            return int(float(handle.read().strip()))
                    return 0

                raw_reads = read_count(raw_count_r1) + read_count(raw_count_r2)
                trimmed_reads = read_count(trimmed_count_r1) + read_count(trimmed_count_r2)
                trim_metrics_path = os.path.join(local_workdir, f"{sample}.trim_metrics.tsv")
                with open(trim_metrics_path, "w") as metrics_handle:
                    metrics_handle.write("metric\tvalue\n")
                    metrics_handle.write(f"raw_reads\t{raw_reads}\n")
                    metrics_handle.write(f"trimmed_reads\t{trimmed_reads}\n")
                shell(f"cp {quote(trim_metrics_path)} {quote(trim_metrics_output)}")
                finish_step_sample(master_config['trim_rule_num'], sample, "run_skewer", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['trim_rule_num'], sample, "run_skewer", tracking["start_time"], "FAIL")
                raise
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
            output.trim_metrics,
        )
