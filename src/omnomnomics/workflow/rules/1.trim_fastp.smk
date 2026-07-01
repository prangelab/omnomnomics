# Rule 1 option fastp

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import gzip
import json
import shlex
import shutil
import subprocess
import tempfile
import time

rule run_fastp:
    wildcard_constraints:
        sample=lane_sample_wildcard_pattern
    input:
        fastq1=lambda wildcards: resolve_fastq_input(
            wildcards.sample,
            "R1" if config["PAIRED"] else "SE",
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
        inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['trim_rule_num']-1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}",
        adapter_mode=config.get("FASTP_ADAPTER_MODE", "overlap"),
        adapter_sequence=config.get("FASTP_ADAPTER_SEQUENCE", ""),
        adapter_sequence_r2=config.get("FASTP_ADAPTER_SEQUENCE_R2", "")
    threads:
        Threads_Per_Rule['1']
    resources:
        mem_mb=Memory_Per_Rule['1'],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule['1']
    run:
        def run_fastp(logfile, trim_tool, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample, trimmed_fastq1, trimmed_fastq2, trim_metrics_output, adapter_mode, adapter_sequence, adapter_sequence_r2):
            log_once(logfile, "step1.header", "Trimming reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_once(logfile, "step1.inputfolder", f"Input folder: {inputfolder}")
            log_once(logfile, "step1.outputfolder", f"Output folder: {outputfolder}")
            log_once(logfile, "step1.trimtool", f"Trim Tool: {trim_tool}")
            log_once(logfile, "step1.fastp_adapter_mode", f"fastp adapter mode: {adapter_mode}")
            tracking = begin_step_sample(master_config['trim_rule_num'], sample, "run_fastp")

            fastp_version = subprocess.check_output(["fastp", "--version"], stderr=subprocess.STDOUT)
            log_once(logfile, "step1.fastp_version", "\n" + fastp_version.decode("utf-8"), "FASTP VERSION")

            sanity_check_dir(logfile, inputfolder, master_config['input_file_types'][master_config['trim_rule_num'] - 1], "step1.sanity")

            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample}.", dir=tmpdir_root if tmpdir_root else None)
            record_step_note(master_config['trim_rule_num'], sample, f"scratch_dir={local_workdir}")

            def quote(path):
                return shlex.quote(path)

            def format_bytes(path):
                try:
                    return f"{os.path.getsize(path)} bytes"
                except OSError as exc:
                    return f"unavailable ({exc})"

            def first_fastq_read_id(path):
                opener = gzip.open if str(path).endswith(".gz") else open
                try:
                    with opener(path, "rt") as handle:
                        first_line = handle.readline().strip()
                    return first_line.split()[0] if first_line else "empty_file"
                except Exception as exc:
                    return f"unavailable ({exc})"

            def stage_input(path):
                local_path = os.path.join(local_workdir, os.path.basename(path))
                stage_command = f"cp {quote(path)} {quote(local_path)}"
                record_step_note(master_config['trim_rule_num'], sample, f"staging {os.path.basename(path)} size={format_bytes(path)}")
                shell(stage_command)
                record_step_note(master_config['trim_rule_num'], sample, f"staged {os.path.basename(local_path)} size={format_bytes(local_path)} first_read={first_fastq_read_id(local_path)}")
                return local_path

            def adapter_args(paired):
                mode = str(adapter_mode or "overlap").strip().lower()
                if mode == "off":
                    return "--disable_adapter_trimming"
                if paired and mode == "auto_detect":
                    return "--detect_adapter_for_pe"
                if mode in {"nextera", "truseq", "explicit"}:
                    args = []
                    if adapter_sequence:
                        args.append(f"--adapter_sequence {quote(adapter_sequence)}")
                    if paired and adapter_sequence_r2:
                        args.append(f"--adapter_sequence_r2 {quote(adapter_sequence_r2)}")
                    return " ".join(args)
                return ""

            def run_fastp_command(command_args):
                command_string = shlex.join(command_args)
                record_step_command(master_config['trim_rule_num'], sample, command_string)
                record_step_note(master_config['trim_rule_num'], sample, "fastp_start")
                start_time = time.time()
                process = subprocess.Popen(
                    command_args,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                for line in process.stdout:
                    print(f"[fastp:{sample}] {line}", end="", flush=True)
                return_code = process.wait()
                elapsed_seconds = time.time() - start_time
                record_step_note(master_config['trim_rule_num'], sample, f"fastp_exit_code={return_code} elapsed_seconds={elapsed_seconds:.1f}")
                if return_code != 0:
                    raise subprocess.CalledProcessError(return_code, command_args)

            def require_output(path, label):
                if not path or not os.path.exists(path):
                    raise RuntimeError(f"fastp did not create expected {label}: {path}")
                record_step_note(master_config['trim_rule_num'], sample, f"{label}_created size={format_bytes(path)}")

            try:
                local_fastq1 = stage_input(fastq1)
                local_fastq2 = stage_input(fastq2) if fastq2 else ""

                local_trimmed_fastq1 = os.path.join(local_workdir, os.path.basename(trimmed_fastq1))
                local_trimmed_fastq2 = os.path.join(local_workdir, os.path.basename(trimmed_fastq2)) if fastq2 else ""
                html_report = os.path.join(local_workdir, f"{sample}.fastp.html")
                json_report = os.path.join(local_workdir, f"{sample}.fastp.json")

                if fastq2:
                    record_step_note(master_config['trim_rule_num'], sample, "running_fastp_paired_end")
                    adapter_option = adapter_args(True)
                    record_step_note(master_config['trim_rule_num'], sample, f"fastp_adapter_option={adapter_option or 'overlap_default'}")
                    fastp_command = [
                        "fastp",
                        "--in1", local_fastq1,
                        "--in2", local_fastq2,
                        "--out1", local_trimmed_fastq1,
                        "--out2", local_trimmed_fastq2,
                        "--thread", str(threads),
                        "--html", html_report,
                        "--json", json_report,
                    ] + shlex.split(adapter_option)
                else:
                    record_step_note(master_config['trim_rule_num'], sample, "running_fastp_single_end")
                    adapter_option = adapter_args(False)
                    record_step_note(master_config['trim_rule_num'], sample, f"fastp_adapter_option={adapter_option or 'single_end_default'}")
                    fastp_command = [
                        "fastp",
                        "--in1", local_fastq1,
                        "--out1", local_trimmed_fastq1,
                        "--thread", str(threads),
                        "--html", html_report,
                        "--json", json_report,
                    ] + shlex.split(adapter_option)

                run_fastp_command(fastp_command)
                require_output(local_trimmed_fastq1, "trimmed_fastq1")
                if local_trimmed_fastq2:
                    require_output(local_trimmed_fastq2, "trimmed_fastq2")
                require_output(html_report, "fastp_html_report")
                require_output(json_report, "fastp_json_report")

                with open(json_report, "r") as handle:
                    fastp_metrics = json.load(handle)
                before_filtering = fastp_metrics.get("summary", {}).get("before_filtering", {})
                after_filtering = fastp_metrics.get("summary", {}).get("after_filtering", {})
                raw_reads = before_filtering.get("total_reads")
                trimmed_reads = after_filtering.get("total_reads")
                trim_metrics_path = os.path.join(local_workdir, f"{sample}.trim_metrics.tsv")
                with open(trim_metrics_path, "w") as metrics_handle:
                    metrics_handle.write("metric\tvalue\n")
                    metrics_handle.write(f"raw_reads\t{raw_reads if raw_reads is not None else 'NA'}\n")
                    metrics_handle.write(f"trimmed_reads\t{trimmed_reads if trimmed_reads is not None else 'NA'}\n")

                copy_fastq1_command = f"cp {quote(local_trimmed_fastq1)} {quote(trimmed_fastq1)}"
                record_step_note(master_config['trim_rule_num'], sample, "copying_trimmed_r1_back")
                shell(copy_fastq1_command)

                if fastq2 and trimmed_fastq2:
                    copy_fastq2_command = f"cp {quote(local_trimmed_fastq2)} {quote(trimmed_fastq2)}"
                    record_step_note(master_config['trim_rule_num'], sample, "copying_trimmed_r2_back")
                    shell(copy_fastq2_command)
                shell(f"cp {quote(trim_metrics_path)} {quote(trim_metrics_output)}")
                finish_step_sample(master_config['trim_rule_num'], sample, "run_fastp", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['trim_rule_num'], sample, "run_fastp", tracking["start_time"], "FAIL")
                raise
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
            output.trim_metrics,
            params.adapter_mode,
            params.adapter_sequence,
            params.adapter_sequence_r2,
        )
