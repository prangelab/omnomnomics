# Rule 1 option trimmomatic

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

rule run_trimmomatic:
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
        trim_tool=config["THETRIMTOOL"],
        seq_type=config["THETYPE"],
        inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['trim_rule_num']-1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}"
    threads:
        Threads_Per_Rule['1']
    resources:
        mem_mb=Memory_Per_Rule['1'],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule['1']
    run:
        def run_trimmomatic(logfile, trim_tool, trim_mem_mb, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample, trimmed_fastq1, trimmed_fastq2):
            trim_mem = f"{max(2048, int(trim_mem_mb * 0.8))}M"
            trim_heap = f"{max(1024, int(trim_mem_mb * 0.4))}M"
            log_once(logfile, "step1.header", "Trimming reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_once(logfile, "step1.inputfolder", f"Input folder: {inputfolder}")
            log_once(logfile, "step1.outputfolder", f"Output folder: {outputfolder}")
            log_once(logfile, "step1.trimtool", f"Trim Tool: {trim_tool}")
            log_once(logfile, "step1.trimheap", f"Trim Heap: {trim_heap}")
            log_once(logfile, "step1.trimmem", f"Trim Mem: {trim_mem}")
            log_it(logfile, f"Sample {sample}: trimming reads...")

            path = os.path.join(config['WORKFLOW_ROOT'], "bin", "Trimmomatic-0.39", "trimmomatic-0.39.jar")
            trimmomatic_version = subprocess.check_output(["java", "-jar", path, "-version"], stderr=subprocess.STDOUT)
            log_once(logfile, "step1.trimmomatic_version", "\n" + trimmomatic_version.decode("utf-8"), "TRIMMOMATIC VERSION")

            sanity_check_dir(logfile, inputfolder, master_config['input_file_types'][master_config['trim_rule_num'] - 1], "step1.sanity")

            if seq_type == "ATAC":
                adapter_file = "NexteraPE-PE.fa" if config["PAIRED"] else "Nextera-SE.fa"
            else:
                adapter_file = "TruSeq3-PE.fa" if config["PAIRED"] else "TruSeq3-SE.fa"

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
                return local_path

            try:
                local_fastq1 = stage_input(fastq1)
                local_fastq2 = stage_input(fastq2) if fastq2 else ""
                local_out_base = os.path.join(local_workdir, f"{sample}.trimmed.fastq.gz")
                adapter_path = os.path.join(config['WORKFLOW_ROOT'], "bin", "Trimmomatic-0.39", "adapters", adapter_file)

                if fastq2:
                    log_it(logfile, f"Running {trim_tool} in Paired End mode.")
                    trimmomatic_command = f"""
                        java -Xms{trim_heap} -Xmx{trim_mem} -jar {quote(path)} PE \
                        -threads {threads} -baseout {quote(local_out_base)} {quote(local_fastq1)} {quote(local_fastq2)} \
                        ILLUMINACLIP:{quote(adapter_path)}:2:30:10:2:True \
                        LEADING:3 TRAILING:3 MINLEN:36
                    """
                else:
                    log_it(logfile, f"Running {trim_tool} in Single End mode.")
                    trimmomatic_command = f"""
                        java -Xms{trim_heap} -Xmx{trim_mem} -jar {quote(path)} SE \
                        -threads {threads} {quote(local_fastq1)} {quote(local_out_base)} \
                        ILLUMINACLIP:{quote(adapter_path)}:2:30:10 \
                        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
                    """

                trimmomatic_command = " ".join(trimmomatic_command.split())
                log_it(logfile, trimmomatic_command, "TRIMMOMATIC COMMAND")
                shell(trimmomatic_command)
                log_it(logfile, f"Trimmomatic completed for {sample}")

                if fastq2:
                    local_pair1 = local_out_base.replace(".trimmed.fastq.gz", ".trimmed_1P.fastq.gz")
                    local_pair2 = local_out_base.replace(".trimmed.fastq.gz", ".trimmed_2P.fastq.gz")
                    if not os.path.exists(local_pair1):
                        raise FileNotFoundError(f"Expected paired-end R1 output for {sample} at {local_pair1}")
                    if not os.path.exists(local_pair2):
                        raise FileNotFoundError(f"Expected paired-end R2 output for {sample} at {local_pair2}")

                    copy_fastq1_command = f"cp {quote(local_pair1)} {quote(trimmed_fastq1)}"
                    log_it(logfile, f"Copying trimmed R1 for {sample} back to project space...")
                    shell(copy_fastq1_command)

                    copy_fastq2_command = f"cp {quote(local_pair2)} {quote(trimmed_fastq2)}"
                    log_it(logfile, f"Copying trimmed R2 for {sample} back to project space...")
                    shell(copy_fastq2_command)
                else:
                    if not os.path.exists(local_out_base):
                        raise FileNotFoundError(f"Expected single-end output for {sample} at {local_out_base}")
                    copy_fastq_command = f"cp {quote(local_out_base)} {quote(trimmed_fastq1)}"
                    log_it(logfile, f"Copying trimmed FASTQ for {sample} back to project space...")
                    shell(copy_fastq_command)
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

        run_trimmomatic(
            logfile,
            params.trim_tool,
            resources.mem_mb,
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
