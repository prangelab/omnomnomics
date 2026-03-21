#Rule 3 option STAR TE

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

rule run_star_te:
    input:
        trimmed_fastq1=f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample3}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample3}}.trimmed.fastq.gz",
        trimmed_fastq2=f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample3}}_R2.trimmed.fastq.gz" if config['PAIRED'] else []
    output:
        bam=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.bam",
        stats=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.STAR_TE_stats.txt",
        extra=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.extra_3.tmp"
    params:
        genome_path=config['STAR_GENOME_DIR'],
        inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}",
        paired=config['PAIRED']
    threads:
        Threads_Per_Rule['3']
    resources:
        mem_mb=Memory_Per_Rule['3'],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule['3']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/benchmarks/{{sample3}}_star_te_benchmark.tsv"
    run:
        log_once(logfile, "step3.header", "Mapping reads...", f"EXECUTING STEP {master_config['map_rule_num']}")
        log_once(logfile, "step3.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step3.outputfolder", f"Output folder: {params.outputfolder}")

        star_version = subprocess.check_output(["STAR", "--version"])
        log_once(logfile, "step3.star_version", "\n" + star_version.decode("utf-8"), "STAR VERSION")
        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['map_rule_num'] - 1], "step3.sanity")

        def run_star_te(genome_path, fastq1, fastq2, paired, threads, outputfolder, sample):
            log_it(logfile, f"Sample {sample}: mapping reads...")
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
                log_it(logfile, f"Staged {path} to {local_path}")
                return local_path

            def reheader_bam(local_bam_path):
                temp_bam_path = local_bam_path + ".tmp"
                header_command = ["samtools", "view", "-H", local_bam_path]
                sed_command = ["sed", "-e", r"s/SN:\([0-9XY]\)/SN:chr\1/", "-e", r"s/SN:MT/SN:chrM/"]

                header_process = subprocess.Popen(header_command, stdout=subprocess.PIPE)
                sed_process = subprocess.Popen(sed_command, stdin=header_process.stdout, stdout=subprocess.PIPE)
                header_process.stdout.close()

                with open(temp_bam_path, "wb") as temp_bam:
                    reheader_process = subprocess.Popen(["samtools", "reheader", "-", local_bam_path], stdin=sed_process.stdout, stdout=temp_bam)
                    sed_process.stdout.close()
                    reheader_process.communicate()

                if not os.path.exists(temp_bam_path):
                    raise FileNotFoundError(f"Failed to create temporary BAM file: {temp_bam_path}")

                os.replace(temp_bam_path, local_bam_path)

            try:
                local_fastq1 = stage_input(fastq1)
                local_fastq2 = stage_input(fastq2) if fastq2 else ""
                local_prefix = os.path.join(local_workdir, f"{sample}.")

                if paired:
                    log_it(logfile, "Running STAR in Paired End mode")
                    star_command = f"""
                        STAR --runThreadN {threads} --genomeDir {quote(genome_path)} \
                        --readFilesIn {quote(local_fastq1)} {quote(local_fastq2)} --readFilesCommand zcat \
                        --outFileNamePrefix {quote(local_prefix)} --outSAMtype BAM Unsorted \
                        --outBAMcompression -1 --genomeLoad NoSharedMemory --twopassMode Basic \
                        --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
                    """
                else:
                    log_it(logfile, "Running STAR in Single End mode")
                    star_command = f"""
                        STAR --runThreadN {threads} --genomeDir {quote(genome_path)} \
                        --readFilesIn {quote(local_fastq1)} --readFilesCommand zcat \
                        --outFileNamePrefix {quote(local_prefix)} --outSAMtype BAM Unsorted \
                        --outBAMcompression -1 --genomeLoad NoSharedMemory --twopassMode Basic --outSAMstrandField intronMotif \
                        --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
                    """

                star_command = " ".join(star_command.split())
                log_it(logfile, star_command, "STAR COMMAND")
                shell(star_command, bench_record=bench_record)

                local_aligned_bam = f"{local_prefix}Aligned.out.bam"
                local_final_bam = os.path.join(local_workdir, f"{sample}.bam")
                local_star_stats = f"{local_prefix}Log.final.out"
                local_final_stats = os.path.join(local_workdir, f"{sample}.STAR_TE_stats.txt")

                if not os.path.exists(local_aligned_bam):
                    raise FileNotFoundError(f"Expected STAR TE BAM output for {sample} at {local_aligned_bam}")
                if not os.path.exists(local_star_stats):
                    raise FileNotFoundError(f"Expected STAR TE stats output for {sample} at {local_star_stats}")

                reheader_bam(local_aligned_bam)
                os.replace(local_aligned_bam, local_final_bam)
                os.replace(local_star_stats, local_final_stats)

                copy_bam_command = f"cp {quote(local_final_bam)} {quote(os.path.join(outputfolder, f'{sample}.bam'))}"
                log_it(logfile, f"Copying BAM for {sample} back to project space...")
                shell(copy_bam_command)

                copy_stats_command = f"cp {quote(local_final_stats)} {quote(os.path.join(outputfolder, f'{sample}.STAR_TE_stats.txt'))}"
                log_it(logfile, f"Copying STAR TE stats for {sample} back to project space...")
                shell(copy_stats_command)
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

            shell(f"""echo "necessity file for aligners. can delete this." > {outputfolder}/{sample}.extra_3.tmp""")

        run_star_te(params.genome_path, input.trimmed_fastq1, input.trimmed_fastq2, params.paired, threads, params.outputfolder, wildcards.sample3)
