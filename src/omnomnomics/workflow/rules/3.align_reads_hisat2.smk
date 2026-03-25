#Rule 3 option HISAT2

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

rule run_hisat2:
    wildcard_constraints:
        sample3=lane_sample_wildcard_pattern
    input:
        trimmed_fastq1= f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample3}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample3}}.trimmed.fastq.gz",
        trimmed_fastq2= f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample3}}_R2.trimmed.fastq.gz" if config['PAIRED'] else []
    output:
        bam=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.bam",
        stats=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.HISAT2_stats.txt",
        extra=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.extra_3.tmp"
        #two extra  output files are create if keepunpaired = 1, but not necesarry to specify because will be made then automatically
    params:
        genome_path=os.path.join(config['HISAT2_GENOME_DIR'], f"{config['THEGENOME']}"),
        keepunpaired=config.get("KEEPUNPAIRED", "0"),
        seq_type=config["THETYPE"],
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}"
    threads:
        Threads_Per_Rule['3']
    resources:
        mem_mb = Memory_Per_Rule['3'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['3']
    run:
        log_once(logfile, "step3.header", "Mapping reads...", f"EXECUTING STEP {master_config['map_rule_num']}")
        log_once(logfile, "step3.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step3.outputfolder", f"Output folder: {params.outputfolder}")

        hisat2_version = subprocess.check_output(["hisat2", "--version"]).decode("utf-8").splitlines()[0]
        log_once(logfile, "step3.hisat2_version", f"\n{hisat2_version}", "HISAT2 VERSION")
        
        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['map_rule_num']-1], "step3.sanity")

        def run_hisat2(seq_type, threads, genome_path, fastq1, fastq2, keepunpaired, inputfolder, outputfolder, sample):
            tracking = begin_step_sample(master_config['map_rule_num'], sample, "run_hisat2")
            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample}.", dir=tmpdir_root if tmpdir_root else None)
            record_step_note(master_config['map_rule_num'], sample, f"scratch_dir={local_workdir}")

            def quote(path):
                return shlex.quote(path)

            def stage_input(path):
                local_path = os.path.join(local_workdir, os.path.basename(path))
                command = f'cp {quote(path)} {quote(local_path)}'
                record_step_note(master_config['map_rule_num'], sample, f"staging {os.path.basename(path)}")
                shell(command)
                return local_path

            try:
                if seq_type == "RNA":
                    if config["PAIRED"]:
                        record_step_note(master_config['map_rule_num'], sample, "running_hisat2_paired_end_rna")
                        MYNAME = os.path.basename(fastq1)
                        record_step_note(master_config['map_rule_num'], sample, f"launching {MYNAME}")
                        local_fastq1 = stage_input(fastq1)
                        local_fastq2 = stage_input(fastq2)
                        local_bam = os.path.join(local_workdir, f"{sample}.bam")
                        local_stats = os.path.join(local_workdir, f"{sample}.HISAT2_stats.txt")
                        if keepunpaired:
                            hisat2_command = f"""
                            hisat2 -p {threads} -x {genome_path} \
                            -1 {local_fastq1} -2 {local_fastq2} --mm --add-chrname --new-summary --dta-cufflinks \
                            --summary-file "{local_stats}" \
                            --un-gz "{outputfolder}/{sample}.unpaired.unaligned.bam" --al-gz "{outputfolder}/{sample}.unpaired.aligned.bam" \
                            | samtools view -b - 1> "{local_bam}"
                            """
                        else:
                            hisat2_command = f"""
                            hisat2 -p {threads} -x {genome_path} \
                            -1 {local_fastq1} -2 {local_fastq2} --mm --add-chrname --new-summary --dta-cufflinks \
                            --summary-file "{local_stats}" \
                            | samtools view -b - 1> "{local_bam}"
                            """
                        hisat2_command = " ".join(hisat2_command.split())
                        record_step_command(master_config['map_rule_num'], sample, hisat2_command)
                        shell(hisat2_command)
                        shell(f'cp {quote(local_bam)} {quote(os.path.join(outputfolder, f"{sample}.bam"))}')
                        shell(f'cp {quote(local_stats)} {quote(os.path.join(outputfolder, f"{sample}.HISAT2_stats.txt"))}')
                    else:
                        record_step_note(master_config['map_rule_num'], sample, "running_hisat2_single_end_rna")
                        local_fastq1 = stage_input(fastq1)
                        local_bam = os.path.join(local_workdir, f"{sample}.bam")
                        local_stats = os.path.join(local_workdir, f"{sample}.HISAT2_stats.txt")
                        hisat2_command = f"""
                        hisat2 -p {threads} -x {genome_path} \
                        -U {local_fastq1} --mm --add-chrname --new-summary --dta-cufflinks \
                        --summary-file "{local_stats}" \
                        | samtools view -b - 1> "{local_bam}"
                        """
                        hisat2_command = " ".join(hisat2_command.split())
                        record_step_command(master_config['map_rule_num'], sample, hisat2_command)
                        shell(hisat2_command)
                        shell(f'cp {quote(local_bam)} {quote(os.path.join(outputfolder, f"{sample}.bam"))}')
                        shell(f'cp {quote(local_stats)} {quote(os.path.join(outputfolder, f"{sample}.HISAT2_stats.txt"))}')
                else: #if ChIP- or ATAC data
                    if config["PAIRED"]:
                        record_step_note(master_config['map_rule_num'], sample, "running_hisat2_paired_end_chip_or_atac")
                        MYNAME = os.path.basename(fastq1)
                        record_step_note(master_config['map_rule_num'], sample, f"launching {MYNAME}")
                        local_fastq1 = stage_input(fastq1)
                        local_fastq2 = stage_input(fastq2)
                        local_bam = os.path.join(local_workdir, f"{sample}.bam")
                        local_stats = os.path.join(local_workdir, f"{sample}.HISAT2_stats.txt")
                        if keepunpaired:
                            hisat2_command = f"""
                            hisat2 -p {threads} -x {genome_path} \
                            -1 {local_fastq1} -2 {local_fastq2} --mm --add-chrname --new-summary --no-spliced-alignment \
                            --un-gz "BAM/{sample}.unpaired.unaligned.bam" --al-gz "BAM/{sample}.unpaired.aligned.bam"\
                            --summary-file "{local_stats}" \
                            | samtools view -b - 1> "{local_bam}"
                            """
                        else:
                            hisat2_command = f"""
                            hisat2 -p {threads} -x {genome_path} \
                            -1 {local_fastq1} -2 {local_fastq2} --mm --add-chrname --new-summary --no-spliced-alignment \
                            --summary-file "{local_stats}" \
                            | samtools view -b - 1> "{local_bam}"
                            """
                        hisat2_command = " ".join(hisat2_command.split())
                        record_step_command(master_config['map_rule_num'], sample, hisat2_command)
                        shell(hisat2_command)
                        shell(f'cp {quote(local_bam)} {quote(os.path.join(outputfolder, f"{sample}.bam"))}')
                        shell(f'cp {quote(local_stats)} {quote(os.path.join(outputfolder, f"{sample}.HISAT2_stats.txt"))}')
                    else:
                        record_step_note(master_config['map_rule_num'], sample, "running_hisat2_single_end_chip_or_atac")
                        local_fastq1 = stage_input(fastq1)
                        local_bam = os.path.join(local_workdir, f"{sample}.bam")
                        local_stats = os.path.join(local_workdir, f"{sample}.HISAT2_stats.txt")
                        hisat2_command = f"""
                        hisat2 -p {threads} -x {genome_path} \
                        -U {local_fastq1} --mm --add-chrname --new-summary --no-spliced-alignment \
                        --summary-file "{local_stats}" \
                        | samtools view -b - 1> "{local_bam}"
                        """
                        hisat2_command = " ".join(hisat2_command.split())
                        record_step_command(master_config['map_rule_num'], sample, hisat2_command)
                        shell(hisat2_command)
                        shell(f'cp {quote(local_bam)} {quote(os.path.join(outputfolder, f"{sample}.bam"))}')
                        shell(f'cp {quote(local_stats)} {quote(os.path.join(outputfolder, f"{sample}.HISAT2_stats.txt"))}')
                finish_step_sample(master_config['map_rule_num'], sample, "run_hisat2", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['map_rule_num'], sample, "run_hisat2", tracking["start_time"], "FAIL")
                raise
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)
            shell(f"""echo "necessity file for aligners. can delete this." > {outputfolder}/{sample}.extra_3.tmp""")
        # Call the function with parameters
        run_hisat2(params.seq_type, threads, params.genome_path, input.trimmed_fastq1, input.trimmed_fastq2, params.keepunpaired, params.inputfolder, params.outputfolder, wildcards.sample3)
