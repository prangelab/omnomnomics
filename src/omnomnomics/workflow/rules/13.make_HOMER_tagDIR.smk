# Optional HOMER tag directory export

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob
import pysam
import shlex
import shutil
import subprocess
import tempfile

rule create_homer_tagDir:
    input:
        filtered_BAM=f"{experiment_dir}/{master_config['input_folders'][master_config['homer_tagdir_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['homer_tagdir_rule_num']-1]}/{{sample}}.filtered.bam",
        bai_BAM = f"{experiment_dir}/{master_config['input_folders'][master_config['homer_tagdir_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.bai" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['homer_tagdir_rule_num']-1]}/{{sample}}.filtered.bam.bai"
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['homer_tagdir_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['homer_tagdir_rule_num']-1]}/{{sample}}.filtered.HOMER_tagDir.tar.gz",
        f"{experiment_dir}/{master_config['output_folders'][master_config['homer_tagdir_rule_num']-1]}/{{sample}}.extra_13.tmp"
    params:
        genome = config['THEGENOME'],
        thetype = config['THETYPE'],
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['homer_tagdir_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['homer_tagdir_rule_num']-1]}"
    threads:
        Threads_Per_Rule['13']
    resources:
        mem_mb = Memory_Per_Rule['13'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['13']
    run:
        log_once(logfile, "step13.header", "Creating optional HOMER tag directory exports...", f"EXECUTING STEP {master_config['homer_tagdir_rule_num']}")
        log_once(logfile, "step13.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step13.outputfolder", f"Output folder: {params.outputfolder}")

        samtools_version = subprocess.check_output("samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_once(logfile, "step13.samtools_version", "\n"+samtools_version.decode("utf-8"), "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['homer_tagdir_rule_num']-1], "step13.homer_sanity")

        # Function to create HOMER tag directories
        def create_homer_tagDir(filepath, outputfolder, genome, thetype):
            tracking = begin_step_sample(master_config['homer_tagdir_rule_num'], wildcards.sample, "create_homer_tagDir")
            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{wildcards.sample}.", dir=tmpdir_root if tmpdir_root else None)
            record_step_note(master_config['homer_tagdir_rule_num'], wildcards.sample, f"scratch_dir={local_workdir}")

            def quote(path):
                return shlex.quote(path)

            basename = os.path.basename(filepath)
            sample_name = basename.replace('.bam', '')
            local_bam = os.path.join(local_workdir, basename)
            local_bai = os.path.join(local_workdir, os.path.basename(input.bai_BAM))
            stage_bam_command = f'cp {quote(filepath)} {quote(local_bam)}'
            stage_bai_command = f'cp {quote(input.bai_BAM)} {quote(local_bai)}'
            record_step_note(master_config['homer_tagdir_rule_num'], wildcards.sample, "staging_bam_and_index")
            shell(stage_bam_command)
            shell(stage_bai_command)

            tag_dir = os.path.join(local_workdir, f"{sample_name}.HOMER_tagDir")
            tar_gz_path = f"{tag_dir}.tar.gz"

            try:
                # Remove existing output files if they exist
                if os.path.exists(tar_gz_path):
                    log_it(logfile, f"Removing existing compressed tag directory {tar_gz_path}")
                    os.remove(tar_gz_path)
                if os.path.exists(tag_dir):
                    log_it(logfile, f"Removing existing uncompressed tag directory {tag_dir}")
                    shell(f"rm -r {quote(tag_dir)}")

                if thetype == "RNA":
                    with pysam.AlignmentFile(local_bam, "rb") as bamfile:
                        paired_end_count = 0
                        for read in bamfile.fetch():
                            if read.is_paired:
                                paired_end_count += 1
                                break

                    if paired_end_count > 0:
                        record_step_note(master_config['homer_tagdir_rule_num'], wildcards.sample, "running_makeTagDirectory_paired_rna")
                        command = f"makeTagDirectory {quote(tag_dir)} {quote(local_bam)} -genome {genome} -sspe"
                    else:
                        record_step_note(master_config['homer_tagdir_rule_num'], wildcards.sample, "running_makeTagDirectory_single_rna")
                        command = f"makeTagDirectory {quote(tag_dir)} {quote(local_bam)} -genome {genome} -single"
                else:
                    record_step_note(master_config['homer_tagdir_rule_num'], wildcards.sample, "running_makeTagDirectory")
                    command = f"makeTagDirectory {quote(tag_dir)} {quote(local_bam)} -genome {genome} -single"

                record_step_command(master_config['homer_tagdir_rule_num'], wildcards.sample, command)
                shell(command)

                # Compress the tag directory
                record_step_note(master_config['homer_tagdir_rule_num'], wildcards.sample, "compressing_homer_tagdir")
                shell(f"cd {quote(local_workdir)} && tar czf {quote(os.path.basename(tar_gz_path))} {quote(os.path.basename(tag_dir))}")

                shell(f'cp {quote(tar_gz_path)} {quote(os.path.join(outputfolder, os.path.basename(tar_gz_path)))}')
                shell(f"""echo "necessity file for homer_tagdir export. can delete this." > {quote(os.path.join(outputfolder, f"{wildcards.sample}.extra_13.tmp"))}""")
                finish_step_sample(master_config['homer_tagdir_rule_num'], wildcards.sample, "create_homer_tagDir", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['homer_tagdir_rule_num'], wildcards.sample, "create_homer_tagDir", tracking["start_time"], "FAIL")
                raise
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

        create_homer_tagDir(input.filtered_BAM, params.outputfolder, params.genome, params.thetype)
