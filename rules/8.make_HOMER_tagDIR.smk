# Rule 8: Create HOMER tag directories

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob
import pysam

rule create_homer_tagDir:
    input:
        filtered_BAM=f"{experiment_dir}/{master_config['input_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.filtered.bam",
        bai_BAM = f"{experiment_dir}/{master_config['input_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.bai" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.filtered.bam.bai"
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.filtered.HOMER_tagDir.tar.gz",
        f"{experiment_dir}/{master_config['output_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.extra_8.tmp"
    params:
        genome = config['THEGENOME'],
        thetype = config['THETYPE'],
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['tagdir_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['tagdir_rule_num']-1]}"
    threads:
        Threads_Per_Rule['8']
    resources:
        mem_mb = Memory_Per_Rule['8']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['tagdir_rule_num']-1]}/benchmarks/{{sample}}_make_homer_tagdirs_benchmark.tsv"
    run:
        log_it(logfile, f"Creating HOMER tag directories...", f"EXECUTING STEP {master_config['tagdir']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        # path = os.path.join(OMNOM_HOME, "bin", "homer", "configureHomer.pl")
        # version = subprocess.check_output("perl {path} -list 2> /dev/null | grep homer",  shell=True, executable='/bin/bash')
        # log_it(logfile, "\n"+version.decode("utf-8"), "VERSION")

        samtools_version = subprocess.check_output("module load samtools && samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+samtools_version.decode("utf-8"), "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['tagdir_rule_num']-1])

        # Function to create HOMER tag directories
        def create_homer_tagDir(filepath, outputfolder, genome, thetype):
            basename = os.path.basename(filepath)
            sample_name = basename.replace('.bam', '')
            tag_dir = os.path.join(outputfolder, f"{sample_name}.HOMER_tagDir")
            tar_gz_path = f"{tag_dir}.tar.gz"

            # Remove existing output files if they exist
            if os.path.exists(tar_gz_path):
                log_it(logfile, f"Removing existing compressed tag directory {tar_gz_path}")
                os.remove(tar_gz_path)
            if os.path.exists(tag_dir):
                log_it(logfile, f"Removing existing uncompressed tag directory {tag_dir}")
                shell(f"rm -r {tag_dir}")

            if thetype == "RNA":
                with pysam.AlignmentFile(filepath, "rb") as bamfile:
                    paired_end_count = 0
                    for read in bamfile.fetch():
                        if read.is_paired:
                            paired_end_count += 1
                            break

                if paired_end_count > 0:
                    log_it(logfile, f"Running makeTagDirectory on paired-end RNA sample {filepath}")
                    shell(f"""
                        module load samtools && \
                        makeTagDirectory {tag_dir} {filepath} -genome {genome} -sspe -single """)
                else:
                    log_it(logfile, f"Running makeTagDirectory on {filepath}")
                    shell(f"""
                        module load samtools && \
                        makeTagDirectory {tag_dir} {filepath} -genome {genome} -single """)
            else:
                log_it(logfile, f"Running makeTagDirectory on {filepath}")
                shell(f"""
                    module load samtools && \
                    makeTagDirectory {tag_dir} {filepath} -genome {genome} -single""")

            # Compress the tag directory
            log_it(logfile, f"Compressing {os.path.basename(tag_dir)} into {tag_dir}.tar.gz")
            shell(f"cd {outputfolder} && tar czf {sample_name}.HOMER_tagDir.tar.gz  {sample_name}.HOMER_tagDir")
            
            # Remove the uncompressed tag directory
            log_it(logfile, f"Removing uncompressed tag directory {tag_dir}")
            shell(f"rm -r {tag_dir}")

            shell(f"""echo "necessity file for touchup_bam. can delete this." > {outputfolder}/{wildcards.sample}.extra_8.tmp""")

        create_homer_tagDir(input.filtered_BAM, params.outputfolder, params.genome, params.thetype)