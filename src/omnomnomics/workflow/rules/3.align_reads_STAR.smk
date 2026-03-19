#Rule 3 option STAR

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import subprocess
import shutil
rule run_star:
    input:
        trimmed_fastq1= f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample3}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample3}}.trimmed.fastq.gz",
        trimmed_fastq2= f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample3}}_R2.trimmed.fastq.gz" if config['PAIRED'] else []
    output:
        bam=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.bam",
        stats=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.STAR_stats.txt",
        extra=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.extra_3.tmp"
    params:
        genome_path=os.path.join(config['STAR_GENOME_DIR'], f"{config['THEGENOME']}"),
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}",
        paired = config['PAIRED']
    threads:
        Threads_Per_Rule['3']
    resources:
        mem_mb = Memory_Per_Rule['3'],
        partition = master_config['partition']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/benchmarks/{{sample3}}_star_benchmark.tsv"
    run:
        log_it(logfile, "Mapping reads...", f"EXECUTING STEP {master_config['map_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        STAR_version = subprocess.check_output("module load STAR && STAR --version", shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+STAR_version.decode("utf-8"), "STAR VERSION")
        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['map_rule_num']-1])

        def run_star(genome_path, fastq1, fastq2, paired, threads, inputfolder, outputfolder, sample3):
        # Run STAR
            if config['PAIRED']:
                log_it(logfile, "Running STAR in Paired End mode")
                MYNAME = os.path.basename(fastq1).replace("_R1", '')
                log_it(logfile, f"Launching: {MYNAME}...")
                shell(
                    f"""
                    module load STAR && \
                    STAR --runThreadN {threads} --genomeDir {params.genome_path} \
                    --readFilesIn {fastq1} {fastq2} --readFilesCommand zcat \
                    --outFileNamePrefix "{outputfolder}/{sample3}." --outSAMtype BAM Unsorted \
                    --outBAMcompression -1 --genomeLoad NoSharedMemory --twopassMode Basic
                    """, bench_record=bench_record
                )
            else:
                log_it(logfile, "Running STAR in Single End mode")
                MYNAME = os.path.basename(fastq1).replace("_R1", '')
                log_it(logfile, f"Launching: {MYNAME}...")
                shell(
                    f"""
                    module load STAR && \
                    STAR --runThreadN {threads} --genomeDir "{params.genome_path}" \
                    --readFilesIn "{fastq1}" --readFilesCommand zcat \
                    --outFileNamePrefix "{outputfolder}/{sample3}." --outSAMtype BAM Unsorted \
                    --outBAMcompression -1 --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --twopassMode Basic\
                    """, bench_record=bench_record
                )


            bam_files = [f for f in os.listdir(f"{outputfolder}") if f.endswith("Aligned.out.bam")]
            for bam in bam_files:
                bam_path = os.path.join(f"{outputfolder}", bam)
                temp_bam_path = bam_path + ".tmp"
                # Modify header
                header_command = f"module load samtools && samtools view -H {bam_path}"
                sed_command = "sed -e 's/SN:\\([0-9XY]\\)/SN:chr\\1/' -e 's/SN:MT/SN:chrM/'"

                header_process = subprocess.Popen(header_command, shell=True, executable='/bin/bash', stdout=subprocess.PIPE)
                sed_process = subprocess.Popen(sed_command, shell=True, executable='/bin/bash', stdin=header_process.stdout, stdout=subprocess.PIPE)
                header_process.stdout.close()

                with open(temp_bam_path, "wb") as temp_bam:
                    reheader_command = f"module load samtools && samtools reheader - {bam_path}"
                    reheader_process = subprocess.Popen(reheader_command, shell=True, executable='/bin/bash', stdin=sed_process.stdout, stdout=temp_bam)
                    sed_process.stdout.close()
                    reheader_process.communicate()

                # Check if the temporary BAM file was created
                if not os.path.exists(temp_bam_path):
                    raise Exception(f"Failed to create temporary BAM file: {temp_bam_path}")

                # Rename temporary BAM to final BAM
                os.rename(temp_bam_path, bam_path)

                # Rename *.Aligned.out.bam files to *_STAR.bam
                old_path = os.path.join(outputfolder, bam)
                new_path = os.path.join(outputfolder, bam.replace(".Aligned.out.bam", ".bam"))
                os.rename(old_path, new_path)

            # Rename *.Log.final.out files to *.STAR_stats.txt
            log_files = [f for f in os.listdir(f"{outputfolder}") if f.endswith("Log.final.out")]
            for log_file in log_files:
                old_path = os.path.join(outputfolder, log_file)
                new_path = os.path.join(outputfolder, log_file.replace(".Log.final.out", ".STAR_stats.txt"))
                os.rename(old_path, new_path)
            shutil.rmtree("BAM/" + sample3 + "._STARgenome")
            shutil.rmtree("BAM/" + sample3 + "._STARpass1")
            os.remove(f"BAM/{sample3}.Log.out")
            os.remove(f"BAM/{sample3}.Log.progress.out")
            os.remove(f"BAM/{sample3}.SJ.out.tab")
            # os.remove(f"BAM/{sample3}.Aligned.out.bam")
            shell(f"""echo "necessity file for aligners. can delete this." > {outputfolder}/{sample3}.extra_3.tmp""")

        run_star(params.genome_path, input.trimmed_fastq1, input.trimmed_fastq2, params.paired, threads, params.inputfolder, params.outputfolder, wildcards.sample3)
