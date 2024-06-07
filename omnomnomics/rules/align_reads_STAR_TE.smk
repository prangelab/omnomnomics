#Rule 3 option STAR TE

## Omnomnomics Snake Rule  ##
import os
import subprocess
rule run_star_te:
    input:
        trimmed_fastqc1= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R1{'_Skewer' if config['THETRIMTOOL'] == 'skewer' else ('_Trimmomatic' if config['THETRIMTOOL'] == 'trimmomatic' else '')}.trimmed.fastq.gz",
        trimmed_fastqc2= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R2{'_Skewer' if config['THETRIMTOOL'] == 'skewer' else ('_Trimmomatic' if config['THETRIMTOOL'] == 'trimmomatic' else '')}.trimmed.fastq.gz" if config['PAIRED'] else None
    output:
        bam=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample}}_STAR_TE.bam",
        stats=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample}}.STAR_TE_stats.txt"
    params:
        genome_path=os.path.join(f"{config['OMNOM_HOME']}", "genomes", "STAR", f"{config['THEGENOME']}"),
        inputfolder = master_config['input_folders'][master_config['map_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['map_rule_num']-1]
    threads:
        Threads_Per_Rule['3']
    resources:
        mem_mb = Memory_Per_Rule['3']
    run:
        logIt(logfile, "Mapping reads...", f"EXECUTING STEP {master_config['map_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")


        def run_star(genome_path, fastq1, fastq2, paired, threads,inputfolder, outputfolder, sample):
            if config['PAIRED']:
                logIt(logfile, "Running STAR in Paired End mode")
                MYNAME = os.path.basename(input.fastq1)
                logIt(logfile, f"Launching: {MYNAME}...")
                shell(
                    f"""
                    STAR --runThreadN {threads} --genomeDir "{params.genome_path}" \
                    --readFilesIn "{input.fastq1}" "{input.fastq2}" --readFilesCommand zcat \
                    --outFileNamePrefix "{outputfolder}/{sample}." --outSAMtype BAM Unsorted \
                    --outBAMcompression -1 --genomeLoad LoadAndRemove --twopassMode Basic\
                    --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
                    """
                )
            else:
                logIt(logfile, "Running STAR in Single End mode")
                MYNAME = os.path.basename(input.fastq1)
                logIt(logfile, f"Launching: {MYNAME}...")
                shell(
                    f"""
                    STAR --runThreadN {threads} --genomeDir "{params.genome_path}" \
                    --readFilesIn "{input.fastq1}" --readFilesCommand zcat \
                    --outFileNamePrefix "{outputfolder}/{sample}." --outSAMtype BAM Unsorted \
                    --outBAMcompression -1 --genomeLoad LoadAndRemove --twopassMode Basic --outSAMstrandField intronMotif \
                    --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
                    """
                )

            bam_files = [f for f in os.listdir(f"{outputfolder}") if f.endswith("Aligned.out.bam")]
            for bam in bam_files:
                bam_path = os.path.join(f"{outputfolder}", "bam")
                temp_bam_path = bam_path + ".tmp"
                # Modify header
                header_process = subprocess.Popen(["samtools", "view", "-H", bam_path], stdout=subprocess.PIPE)
                sed_process = subprocess.Popen(
                    ["sed", "-e", "s/SN:\\([0-9XY]\\)/SN:chr\\1/", "-e", "s/SN:MT/SN:chrM/"], stdin=header_process.stdout, stdout=subprocess.PIPE)
                header_process.stdout.close()
                with open(temp_bam_path, "wb") as temp_bam:
                    reheader_process = subprocess.Popen(["samtools", "reheader", "-", bam_path], stdin=sed_process.stdout, stdout=temp_bam)
                    sed_process.stdout.close()
                    reheader_process.communicate()
                # Rename temporary BAM to final BAM
                os.rename(temp_bam_path, bam_path)

                # Rename *.Aligned.out.bam files to *_STAR.bam
                old_path = os.path.join(outputfolder, bam)
                new_path = os.path.join(outputfolder, bam.replace("Aligned.out.bam", "_STAR_TE.bam"))
                os.rename(old_path, new_path)

            # Rename *.Log.final.out files to *.STAR_stats.txt
            log_files = [f for f in os.listdir(f"{outputfolder}") if f.endswith("Log.final.out")]
            for log_file in log_files:
                old_path = os.path.join(outputfolder, log_file)
                new_path = os.path.join(outputfolder, log_file.replace(".Log.final.out", ".STAR_TE_stats.txt"))
                os.rename(old_path, new_path)


        run_star(params.genome_path, input.fastq1, input.fastq2, params.paired, threads, params.inputfolder, params.outputfolder, wildcards.sample)





