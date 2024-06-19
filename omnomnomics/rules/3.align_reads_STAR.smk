#Rule 3 option STAR

## Omnomnomics Snake Rule  ##
import os
import subprocess
import shutil
rule run_star:
    input:
        # trimmed_fastq1= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R1{'_Skewer' if config['THETRIMTOOL'] == 'skewer' else ('_Trimmomatic' if config['THETRIMTOOL'] == 'trimmomatic' else '')}.trimmed.fastq.gz",
        # trimmed_fastq2= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R2{'_Skewer' if config['THETRIMTOOL'] == 'skewer' else ('_Trimmomatic' if config['THETRIMTOOL'] == 'trimmomatic' else '')}.trimmed.fastq.gz" if config['PAIRED'] else None
        trimmed_fastq1= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R2.trimmed.fastq.gz" if config['PAIRED'] else []
    output:
        bam=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample}}.bam",
        stats=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample}}.STAR_stats.txt",
        extra=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample}}.extra_3.tmp"
    params:
        genome_path=os.path.join(f"{config['OMNOM_HOME']}", "genomes", "STAR", f"{config['THEGENOME']}"),
        inputfolder = master_config['input_folders'][master_config['map_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['map_rule_num']-1],
        paired = config['PAIRED']
    threads:
        Threads_Per_Rule['3']
    resources:
        #mem_mb = Memory_Per_Rule['3']
        mem_mb = 128000
    benchmark:
        f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample}}_star_benchmark.tsv"
    run:
        log_it(logfile, "Mapping reads...", f"EXECUTING STEP {master_config['map_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        STAR_version = subprocess.check_output("module load STAR && STAR --version", shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+STAR_version.decode("utf-8"), "STAR VERSION")
        print(STAR_version.decode("utf-8"))
        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['map_rule_num']-1])

        def run_star(genome_path, fastq1, fastq2, paired, threads, inputfolder, outputfolder, sample):
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
                    --outFileNamePrefix "{outputfolder}/{sample}." --outSAMtype BAM Unsorted \
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
                    --outFileNamePrefix "{outputfolder}/{sample}." --outSAMtype BAM Unsorted \
                    --outBAMcompression -1 --genomeLoad NoSharedMemory --outSAMstrandField intronMotif --twopassMode Basic\
                    """, bench_record=bench_record
                )


            bam_files = [f for f in os.listdir(f"{outputfolder}") if f.endswith("Aligned.out.bam")]
            print(bam_files)
            for bam in bam_files:
                print(bam)
                bam_path = os.path.join(f"{outputfolder}", bam)
                print(bam_path)
                temp_bam_path = bam_path + ".tmp"
                print(temp_bam_path)
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
                print(old_path)
                new_path = os.path.join(outputfolder, bam.replace(".Aligned.out.bam", ".bam"))
                print(new_path)
                os.rename(old_path, new_path)

            # Rename *.Log.final.out files to *.STAR_stats.txt
            log_files = [f for f in os.listdir(f"{outputfolder}") if f.endswith("Log.final.out")]
            for log_file in log_files:
                old_path = os.path.join(outputfolder, log_file)
                new_path = os.path.join(outputfolder, log_file.replace(".Log.final.out", ".STAR_stats.txt"))
                os.rename(old_path, new_path)
            shutil.rmtree("BAM/" + sample + "._STARgenome")
            shutil.rmtree("BAM/" + sample + "._STARpass1")
            os.remove(f"BAM/{sample}.Log.out")
            os.remove(f"BAM/{sample}.Log.progress.out")
            os.remove(f"BAM/{sample}.SJ.out.tab")
            # os.remove(f"BAM/{sample}.Aligned.out.bam")
            shell(f"""echo "necessity file for aligners. can delete this." > {outputfolder}/{sample}.extra_3.tmp""")

        run_star(params.genome_path, input.trimmed_fastq1, input.trimmed_fastq2, params.paired, threads, params.inputfolder, params.outputfolder, wildcards.sample)