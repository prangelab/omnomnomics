# Rule 2 fastqc

## Omnomnomics Snake Rule  ##
import os

rule run_fastqc:
   input:
        # fastq1="trimmed_FASTQ/{sample}_R1.trimmed.fastq.gz",
        # fastq2="trimmed_FASTQ/{sample}_R2.trimmed.fastq.gz" if config["PAIRED"] else None
        fastq1=f"{master_config['input_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.fastq.gz",
        fastq2=f"{master_config['input_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.fastq.gz" if config["PAIRED"] else None
   output:
    #    report1="fastqc_reports/{sample}_R1_fastqc.html",
    #    report2="fastqc_reports/{sample}_R2_fastqc.html" if config["PAIRED"] else None
        report1=f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.trimmed_fastq.html",
        report2=f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.trimmed_fastq.html" if config["PAIRED"] else None
        report3=f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.trimmed_fastq.zip",
        report4=f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.trimmed_fastq.zip" if config["PAIRED"] else None
   params:
        outputfolder = master_config['output_folders'][master_config['qc_rule_num']-1]
   threads:
       Threads_Per_Rule['2']
   resources:
       mem_mb = Memory_Per_Rule['2']
   run:
       def run_fastqc(threads, fastq1, fastq2, outputfolder):
           log_it(logfile, "Generating FastQC reports...", f"EXECUTING STEP {master_config['qc']}")
           log_it(logfile, "Input folder: trimmed_FASTQ")
           log_it(logfile, "Output folder: fastqc_reports")


           if fastq2:
               log_it(logfile, "Running FastQC in paired end mode...")
               fastqc_command = f"""
                    module load fastqc && \
                    fastqc -t {threads} -o {outputfolder} {fastq1} {fastq2}
               """
           else:
               log_it(logfile, "Running FastQC in single end mode...")
               fastqc_command = f"""
                    module load fastqc && \
                    fastqc -t {threads} -o {outputfolder} {fastq1}
               """

           # Run the FastQC command
           shell(fastqc_command)

       # Call the function with parameters
       run_fastqc(threads, input.fastq1, input.fastq2, params.outputfolder)