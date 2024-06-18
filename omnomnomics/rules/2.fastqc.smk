# Rule 2 fastqc

## Omnomnomics Snake Rule  ##
import os
# def get_fastqc_input(wildcards):
#     input_folder = master_config['input_folders'][master_config['qc_rule_num']-1]
#     #sample = wildcards.sample
#     semaphore_path = os.path.join(experiment_dir, "omnomnomics.semaphore")

#     suffix = ''
#     with open(semaphore_path, 'r') as file:
#         lines = file.readlines()
#         if len(lines) >= 2 and lines[1].strip() == "1":  
#             if config['THEMAPTOOL'] == 'skewer':
#                 suffix = '_Skewer'
#             elif config['THEMAPTOOL'] == 'trimmomatic':
#                 suffix = '_Trimmomatic'

#     if config['PAIRED']:
#         fastq1 = f"{input_folder}/{{sample}}_R1{suffix}.trimmed.fastq.gz"
#         fastq2 = f"{input_folder}/{{sample}}_R2{suffix}.trimmed.fastq.gz"
#         print(f"DEBUG: Paired input for {sample}: {fastq1}, {fastq2}")
#         return [fastq1, fastq2]
#     else:
#         fastq = f"{input_folder}/{{sample}}{suffix}.trimmed.fastq.gz"
#         print(f"DEBUG: Single-end input for {sample}: {fastq}")
#         return [fastq, ]
rule run_fastqc:
    input:
        # fastq1="trimmed_FASTQ/{sample}_R1.trimmed.fastq.gz",
        # fastq2="trimmed_FASTQ/{sample}_R2.trimmed.fastq.gz" if config["PAIRED"] else None
        # fastq1=f"{master_config['input_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.fastq.gz",
        # fastq2=f"{master_config['input_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.fastq.gz" if config["PAIRED"] else None
        # get_fastqc_input
        trimmed_fastq1= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R2.trimmed.fastq.gz" if config['PAIRED'] else None
    output:
    #    report1="fastqc_reports/{sample}_R1_fastqc.html",
    #    report2="fastqc_reports/{sample}_R2_fastqc.html" if config["PAIRED"] else None
        report1=f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.trimmed_fastqc.html" if config["PAIRED"] else f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}.trimmed_fastqc.html",
        report2=f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.trimmed_fastqc.html" if config["PAIRED"] else None,
        report3=f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R1.trimmed_fastqc.zip" if config["PAIRED"] else f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}.trimmed_fastqc.zip",
        report4=f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_R2.trimmed_fastqc.zip" if config["PAIRED"] else None

    params:
        inputfolder = master_config['input_folders'][master_config['qc_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['qc_rule_num']-1]
    threads:
        Threads_Per_Rule['2']
    resources:
        mem_mb = Memory_Per_Rule['2']
    benchmark:
        f"{master_config['output_folders'][master_config['qc_rule_num']-1]}/{{sample}}_fastqc_benchmark.tsv"
    run:
        def run_fastqc(threads, input, outputfolder):
            log_it(logfile, "Generating FastQC reports...", f"EXECUTING STEP {master_config['qc']}")
            log_it(logfile, "Input folder: trimmed_FASTQ")
            log_it(logfile, "Output folder: fastqc_reports")

            fastqc_version = subprocess.check_output("module load fastqc && fastqc --version", shell=True, executable='/bin/bash')
            log_it(logfile, "\n"+fastqc_version.decode("utf-8"), "FASTQC VERSION")
            print(fastqc_version.decode("utf-8"))
            sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['qc_rule_num']-1])

            if config['PAIRED']:
                log_it(logfile, "Running FastQC in paired end mode...")
                fastqc_command = f"""
                    module load fastqc && \
                    fastqc -t {threads} -o {outputfolder} {input.trimmed_fastq1} {input.trimmed_fastq2}
                """
            else:
                log_it(logfile, "Running FastQC in single end mode...")
                fastqc_command = f"""
                    module load fastqc && \
                    fastqc -t {threads} -o {outputfolder} {input.trimmed_fastq1}
                """

            # Run the FastQC command
            shell(fastqc_command)

        # Call the function with parameters
        run_fastqc(threads, input, params.outputfolder)