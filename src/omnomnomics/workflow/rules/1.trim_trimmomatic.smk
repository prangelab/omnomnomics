# Rule 1 option trimmomatic

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob

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
        trim_heap=config["THEHEAPINIT"],
        trim_mem=config["THEMEM"],
        seq_type=config["THETYPE"],
        inputfolder =f"{experiment_dir}/{master_config['input_folders'][master_config['trim_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}"
    threads:
        Threads_Per_Rule['1']
    resources:
        mem_mb = Memory_Per_Rule['1'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['1']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/benchmarks/{{sample}}_trimmomatic_benchmark.tsv"
    run:
        def run_trimmomatic(trim_tool, trim_heap, trim_mem, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample):
            log_it(logfile, "Trimming Reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_it(logfile, f"Input folder: {inputfolder}")
            log_it(logfile, f"Output folder: {outputfolder}")
            log_it(logfile, f"Trim Tool: {trim_tool}")
            log_it(logfile, f"Trim Heap: {trim_heap}")
            log_it(logfile, f"Trim Mem: {trim_mem}")

            path = os.path.join(config['WORKFLOW_ROOT'], "bin", "Trimmomatic-0.39", "trimmomatic-0.39.jar")
            trimmomatic_version = subprocess.check_output(f"module load java && java -jar {path} -version",  shell=True, executable='/bin/bash')
            log_it(logfile, "\n"+trimmomatic_version.decode("utf-8"), "trimmomatic VERSION")

            sanity_check_dir(logfile, inputfolder, master_config['input_file_types'][master_config['trim_rule_num']-1])

            if seq_type == "ATAC":
                adapter_file = "NexteraPE-PE.fa" if config["PAIRED"] else "Nextera-SE.fa"
            else:
                adapter_file = "TruSeq3-PE.fa" if config["PAIRED"] else "TruSeq3-SE.fa"

            if fastq2:
                log_it(logfile, f"Running {trim_tool} in Paired End mode.")
                out_base = f"{outputfolder}/{sample}.trimmed.fastq.gz"
                java_command = f"""
                    module load java && \
                    java -Xms{trim_heap} -Xmx{trim_mem} -jar "{config['WORKFLOW_ROOT']}/bin/Trimmomatic-0.39/trimmomatic-0.39.jar" PE \
                    -threads {threads} -baseout {out_base} {fastq1} {fastq2} \
                    ILLUMINACLIP:{config['WORKFLOW_ROOT']}/bin/Trimmomatic-0.39/adapters/{adapter_file}:2:30:10:2:True \
                    LEADING:3 TRAILING:3 MINLEN:36
                """
            else:
                log_it(logfile, f"Running {trim_tool} in Single End mode.")
                out_base = f"{outputfolder}/{sample}.trimmed.fastq.gz"
                java_command = f"""
                    module load java && \
                    java -Xms{trim_heap} -Xmx{trim_mem} -jar "{config['WORKFLOW_ROOT']}/bin/Trimmomatic-0.39/trimmomatic-0.39.jar" SE \
                    -threads {threads} {fastq1} {out_base} \
                    ILLUMINACLIP:{config['WORKFLOW_ROOT']}/bin/Trimmomatic-0.39/adapters/{adapter_file}:2:30:10 \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
                """
            
            # Run the trimmomatic command
            shell(java_command)

            # Rename R1 trimmed files
            log_it(logfile, "Renaming trimmed results ...")
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*1P.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('.trimmed_1P.fastq.gz', '_R1.trimmed.fastq.gz'))
                os.rename(file_path, new_name)

            # Rename R2 trimmed files
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*2P.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('.trimmed_2P.fastq.gz', '_R2.trimmed.fastq.gz'))
                os.rename(file_path, new_name)

            # Remove unpaired files if PAIRED is 1
            if config['PAIRED']:
                log_it(logfile, "Removing unpaired files")
                for file_path in glob.glob(os.path.join(f"{outputfolder}", '*U.fastq.gz')):
                    os.remove(file_path)

        # Call the function with parameters
        run_trimmomatic(params.trim_tool, (str(int(params.trim_heap[:-1])//2)+'M'), (str(int(params.trim_mem[:-1])//2)+'M'), params.seq_type, threads, input.fastq1, input.fastq2, params.inputfolder, params.outputfolder, wildcards.sample)
