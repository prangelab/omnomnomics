# Rule 1 option trimmomatic

## Omnomnomics Snake Rule  ##
import os
import glob

rule run_trimmomatic:
    input:
        fastq1=f"{master_config['input_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R1.fastq.gz",
        fastq2=f"{master_config['input_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R2.fastq.gz" if config["PAIRED"] else None
        #still account for unpaired scenario
    output:
        trimmed_fastq1=f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R1_Trimmomatic.trimmed.fastq.gz",
        trimmed_fastq2=f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R2_Trimmomatic.trimmed.fastq.gz" if config["PAIRED"] else None
        #still account for unpaired scenario
    params:
        trim_tool=config["THETRIMTOOL"],
        trim_heap=config["THEHEAPINIT"],
        trim_mem=config["THEMEM"],
        seq_type=config["THETYPE"],
        inputfolder = master_config['input_folders'][master_config['trim_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['trim_rule_num']-1]
    threads:
        Threads_Per_Rule['1']
    resources:
        #mem_mb = 512000
        mem_mb = Memory_Per_Rule['1']
    run:
        def run_trimmomatic(trim_tool, trim_heap, trim_mem, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample):
            log_it(logfile, "Trimming Reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_it(logfile, f"Input folder: {inputfolder}")
            log_it(logfile, f"Output folder: {outputfolder}")
            log_it(logfile, f"Trim Tool: {trim_tool}")
            log_it(logfile, f"Trim Heap: {trim_heap}")
            log_it(logfile, f"Trim Mem: {trim_mem}")

            if seq_type == "ATAC":
                adapter_file = "NexteraPE-PE.fa" if config["PAIRED"] else "Nextera-SE.fa"
            else:
                adapter_file = "TruSeq3-PE.fa" if config["PAIRED"] else "TruSeq3-SE.fa"

            if fastq2:
                log_it(logfile, f"Running {trim_tool} in Paired End mode.")
                out_base = f"{outputfolder}/{sample}.trimmed.fastq.gz"
                java_command = f"""
                    module load java && \
                    java -Xms{trim_heap} -Xmx{trim_mem} -jar $OMNOM_HOME/bin/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                    -threads {threads} -baseout {out_base} {fastq1} {fastq2} \
                    ILLUMINACLIP:$OMNOM_HOME/bin/Trimmomatic-0.39/adapters/{adapter_file}:2:30:10:2:True \
                    LEADING:3 TRAILING:3 MINLEN:36
                """
            else:
                log_it(logfile, f"Running {trim_tool} in Single End mode.")
                out_base = f"{outputfolder}/{sample}.trimmed.fastq.gz"
                java_command = f"""
                    module load java && \
                    java -Xms{trim_heap} -Xmx{trim_mem} -jar $OMNOM_HOME/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
                    -threads {threads} {fastq1} {out_base} \
                    ILLUMINACLIP:$OMNOM_HOME/bin/Trimmomatic-0.39/adapters/{adapter_file}:2:30:10 \
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
                """
            
            # Run the trimmomatic command
            shell(java_command)

            # Rename R1 trimmed files
            log_it(logfile, "Renaming trimmed results ...")
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*1P.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('.trimmed_1P.fastq.gz', '_R1_Trimmomatic.trimmed.fastq.gz'))
                os.rename(file_path, new_name)

            # Rename R2 trimmed files
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*2P.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('.trimmed_2P.fastq.gz', '_R2_Trimmomatic.trimmed.fastq.gz'))
                os.rename(file_path, new_name)

            # Remove unpaired files if PAIRED is 1
            if config['PAIRED']:
                log_it(logfile, "Removing unpaired files")
                for file_path in glob.glob(os.path.join(f"{outputfolder}", '*U.fastq.gz')):
                    os.remove(file_path)
            
            # Rename trimmed files for Trimmomatic
            # for filename in os.listdir('trimmed_FASTQ'):
            #     file_path = os.path.join(directory, filename)
            #     if os.path.isfile(file_path):
            #         base, ext = os.path.splitext(filename)
            #         new_filename = f"{base}_Trimmomatic{ext}"
            #         new_file_path = os.path.join(directory, new_filename)
            #         os.rename(file_path, new_file_path)

        # Call the function with parameters
        run_trimmomatic(params.trim_tool, (str(int(params.trim_heap[:-1])//2)+'M'), (str(int(params.trim_mem[:-1])//2)+'M'), params.seq_type, threads, input.fastq1, input.fastq2, params.inputfolder, params.outputfolder, wildcards.sample)