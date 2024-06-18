# Rule 1 option skewer

## Omnomnomics Snake Rule  ##
import os
import glob
import subprocess

rule run_skewer:
    input:
        fastq1=f"{master_config['input_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R1.fastq.gz",
        fastq2=f"{master_config['input_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R2.fastq.gz" if config["PAIRED"] else None
    output:
        trimmed_fastq1=f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2=f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_R2.trimmed.fastq.gz" if config["PAIRED"] else None
    params:
        seq_type=config["THETYPE"],
        inputfolder = master_config['input_folders'][master_config['trim_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['trim_rule_num']-1]
    threads: 
        Threads_Per_Rule['1']
        # max((master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1]
        #  if 'mincores_single_sample_step1_9' in master_config
        #  and isinstance(master_config['mincores_single_sample_step1_9'], list)
        #  and len(master_config['mincores_single_sample_step1_9']) > (master_config['trim_rule_num'] - 1)
        #  and master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] is not None
        #  and isinstance(master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1], int)
        #  and master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] > master_config['min_slice_cores']
        #  else master_config['min_slice_cores']),
        #  min((master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num'] - 1]
        #  if 'maxcores_single_sample_step1_9' in master_config
        #  and isinstance(master_config['maxcores_single_sample_step1_9'], list)
        #  and len(master_config['maxcores_single_sample_step1_9']) > (master_config['trim_rule_num'] - 1)
        #  and master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] is not None
        #  and isinstance(master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num'] - 1], int)
        #  and master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] < master_config['cores_per_node']
        #  else master_config['cores_per_node']), ((master_config['nodes_in_partition']*master_config['cores_per_node'])/num_samples)) )

        #max(master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num']-1], min(master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num']-1],((master_config['nodes_in_partition']*master_config['cores_per_node'])/num_samples)))
    resources:
        # mem_mb = (master_config['min_mem_mb'][master_config['trim_rule_num'] - 1]
        #     if 'min_mem_mb' in master_config
        #     and isinstance(master_config['min_mem_mb'], list)
        #     and len(master_config['min_mem_mb']) > (master_config['trim_rule_num'] - 1)
        #     and master_config['min_mem_mb'][master_config['trim_rule_num'] - 1] is not None
        #     and isinstance(master_config['min_mem_mb'][master_config['trim_rule_num'] - 1], int)
        #     and master_config['min_mem_mb'][master_config['trim_rule_num'] - 1] > master_config['min_slice_mem']
        #     else (master_config['min_slice_mem'] 
        #         if('min_mem_mb' in master_config
        #             and isinstance(master_config['min_mem_mb'], list)
        #             and len(master_config['min_mem_mb']) > (master_config['trim_rule_num'] - 1)
        #             and master_config['min_mem_mb'][master_config['trim_rule_num'] - 1] is not None
        #             and isinstance(master_config['min_mem_mb'][master_config['trim_rule_num'] - 1], int)
        #             and master_config['min_mem_mb'][master_config['trim_rule_num'] - 1] <= master_config['min_slice_mem'])
        #         else (master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1]
        #         if 'mincores_single_sample_step1_9' in master_config
        #         and isinstance(master_config['mincores_single_sample_step1_9'], list)
        #         and len(master_config['mincores_single_sample_step1_9']) > (master_config['trim_rule_num'] - 1)
        #         and master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] is not None
        #         and isinstance(master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1], int)
        #         and master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] > master_config['min_slice_cores']
        #         else master_config['min_slice_cores'])* master_config['max_mem_per_core_mb']))
        mem_mb = Memory_Per_Rule['1']
    benchmark:
        f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}_skewer_benchmark.tsv"
    run:
        def run_skewer(logfile, trim_tool, seq_type, threads, fastq1, fastq2, inputfolder, outputfolder, sample):
            log_it(logfile, "Trimming reads...", f"EXECUTING STEP {master_config['trim_rule_num']}")
            log_it(logfile, f"Input folder: {inputfolder}")
            log_it(logfile, f"Output folder: {outputfolder}")
            log_it(logfile, f"Trim Tool: {trim_tool}")

            skewer_version = subprocess.check_output(["skewer", "--version"])
            log_it(logfile, "\n"+skewer_version.decode("utf-8"), "SKEWER VERSION")
            print(skewer_version.decode("utf-8"))


            sanity_check_dir(logfile, inputfolder,  master_config['input_file_types'][master_config['trim_rule_num']-1])

            if seq_type == "ATAC":
                adapter_option = "-x CTGTCTCTTATACACATCT -y AGATGTGTATAAGAGACAG" if config["PAIRED"] else "-x CTGTCTCTTATACACATCT"
            else:
                adapter_option = ""

            if fastq2:
                log_it(logfile, "Running skewer in Paired End mode.")
                skewer_command = f"""
                    skewer --quiet {adapter_option} -m pe -q 15 -Q 15 -z -t {threads} -o "{outputfolder}/{sample}" {fastq1} {fastq2}
                """
            else:
                log_it(logfile, "Running skewer in Single End mode.")
                skewer_command = f"""
                    skewer --quiet {adapter_option} -m any -q 15 -Q 15 -z -t {threads} -o "{outputfolder}/{sample}" {fastq1}
                """

            # Run the skewer command
            shell(skewer_command, bench_record=bench_record)

            # Rename R1 trimmed files
            log_it(logfile, "Renaming trimmed results...")
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*pair1.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('-trimmed-pair1.fastq.gz', '_R1.trimmed.fastq.gz')) ######check what it outputs with unpaired
                os.rename(file_path, new_name)

            # Rename R2 trimmed files
            for file_path in glob.glob(os.path.join(f"{outputfolder}", '*pair2.fastq.gz')):
                base_name = os.path.basename(file_path)
                new_name = os.path.join(f"{outputfolder}", base_name.replace('-trimmed-pair2.fastq.gz', '_R2.trimmed.fastq.gz'))
                os.rename(file_path, new_name)
            os.remove((os.path.join(f"{outputfolder}", f"{sample}" + "-trimmed.log" )))

            # # Rename trimmed files for Skewer
            # for filename in os.listdir('trimmed_FASTQ'):
            #     file_path = os.path.join(directory, filename)
            #     if os.path.isfile(file_path):
            #         base, ext = os.path.splitext(filename)
            #         new_filename = f"{base}_Skewer{ext}"
            #         new_file_path = os.path.join(directory, new_filename)
            #         os.rename(file_path, new_file_path)

        # Call the function with parameters
        run_skewer(logfile, config["THETRIMTOOL"], params.seq_type, threads, input.fastq1, input.fastq2, params.inputfolder, params.outputfolder, wildcards.sample)