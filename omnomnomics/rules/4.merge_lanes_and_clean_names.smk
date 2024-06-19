#Rule 4 merge lanes and clean names

## Omnomnomics Snake Rule  ##
import os
import glob

def input_function(wildcards):
    input_folder = master_config['input_folders'][master_config['merge_rule_num']-1]
    
    #input_files = [f"{input_folder}/{sample}_L0{i:02}.bam" for i in range(1, 100) if os.path.exists(f"{input_folder}/{sample}_L0{i:02}.bam")]

    # Use glob.glob to list all files matching the pattern
    input_files = glob.glob(f"{input_folder}/{wildcards.sample}_L*.bam")
    
    return input_files

rule merge_bam:
    input:
        # bam_files = glob.glob(f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}_L00{{num}}.bam")
        bam_files = input_function,
        extra_files = expand(f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}.extra_3.tmp", sample = samples) if 3 in themode else []
        #extra_files = expand(f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}.extra_3.tmp", sample = samples) if 3 in themode else None
        #extra_file = f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}.HISAT2_stats.txt"
        #extra_files = f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}.HISAT2_stats.txt" if config['THEMAPTOOL'] == "HISAT2" else f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}.STAR_stats.txt"
        #bam_files =  glob.glob(f"{input_folder}/{{sample}}_L00*.bam")
        #bam_files = f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}_L004.bam"
        #bam_files = [f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}_L0{i:02}.bam" for i in range(1, 100) if os.path.exists(f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{{sample}}_L00{i:02}.bam")]
    output:
        f"{master_config['output_folders'][master_config['merge_rule_num']-1]}/{{sample}}.bam",
        f"{master_config['output_folders'][master_config['merge_rule_num']-1]}/{{sample}}.extra_4.tmp"
    params:
        inputfolder = master_config['input_folders'][master_config['merge_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['merge_rule_num']-1]
    threads:
        Threads_Per_Rule['4']
    resources:
        mem_mb = Memory_Per_Rule['4']
    benchmark:
        f"{master_config['output_folders'][master_config['merge_rule_num']-1]}/{{sample}}_mergebam_benchmark.tsv"
    run:
        log_it(logfile, "Merging lanes...", f"EXECUTING STEP {master_config['merge_rule_num']}")
        log_it(logfile, "Input folder: BAM")
        log_it(logfile, "Output folder: BAM")
        # mouse_F_ldlr_PlaqueMacrophages_RNA_2301_DMSO_poolA_RK_I20231020_L004.bam
        # print(input.bam_files)
        # input_files = [f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{wildcards.sample}_L0{i:02}.bam" for i in range(1, 100) if os.path.exists(f"{master_config['input_folders'][master_config['merge_rule_num']-1]}/{wildcards.sample}_L0{i:02}.bam")]

        # for file in input_files:
        #     print(file)

        # input_files = glob.glob(f"{input_folder}/{wildcards.sample}_L*.bam")
        # print(input_files)
            
        merging_version = subprocess.check_output("module load samtools && samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+merging_version.decode("utf-8"), "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['merge_rule_num']-1])

        def merge_bam_files(input_bamfiles, threads):
            # Check if we have lane info
            print(f"INPUT BAM MERGE FILES {input_bamfiles}")
            if input_bamfiles:
                # Initialize a dictionary to store sample names and their associated lane numbers with aligner extensions
                sample_lane_dict = {}

                # Iterate through each BAM file in the list of bam_files
                for bam_file in input_bamfiles:
                    # Extract the base name of the BAM file (i.e., the file name without the directory path)
                    basename = os.path.basename(bam_file)
                    
                    # Extract sample name, lane information, and aligner extension
                    parts = basename.split('_L00')
                    sample = parts[0]
                    lane_and_extension = parts[1].split('.', 1)
                    lane_number = 'L00' + lane_and_extension[0]  # Get the full lane number (e.g., L001, L002)
                    ext = '.' + lane_and_extension[1]  # Get the aligner extension or empty string

                    # Create a unique key for the sample and aligner extension
                    sample_key = f"{sample}{ext}"

                    # If the sample_key is not already in the dictionary, add it with an empty set as its value
                    if sample_key not in sample_lane_dict:
                        sample_lane_dict[sample_key] = set()
                    
                    # Add the lane number to the set of lane numbers for this sample_key
                    sample_lane_dict[sample_key].add(lane_number)

                # Check for samples with multiple lane numbers
                #sample_keys_with_multiple_lanes = [sample_key for sample_key, lanes in sample_lane_dict.items() if len(lanes) > 1]

                # Iterate over the samples with multiple lanes
                for sample_key, lanes in sample_lane_dict.items():
                    # Get number of lanes (L00n)
                    num_lanes = len(lanes)
                    # Merge only if more than one lane was run
                    if num_lanes > 1:
                        # Collect the set of files per sample
                        parts = sample_key.split('.', 1)
                        sample = parts[0]
                        ext = '.'+parts[1]
                        #All files for this samples with different lane numbers but same aligners extension
                        bam_list = glob.glob(f"BAM/{sample}_L00*{ext}")

                        # Set clean names
                        myname = parts[0]   # the sample name without lane info or aligner extension

                        # Run samtools merge
                        log_it(logfile, f"Merging {myname} lanes...")
                        shell(f"""
                            module load samtools && \
                            samtools merge -@ {threads} -o BAM/{myname}{ext} {' '.join(bam_list)}""")
                        shell(f"""echo "necessity file for merge bams. can delete this." > BAM/{wildcards.sample}.extra.tmp""")
                    else:
                        # If there is only one lane we don't need to bother with merging and can simply clean the name
                        parts = sample_key.split('.', 1)
                        from_file = glob.glob(f"BAM/{parts[0]}_L00*.{parts[1]}")[0]
                        to_file = os.path.join('BAM', parts[0] + ('.'+parts[1]))
                        log_it(logfile, f"Renaming {from_file} to {to_file}...")
                        shell(f"mv {from_file} {to_file}")
                        shell(f"""echo "necessity file for merge bams. can delete this." > BAM/{wildcards.sample}.extra_4.tmp""")

                # # Remove old BAM files with lane info
                # for bam_file in glob.glob(f"BAM/*_L00*.bam"):
                #     shell(f"rm {bam_file}")
            else:
                log_it(logfile, "No lane info found!")
                shell(f"""echo "necessity file for merge bams. can delete this." > BAM/{wildcards.sample}.extra_4.tmp""")

        merge_bam_files(input.bam_files, threads)