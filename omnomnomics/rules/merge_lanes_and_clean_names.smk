#Rule 4 merge lanes and clean names

## Omnomnomics Snake Rule  ##
import os
import glob

rule merge_bam:
    input:
        bam_files = glob.glob(f"BAM/{{sample}}_L00*.bam")
    output:
        bam="merged_BAM/{sample}.bam"
    params:
        infolder="trimmed_FASTQ",
        outfolder="merged_BAM",
        maxcores=config["MAXCORES"],
        step_merge=4,
        inputfiletype=".bam"
    threads:
        Threads_Per_Rule['4']
    resources:
        mem_mb = Memory_Per_Rule['4']
    run:
        log_it(logfile, "Merging lanes...", f"EXECUTING STEP {master_config['merge_rule_num']}")
        log_it(logfile, "Input folder: BAM")
        log_it(logfile, "Output folder: BAM")

        def merge_bam_files(input_bamfiles, threads):
            # Check if we have lane info
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
                    aligner_extension = '.' + lane_and_extension[1]  # Get the aligner extension or empty string

                    # Create a unique key for the sample and aligner extension
                    sample_key = f"{sample}{aligner_extension}"

                    # If the sample_key is not already in the dictionary, add it with an empty set as its value
                    if sample_key not in sample_lane_dict:
                        sample_lane_dict[sample_key] = set()
                    
                    # Add the lane number to the set of lane numbers for this sample_key
                    sample_lane_dict[sample_key].add(lane_number)

                # Check for samples with multiple lane numbers
                #sample_keys_with_multiple_lanes = [sample_key for sample_key, lanes in sample_lane_dict.items() if len(lanes) > 1]

                # Iterate over the samples with multiple lanes
                for sample_key in sample_lane_dict.items():
                    # Get number of lanes (L00n)
                    num_lanes = len(sample_lane_dict[sample_key])
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
                        shell(f"samtools merge -@ {threads} -o BAM/{myname}{ext} {' '.join(bam_list)}")
                    else:
                        # If there is only one lane we don't need to bother with merging and can simply clean the name
                        parts = sample_key.split('.', 1)
                        from_file = glob.glob(f"BAM/{parts[0]}_L00*.{parts[1]}")
                        to_file = os.path.join('BAM', parts[0], ('.'+parts[1]))
                        log_it(logfile, f"Renaming {from_file} to {to_file}...")
                        shell(f"mv {from_file} {to_file}")

                # Remove old BAM files with lane info
                for bam_file in glob.glob(f"BAM/*_L00*.bam"):
                    shell(f"rm {bam_file}")
            else:
                log_it(logfile, "No lane info found!")

        merge_bam_files(input.bam_files, threads)