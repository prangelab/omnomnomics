#Rule 4 merge lanes and clean names

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob
import shlex
import shutil
import subprocess
import tempfile

rule merge_bam:
    wildcard_constraints:
        sample4=merged_sample_wildcard_pattern
    input:
        extra_files=lambda wildcards: [
            f"{experiment_dir}/{master_config['input_folders'][master_config['merge_rule_num']-1]}/{lane_sample}.extra_3.tmp"
            for lane_sample in lane_samples_for_merged_sample(wildcards.sample4)
        ] if 3 in themode else []
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}/{{sample4}}.bam",
        f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}/{{sample4}}.extra_4.tmp"
    params:
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['merge_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}"
    threads:
        Threads_Per_Rule['4']
    resources:
        mem_mb = Memory_Per_Rule['4'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['4']
    run:
        log_once(logfile, "step4.header", "Merging lanes...", f"EXECUTING STEP {master_config['merge_rule_num']}")
        log_once(logfile, "step4.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step4.outputfolder", f"Output folder: {params.outputfolder}")
            
        merging_version = subprocess.check_output("samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_once(logfile, "step4.samtools_version", "\n"+merging_version.decode("utf-8"), "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['merge_rule_num']-1])

        def merge_bam_files(threads, inputfolder, sample4):
            tracking = begin_step_sample(master_config['merge_rule_num'], sample4, "merge_bam")
            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample4}.", dir=tmpdir_root if tmpdir_root else None)
            record_step_note(master_config['merge_rule_num'], sample4, f"scratch_dir={local_workdir}")

            def quote(path):
                return shlex.quote(path)

            # Check if we have lane info
            input_bamfiles = glob.glob(f"{inputfolder}/{sample4}_L*.bam")
            try:
                if input_bamfiles:
                    # Initialize a dictionary to store sample4 names and their associated lane numbers with aligner extensions
                    sample_lane_dict = {}

                    # Iterate through each BAM file in the list of bam_files
                    for bam_file in input_bamfiles:
                        # Extract the base name of the BAM file (i.e., the file name without the directory path)
                        basename = os.path.basename(bam_file)
                        
                        # Extract sample4 name, lane information, and aligner extension
                        parts = basename.split('_L00')
                        sample2 = parts[0]
                        lane_and_extension = parts[1].split('.', 1)
                        lane_number = 'L00' + lane_and_extension[0]  # Get the full lane number (e.g., L001, L002)
                        ext = '.' + lane_and_extension[1]  # Get the aligner extension or empty string

                        # Create a unique key for the sample4 and aligner extension
                        sample_key = f"{sample2}{ext}"

                        # If the sample_key is not already in the dictionary, add it with an empty set as its value
                        if sample_key not in sample_lane_dict:
                            sample_lane_dict[sample_key] = set()
                        
                        # Add the lane number to the set of lane numbers for this sample_key
                        sample_lane_dict[sample_key].add(lane_number)

                    # Iterate over the samples with multiple lanes
                    for sample_key, lanes in sample_lane_dict.items():
                        # Get number of lanes (L00n)
                        num_lanes = len(lanes)
                        # Merge only if more than one lane was run
                        if num_lanes > 1:
                            # Collect the set of files per sample
                            parts = sample_key.split('.', 1)
                            sample3 = parts[0]
                            ext = '.'+parts[1]
                            #All files for this samples with different lane numbers but same aligners extension
                            bam_list = glob.glob(f"{inputfolder}/{sample3}_L0*{ext}")

                            # Set clean names
                            myname = parts[0]   # the sample name without lane info or aligner extension

                            # Run samtools merge
                            record_step_note(master_config['merge_rule_num'], myname, "merging_lanes")
                            local_bam_list = []
                            for bamfile in bam_list:
                                local_bam = os.path.join(local_workdir, os.path.basename(bamfile))
                                shell(f'cp {quote(bamfile)} {quote(local_bam)}')
                                local_bam_list.append(local_bam)
                            local_output = os.path.join(local_workdir, f"{myname}{ext}")
                            merge_command = f"samtools merge -@ {threads} -o {quote(local_output)} {' '.join(quote(bam) for bam in local_bam_list)}"
                            record_step_command(master_config['merge_rule_num'], myname, merge_command)
                            shell(merge_command)
                            shell(f'cp {quote(local_output)} {quote(os.path.join(inputfolder, f"{myname}{ext}"))}')
                            shell(f"""echo "necessity file for merge bams. can delete this." > {inputfolder}/{wildcards.sample4}.extra_4.tmp""")
                        else:
                            # If there is only one lane we don't need to bother with merging and can simply clean the name
                            parts = sample_key.split('.', 1)
                            from_file = glob.glob(f"{inputfolder}/{parts[0]}_L00*.{parts[1]}")[0]
                            to_file = os.path.join(inputfolder, parts[0] + ('.'+parts[1]))
                            record_step_note(master_config['merge_rule_num'], parts[0], "renaming_single_lane_bam")
                            shell(f"cp {quote(from_file)} {quote(to_file)}")
                            shell(f"""echo "necessity file for merge bams. can delete this." > {inputfolder}/{wildcards.sample4}.extra_4.tmp""")

                else:
                    record_step_note(master_config['merge_rule_num'], sample4, "no_lane_info_found")
                    shell(f"""echo "necessity file for merge bams. can delete this." > {inputfolder}/{wildcards.sample4}.extra_4.tmp""")
                finish_step_sample(master_config['merge_rule_num'], sample4, "merge_bam", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['merge_rule_num'], sample4, "merge_bam", tracking["start_time"], "FAIL")
                raise
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

        merge_bam_files(threads, params.inputfolder, wildcards.sample4)
