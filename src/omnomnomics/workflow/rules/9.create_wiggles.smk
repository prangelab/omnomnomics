# Rule 9 Create Wiggles

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob

rule create_wiggles:
    input:
        input_tar_gz_file = f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.filtered.HOMER_tagDir.tar.gz"
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num']-1]}/{{sample}}.extra_9.tmp"
    params:
        thetype = config['THETYPE'],
        genome = config['THEGENOME'],
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num']-1]}"
    threads:
        Threads_Per_Rule['9']
    resources:
        mem_mb = Memory_Per_Rule['9'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['9']
    run:
        def create_wig(input_tar_gz_file, inputfolder, outputfolder, thetype, genome):
            log_it(logfile, "Creating Wiggles and TrackHubs...", f"EXECUTING STEP {master_config['wig']}")
            log_it(logfile, f"Input folder: {params.inputfolder}")
            log_it(logfile, f"Output folder: {params.outputfolder}")

            sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['wig_rule_num']-1])
            basename = os.path.basename(input_tar_gz_file)
            tag_dir = os.path.join(inputfolder, basename.replace(".tar.gz", ""))
            tag_dir_short = basename.replace(".tar.gz", "")
            
            # unpack the tar.gz file
            log_it(logfile, f"Unpacking {basename} in {inputfolder}")
            shell(f"""
                    mkdir -p {tag_dir} && \
                    cd {inputfolder} && \
                    tar --strip-components=1 -xzf {basename} -C {tag_dir_short}
                    
            """)
            
            if os.path.exists(tag_dir):
                log_it(logfile, f"Unpacked {basename} successfully")
            else:
                log_it(logfile, f"Failed to unpack {basename}")

            if thetype == "RNA": 
                log_it(logfile, f"Creating trackhub from {tag_dir}")
                base_name = os.path.basename(tag_dir).replace(".HOMER_tagDir", '')
                log_it(logfile, f"makeMultiWigHub.pl {master_config['output_folders'][master_config['wig_rule_num']-1]}/{base_name}.hub {genome} -d {master_config['input_folders'][master_config['wig_rule_num']-1]}/{tag_dir_short} -force -fsize 1e8 -strand -webdir . -url \"https:/www.macrophages.eu/UCSCtracks/\"")
                shell(f"makeMultiWigHub.pl {master_config['output_folders'][master_config['wig_rule_num']-1]}/{base_name}.hub {genome} -d {master_config['input_folders'][master_config['wig_rule_num']-1]}/{tag_dir_short} -force -fsize 1e8 -strand -webdir . -url \"https:/www.macrophages.eu/UCSCtracks/\"")
                
            else:
                log_it(logfile, f"Creating bigwig from {tag_dir}")
                basename = os.path.basename(tag_dir).replace(".HOMER_tagDir", '')
                log_it(logfile, f"makeUCSCfile {tag_dir} -fsize 1e8 > {outputfolder}/{basename}.bw")
                shell(f"makeUCSCfile {tag_dir} -fsize 1e8 > {outputfolder}/{basename}.bw")

            # Fix results by converting bedGraph to bigWig if needed
            for hub_file in glob.glob(os.path.join(outputfolder, "*.hub")):
                log_it(logfile, f"Fixing trackhub {hub_file}...")
                log_it(logfile, f"fix_HOMER_trackHub.sh -i {hub_file} -g {genome}")
                path = os.path.join(config['SCRIPT_DIR'], "fix_HOMER_trackHub.sh")
                shell(f"{path} -i {hub_file} -g {genome} -c {config['GENOME_AUX_DIR']}")

            for bw_file in glob.glob(os.path.join(outputfolder, "*.bw")):
                log_it(logfile, f"Fixing bigwig {bw_file}...")
                log_it(logfile, f"fix_HOMER_bigwig.sh -i {bw_file} -g {genome}")
                path = os.path.join(config['SCRIPT_DIR'], "fix_HOMER_bigwig.sh")
                shell(f"{path} -i {bw_file} -g {genome} -c {config['GENOME_AUX_DIR']}")

            # Remove the uncompressed tag directory
            log_it(logfile, f"Removing uncompressed tag directory {tag_dir}")
            shell(f"rm -r {tag_dir}")

            shell(f"""echo "necessity file for create wiggs. can delete this." > {outputfolder}/{wildcards.sample}.extra_9.tmp""")
            
        create_wig(input.input_tar_gz_file, params.inputfolder, params.outputfolder, params.thetype, params.genome)
