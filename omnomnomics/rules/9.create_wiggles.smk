# Rule 9 Create Wiggles

## Omnomnomics Snake Rule  ##
import os
import glob

rule create_wiggles:
    input:
        input_tar_gz_file = f"{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz" if config['THETYPE'] != "CHIP" else f"{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.filtered.HOMER_tagDir.tar.gz"
        # input2 = glob.glob(f"{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.HOMER_tagDir/*")
    output:
        #f"{master_config['output_folders'][master_config['wig_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bw" if config['THETYPE'] != "CHIP" else f"{master_config['output_folders'][master_config['wig_rule_num']-1]}/{{sample}}.filtered.bw"
        f"{master_config['output_folders'][master_config['wig_rule_num']-1]}/{{sample}}.extra_9.tmp"
    params:
        thetype = config['THETYPE'],
        genome = config['THEGENOME'],
        inputfolder =master_config['input_folders'][master_config['wig_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['wig_rule_num']-1]
    threads:
        Threads_Per_Rule['9']
    resources:
        mem_mb = Memory_Per_Rule['9']
    benchmark:
        f"{master_config['output_folders'][master_config['wig_rule_num']-1]}/{{sample}}_create_wiggles_benchmark.tsv"
    run:
        def create_wig(input_tar_gz_file, inputfolder, outputfolder, thetype, genome):
            log_it(logfile, "Creating Wiggles and TrackHubs...", f"EXECUTING STEP {master_config['wig']}")
            log_it(logfile, f"Input folder: HOMER_tagDirs")
            log_it(logfile, f"Output folder: BigWigs")

            # path = os.path.join(OMNOM_HOME, "bin", "homer", "configureHomer.pl")
            # version = subprocess.check_output("perl {path} -list 2> /dev/null | grep homer",  shell=True, executable='/bin/bash')
            # log_it(logfile, "\n"+version.decode("utf-8"), "VERSION")

            sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['wig_rule_num']-1])
            basename = os.path.basename(input_tar_gz_file)
            tag_dir = os.path.join(inputfolder, basename.replace(".tar.gz", ""))
            
            # unpack the tar.gz file
            log_it(logfile, f"Unpacking {basename} into {inputfolder}")
            shell(f"""
                cd {inputfolder} && \
                tar --strip-components=1 -xzf {basename}
            """)
            
            if os.path.exists(tag_dir):
                print("TARGZ UNPACKED SUCCESSFULLY")
                log_it(logfile, f"Unpacked {basename} successfully")
            else:
                print("DIDN'T UNPACK SUCCESSFULLY")
                log_it(logfile, f"Failed to unpack {basename}")
            print(f"TAGDIR IS {tag_dir}")

            if thetype == "RNA":
                log_it(logfile, f"Creating trackhub from {tag_dir}")
                base_name = os.path.basename(tag_dir).replace(".HOMER_tagDir", '')
                log_it(logfile, f"makeMultiWigHub.pl {outputfolder}/{base_name}.hub {genome} -d {tag_dir} -fsize 1e8 -strand -webdir . -url \"https:/www.macrophages.eu/UCSCtracks/\"")
                shell(f"makeMultiWigHub.pl {outputfolder}/{base_name}.hub {genome} -d {tag_dir} -fsize 1e8 -strand -webdir . -url \"https:/www.macrophages.eu/UCSCtracks/\"")
            else:
                log_it(logfile, f"Creating bigwig from {tag_dir}")
                basename = os.path.basename(tag_dir).replace(".HOMER_tagDir", '')
                log_it(logfile, f"makeUCSCfile {tag_dir} -fsize 1e8 > {outputfolder}/{basename}.bw")
                shell(f"makeUCSCfile {tag_dir} -fsize 1e8 > {outputfolder}/{basename}.bw")

            # Fix results by converting bedGraph to bigWig if needed
            for hub_file in glob.glob(os.path.join(outputfolder, "*.hub")):
                log_it(logfile, f"Fixing trackhub {hub_file}...")
                log_it(logfile, f"fix_HOMER_trackHub.sh -i {hub_file} -g {genome}")
                path = os.path.join(OMNOM_HOME, "bin", "scripts", "fix_HOMER_trackHub.sh") 
                shell(f"{path} -i {hub_file} -g {genome}")

            for bw_file in glob.glob(os.path.join(outputfolder, "*.bw")):
                log_it(logfile, f"Fixing bigwig {bw_file}...")
                log_it(logfile, f"fix_HOMER_bigwig.sh -i {bw_file} -g {genome}")
                path = os.path.join(OMNOM_HOME, "bin", "scripts", "fix_HOMER_bigwig.sh") 
                shell(f"{path} -i {bw_file} -g {genome}")

            shell(f"""echo "necessity file for create wiggs. can delete this." > {outputfolder}/{wildcards.sample}.extra_9.tmp""")
            
        create_wig(input.input_tar_gz_file, params.inputfolder, params.outputfolder, params.thetype, params.genome)


    ##can do it like this if have problem with the nullglob
    # def create_wig(input_tar_gz_file, inputfolder, outputfolder, thetype, genome):
    #            log_it(logfile, "Creating Wiggles and TrackHubs...", f"EXECUTING STEP {master_config['wig']}")
    #            log_it(logfile, f"Input folder: HOMER_tagDirs")
    #            log_it(logfile, f"Output folder: BigWigs")

    #            # Extract the tar.gz file
    #            os.system(f"tar --strip-components=1 -xzf {input_tar_gz_file}")

    #            # Using a shell script to handle nullglob and subsequent commands
    #            shell_script = f"""
    #            #!/bin/bash
    #            shopt -s nullglob
    #            inputfolder="{inputfolder}"
    #            outputfolder="{outputfolder}"
    #            thetype="{thetype}"
    #            genome="{genome}"

    #            for i in "$inputfolder"/*.HOMER_tagDir/*; do
    #                if [[ "$thetype" == "RNA" ]]; then
    #                    echo "Creating trackhub from $i"
    #                    makeMultiWigHub.pl "$outputfolder"/$(basename "$i").hub "$genome" -d "$i" -fsize 1e8 -strand -webdir . -url "https:/www.macrophages.eu/UCSCtracks/"
    #                else
    #                    echo "Creating bigwig from $i"
    #                    makeUCSCfile "$i" -fsize 1e8 > "$outputfolder"/$(basename "$i").bw
    #                fi
    #            done

    #            for hub_file in "$outputfolder"/*.hub; do
    #                echo "Fixing trackhub $hub_file..."
    #                fix_HOMER_trackHub.sh -i "$hub_file" -g "$genome"
    #            done

    #            for bw_file in "$outputfolder"/*.bw; do
    #                echo "Fixing bigwig $bw_file..."
    #                fix_HOMER_bigwig.sh -i "$bw_file" -g "$genome"
    #            done

    #            shopt -u nullglob
    #            """

    #            # Write the script to a temporary file and execute it
    #            script_path = "temp_script.sh"
    #            with open(script_path, "w") as script_file:
    #                script_file.write(shell_script)
            
    #            os.system(f"bash {script_path}")
    #            os.remove(script_path)
        
    #        create_wig(input.input_tar_gz_file, params.inputfolder, params.outputfolder, params.thetype, params.genome)