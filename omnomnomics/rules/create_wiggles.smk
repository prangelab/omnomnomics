# Rule 9 Create Wiggles

## Omnomnomics Snake Rule  ##
import os
import glob


rule create_wiggles:
    input:
        input_tar_gz_file=glob.glob(f"{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.HOMER_tagDir.tar.gz")
        # input2 = glob.glob(f"{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.HOMER_tagDir/*")
    output:
        f"{master_config['output_folders'][master_config['wig_rule_num']-1]}/{{sample}}.bw"
    params:
        thetype = config['THETYPE']
        genome = config['THEGENOME']
        inputfolder =master_config['input_folders'][master_config['wig_rule_num']-1]
        outputfolder = master_config['output_folders'][master_config['wig_rule_num']-1]
    threads:
        Threads_Per_Rule['9']
    resources:
        mem_mb = Memory_Per_Rule['9']
    run:
        def create_wig(input_tar_gz_file, inputfolder, outputfolder, thetype, genome):
            log_it(logfile, "Creating Wiggles and TrackHubs...", f"EXECUTING STEP {master_config['wig']}")
            log_it(logfile, f"Input folder: HOMER_tagDirs")
            log_it(logfile, f"Output folder: BigWigs")

            version = subprocess.check_output(["perl", os.path.join(OMNOM_HOME, "bin", "homer", "configureHomer.pl"), "-list", "2>", "/dev/null", "|", "grep", "homer"])
            log_it(logfile, "\n"+version.decode("utf-8"), "VERSION")
            print(version.decode("utf-8"))
            sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['wig_rule_num']-1])

            os.system(f"tar --strip-components=1 -xzf {input_tar_gz_file}")

            tag_dir = f"{inputfolder}/{{sample}}.HOMER_tagDir"
            if thetype == "RNA":
                log_it(logfile, f"Creating trackhub from {tag_dir}")
                #log_it(f"makeMultiWigHub.pl {outfolder}/{os.path.basename(i)}.hub {genome} -d {i} -fsize 1e8 -strand -webdir . -url \"https:/www.macrophages.eu/UCSCtracks/\"")
                shell(f"makeMultiWigHub.pl {outputfolder}/{os.path.basename(tag_dir).replace(".HOMER_tagDir", '')}.hub {genome} -d {tag_dir} -fsize 1e8 -strand -webdir . -url \"https:/www.macrophages.eu/UCSCtracks/\"")
            else:
                log_it(logfile, f"Creating bigwig from tag_diri}")
                #log_it(f"makeUCSCfile {i} -fsize 1e8 > {outfolder}/{os.path.basename(i)}.bw")
                shell(f"makeUCSCfile {tag_dir} -fsize 1e8 > {outputfolder}/{os.path.basename(tag_dir).replace(".HOMER_tagDir", '')}.bw")

            # Set nullglob to avoid running on the literal wildcard string if no files exist
            os.system("shopt -s nullglob")

            # Fix results by converting bedGraph to bigWig if needed
            for hub_file in glob.glob(os.path.join(outputfolder, "*.hub")):
                log_it(logfile, f"Fixing trackhub {hub_file}...")
                #log_it(f"fix_HOMER_trackHub.sh -i {hub_file} -g {genome}")
                shell(f"fix_HOMER_trackHub.sh -i {hub_file} -g {genome}")


            for bw_file in glob.glob(os.path.join(THEOUTFOLDER, "*.bw")):
                log_it(logfile, f"Fixing bigwig {bw_file}...")
                #log_it(f"fix_HOMER_bigwig.sh -i {bw_file} -g {genome}")
                shell(f"fix_HOMER_bigwig.sh -i {bw_file} -g {genome}")

            # Unset nullglob for safety
            os.system("shopt -u nullglob")
        
        create_wig(input.input_tar_gz_file, params.inputfolder, params.outfolder, params.thetype, params.genome)


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