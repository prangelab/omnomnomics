# Rule 11: Call Peaks

## Omnomnomics Snake Rule ##
import os
import re
import glob
import subprocess


rule call_peaks:
    input:
        bam_files = glob.glob(f"{master_config['input_folders'][master_config['callpeaks_rule_num']-1]}/*.bam")
    output:
        all_merged_peaks = f"{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/all_groups.merged_peaks.bed"
    params:
        thetype=config['THETYPE'],  
        broad=config['BROAD'],
        input_sample=config['INPUT'],
        homer_input=config['HOMERINPUT'],
        genome=config['THEGENOME'],  
        experiment_dir=config['EXPERIMENT_DIR'], 
        name_fields=config['NAMEFIELDS'],  
        separator=config['THESEPARATOR'],
        thetype_field = ... ,
        style=config['THESTYLE'],
        size=config['HOMERSIZE'],
        mindist=config['HOMERMINDIST'],
        inputfolder=master_config['input_folders'][master_config['callpeaks_rule_num']-1],
        outputfolder=master_config['output_folders'][master_config['callpeaks_rule_num']-1]
    threads:
        Threads_Per_Rule['11']
    resources:
        mem_mb = Memory_Per_Rule['11']
    run: 
        log_it(logfile, "Calling Peaks...", f"EXECUTING STEP {master_config['callpeaks_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")
        def call_peaks(logfile, thetype):
            if thetype == "RNA":
                log_it("Not a ChIP- or ATAC-seq experiment, skipping this step...")
                #might have to make the output files, or condition in output on if type = RNA, and then already not execute this step. find a handy way to account for this.
                return
            
            log_it(logfile, "Calling peaks...", f"EXECUTING STEP {step_callpeaks}")
            log_it(logfile, f"The type = {thetype}")

            # Report version
            log_it(logfile, subprocess.getoutput("macs3 --version"))
            log_it(logfile, subprocess.getoutput(f"perl {os.getenv('OMNOM_HOME')}/bin/homer/configureHomer.pl -list 2> /dev/null | grep homer"))
            print(version.decode("utf-8"))

            version = subprocess.check_output(["perl", os.path.join(OMNOM_HOME, "bin", "homer", "configureHomer.pl"), "-list", "2>", "/dev/null", "|", "grep", "homer"])
            log_it(logfile, "\n"+version.decode("utf-8"), "VERSION")
            print(version.decode("utf-8"))

            if thetype == "CHIP":
                #If ChIP, call peaks
                log_it(logfile, "Finding ChIP enriched regions...")

                log_it(logfile, "Calling peaks with MACS3...")
                
                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder,  master_config['input_file_types'][master_config['callpeaks_rule_num']-1]) ###############might have to adjust this based on what I do with input and output folders
                
                # Distribute samples over groups based on the file name pattern
                bams = [os.path.basename(bam).split(separator)[thetype_field - 1] for bam in glob.glob(f"{inputfolder}/*.bam")]
                chip_groups = sorted(set(bams))
                
                for group in chip_groups:
                    # BAMS=( $(ls "$THEINFOLDER" | grep $(echo $THEGROUP | sed "s|$THESEPARATOR|.*|g") - | grep ".bam$") )
                    bams = [bam for bam in os.listdir(inputfolder) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam)] 
                    bams = [os.path.join(inputfolder, bam) for bam in bams]

                    # TRACKS=($(ls "$THEINFOLDER" | grep $(echo $SUPERHUB | sed "s|$THESEPARATOR|.*|g") - | grep ".bw$"))
                    # tracks = [track for track in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", track)]
                    
                    log_it(logfile, f"Calling peaks for: {group}...")
                    log_it(logfile, f"Files in group: {', '.join(bams)}")

                    the_name = os.path.basename(bams[0]).split(separator)[name_fields - 1]

                    if input_sample == "NA":#No input sample provided
                        if broad == "1": # Check if we should call broad peaks
                            log_it(logfile, "Calling peaks with MACS3, using p and q values 1e-6 and 1e-9...")
                            for ext in ["q-9.broad", "q-6.broad", "p-9.broad", "p-6.broad"]:
                                shell(f"""  micromamba activate macs3 && \
                                macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -{'q' if ext.startswith('q') else 'p'} {'1e-9' if ext.endswith('9.broad') else '1e-6'} --broad --verbose 0 """)
                        else:
                            log_it(logfile, "Calling peaks with MACS3, using p and q values 1e-6 and 1e-9...")
                            for ext in ["q-9", "q-6", "p-9", "p-6"]:
                                shell(f""" micromamba activate macs3 && \
                                macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -{'q' if ext.startswith('q') else 'p'} {'1e-9' if ext.endswith('9') else '1e-6'} --verbose 0""")
                    else:# We have input!
                        if broad == "1": # Check if we should call broad peaks
                            log_it(logfile, "Calling peaks with MACS3, using p and q values 1e-6 and 1e-9...")
                            for ext in ["q-9.broad", "q-6.broad", "p-9.broad", "p-6.broad"]:
                                shell(f"""micromamba activate macs3 && \
                                macs3 callpeak -t {' '.join(bams)} -c {input_sample} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -{'q' if ext.startswith('q') else 'p'} {'1e-9' if ext.endswith('9.broad') else '1e-6'} --broad --verbose 0""") #################verbose 2?
                        else:
                            log_it(logfile, "Calling peaks with MACS3, using p and q values 1e-6 and 1e-9...")
                            for ext in ["q-9", "q-6", "p-9", "p-6"]:
                                shell(f""" micromamba activate macs3 && \
                                macs3 callpeak -t {' '.join(bams)} -c {input_sample} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -{'q' if ext.startswith('q') else 'p'} {'1e-9' if ext.endswith('9') else '1e-6'} --verbose 0""")
                # Clean up the output
                for group in chip_groups: 
                    the_name = os.path.basename(
                        next(bam for bam in os.listdir(inputfolder) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam))
                    ).split(separator)[name_fields - 1] ##### why -1? # Set group name
                    log_it(logfile, the_name)
                    for file in glob.glob(f"{outputfolder}/{the_name}.MACS3*"):
                        if not re.search(r"(broadPeak|narrowPeak)$", file): ############correct translation?
                            os.remove(file) # Delete all redundant MACS3 output
                    #can also do: shell(f"find {outputfolder} -type f -name {$the_name}.MACS3* ! \( -name *broadPeak -o -name *narrowPeak \) -delete &")

                log_it(logfile, "Cleaning up MACS3 output (keeping only narrowPeak and broadPeak files)...")

                log_it(logfile, "Converting to a clean 3 column BED format...")
                for file in glob.glob(f"{outputfolder}/*{'narrowPeak' if broad != '1' else 'broadPeak'}"):
                    sorted_bed = f"{outputfolder}/{os.path.basename(file).replace('narrowPeak', '').replace('broadPeak', '')}.bed"
                    shell(f"cut -f1-3 {file} | sort -k1,1 -k2,2n -k3,3n > {sorted_bed}")
                    os.remove(file)

                # Call peaks using HOMER
                inputfolder = "HOMER_tagDirs"
                log_it(logfile, "Calling peaks with HOMER...")
                log_it(logfile, f"Input folder: {inputfolder}")
                
                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder,  master_config['input_file_types'][master_config['callpeaks_rule_num']-1]) ###############might have to adjust this based on what I do with input and output folders
                
                # See if we need to unpack the tag dirs
                if glob.glob(f"{inputfolder}/*tagDir.tar.gz"):
                    for tagdir in glob.glob(f"{inputfolder}/*tagDir.tar.gz"):
                        tagdir_basename = os.path.basename(tagdir)
                        shell(f"cd {inputfolder} && tar --strip-components=1 -xzf {tagdir_basename}")
                # Fetch tagdirs
                tagdirs = glob.glob(f"{inputfolder}/*.HOMER_tagDir")
                
                if homer_input == "NA": # No input sample was provided
                    if broad == "1": # Check if we should call broad peaks
                        log_it(logfile, "Calling broad peaks with HOMER...") 
                        for tagdir in tagdirs:
                            log_it(logfile, "findPeaks {tagdir} -style {thestyle} -size {homer_size} -minDist {homer_mindist} -region -o auto")
                            shell(f"findPeaks {tagdir} -style {thestyle} -size {homer_size} -minDist {homer_mindist} -region -o auto")
                    else:
                        log_it(logfile, "Calling peaks with HOMER...")
                        for tagdir in tagdirs:
                            log_it(logfile, "findPeaks {tagdir} -style {thestyle} -o auto")
                            shell(f"findPeaks {tagdir} -style {thestyle} -o auto")
                else:  # We have input!
                    if broad == "1": # Check if we should call broad peaks
                        log_it(logfile, "Calling broad peaks with HOMER...")
                        for tagdir in tagdirs:
                            log_it(logfile, "findPeaks {tagdir} -i {os.getenv('HOMERINPUT')} -style {thestyle} -size {homer_size} -minDist {homer_mindist} -region -o auto")
                            shell(f"findPeaks {tagdir} -i {os.getenv('HOMERINPUT')} -style {thestyle} -size {homer_size} -minDist {homer_mindist} -region -o auto")
                    else:
                        log_it(logfile, "Calling peaks with HOMER...")
                        for tagdir in tagdirs:
                            log_it(logfile, "findPeaks {tagdir} -i {os.getenv('HOMERINPUT')} -style {thestyle} -o auto")
                            shell(f"findPeaks {tagdir} -i {os.getenv('HOMERINPUT')} -style {thestyle} -o auto")
                
                # Convert peaks to BED format and clean up unwanted contigs (chrUn | alt | random)
                log_it(logfile, "Converting peaks to BED format and cleaning up unwanted contigs (chrUn | alt | random)...")
                for tagdir in tagdirs:
                    log_it(logfile, f"Testpp {tagdir}")
                    if thestyle == "factor":
                        shell(f"pos2bed.pl {tagdir}/peaks.txt | grep -v '#\\|alt\\|Un\\|random' | sort -k1,1 -k2,2n -k3,3n | cut -f-3 > {tagdir}/peaks.bed &")
                    elif thestyle == "histone":
                        shell(f"pos2bed.pl {tagdir}/regions.txt | grep -v '#\\|alt\\|Un\\|random' | sort -k1,1 -k2,2n -k3,3n | cut -f-3 > {tagdir}/regions.bed &")

                ## Distribute samples over groups based on the file name pattern
                # Build group list
                # for tagdir in tagdirs:
                #     shell(f"basename {tagdir} | cut -f {thetype_field} -d {separator} >> groups.tmp.txt")
                # chip_groups = subprocess.getoutput("sort groups.tmp.txt | uniq").split() ######## different way to do this is way I did it at the top and in merge wiggles doesnt work. 

                tagdirs2 = [os.path.basename(tagdir).split(separator)[thetype_field - 1] for tagdir in tagdirs]
                chip_groups = sorted(set(tagdirs2))

                # Iterate over the groups to merge the BED files
                for group in chip_groups:
                    log_it(logfile, f"Merging peaks for: {group}...")
                    log_it(logfile, f"Samples in group: {', '.join(beds)}")
                    # Fetch samples
                    beds = [os.path.join(inputfolder, bed) for bed in os.listdir(inputfolder) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.HOMER_tagDir$", bed)]
                    if thestyle == "factor":
                        beds = [f"{bed}/peaks.bed" for bed in beds]
                    elif thestyle == "histone":
                        beds = [f"{bed}/regions.bed" for bed in beds]
                    
                    the_name = os.path.basename(
                        next(bed for bed in os.listdir(inputfolder) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.HOMER_tagDir$", bed))
                    ).split(separator)[name_fields - 1]

                    sorted_bed = f"{outputfolder}/{the_name}.HOMER.merged_peaks.bed"
                    shell(f"cat {' '.join(beds)} | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {sorted_bed}")
            else:
                # If ATAC, call open regions
                log_it("Calling open regions...", f"EXECUTING STEP {step_callpeaks}")
                
                inputfolder = sys.argv[1]
                outputfolder = sys.argv[2]

                log_it(f"Input folder: {inputfolder}")
                log_it(f"Output folder: {outputfolder}")

                os.makedirs(outputfolder, exist_ok=True)
                
                # Sanity check the working dir
                sanity_check_dir(inputfolder, ".filtered.bam", step_callpeaks)
                
                # Distribute samples over groups based on the file name pattern
                bams = [os.path.basename(bam).split(separator)[thetype_field - 1] for bam in glob.glob(f"{inputfolder}/*.bam")]
                atac_groups = sorted(set(bams))

                for group in atac_groups:
                    bams = [bam for bam in os.listdir(inputfolder) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam)]
                    bams = [os.path.join(inputfolder, bam) for bam in bams]
                    
                    log_it(f"Calling open regions for: {group}...")
                    log_it(f"Files in group: {', '.join(bams)}")

                    the_name = os.path.basename(
                        next(bam for bam in os.listdir(inputfolder) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam))
                    ).split(separator)[name_fields - 1]

                    log_it("Calling peaks with MACS3, q value 1e-9...")
                    shell(f"macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.q-9 -q 1e-9 --verbose 0 &")
                
                log_it("Waiting for MACS to finish calling open regions...")
                shell("wait")

                log_it("Cleaning up MACS3 output (keeping only narrowPeak files)...")
                for group in atac_groups:
                    the_name = os.path.basename(
                        next(bam for bam in os.listdir(inputfolder) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam))
                    ).split(separator)[name_fields - 1]
                    log_it(the_name)
                    shell(f"find {outputfolder} -type f -name {the_name}.MACS3* ! -name *narrowPeak -delete &")

                shell("wait")
                
                log_it("Converting to a clean 3 column BED format...")
                for file in glob.glob(f"{outputfolder}/*narrowPeak"):
                    sorted_bed = f"{outputfolder}/{os.path.basename(file).replace('narrowPeak', '')}.bed"
                    shell(f"cut -f1-3 {file} | sort -k1,1 -k2,2n -k3,3n > {sorted_bed}")
                    os.remove(file)

                # Concatenate the peak files
                shell(f"cat {outputfolder}/*.bed | sort -k1,1 -k2,2n -k3,3n | mergeBed -i - > {outputfolder}/all_groups.merged_peaks.bed")
                log_it(f"Total merged peaks: {subprocess.getoutput(f'wc -l {outputfolder}/all_groups.merged_peaks.bed').strip().split()[0]}")

            # Clean up
            os.remove("groups.tmp.txt")
            shell("micromamba deactivate")


        call_peaks(
            logfile,
            thetype="CHIP",
            inputfolder="input_folder",
            outputfolder="output_folder",
            separator="_",
            thetype_field=,
            name_fields=0,  # Adjust index based on separator
            broad=1,
            input_sample="NA",
            homer_input="input_file",
            homer_size=200,
            homer_mindist=100,
            thestyle="histone"
        )