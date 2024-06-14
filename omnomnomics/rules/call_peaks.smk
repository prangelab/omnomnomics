# Rule 11: Call Peaks

## Omnomnomics Snake Rule  ##
import os
import re
import glob
import subprocess

def input_function(wildcards):
    input_folder1 = master_config['input_folders'][master_config['callpeaks_rule_num']-1][0]
    input_folder2 = master_config['input_folders'][master_config['callpeaks_rule_num']-1][1]
    if config['THETYPE'] == "CHIP":
        input_files = glob.glob(f"{input_folder1}/*.bam")
        input_files2 = glob.glob(f"{input_folder2}/*.tar.gz")
        for file in input_files2:
            input_files.append(file)
    else:
        input_files = glob.glob(f"{input_folder1}/*.bam")
    return input_files

rule call_peaks:
    input:
        input_function
        #bam_files = glob.glob(f"{master_config['input_folders'][master_config['callpeaks_rule_num']-1]}/*.bam")
        #filtered_bam = rules.touchup_bam.output,
        #homer_tagDirs = rules.create_homer_tagDir.output if config['THETYPE'] == "CHIP" else None
        #output of rule 5 and 8? of runt hij dan altijd rule 5 en 8 ook als ze niet specified zijn en output al aanwezig is? ik denk het niet dus dan prima maar wel even testen.
    output:
        # f"{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/{{THENAME}}.bed", ###, just add a random file for this rule to rule all and create that file here and then delete it in final steps
        # f"{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/all_groups.merged_peaks.bed" if config['THETYPE'] == 'ATAC' else None,
        # f"{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/all_groups.merged_peaks.bed" if config['THETYPE'] == 'CHIP' else None
        f"{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/extra.tmp"
    params:
        thetype= lambda wildcards: config['THETYPE'],  
        broad= lambda wildcards: config['BROAD'],
        input_sample= lambda wildcards: config['INPUT'],
        homer_input= lambda wildcards: config['HOMERINPUT'],
        genome= lambda wildcards: config['THEGENOME'],  
        experiment_dir= lambda wildcards: config['EXPERIMENT_DIR'], 
        name_fields= lambda wildcards: config['NAMEFIELDS'],  
        separator= lambda wildcards: config['THESEPARATOR'],
        thetype_field = lambda wildcards: config['THETYPEFIELD'] ,
        style= lambda wildcards: config['THESTYLE'],
        homersize= lambda wildcards: config['HOMERSIZE'],
        homermindist= lambda wildcards: config['HOMERMINDIST'],
        inputfolder1=lambda wildcards: master_config['input_folders'][master_config['callpeaks_rule_num']-1][0],
        inputfolder2= lambda wildcards: master_config['input_folders'][master_config['callpeaks_rule_num']-1][1],
        outputfolder= lambda wildcards: master_config['output_folders'][master_config['callpeaks_rule_num']-1]
    threads:
        lambda wildcards: Threads_Per_Rule['11']
    resources:
        #mem_mb = lambda wildcards: Memory_Per_Rule['11']
    run: 
        log_it(logfile, "Calling Peaks...", f"EXECUTING STEP {master_config['callpeaks_rule_num']}")
        log_it(logfile, f"Input folders: {params.inputfolder1} and {params.inputfolder2}")
        log_it(logfile, f"Output folder: {params.outputfolder}")
        def call_peaks(logfile, thetype, inputfolder1, inputfolder2, outputfolder, separator, thetype_field, name_fields, broad, input_sample, homer_input, homer_size, homer_mindist, thestyle):
            if thetype == "RNA":
                log_it(logfile, "Not a ChIP- or ATAC-seq experiment, skipping this step...")
                return
            
            log_it(logfile, "Calling peaks...", f"EXECUTING STEP {step_callpeaks}") #place this down below only if type is CHiP, if type is ATAC, call open regions?
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
                sanity_check_dir(logfile, inputfolder1,  master_config['input_file_types'][master_config['callpeaks_rule_num']-1]) ###############might have to adjust this based on what I do with input and output folders
                
                # Distribute samples over groups based on the file name pattern
                bams = [os.path.basename(bam).split(separator)[thetype_field - 1] for bam in glob.glob(f"{inputfolder1}/*.bam")]
                chip_groups = sorted(set(bams))
                
                for group in chip_groups:
                    # BAMS=( $(ls "$THEINFOLDER" | grep $(echo $THEGROUP | sed "s|$THESEPARATOR|.*|g") - | grep ".bam$") )
                    bams = [bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam)] 
                    bams = [os.path.join(inputfolder1, bam) for bam in bams]

                    # TRACKS=($(ls "$THEINFOLDER" | grep $(echo $SUPERHUB | sed "s|$THESEPARATOR|.*|g") - | grep ".bw$"))
                    # tracks = [track for track in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", track)]
                    
                    log_it(logfile, f"Calling peaks for: {group}...")
                    log_it(logfile, f"Files in group: {', '.join(bams)}")

                    bam_basename = os.path.basename(bams[0])
                    fields = bam_basename.split(separator)

                    # Function to parse the NAMEFIELDS and return the specified field indices
                    def parse_namefields(namefields):
                        indices = []
                        for part in namefields.split(','):
                            if '-' in part:
                                start, end = map(int, part.split('-'))
                                indices.extend(range(start, end + 1))
                            else:
                                indices.append(int(part))
                        return sorted(set(indices))  # Remove duplicates and sort the indices
                    field_indices = parse_namefields(name_fields)

                    the_name = separator.join(fields[i-1] for i in field_indices if i-1 < len(fields))

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
                                macs3 callpeak -t {' '.join(bams)} -c {input_sample} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -{'q' if ext.startswith('q') else 'p'} {'1e-9' if ext.endswith('9.broad') else '1e-6'} --broad --verbose 0""")
                        else:
                            log_it(logfile, "Calling peaks with MACS3, using p and q values 1e-6 and 1e-9...")
                            for ext in ["q-9", "q-6", "p-9", "p-6"]:
                                shell(f""" micromamba activate macs3 && \
                                macs3 callpeak -t {' '.join(bams)} -c {input_sample} --outdir {outputfolder} -n {the_name}.MACS3.{ext} -{'q' if ext.startswith('q') else 'p'} {'1e-9' if ext.endswith('9') else '1e-6'} --verbose 0""")
                # Clean up the output
                for group in chip_groups: 
                    the_name = os.path.basename(
                        next(bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam))
                    ).split(separator)[name_fields - 1] ##### why -1? # Set group name
                    # for file in glob.glob(f"{outputfolder}/{the_name}.MACS3*"): #################alternate (possible) translation if code under this somehow doesnt work
                    #     if not re.search(r"(broadPeak|narrowPeak)$", file): 
                    #         os.remove(file) 
                    shell(f"find {outputfolder} -type f -name {the_name}.MACS3* ! \( -name *broadPeak -o -name *narrowPeak \) -delete ") # Delete all redundant MACS3 output


                log_it(logfile, "Cleaning up MACS3 output (keeping only narrowPeak and broadPeak files)...")

                log_it(logfile, "Converting to a clean 3 column BED format...")
                for file in glob.glob(f"{outputfolder}/*{'narrowPeak' if broad != '1' else 'broadPeak'}"):
                    sorted_bed = f"{outputfolder}/{os.path.basename(file).replace('narrowPeak', '').replace('broadPeak', '')}.bed"
                    shell(f"cut -f1-3 {file} | sort -k1,1 -k2,2n -k3,3n > {sorted_bed}")
                    os.remove(file)

                # Call peaks using HOMER
                #inputfolder = "HOMER_tagDirs"
                log_it(logfile, "Calling peaks with HOMER...")
                log_it(logfile, f"Input folder: {inputfolder2}")
                
                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder2,  master_config['input_file_types'][master_config['callpeaks_rule_num']-1]) ###############might have to adjust this based on what I do with input and output folders
                
                # See if we need to unpack the tag dirs
                if glob.glob(f"{inputfolder2}/*tagDir.tar.gz"):
                    for tagdir in glob.glob(f"{inputfolder2}/*tagDir.tar.gz"):
                        tagdir_basename = os.path.basename(tagdir)
                        shell(f"cd {inputfolder2} && tar --strip-components=1 -xzf {tagdir_basename}")
                # Fetch tagdirs
                tagdirs = glob.glob(f"{inputfolder2}/*.HOMER_tagDir")
                
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
                        shell(f"pos2bed.pl {tagdir}/peaks.txt | grep -v '#\\|alt\\|Un\\|random' | sort -k1,1 -k2,2n -k3,3n | cut -f-3 > {tagdir}/peaks.bed ")
                    elif thestyle == "histone":
                        shell(f"pos2bed.pl {tagdir}/regions.txt | grep -v '#\\|alt\\|Un\\|random' | sort -k1,1 -k2,2n -k3,3n | cut -f-3 > {tagdir}/regions.bed ")

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
                    beds = [os.path.join(inputfolder2, bed) for bed in os.listdir(inputfolder2) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.HOMER_tagDir$", bed)]
                    if thestyle == "factor":
                        beds = [f"{bed}/peaks.bed" for bed in beds]
                    elif thestyle == "histone":
                        beds = [f"{bed}/regions.bed" for bed in beds]
                    
                    the_name = os.path.basename(
                        next(bed for bed in os.listdir(inputfolder2) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.HOMER_tagDir$", bed))
                    ).split(separator)[name_fields - 1]

                    sorted_bed = f"{outputfolder}/{the_name}.HOMER.merged_peaks.bed"
                    shell(f"cat {' '.join(beds)} | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {sorted_bed}")
            else:
                # If ATAC, call open regions
                log_it(logfile, "Calling open regions...")

                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder1,  master_config['input_file_types'][master_config['callpeaks_rule_num']-1]) ###############might have to adjust this based on what I do with input and output folders

                # Distribute samples over groups based on the file name pattern
                bams = [os.path.basename(bam).split(separator)[thetype_field - 1] for bam in glob.glob(f"{inputfolder1}/*.bam")]
                atac_groups = sorted(set(bams))

                for group in atac_groups: # Iterate over the groups
                    bams = [bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam)]
                    bams = [os.path.join(inputfolder1, bam) for bam in bams] # Fetch BAM files
                    
                    log_it(logfile, f"Calling open regions for: {group}...")
                    log_it(logfile, f"Files in group: {', '.join(bams)}")
                    # Set group name
                    the_name = os.path.basename(
                        next(bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam))
                    ).split(separator)[name_fields - 1]

                    # Call open regions with MACS3 hmmr
                    # Turns out this MACS subroutine is till buggy and not functining prooerply, check back with it later but for now just call as ChIP peaks
                    #logIt "Calling open regions with MACS3 hmmratac..."
                    #logIt "	macs3 hmmratac -b $BAMS --outdir \"$THEOUTFOLDER\" -n $NAMEFIELDS --verbose 0 &"
                    #macs3 hmmratac -b ${BAMS[*]} --outdir "$THEOUTFOLDER" -n $THENAME --verbose 2 &

                    log_it(logfile, "Calling peaks with MACS3, q value 1e-9...")
                    log_it(logfile, f"macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.q-9 -q 1e-9 --verbose 0" )
                    shell(f""" micromamba activate macs3 && \
                        macs3 callpeak -t {' '.join(bams)} --outdir {outputfolder} -n {the_name}.MACS3.q-9 -q 1e-9 --verbose 0 """)

                log_it(logfile, "Cleaning up MACS3 output (keeping only narrowPeak files)...")
                for group in atac_groups: ## Clean up the output
                    # Set group name
                    the_name = os.path.basename(
                        next(bam for bam in os.listdir(inputfolder1) if re.match(re.escape(group).replace(re.escape(separator), ".*") + ".*\\.bam$", bam))
                    ).split(separator)[name_fields - 1]
                    # Delete all redundant MACS3 output
                    shell(f"find {outputfolder} -type f -name {the_name}.MACS3* ! -name *narrowPeak -delete")
                
                # Create clean BED files
                log_it(logfile, "Converting to a clean 3 column BED format...")
                for file in glob.glob(f"{outputfolder}/*narrowPeak"):
                    sorted_bed = f"{outputfolder}/{os.path.basename(file).replace('narrowPeak', '')}.bed"
                    shell(f"cut -f1-3 {file} | sort -k1,1 -k2,2n -k3,3n > {sorted_bed}")
                    os.remove(file)

                # Concatenate the peak files
                shell(f"cat {outputfolder}/*.bed | sort -k1,1 -k2,2n -k3,3n | mergeBed -i - > {outputfolder}/all_groups.merged_peaks.bed")
                log_it(logfile, f"Total merged peaks: {subprocess.getoutput(f'wc -l {outputfolder}/all_groups.merged_peaks.bed').strip().split()[0]}")
            
            
            # Clean up
            # os.remove("groups.tmp.txt")
            # shell("micromamba deactivate")

        call_peaks(
            logfile,
            thetype=params.thetype,
            inputfolder1=params.inputfolder1,
            inputfolder2=params.inputfolder2,
            outputfolder=params.outputfolder,
            separator=params.separator,
            thetype_field=int(params.thetype_field),
            name_fields=params.name_fields,  
            broad=params.broad,
            input_sample=params.input_sample,
            homer_input=int(params.homer_input),
            homer_size=int(params.homersize),
            homer_mindist=params.homermindist,
            thestyle=params.style
        )
        shell(f"""echo "necessity file for callpeaks. can delete this." > peak_calling/extra.tmp""")


#### for the 4 occasions of namefields with re
import os
import re

def parse_namefields(namefields_str):
    fields = []
    for part in namefields_str.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            fields.extend(range(start, end + 1))
        else:
            fields.append(int(part))
    return fields

# Define the variables
THEINFOLDER = "/data/input"
THEGROUP = "sample_group"
THESEPARATOR = "_"
NAMEFIELDS = "1,2,4-6"

# Step 1: List the contents of the input folder
all_files = os.listdir(THEINFOLDER)

# Step 2: Filter the list using a regex based on THEGROUP and THESEPARATOR
group_pattern = re.escape(THEGROUP).replace(re.escape(THESEPARATOR), ".*")
filtered_files = [f for f in all_files if re.match(group_pattern + r".*", f)]

# Step 3: Further filter the list to only include .bam files
bam_files = [f for f in filtered_files if f.endswith(".bam")]

# Step 4: Select the first match from the filtered list
if bam_files:
    first_bam = bam_files[0]
else:
    raise FileNotFoundError("No matching BAM files found.")

# Step 5: Extract the basename of the selected file
basename = os.path.basename(first_bam)

# Step 6: Split the basename using THESEPARATOR and select the fields specified by NAMEFIELDS
fields = basename.split(THESEPARATOR)
namefields_indices = parse_namefields(NAMEFIELDS)
THENAME = THESEPARATOR.join(fields[i-1] for i in namefields_indices if i-1 < len(fields))

print(THENAME)
