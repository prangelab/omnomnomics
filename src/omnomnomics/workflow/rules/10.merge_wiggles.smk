# Rule 10 merge wiggles

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import re
import shutil
import subprocess

rule merge_wiggles:
    input:
        extra_file = expand(f"{experiment_dir}/{master_config['input_folders'][master_config['mergewig_rule_num']-1]}/{{sample}}.extra_9.tmp", sample = samples2) if 9 in themode else []
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['mergewig_rule_num']-1]}/extra_10.tmp"
    params:
        thegenome=lambda wildcards: config['THEGENOME'],
        coltable=lambda wildcards: config["THECOLTABLE"],
        theappendix=lambda wildcards: config['THEAPPENDIX'],
        thetype=lambda wildcards: config['THETYPE'],
        theseparator=lambda wildcards: config['THESEPARATOR'],
        thecolfield=lambda wildcards: config['THECOLFIELD'],
        thenamefields=lambda wildcards: config['NAMEFIELDS'],
        thetypefield=lambda wildcards: config['THETYPEFIELD'],
        thehubmail=lambda wildcards: config['THEHUBMAIL'],
        theoverlay=lambda wildcards: config['THEOVERLAY'],
        inputfolder=lambda wildcards: f"{experiment_dir}/{master_config['input_folders'][master_config['mergewig_rule_num']-1]}",
        outputfolder=lambda wildcards: f"{experiment_dir}/{master_config['output_folders'][master_config['mergewig_rule_num']-1]}"
    threads:
        lambda wildcards: Threads_Per_Rule['10']
    resources:
      mem_mb=lambda wildcards: Memory_Per_Rule['10'],
      partition=lambda wildcards: master_config['partition'],
      runtime=lambda wildcards: Runtime_Per_Rule['10']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['mergewig_rule_num']-1]}/benchmarks/merge_wiggles_benchmark.tsv"
    run:
        log_it(logfile, "Merging Wiggles and TrackHubs...", f"EXECUTING STEP {master_config['mergewig_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        # Function to parse the NAMEFIELDS and return the specified field indices
        def parse_namefields(namefields):
            indices = []
            for part in namefields.split(','): #Split all the fields
                if '-' in part:
                    start, end = map(int, part.split('-')) 
                    indices.extend(range(start, end + 1)) #Add the wanted namefields
                else:
                    indices.append(int(part)) #Add the wanted namefields
            return sorted(set(indices))  # Remove duplicates and sort the indices

        def merge_wig(input_folder, output_folder, type_, col_field, separator, type_field, name_fields, col_table, appendix, genome, hub_mail, overlay):
            log_it(logfile, f"Merging Wiggles and TrackHubs...\nInput folder: {input_folder}\nOutput folder: {output_folder}")

            # RNA-seq will already be in hubs so we need to copy files over into a fresh 'group' hub
            if type_ == "RNA":
                # Sanity check the working dir
                sanity_check_dir(logfile, input_folder,  master_config['input_file_types'][master_config['mergewig_rule_num']-1][0])
                
                # Build coltype list
                coltypes = sorted(set(
                    os.path.basename(hub).split(separator)[col_field - 1] for hub in os.listdir(input_folder) if hub.endswith(".hub")
                ))

                # Build merged hub list
                hubtypes = sorted(set(
                    os.path.basename(hub).split(separator)[type_field - 1] for hub in os.listdir(input_folder) if hub.endswith(".hub")
                ))

                # Set colors
                with open(os.path.join(config['COLOR_DATA_FOLDER_DEFAULT'], "poscols.hub")) as f:
                    poscols = f.read().splitlines()

                with open(os.path.join(config['COLOR_DATA_FOLDER_DEFAULT'], "negcols.hub")) as f:
                    negcols = f.read().splitlines()

                # Iterate over the hub types
                for superhub in hubtypes:
                    hubs = [
                        hub for hub in os.listdir(input_folder)
                        if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.hub$", hub)
                    ]
                    # Iterate over individual hubs
                    for hub in hubs:
                        # Set name, color and type (desired hub) from the folder name
                        hub_path = os.path.join(input_folder, hub)
                        hub_basename = os.path.basename(hub_path)
                        fields = hub_basename.split(separator)

                        field_indices = parse_namefields(name_fields)

                        name = separator.join(fields[i-1] for i in field_indices if i-1 < len(fields))

                        htype = hub_basename.split(separator)[type_field - 1]
                        coltype = hub_basename.split(separator)[col_field - 1]

                        log_it(logfile, f"Processing {superhub} {name}...")

                        # Set track colors based on color field
                        poscol = next(poscols[i] for i, ct in enumerate(coltypes) if ct == coltype)
                        negcol = next(negcols[i] for i, ct in enumerate(coltypes) if ct == coltype)

                        merged_hub_folder = os.path.join(input_folder, f"{htype}.{appendix}")
                        genome_folder = os.path.join(merged_hub_folder, genome)

                        # Check if the merged hub folder already exists. If not, make it and populate it.
                        if not os.path.exists(merged_hub_folder):

                            # Make dir structure
                            os.makedirs(genome_folder)

                            # Copy the hub files
                            shutil.copyfile(os.path.join(hub_path, "hub.txt"), os.path.join(merged_hub_folder, "hub.txt"))
                            shutil.copyfile(os.path.join(hub_path, "genomes.txt"), os.path.join(merged_hub_folder, "genomes.txt"))
                            shutil.copyfile(os.path.join(hub_path, f"{genome}/trackDb.txt"), os.path.join(genome_folder, "trackDb.txt"))

                            # Edit the hub file and trackDb
                            subprocess.run(["sed", "-i", f"s|{hub_basename}|{htype}.{appendix}|", os.path.join(merged_hub_folder, "hub.txt")])
                            subprocess.run(["sed", "-i", f"s|maxHeightPixels [0-9]\+:[0-9]\+:[0-9]\+|maxHeightPixels 64:32:8|;s|{hub_basename}|{htype}.{appendix}|", os.path.join(genome_folder, "trackDb.txt")])
                        else:  # If it does already exist, append track info to the trackDb.
                            subprocess.run(['sed', '-n', '-e', '/^$/,$p', os.path.join(input_folder, hub, genome, 'trackDb.txt'), '>>', os.path.join(input_folder, htype + '.' + appendix, genome, 'trackDb.txt')])

                        # Edit the trackDb
                        subprocess.run(['sed', '-i', f's|{hub_basename}|{htype}.{appendix}|;s|{os.path.splitext(hub_basename)[0]}.HOMER_tagDir+|{name}.{appendix}+|;s|{os.path.splitext(hub_basename)[0]}.HOMER_tagDir-|{name}.{appendix}-|', os.path.join(input_folder, htype + '.' + appendix, genome, 'trackDb.txt')])
                        path = os.path.join(input_folder, htype + '.' + appendix, genome, 'trackDb.txt')
                        homerposcol = master_config['homerposcol']
                        altposcol = master_config['altposcol']
                        homernegcol = master_config['homernegcol']
                        altnegcol = master_config['altnegcol']
                        subprocess.run(['sed', '-i', '-E', f's|{homerposcol}\|{altposcol}|{poscol}|;s|{homernegcol}\|{altnegcol}|{negcol}|;', f"{path}" ])
                        subprocess.run(['sed', '-i', 's/\"//g', os.path.join(input_folder, htype + '.' + appendix, genome, 'trackDb.txt')])

                        # Copy the bigWig files
                        shutil.copyfile(os.path.join(hub_path, f"{genome}/{hub_basename.replace('.hub', '.HOMER_tagDirneg.ucsc.bigWig')}"), os.path.join(genome_folder, f"{hub_basename.replace('.hub', '.HOMER_tagDirneg.ucsc.bigWig')}"))
                        shutil.copyfile(os.path.join(hub_path, f"{genome}/{hub_basename.replace('.hub', '.HOMER_tagDirpos.ucsc.bigWig')}"), os.path.join(genome_folder, f"{hub_basename.replace('.hub', '.HOMER_tagDirpos.ucsc.bigWig')}"))

                    # Copy merged hub to output folder if needed
                    if input_folder != output_folder:
                        log_it(logfile, f"Moving merged hub {htype}.{appendix} to {output_folder}...")
                        # If a merged hub already exists in the output folder, remove it first.
                        if os.path.isdir(os.path.join(output_folder, f'{htype}.{appendix}')):
                            shutil.rmtree(os.path.join(output_folder, f'{htype}.{appendix}'))
                        shutil.move(merged_hub_folder, output_folder)

            else:
                ## If ChIP or ATAC, we start from BigWigs and need to create the 'group' hubs from scratch	
                # Sanity check the working dir
                sanity_check_dir(logfile, input_folder, master_config['input_file_types'][master_config['mergewig_rule_num']-1][1])
                
                # Set color default
                thecol = "0,255,0"

                ## Set color table(s)
                # Check if we have a .txt file
                col_array = []
                if col_table.endswith(".txt"):
                    # We have a list file of color tables, let's put them in the array
                    log_it(logfile, f"Color table file: {col_table}")
                    with open(col_table) as f:
                        col_array = f.read().splitlines()
                else:
                    # Not a list file, so it's a color table. Let's put it in the array
                    col_array.append(col_table)

                # Build hub list
                hubtypes = sorted(set(
                    os.path.basename(bw).split(separator)[type_field - 1] for bw in os.listdir(input_folder) if bw.endswith(".bw")
                ))

                # Populate the hubs
                ctabletracker = 0
                for superhub in hubtypes:
                    # Build coltype list
                    coltypes = sorted(set(
                        os.path.basename(bw).split(separator)[col_field - 1] for bw in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", bw)
                    ))

                    # Loop around to reuse the first element if we have less color tables than hubs)
                    if ctabletracker == len(col_array):
                        ctabletracker = 0

                    # Get total number of colors required
                    ncols = len([bw for bw in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", bw)])
                    # Get total number of colors in the color table.
                    tcols = sum(1 for _ in open(col_array[ctabletracker]))

                    # Make sure there are enough colors in the table
                    if tcols <= ncols:
                        diffcols = ncols - tcols
                        col_int = 1
                    else:
                        # Get the interval we need to spread our colors through the table. (We do -1 on total colors to avoid taking the last near black or white row)
                        col_int = (tcols - 1) // ncols

                    # Get the corresponding lines from the color tabl
                    mycols = [line.strip() for idx, line in enumerate(open(col_array[ctabletracker])) if idx % col_int == 0]

                    if tcols <= ncols:
                        mycols.extend([line.strip() for idx, line in enumerate(open(col_array[ctabletracker])) if idx % col_int == 0][:diffcols])

                    # Iterate over individual hubs
                    for bw in [f for f in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", f)]:
                        bw_path = os.path.join(input_folder, bw)
                        bw_basename = os.path.basename(bw_path)
                        fields = bw_basename.split(separator)
                        field_indices = parse_namefields(name_fields)

                        #Set name,color and type (desired hub) from the folder name
                        name = separator.join(fields[i-1] for i in field_indices if i-1 < len(fields))
                        htype = bw_basename.split(separator)[type_field - 1]
                        coltype = bw_basename.split(separator)[col_field - 1]

                        # Set track color based on track number
                        tracks = [track for track in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", track)]
                        for i, track_file in enumerate(tracks):
                            if track_file == bw:
                                break
                        thecol  = mycols[i]

                        log_it(logfile, f"Processing {name}...")

                        merged_hub_folder = os.path.join(input_folder, f"{htype}.{appendix}")
                        genome_folder = os.path.join(merged_hub_folder, genome)
                        
                        # Check if the merged hub folder already exists. If not, make it and populate it.
                        if not os.path.exists(merged_hub_folder):
                            os.makedirs(genome_folder)

                            with open(os.path.join(merged_hub_folder, "genomes.txt"), "w") as genomes_file:
                                genomes_file.write(f"genome {genome}\ntrackDb {genome}/trackDb.txt\n")

                            with open(os.path.join(merged_hub_folder, "hub.txt"), "w") as hub_file:
                                hub_file.write(f"hub {htype}.{appendix}\nshortLabel {htype}.{appendix}\nlongLabel {htype}.{appendix}\ngenomesFile genomes.txt\nemail {hub_mail}\n")

                            with open(os.path.join(genome_folder, "trackDb.txt"), "w") as trackdb_file:
                                trackdb_file.write(f"track {htype}.{appendix}\ncontainer multiWig\nnoInherit on\nshortLabel {htype}.{appendix}\nlongLabel {htype}.{appendix}\ntype bigWig\nconfigurable on\nvisibility full\naggregate {overlay}\nshowSubtrackColorOnUi on\nautoScale on\nwindowingFunction maximum\nsmoothingWindow 2\npriority 1.4\nyLineMark 0\nyLineOnOff on\nmaxHeightPixels 64:32:8\n")

                        with open(os.path.join(genome_folder, "trackDb.txt"), "a") as trackdb_file:
                            trackdb_file.write(f"\ntrack {name}.{appendix}\nbigDataUrl {bw_basename}\nshortLabel {name}.{appendix}\nlongLabel {name}.{appendix}\ntype bigWig\nparent {htype}.{appendix}\ncolor {thecol}\n")
                        #Copy the bigWig file
                        shutil.copyfile(bw_path, os.path.join(genome_folder, bw_basename))

                    ctabletracker += 1

                # Copy merged hubs to output folder if needed
                if input_folder != output_folder:
                    log_it(logfile, f"Moving merged hubs from {input_folder} to {output_folder}...")
                    for folder in os.listdir(input_folder):
                        folder_path = os.path.join(input_folder, folder)
                        if os.path.isdir(folder_path) and folder.endswith(appendix):
                            # If a merged hub already exists in the output folder, remove it first.
                            if os.path.isdir(os.path.join(output_folder, folder)):
                                shutil.rmtree(os.path.join(output_folder, folder))
                            shutil.move(folder_path, output_folder)

            shell(f"""echo "necessity file for merge wiggle. can delete this." > {output_folder}/extra_10.tmp""")
            log_it(logfile, "Merge Wiggles and TrackHubs completed!")

        merge_wig(
            input_folder=params.inputfolder, output_folder=params.outputfolder,
            type_=params.thetype, col_field=int(params.thecolfield), separator=params.theseparator,
            type_field=int(params.thetypefield), name_fields=params.thenamefields,
            col_table=params.coltable, appendix=params.theappendix,
            genome=params.thegenome, hub_mail=params.thehubmail,
            overlay=params.theoverlay
        )
