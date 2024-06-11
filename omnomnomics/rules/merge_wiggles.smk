# Rule 10 merge wiggles

## Omnomnomics Snake Rule  ##
import os
import re
import shutil
import subprocess

rule merge_wig:
    input:
        theinfolder=lambda wildcards: master_config['input_folders'][master_config['mergewig_rule_num']-1],
        coltable=config["COLTABLE"]
        #use rules."rule_name".output
    output:
        theoutfolder=lambda wildcards: master_config['output_folders'][master_config['mergewig_rule_num']-1]
        #directory()
    params:
        thegenome=config['THEGENOME'],
        theappendix=config['THEAPPENDIX'],
        thetype= config['THETYPE']
        theseparator=config['THESEPARATOR'],
        thecolfield=config['THECOLFIELD'],
        thenamefields=config['THENAMEFIELDS'],
        thetypefield=config['THETYPEFIELD'],
        thehubmail=config['THEHUBMAIL'],
        theoverlay=config['THEOVERLAY'], 
        inputfolder = master_config['input_folders'][master_config['mergewig_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['mergewig_rule_num']-1]
    threads:
        Threads_Per_Rule['10']
    resources:
        mem_mb = Memory_Per_Rule['10']
    run:
        logIt(logfile, "Merging Wiggles and TrackHubs...", f"EXECUTING STEP {master_config['mergewig_rule_num']}")
        logIt(logfile, f"Input folder: {input.theinfolder}")
        logIt(logfile, f"Output folder: {output.theoutfolder}")

        def merge_wig(input_folder, output_folder, type_, col_field, separator, type_field, name_fields, col_table, appendix, genome, hub_mail, overlay):
            log_it(f"Merging Wiggles and TrackHubs...\nInput folder: {input_folder}\nOutput folder: {output_folder}")

            if type_ == "RNA":
                sanity_check_dir(logfile, input_folder,  master_config['input_file_types'][master_config['mergewig_rule_num']-1][0])
                coltypes = sorted(set(
                    os.path.basename(hub).split(separator)[col_field - 1] for hub in os.listdir(input_folder) if hub.endswith(".hub")
                ))

                hubtypes = sorted(set(
                    os.path.basename(hub).split(separator)[type_field - 1] for hub in os.listdir(input_folder) if hub.endswith(".hub")
                ))

                with open(os.path.join(os.environ["OMNOM_HOME"], "bin/color_data_for_hubs/poscols.hub")) as f:
                    poscols = f.read().splitlines()

                with open(os.path.join(os.environ["OMNOM_HOME"], "bin/color_data_for_hubs/negcols.hub")) as f:
                    negcols = f.read().splitlines()

                for superhub in hubtypes:
                    hubs = [
                        hub for hub in os.listdir(input_folder)
                        if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.hub$", hub)
                    ]
                    for hub in hubs:
                        hub_path = os.path.join(input_folder, hub)
                        hub_basename = os.path.basename(hub_path)
                        name = hub_basename.split(separator)[name_fields - 1]
                        htype = hub_basename.split(separator)[type_field - 1]
                        coltype = hub_basename.split(separator)[col_field - 1]

                        log_it(f"Processing {superhub} {name}...")

                        poscol = next(poscols[i] for i, ct in enumerate(coltypes) if ct == coltype)
                        negcol = next(negcols[i] for i, ct in enumerate(coltypes) if ct == coltype)

                        merged_hub_folder = os.path.join(input_folder, f"{htype}.{appendix}")
                        genome_folder = os.path.join(merged_hub_folder, genome)

                        if not os.path.exists(merged_hub_folder):
                            os.makedirs(genome_folder)
                            shutil.copyfile(os.path.join(hub_path, "hub.txt"), os.path.join(merged_hub_folder, "hub.txt"))
                            shutil.copyfile(os.path.join(hub_path, "genomes.txt"), os.path.join(merged_hub_folder, "genomes.txt"))
                            shutil.copyfile(os.path.join(hub_path, f"{genome}/trackDb.txt"), os.path.join(genome_folder, "trackDb.txt"))

                            subprocess.run(["sed", "-i", f"s|{hub_basename}|{htype}.{appendix}|", os.path.join(merged_hub_folder, "hub.txt")])
                            subprocess.run(["sed", "-i", f"s|maxHeightPixels [0-9]+:[0-9]+:[0-9]+|maxHeightPixels 64:32:8|;s|{hub_basename}|{htype}.{appendix}|", os.path.join(genome_folder, "trackDb.txt")])

                            # with open(os.path.join(merged_hub_folder, "hub.txt"), "a") as hub_file:
                            #     hub_file.write(f"\ntrack {htype}.{appendix}\ncontainer multiWig\nnoInherit on\nshortLabel {htype}.{appendix}\nlongLabel {htype}.{appendix}\ntype bigWig\nconfigurable on\nvisibility full\naggregate {overlay}\nshowSubtrackColorOnUi on\nautoScale on\nwindowingFunction maximum\nsmoothingWindow 2\npriority 1.4\nyLineMark 0\nyLineOnOff on\nmaxHeightPixels 64:32:8\n")
                        else:
                            # with open(os.path.join(genome_folder, "trackDb.txt"), "a") as trackdb_file:
                            #     trackdb_file.write(f"\ntrack {name}.{appendix}\nbigDataUrl {hub_basename}\nshortLabel {name}.{appendix}\nlongLabel {name}.{appendix}\ntype bigWig\nparent {htype}.{appendix}\ncolor {poscol if coltype in coltypes else negcol}\n")
                            subprocess.run(['sed', '-n', '-e', '/^$/,$p', os.path.join(input_folder, hub, genome, 'trackDb.txt'), '>>', os.path.join(input_folder, htype + '.' + appendix, genome, 'trackDb.txt')])

                        subprocess.run(['sed', '-i', f's|{hub_basename}|{htype}.{appendix}|;s|{os.path.splitext(hub_basename)[0]}.HOMER_tagDir+|{name}.{appendix}+|;s|{os.path.splitext(hub_basename)[0]}.HOMER_tagDir-|{name}.{appendix}-|', os.path.join(input_folder, htype + '.' + appendix, genome, 'trackDb.txt')])
                        subprocess.run(['sed', '-i', '-E', f's|{master_config['homerposcol']}|{master_config['altposcol']}|{poscol}|;s|{master_config['homernegcol']}|{master_config['altnegcol']}|{negcol}|;', os.path.join(input_folder, htype + '.' + appendix, genome, 'trackDb.txt')])
                        subprocess.run(['sed', '-i', 's/\"//g', os.path.join(input_folder, htype + '.' + appendix, genome, 'trackDb.txt')])

                        shutil.copyfile(os.path.join(hub_path, f"{genome}/{hub_basename.replace('.hub', '.HOMER_tagDirneg.ucsc.bigWig')}"), os.path.join(genome_folder, f"{hub_basename.replace('.hub', '.HOMER_tagDirneg.ucsc.bigWig')}"))
                        shutil.copyfile(os.path.join(hub_path, f"{genome}/{hub_basename.replace('.hub', '.HOMER_tagDirpos.ucsc.bigWig')}"), os.path.join(genome_folder, f"{hub_basename.replace('.hub', '.HOMER_tagDirpos.ucsc.bigWig')}"))

                    if input_folder != output_folder:
                        log_it(f"Moving merged hub {merged_hub_folder} to {output_folder}...")
                        if os.path.isdir(os.path.join(output_folder, merged_hub_folder)):
                            shutil.rmtree(os.path.join(output_folder, merged_hub_folder))
                        shutil.move(merged_hub_folder, output_folder)

            else:
                sanity_check_dir(logfile, input_folder, ".bw")
                thecol = "0,255,0"
                col_array = []
                if col_table.endswith(".txt"):
                    log_it(f"Color table file: {col_table}")
                    with open(col_table) as f:
                        col_array = f.read().splitlines()
                else:
                    col_array.append(col_table)

                hubtypes = sorted(set(
                    os.path.basename(bw).split(separator)[type_field - 1] for bw in os.listdir(input_folder) if bw.endswith(".bw")
                ))

                ctabletracker = 0
                for superhub in hubtypes:
                    coltypes = sorted(set(
                        os.path.basename(bw).split(separator)[col_field - 1] for bw in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", bw)
                    ))

                    if ctabletracker == len(col_array):
                        ctabletracker = 0

                    ncols = len([bw for bw in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", bw)])
                    tcols = sum(1 for _ in open(col_array[ctabletracker]))

                    if tcols <= ncols:
                        diffcols = ncols - tcols
                        col_int = 1
                    else:
                        col_int = (tcols - 1) // ncols

                    mycols = [line.strip() for idx, line in enumerate(open(col_array[ctabletracker])) if idx % col_int == 0]

                    if tcols <= ncols:
                        mycols.extend([line.strip() for idx, line in enumerate(open(col_array[ctabletracker])) if idx % col_int == 0][:diffcols])

                    for bw in [f for f in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", f)]:
                        bw_path = os.path.join(input_folder, bw)
                        bw_basename = os.path.basename(bw_path)
                        name = bw_basename.split(separator)[name_fields - 1] ### why -1? in bash code not
                        htype = bw_basename.split(separator)[type_field - 1]
                        coltype = bw_basename.split(separator)[col_field - 1]

                        # Set track color based on track number
                        tracks = [track for track in os.listdir(input_folder) if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", track)]
                        for i, track_file in enumerate(TRACKS):
                            if track_file == bw:
                                break
                        thecol  = mycols[i]

                        log_it(f"Processing {name}...")

                        merged_hub_folder = os.path.join(input_folder, f"{htype}.{appendix}")
                        genome_folder = os.path.join(merged_hub_folder, genome)

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

                        shutil.copyfile(bw_path, os.path.join(genome_folder, bw_basename))

                    ctabletracker += 1

                # Copy merged hubs to output folder if needed
                if input_folder != output_folder:
                    log_it(f"Moving merged hubs from {input_folder} to {output_folder}...")
                    for folder in os.listdir(input_folder):
                        folder_path = os.path.join(input_folder, folder)
                        if os.path.isdir(folder_path) and folder.endswith(appendix):
                            # If a merged hub already exists in the output folder, remove it first.
                            if os.path.isdir(os.path.join(output_folder, folder)):
                                shutil.rmtree(os.path.join(output_folder, folder))
                            shutil.move(folder_path, output_folder)

                log_it("Merge Wiggles and TrackHubs completed!")

        merge_wig(
            input_folder=params.inputfolder, output_folder=params.outputfolder,
            type_=params.thetype, col_field=params.thecolfield, separator=params.theseparator,
            type_field=params.thetypefield, name_fields=params.thenamefields,
            col_table=params.coltable, appendix=params.theappendix,
            genome=params.thegenome, hub_mail=params.thehubmail,
            overlay=params.theoverlay
        )
