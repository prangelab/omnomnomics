# Rule 9 merge wiggles

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import re
import shutil

rule merge_wiggles:
    input:
        extra_file=expand(f"{experiment_dir}/{master_config['input_folders'][master_config['mergewig_rule_num']-1]}/{{sample}}.extra_8.tmp", sample=samples2) if 8 in themode else []
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['mergewig_rule_num']-1]}/extra_9.tmp"
    params:
        thegenome=lambda wildcards: config['THEGENOME'],
        coltable=lambda wildcards: config["THECOLTABLE"],
        theappendix=lambda wildcards: config['THEAPPENDIX'],
        theseparator=lambda wildcards: config['THESEPARATOR'],
        thecolfield=lambda wildcards: config['THECOLFIELD'],
        thenamefields=lambda wildcards: config['NAMEFIELDS'],
        thetypefield=lambda wildcards: config['THETYPEFIELD'],
        thehubmail=lambda wildcards: config['THEHUBMAIL'],
        theoverlay=lambda wildcards: config['THEOVERLAY'],
        inputfolder=lambda wildcards: f"{experiment_dir}/{master_config['input_folders'][master_config['mergewig_rule_num']-1]}",
        outputfolder=lambda wildcards: f"{experiment_dir}/{master_config['output_folders'][master_config['mergewig_rule_num']-1]}"
    threads:
        lambda wildcards: Threads_Per_Rule['9']
    resources:
        mem_mb=lambda wildcards: Memory_Per_Rule['9'],
        partition=lambda wildcards: master_config['partition'],
        runtime=lambda wildcards: Runtime_Per_Rule['9']
    run:
        log_once(logfile, "step9.header", "Merging BigWigs and TrackHubs...", f"EXECUTING STEP {master_config['mergewig_rule_num']}")
        log_once(logfile, "step9.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step9.outputfolder", f"Output folder: {params.outputfolder}")

        def parse_namefields(namefields):
            indices = []
            for part in namefields.split(','):
                if '-' in part:
                    start, end = map(int, part.split('-'))
                    indices.extend(range(start, end + 1))
                else:
                    indices.append(int(part))
            return sorted(set(indices))

        def ensure_hub_structure(hub_folder, genome_folder, hub_name, genome, hub_mail, overlay):
            if os.path.exists(hub_folder):
                return

            os.makedirs(genome_folder)

            with open(os.path.join(hub_folder, "genomes.txt"), "w") as genomes_file:
                genomes_file.write(f"genome {genome}\ntrackDb {genome}/trackDb.txt\n")

            with open(os.path.join(hub_folder, "hub.txt"), "w") as hub_file:
                hub_file.write(
                    f"hub {hub_name}\n"
                    f"shortLabel {hub_name}\n"
                    f"longLabel {hub_name}\n"
                    f"genomesFile genomes.txt\n"
                    f"email {hub_mail}\n"
                )

            with open(os.path.join(genome_folder, "trackDb.txt"), "w") as trackdb_file:
                trackdb_file.write(
                    f"track {hub_name}\n"
                    "container multiWig\n"
                    "noInherit on\n"
                    f"shortLabel {hub_name}\n"
                    f"longLabel {hub_name}\n"
                    "type bigWig\n"
                    "configurable on\n"
                    "visibility full\n"
                    f"aggregate {overlay}\n"
                    "showSubtrackColorOnUi on\n"
                    "autoScale on\n"
                    "windowingFunction maximum\n"
                    "smoothingWindow 2\n"
                    "priority 1.4\n"
                    "yLineMark 0\n"
                    "yLineOnOff on\n"
                    "maxHeightPixels 64:32:8\n"
                )

        def merge_wig(input_folder, output_folder, col_field, separator, type_field, name_fields, col_table, appendix, genome, hub_mail, overlay):
            sanity_check_dir(logfile, input_folder, master_config['input_file_types'][master_config['mergewig_rule_num'] - 1], "step9.sanity")

            col_array = []
            if col_table.endswith(".txt"):
                log_it(logfile, f"Color table file: {col_table}")
                with open(col_table) as handle:
                    col_array = handle.read().splitlines()
            else:
                col_array.append(col_table)

            field_indices = parse_namefields(name_fields)
            hubtypes = sorted(set(
                os.path.basename(bw).split(separator)[type_field - 1]
                for bw in os.listdir(input_folder)
                if bw.endswith(".bw")
            ))

            ctabletracker = 0
            for superhub in hubtypes:
                tracks = [
                    bw for bw in os.listdir(input_folder)
                    if re.match(re.escape(superhub).replace(re.escape(separator), ".*") + ".*\\.bw$", bw)
                ]
                tracks.sort()
                if not tracks:
                    continue

                coltypes = sorted(set(
                    os.path.basename(bw).split(separator)[col_field - 1]
                    for bw in tracks
                ))

                if ctabletracker == len(col_array):
                    ctabletracker = 0

                with open(col_array[ctabletracker]) as handle:
                    table_colors = [line.strip() for line in handle if line.strip()]

                if not table_colors:
                    raise ValueError(f"Color table {col_array[ctabletracker]} does not contain any colors")

                if len(table_colors) < len(tracks):
                    repeats = (len(tracks) // len(table_colors)) + 1
                    mycols = (table_colors * repeats)[:len(tracks)]
                else:
                    stride = max(1, (len(table_colors) - 1) // max(1, len(tracks)))
                    mycols = [table_colors[idx] for idx in range(0, len(table_colors), stride)][:len(tracks)]

                hub_name = f"{superhub}.{appendix}"
                hub_folder = os.path.join(output_folder, hub_name)
                genome_folder = os.path.join(hub_folder, genome)

                if os.path.isdir(hub_folder):
                    shutil.rmtree(hub_folder)

                ensure_hub_structure(hub_folder, genome_folder, hub_name, genome, hub_mail, overlay)

                for i, bw in enumerate(tracks):
                    bw_path = os.path.join(input_folder, bw)
                    fields = bw.split(separator)
                    name = separator.join(fields[idx - 1] for idx in field_indices if idx - 1 < len(fields))
                    color = mycols[i]
                    log_it(logfile, f"Processing {name}...")

                    with open(os.path.join(genome_folder, "trackDb.txt"), "a") as trackdb_file:
                        trackdb_file.write(
                            f"\ntrack {name}.{appendix}\n"
                            f"bigDataUrl {bw}\n"
                            f"shortLabel {name}.{appendix}\n"
                            f"longLabel {name}.{appendix}\n"
                            "type bigWig\n"
                            f"parent {hub_name}\n"
                            f"color {color}\n"
                        )

                    shutil.copyfile(bw_path, os.path.join(genome_folder, bw))

                ctabletracker += 1

            shell(f"""echo "necessity file for merge wiggle. can delete this." > {output_folder}/extra_9.tmp""")
            log_it(logfile, "Merge Wiggles and TrackHubs completed!")

        merge_wig(
            input_folder=params.inputfolder,
            output_folder=params.outputfolder,
            col_field=int(params.thecolfield),
            separator=params.theseparator,
            type_field=int(params.thetypefield),
            name_fields=params.thenamefields,
            col_table=params.coltable,
            appendix=params.theappendix,
            genome=params.thegenome,
            hub_mail=params.thehubmail,
            overlay=params.theoverlay,
        )
