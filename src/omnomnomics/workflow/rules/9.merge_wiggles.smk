# Rule 9 merge wiggles

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import shutil

RNA_PLUS_TRACK_COLOR = "255,0,0"
RNA_MINUS_TRACK_COLOR = "0,0,255"

rule merge_wiggles:
    input:
        extra_file=expand(f"{experiment_dir}/{master_config['input_folders'][master_config['mergewig_rule_num']-1]}/{{sample}}.extra_8.tmp", sample=samples2) if 8 in themode else []
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['mergewig_rule_num']-1]}/extra_9.tmp"
    params:
        thetype=lambda wildcards: config['THETYPE'],
        thegenome=lambda wildcards: config['THEGENOME'],
        coltable=lambda wildcards: config["THECOLTABLE"],
        theappendix=lambda wildcards: config['THEAPPENDIX'],
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
        if params.thetype == "RNA":
            log_once(
                logfile,
                "step9.rna_track_colors",
                f"RNA track colors are fixed to plus strand {RNA_PLUS_TRACK_COLOR} and minus strand {RNA_MINUS_TRACK_COLOR}.",
            )

        bw_files_for_step = []
        if os.path.isdir(params.inputfolder):
            bw_files_for_step = [bw_name for bw_name in os.listdir(params.inputfolder) if bw_name.endswith(".bw")]

        estimated_copy_bytes = sum(
            os.path.getsize(os.path.join(params.inputfolder, bw_name))
            for bw_name in bw_files_for_step
        )
        if evaluate_space_heavy_step(
            logfile,
            master_config['mergewig_rule_num'],
            estimated_extra_bytes=estimated_copy_bytes,
            delete_partial_hubs=True,
        ):
            os.makedirs(params.outputfolder, exist_ok=True)
            shell(f"""echo "necessity file for merge wiggle. can delete this." > {params.outputfolder}/extra_9.tmp""")
            log_it(logfile, "Skipping trackhub creation due to max project size constraint.")
            return

        if not bw_files_for_step:
            os.makedirs(params.outputfolder, exist_ok=True)
            log_it(
                logfile,
                (
                    f"No .bw files found in {params.inputfolder}. "
                    "Skipping trackhub creation for step 9."
                ),
                "WARNING",
            )
            shell(f"""echo "necessity file for merge wiggle. can delete this." > {params.outputfolder}/extra_9.tmp""")
            return

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

        def sample_root_from_bw(filename, thetype):
            if thetype == "RNA":
                if filename.endswith(".plus.bw"):
                    return filename[:-len(".plus.bw")]
                if filename.endswith(".minus.bw"):
                    return filename[:-len(".minus.bw")]
            if filename.endswith(".bw"):
                return filename[:-len(".bw")]
            return filename

        def append_track(trackdb_path, track_name, bigwig_name, short_label, long_label, parent_name, color, negate=False):
            with open(trackdb_path, "a") as trackdb_file:
                trackdb_file.write(
                    f"\ntrack {track_name}\n"
                    f"bigDataUrl {bigwig_name}\n"
                    f"shortLabel {short_label}\n"
                    f"longLabel {long_label}\n"
                    "type bigWig\n"
                    f"parent {parent_name}\n"
                    f"color {color}\n"
                    "alwaysZero on\n"
                )
                if negate:
                    trackdb_file.write("negateValues on\n")

        def merge_wig(input_folder, output_folder, thetype, col_table, appendix, genome, hub_mail, overlay):
            col_array = []
            if col_table.endswith(".txt"):
                log_it(logfile, f"Color table file: {col_table}")
                with open(col_table) as handle:
                    col_array = handle.read().splitlines()
            else:
                col_array.append(col_table)

            bw_files = [bw for bw in os.listdir(input_folder) if bw.endswith(".bw")]
            sample_roots = sorted(set(sample_root_from_bw(bw, thetype) for bw in bw_files))
            hubtypes = sorted(set(sample_type_for_sample(sample_root) for sample_root in sample_roots))

            ctabletracker = 0
            for superhub in hubtypes:
                sample_tracks = [
                    sample_root for sample_root in sample_roots
                    if sample_type_for_sample(sample_root) == superhub
                ]
                sample_tracks.sort(key=lambda sample_root: sample_id_for_sample(sample_root))
                if not sample_tracks:
                    continue

                color_categories = []
                for sample_root in sample_tracks:
                    color_category = sample_color_for_sample(sample_root)
                    if color_category not in color_categories:
                        color_categories.append(color_category)

                if ctabletracker == len(col_array):
                    ctabletracker = 0

                with open(col_array[ctabletracker]) as handle:
                    table_colors = [line.strip() for line in handle if line.strip()]

                if not table_colors:
                    raise ValueError(f"Color table {col_array[ctabletracker]} does not contain any colors")

                if len(table_colors) < len(color_categories):
                    repeats = (len(color_categories) // len(table_colors)) + 1
                    category_palette = (table_colors * repeats)[:len(color_categories)]
                else:
                    stride = max(1, (len(table_colors) - 1) // max(1, len(color_categories)))
                    category_palette = [table_colors[idx] for idx in range(0, len(table_colors), stride)][:len(color_categories)]
                category_to_color = dict(zip(color_categories, category_palette))

                hub_name = f"{superhub}.{appendix}"
                hub_folder = os.path.join(output_folder, hub_name)
                genome_folder = os.path.join(hub_folder, genome)

                if os.path.isdir(hub_folder):
                    shutil.rmtree(hub_folder)

                ensure_hub_structure(hub_folder, genome_folder, hub_name, genome, hub_mail, overlay)

                for sample_root in sample_tracks:
                    name = sample_id_for_sample(sample_root)
                    color_category = sample_color_for_sample(sample_root)
                    color = category_to_color[color_category]
                    log_it(logfile, f"Processing {name}...")
                    trackdb_path = os.path.join(genome_folder, "trackDb.txt")

                    if thetype == "RNA":
                        plus_bw = f"{sample_root}.plus.bw"
                        minus_bw = f"{sample_root}.minus.bw"
                        plus_bw_path = os.path.join(input_folder, plus_bw)
                        minus_bw_path = os.path.join(input_folder, minus_bw)
                        if not os.path.exists(plus_bw_path) or not os.path.exists(minus_bw_path):
                            raise FileNotFoundError(
                                f"Expected stranded RNA BigWigs {plus_bw} and {minus_bw} in {input_folder}"
                            )

                        append_track(
                            trackdb_path,
                            f"{name}.{appendix}.plus",
                            plus_bw,
                            f"{name}+",
                            f"{name} plus strand",
                            hub_name,
                            RNA_PLUS_TRACK_COLOR,
                        )
                        append_track(
                            trackdb_path,
                            f"{name}.{appendix}.minus",
                            minus_bw,
                            f"{name}-",
                            f"{name} minus strand",
                            hub_name,
                            RNA_MINUS_TRACK_COLOR,
                            negate=True,
                        )

                        shutil.copyfile(plus_bw_path, os.path.join(genome_folder, plus_bw))
                        shutil.copyfile(minus_bw_path, os.path.join(genome_folder, minus_bw))
                    else:
                        bw = f"{sample_root}.bw"
                        bw_path = os.path.join(input_folder, bw)
                        append_track(
                            trackdb_path,
                            f"{name}.{appendix}",
                            bw,
                            f"{name}.{appendix}",
                            f"{name}.{appendix}",
                            hub_name,
                            color,
                        )
                        shutil.copyfile(bw_path, os.path.join(genome_folder, bw))

                ctabletracker += 1

            shell(f"""echo "necessity file for merge wiggle. can delete this." > {output_folder}/extra_9.tmp""")
            log_it(logfile, "Merge Wiggles and TrackHubs completed!")

        merge_wig(
            input_folder=params.inputfolder,
            output_folder=params.outputfolder,
            thetype=params.thetype,
            col_table=params.coltable,
            appendix=params.theappendix,
            genome=params.thegenome,
            hub_mail=params.thehubmail,
            overlay=params.theoverlay,
        )
