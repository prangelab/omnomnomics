# Rule 14: Analyze Peaks (pre-DE, chromatin-focused)

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import csv
import glob
import itertools
import os
import shlex
import shutil
import subprocess
import tempfile


def analyze_peaks_input(_wildcards):
    input_files = []
    if config["THETYPE"] in {"ATAC", "CHIP"}:
        input_files.append(
            f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/extra_{master_config['peakqc_rule_num']}.tmp"
        )
        if config["THETYPE"] == "ATAC":
            input_files.append(
                f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/peak_qc/peak_annotations/atac.all_groups.merged_peaks.annotated.bed"
            )
            input_files.extend(
                glob.glob(
                    f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/peak_qc/filtered_peaks/*.bed"
                )
            )
            input_files.extend(
                f"{experiment_dir}/{master_config['input_folders'][master_config['analyzepeaks_rule_num'] - 1][0]}/{sample}.sorted.dups_marked.filtered.bam"
                for sample in samples2
            )
        else:
            input_files.append(
                f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/peak_qc/peak_annotations/chip.all_groups.merged_peaks.annotated.bed"
            )
            input_files.extend(
                glob.glob(
                    f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num'] - 1]}/*.bed"
                )
            )
            input_files.extend(
                f"{experiment_dir}/{master_config['input_folders'][master_config['analyzepeaks_rule_num'] - 1][0]}/{sample}.filtered.bam"
                for sample in samples2
            )
    return input_files


rule analyze_peaks:
    input:
        analyze_peaks_input
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['analyzepeaks_rule_num'] - 1]}/extra_{master_config['analyzepeaks_rule_num']}.tmp"
    params:
        thetype=config["THETYPE"],
        broad_mode=str(config.get("BROAD_MODE", "off")).strip().lower(),
        bam_inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['analyzepeaks_rule_num'] - 1][0]}",
        peak_inputfolder=(
            f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/peak_qc/filtered_peaks"
            if config["THETYPE"] == "ATAC"
            else f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num'] - 1]}"
        ),
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['analyzepeaks_rule_num'] - 1]}",
    threads:
        Threads_Per_Rule[str(master_config["analyzepeaks_rule_num"])]
    resources:
        mem_mb=Memory_Per_Rule[str(master_config["analyzepeaks_rule_num"])],
        partition=master_config["partition"],
        runtime=Runtime_Per_Rule[str(master_config["analyzepeaks_rule_num"])]
    run:
        tracking = begin_step_sample(master_config["analyzepeaks_rule_num"], "aggregate", "analyze_peaks")
        log_once(
            logfile,
            "step14.header",
            "Analyzing chromatin peak sets (unions, intersections, signal heatmaps/profiles)...",
            f"EXECUTING STEP {master_config['analyzepeaks_rule_num']}",
        )
        log_once(logfile, "step14.inputfolder", f"Input folders: {params.bam_inputfolder} and {params.peak_inputfolder}")
        log_once(logfile, "step14.outputfolder", f"Output folder: {params.outputfolder}")

        def quote(path):
            return shlex.quote(path)

        def write_tmp_file(outputfolder):
            shell(
                f"""echo "necessity file for analyze peaks. can delete this." > {quote(os.path.join(outputfolder, f"extra_{master_config['analyzepeaks_rule_num']}.tmp"))}"""
            )

        def require_tool(executable_name):
            if shutil.which(executable_name) is None:
                raise RuntimeError(
                    f"Required executable '{executable_name}' not found on PATH for analyze_peaks step."
                )

        def ensure_dir(path):
            os.makedirs(path, exist_ok=True)
            return path

        def load_annotation_rows(annotation_bed):
            with open(annotation_bed, newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames is None:
                    raise ValueError(f"Annotation BED has no header: {annotation_bed}")
                return list(reader)

        def build_group_peak_sets(peak_inputfolder, sets_dir):
            peak_files = sorted(glob.glob(os.path.join(peak_inputfolder, "*.bed")))
            grouped = {}
            for path in peak_files:
                base = os.path.basename(path)
                if base == "all_groups.merged_peaks.bed":
                    continue
                if ".MACS3." in base:
                    group = base.split(".MACS3.", 1)[0]
                else:
                    group = os.path.splitext(base)[0]
                grouped.setdefault(group, []).append(path)

            union_paths = {}
            for group, files in grouped.items():
                union_path = os.path.join(sets_dir, f"{group}.union.bed")
                shell(
                    f"cat {' '.join(quote(p) for p in files)} | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {quote(union_path)}"
                )
                union_paths[group] = union_path
            return union_paths

        def write_set_manifest(union_paths, manifest_path):
            with open(manifest_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["set_name", "set_type", "peak_bed", "peak_count"])
                for group, bed_path in sorted(union_paths.items()):
                    peak_count = int(
                        subprocess.check_output(
                            f"wc -l {quote(bed_path)} | awk '{{print $1}}'",
                            shell=True,
                            executable="/bin/bash",
                            text=True,
                        ).strip()
                    )
                    writer.writerow([group, "group_union", bed_path, peak_count])

        def write_intersections(union_paths, intersections_dir, summary_path):
            ensure_dir(intersections_dir)
            unique_dir = ensure_dir(os.path.join(intersections_dir, "unique"))
            shared_dir = ensure_dir(os.path.join(intersections_dir, "shared"))
            matrix_rows = []
            groups = sorted(union_paths)

            for group_a, group_b in itertools.combinations(groups, 2):
                bed_a = union_paths[group_a]
                bed_b = union_paths[group_b]
                pair_name = f"{group_a}__{group_b}"
                shared_bed = os.path.join(shared_dir, f"{pair_name}.shared.bed")
                unique_a_bed = os.path.join(unique_dir, f"{group_a}.vs.{group_b}.unique.bed")
                unique_b_bed = os.path.join(unique_dir, f"{group_b}.vs.{group_a}.unique.bed")
                shell(f"bedtools intersect -a {quote(bed_a)} -b {quote(bed_b)} -u > {quote(shared_bed)}")
                shell(f"bedtools intersect -a {quote(bed_a)} -b {quote(bed_b)} -v > {quote(unique_a_bed)}")
                shell(f"bedtools intersect -a {quote(bed_b)} -b {quote(bed_a)} -v > {quote(unique_b_bed)}")

                shared_count = int(subprocess.check_output(f"wc -l {quote(shared_bed)} | awk '{{print $1}}'", shell=True, executable="/bin/bash", text=True).strip())
                unique_a_count = int(subprocess.check_output(f"wc -l {quote(unique_a_bed)} | awk '{{print $1}}'", shell=True, executable="/bin/bash", text=True).strip())
                unique_b_count = int(subprocess.check_output(f"wc -l {quote(unique_b_bed)} | awk '{{print $1}}'", shell=True, executable="/bin/bash", text=True).strip())

                jaccard_output = subprocess.check_output(
                    f"bedtools jaccard -a {quote(bed_a)} -b {quote(bed_b)}",
                    shell=True,
                    executable="/bin/bash",
                    text=True,
                ).strip().splitlines()
                jaccard = "NA"
                if len(jaccard_output) > 1:
                    jaccard_fields = jaccard_output[1].split("\t")
                    if len(jaccard_fields) >= 3:
                        jaccard = jaccard_fields[2]
                matrix_rows.append(
                    [group_a, group_b, shared_count, unique_a_count, unique_b_count, jaccard]
                )

            with open(summary_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["group_a", "group_b", "shared_peak_count", "unique_a_peak_count", "unique_b_peak_count", "jaccard"])
                writer.writerows(matrix_rows)

        def run_signal_plots(bam_inputfolder, union_paths, signal_dir, assay_type):
            require_tool("bamCoverage")
            require_tool("computeMatrix")
            require_tool("plotHeatmap")
            require_tool("plotProfile")
            matrices_dir = ensure_dir(os.path.join(signal_dir, "matrices"))
            heatmaps_dir = ensure_dir(os.path.join(signal_dir, "heatmaps"))
            profiles_dir = ensure_dir(os.path.join(signal_dir, "profiles"))

            with tempfile.TemporaryDirectory(prefix="omnomnomics_deeptools_bw_") as tmpdir:
                bigwigs = []
                sample_labels = []
                bam_suffix = ".sorted.dups_marked.filtered.bam" if assay_type == "ATAC" else ".filtered.bam"
                for sample in samples2:
                    bam_path = os.path.join(bam_inputfolder, f"{sample}{bam_suffix}")
                    if not os.path.exists(bam_path):
                        continue
                    bw_path = os.path.join(tmpdir, f"{sample}.bw")
                    shell(
                        f"bamCoverage --bam {quote(bam_path)} --outFileName {quote(bw_path)} --binSize 25 --normalizeUsing CPM --numberOfProcessors {threads}"
                    )
                    bigwigs.append(bw_path)
                    sample_labels.append(sample)

                if not bigwigs:
                    log_it(logfile, "No BAM files found for deepTools signal plotting. Skipping signal plots.")
                    return

                for group, region_bed in sorted(union_paths.items()):
                    group_label = group.replace("_", " ")
                    matrix_path = os.path.join(matrices_dir, f"{group}.matrix.gz")
                    heatmap_path = os.path.join(heatmaps_dir, f"{group}.heatmap.pdf")
                    profile_path = os.path.join(profiles_dir, f"{group}.profile.pdf")
                    shell(
                        f"computeMatrix reference-point --referencePoint center -b 3000 -a 3000 "
                        f"-R {quote(region_bed)} -S {' '.join(quote(x) for x in bigwigs)} "
                        f"--skipZeros --binSize 50 --numberOfProcessors {threads} -o {quote(matrix_path)}"
                    )
                    shell(
                        f"plotHeatmap -m {quote(matrix_path)} -out {quote(heatmap_path)} "
                        f"--whatToShow 'heatmap and colorbar' --sortRegions descend "
                        f"--plotTitle {quote(group_label)} --regionsLabel {quote(group_label)} "
                        f"--samplesLabel {' '.join(quote(label) for label in sample_labels)} "
                        f"--xAxisLabel 'distance from center (bp)' --refPointLabel center"
                    )
                    shell(
                        f"plotProfile -m {quote(matrix_path)} -out {quote(profile_path)} "
                        f"--perGroup --plotTitle {quote(group_label)} "
                        f"--regionsLabel {quote(group_label)} "
                        f"--samplesLabel {' '.join(quote(label) for label in sample_labels)} "
                        f"--refPointLabel center"
                    )

        def build_peak_metadata_for_r(annotation_bed, output_path):
            expected = {
                "#chrom",
                "start",
                "end",
                "underscore",
                "genomic_region",
                "assigned_genes",
                "nearest_gene",
                "distance_to_nearest_gene_bp",
                "width_bp",
            }
            with open(annotation_bed, newline="") as inp:
                reader = csv.DictReader(inp, delimiter="\t")
                if reader.fieldnames is None:
                    raise ValueError(f"Peak annotation BED has no header: {annotation_bed}")
                missing = [field for field in expected if field not in reader.fieldnames]
                if missing:
                    raise ValueError(
                        "Peak annotation BED is missing required columns: "
                        + ", ".join(missing)
                    )
                rows = list(reader)

            with open(output_path, "w", newline="") as out:
                writer = csv.writer(out, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    [
                        "underscore",
                        "chrom",
                        "start",
                        "end",
                        "genomic_region",
                        "assigned_genes",
                        "nearest_gene",
                        "distance_to_nearest_gene_bp",
                        "width_bp",
                    ]
                )
                for row in rows:
                    writer.writerow(
                        [
                            row["underscore"],
                            row["#chrom"],
                            row["start"],
                            row["end"],
                            row["genomic_region"],
                            row["assigned_genes"],
                            row["nearest_gene"],
                            row["distance_to_nearest_gene_bp"],
                            row["width_bp"],
                        ]
                    )

        def write_genebody_feature_metadata(annotation_bed, summary_dir):
            metadata_path = os.path.join(summary_dir, "feature_metadata_for_r.tsv")
            build_peak_metadata_for_r(annotation_bed, metadata_path)
            legacy_metadata_path = os.path.join(summary_dir, "peak_metadata_for_r.tsv")
            shutil.copy2(metadata_path, legacy_metadata_path)
            return metadata_path, legacy_metadata_path

        def write_genebody_summary(annotation_bed, feature_bed, summary_dir):
            rows = load_annotation_rows(annotation_bed)
            summary_path = os.path.join(summary_dir, "gene_body_summary.tsv")
            region_counts = {}
            widths = []
            assigned_gene_rows = 0
            nearest_gene_rows = 0
            for row in rows:
                region = str(row.get("genomic_region", "NA")).strip() or "NA"
                region_counts[region] = region_counts.get(region, 0) + 1
                try:
                    widths.append(int(row.get("width_bp", "0")))
                except ValueError:
                    pass
                if str(row.get("assigned_genes", "NA")).strip() not in {"", "NA"}:
                    assigned_gene_rows += 1
                if str(row.get("nearest_gene", "NA")).strip() not in {"", "NA"}:
                    nearest_gene_rows += 1
            feature_count = int(
                subprocess.check_output(
                    f"wc -l {quote(feature_bed)} | awk '{{print $1}}'",
                    shell=True,
                    executable="/bin/bash",
                    text=True,
                ).strip()
            )
            widths_sorted = sorted(widths)
            median_width = widths_sorted[len(widths_sorted) // 2] if widths_sorted else 0
            with open(summary_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["metric", "value"])
                writer.writerow(["feature_type", "gene_body"])
                writer.writerow(["feature_count", feature_count])
                writer.writerow(["median_width_bp", median_width])
                writer.writerow(["assigned_gene_rows", assigned_gene_rows])
                writer.writerow(["nearest_gene_rows", nearest_gene_rows])
                for region_name in sorted(region_counts):
                    writer.writerow([f"region_count_{region_name}", region_counts[region_name]])
            return summary_path

        def write_feature_summary(annotation_bed, feature_bed, summary_dir, feature_type, summary_name):
            rows = load_annotation_rows(annotation_bed)
            summary_path = os.path.join(summary_dir, summary_name)
            region_counts = {}
            widths = []
            assigned_gene_rows = 0
            nearest_gene_rows = 0
            for row in rows:
                region = str(row.get("genomic_region", "NA")).strip() or "NA"
                region_counts[region] = region_counts.get(region, 0) + 1
                try:
                    widths.append(int(row.get("width_bp", "0")))
                except ValueError:
                    pass
                if str(row.get("assigned_genes", "NA")).strip() not in {"", "NA"}:
                    assigned_gene_rows += 1
                if str(row.get("nearest_gene", "NA")).strip() not in {"", "NA"}:
                    nearest_gene_rows += 1
            feature_count = int(
                subprocess.check_output(
                    f"wc -l {quote(feature_bed)} | awk '{{print $1}}'",
                    shell=True,
                    executable="/bin/bash",
                    text=True,
                ).strip()
            )
            widths_sorted = sorted(widths)
            median_width = widths_sorted[len(widths_sorted) // 2] if widths_sorted else 0
            with open(summary_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["metric", "value"])
                writer.writerow(["feature_type", feature_type])
                writer.writerow(["feature_count", feature_count])
                writer.writerow(["median_width_bp", median_width])
                writer.writerow(["assigned_gene_rows", assigned_gene_rows])
                writer.writerow(["nearest_gene_rows", nearest_gene_rows])
                for region_name in sorted(region_counts):
                    writer.writerow([f"region_count_{region_name}", region_counts[region_name]])
            return summary_path

        def run_genebody_signal_plots(bam_inputfolder, feature_bed, signal_dir):
            require_tool("bamCoverage")
            require_tool("computeMatrix")
            require_tool("plotHeatmap")
            require_tool("plotProfile")
            matrices_dir = ensure_dir(os.path.join(signal_dir, "matrices"))
            heatmaps_dir = ensure_dir(os.path.join(signal_dir, "heatmaps"))
            profiles_dir = ensure_dir(os.path.join(signal_dir, "profiles"))

            with tempfile.TemporaryDirectory(prefix="omnomnomics_genebody_bw_") as tmpdir:
                bigwigs = []
                sample_labels = []
                for sample in samples2:
                    bam_path = os.path.join(bam_inputfolder, f"{sample}.filtered.bam")
                    if not os.path.exists(bam_path):
                        continue
                    bw_path = os.path.join(tmpdir, f"{sample}.bw")
                    shell(
                        f"bamCoverage --bam {quote(bam_path)} --outFileName {quote(bw_path)} "
                        f"--binSize 25 --normalizeUsing CPM --numberOfProcessors {threads}"
                    )
                    bigwigs.append(bw_path)
                    sample_labels.append(sample)

                if not bigwigs:
                    log_it(logfile, "No BAM files found for gene-body signal plotting. Skipping pre-DE gene-body summaries.")
                    return None

                matrix_path = os.path.join(matrices_dir, "all_gene_bodies.scale_regions.matrix.gz")
                heatmap_path = os.path.join(heatmaps_dir, "all_gene_bodies.scale_regions.heatmap.pdf")
                profile_path = os.path.join(profiles_dir, "all_gene_bodies.scale_regions.profile.pdf")
                shell(
                    f"computeMatrix scale-regions "
                    f"-R {quote(feature_bed)} -S {' '.join(quote(x) for x in bigwigs)} "
                    f"--beforeRegionStartLength 1000 --afterRegionStartLength 1000 "
                    f"--regionBodyLength 5000 --skipZeros --binSize 50 "
                    f"--numberOfProcessors {threads} -o {quote(matrix_path)}"
                )
                shell(
                    f"plotHeatmap -m {quote(matrix_path)} -out {quote(heatmap_path)} "
                    f"--whatToShow 'heatmap and colorbar' --sortRegions descend "
                    f"--plotTitle {quote('All gene bodies')} --regionsLabel {quote('gene bodies')} "
                    f"--samplesLabel {' '.join(quote(label) for label in sample_labels)} "
                    f"--startLabel TSS --endLabel TES --xAxisLabel {quote('gene body scaled to 5 kb')} "
                    f"--heatmapWidth 8 --heatmapHeight 14"
                )
                shell(
                    f"plotProfile -m {quote(matrix_path)} -out {quote(profile_path)} "
                    f"--perGroup --plotTitle {quote('All gene bodies')} "
                    f"--regionsLabel {quote('gene bodies')} "
                    f"--samplesLabel {' '.join(quote(label) for label in sample_labels)} "
                    f"--startLabel TSS --endLabel TES"
                )
                return {
                    "matrix": matrix_path,
                    "heatmap": heatmap_path,
                    "profile": profile_path,
                }

        try:
            if params.thetype not in {"ATAC", "CHIP"}:
                log_it(logfile, "Analyze peaks step is supported for ATAC/CHIP only. Skipping.")
                write_tmp_file(params.outputfolder)
                finish_step_sample(master_config["analyzepeaks_rule_num"], "aggregate", "analyze_peaks", tracking["start_time"], "OK")
                return
            if params.thetype == "CHIP" and params.broad_mode == "genebody":
                require_tool("bedtools")
                summary_dir = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks", "summary"))
                signal_dir = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks", "signal"))
                feature_bed = os.path.join(params.peak_inputfolder, "all_groups.merged_peaks.bed")
                annotation_bed = (
                    f"{experiment_dir}/"
                    f"{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}"
                    f"/peak_qc/peak_annotations/chip.all_groups.merged_peaks.annotated.bed"
                )
                if not os.path.isfile(feature_bed):
                    raise FileNotFoundError(f"Expected gene-body feature BED not found: {feature_bed}")
                if not os.path.isfile(annotation_bed):
                    raise FileNotFoundError(f"Expected gene-body annotation BED not found: {annotation_bed}")
                metadata_path, legacy_metadata_path = write_genebody_feature_metadata(annotation_bed, summary_dir)
                gene_body_summary = write_genebody_summary(annotation_bed, feature_bed, summary_dir)
                signal_outputs = run_genebody_signal_plots(params.bam_inputfolder, feature_bed, signal_dir)
                report_path = os.path.join(summary_dir, "analyze_peaks_report.tsv")
                with open(report_path, "w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["metric", "value"])
                    writer.writerow(["analysis_mode", "genebody"])
                    writer.writerow(["feature_bed", feature_bed])
                    writer.writerow(["annotation_bed", annotation_bed])
                    writer.writerow(["feature_metadata_for_r", metadata_path])
                    writer.writerow(["legacy_peak_metadata_for_r", legacy_metadata_path])
                    writer.writerow(["gene_body_summary", gene_body_summary])
                    if signal_outputs is not None:
                        writer.writerow(["signal_matrix", signal_outputs["matrix"]])
                        writer.writerow(["signal_heatmap", signal_outputs["heatmap"]])
                        writer.writerow(["signal_profile", signal_outputs["profile"]])
                log_it(logfile, f"ChIP gene-body pre-DE summary: {gene_body_summary}")
                write_tmp_file(params.outputfolder)
                finish_step_sample(master_config["analyzepeaks_rule_num"], "aggregate", "analyze_peaks", tracking["start_time"], "OK")
                return
            if params.thetype == "CHIP" and params.broad_mode == "diffuse":
                require_tool("bedtools")
                summary_dir = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks", "summary"))
                signal_dir = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks", "signal"))
                feature_bed = os.path.join(params.peak_inputfolder, "all_groups.merged_peaks.bed")
                annotation_bed = (
                    f"{experiment_dir}/"
                    f"{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}"
                    f"/peak_qc/peak_annotations/chip.all_groups.merged_peaks.annotated.bed"
                )
                if not os.path.isfile(feature_bed):
                    raise FileNotFoundError(f"Expected diffuse bin BED not found: {feature_bed}")
                if not os.path.isfile(annotation_bed):
                    raise FileNotFoundError(f"Expected diffuse bin annotation BED not found: {annotation_bed}")
                metadata_path = os.path.join(summary_dir, "feature_metadata_for_r.tsv")
                build_peak_metadata_for_r(annotation_bed, metadata_path)
                legacy_metadata_path = os.path.join(summary_dir, "peak_metadata_for_r.tsv")
                shutil.copy2(metadata_path, legacy_metadata_path)
                diffuse_summary = write_feature_summary(
                    annotation_bed,
                    feature_bed,
                    summary_dir,
                    "diffuse_bin",
                    "diffuse_bin_summary.tsv",
                )
                run_signal_plots(
                    params.bam_inputfolder,
                    {"all_bins": feature_bed},
                    signal_dir,
                    params.thetype,
                )
                report_path = os.path.join(summary_dir, "analyze_peaks_report.tsv")
                with open(report_path, "w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["metric", "value"])
                    writer.writerow(["analysis_mode", "diffuse"])
                    writer.writerow(["feature_bed", feature_bed])
                    writer.writerow(["annotation_bed", annotation_bed])
                    writer.writerow(["feature_metadata_for_r", metadata_path])
                    writer.writerow(["legacy_peak_metadata_for_r", legacy_metadata_path])
                    writer.writerow(["diffuse_bin_summary", diffuse_summary])
                log_it(logfile, f"ChIP diffuse pre-DE summary: {diffuse_summary}")
                write_tmp_file(params.outputfolder)
                finish_step_sample(master_config["analyzepeaks_rule_num"], "aggregate", "analyze_peaks", tracking["start_time"], "OK")
                return

            require_tool("bedtools")
            summary_dir = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks", "summary"))
            sets_dir = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks", "sets", "group_unions"))
            intersections_dir = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks", "intersections"))
            signal_dir = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks", "signal"))

            union_paths = build_group_peak_sets(params.peak_inputfolder, sets_dir)
            if len(union_paths) < 2:
                log_it(logfile, "Analyze peaks requires at least two biological groups with peak BED files. Skipping intersections and signal plots.")
                write_tmp_file(params.outputfolder)
                finish_step_sample(master_config["analyzepeaks_rule_num"], "aggregate", "analyze_peaks", tracking["start_time"], "OK")
                return

            set_manifest_path = os.path.join(summary_dir, "set_manifest.tsv")
            write_set_manifest(union_paths, set_manifest_path)
            overlap_summary_path = os.path.join(summary_dir, "overlap_matrix.tsv")
            write_intersections(union_paths, intersections_dir, overlap_summary_path)
            run_signal_plots(params.bam_inputfolder, union_paths, signal_dir, params.thetype)
            peak_annotation_bed = (
                f"{experiment_dir}/"
                f"{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}"
                f"/peak_qc/peak_annotations/{params.thetype.lower()}.all_groups.merged_peaks.annotated.bed"
            )
            peak_metadata_tsv = os.path.join(summary_dir, "peak_metadata_for_r.tsv")
            build_peak_metadata_for_r(peak_annotation_bed, peak_metadata_tsv)

            report_path = os.path.join(summary_dir, "analyze_peaks_report.tsv")
            with open(report_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["metric", "value"])
                writer.writerow(["group_count", len(union_paths)])
                writer.writerow(["set_manifest", set_manifest_path])
                writer.writerow(["overlap_matrix", overlap_summary_path])
                writer.writerow(["peak_annotation_bed", peak_annotation_bed])
                writer.writerow(["peak_metadata_for_r", peak_metadata_tsv])

            write_tmp_file(params.outputfolder)
            finish_step_sample(master_config["analyzepeaks_rule_num"], "aggregate", "analyze_peaks", tracking["start_time"], "OK")
        except Exception:
            finish_step_sample(master_config["analyzepeaks_rule_num"], "aggregate", "analyze_peaks", tracking["start_time"], "FAIL")
            raise
