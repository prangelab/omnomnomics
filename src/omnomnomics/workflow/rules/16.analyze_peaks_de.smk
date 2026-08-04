# Rule 16: Analyze Peaks DE (post-DE, chromatin-focused)

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import csv
import glob
import gzip
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess

from omnomnomics.peak_annotation import promoter_intervals_by_gene


def analyze_peaks_de_input(_wildcards):
    input_files = []
    if config["THETYPE"] in {"ATAC", "CHIP"}:
        input_files.append(
            ancient(
                f"{experiment_dir}/{master_config['output_folders'][master_config['dechrom_rule_num'] - 1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.chrom.results.zip"
            )
        )
        input_files.append(
            ancient(f"{experiment_dir}/{master_config['output_folders'][master_config['dechrom_rule_num'] - 1]}/peak_metadata.tsv")
        )
        bam_inputfolder = master_config["output_folders"][master_config["touchup_rule_num"] - 1]
        bam_suffix = ".sorted.dups_marked.filtered.bam" if config["THETYPE"] == "ATAC" else ".filtered.bam"
        input_files.extend(
            f"{experiment_dir}/{bam_inputfolder}/{sample}{bam_suffix}"
            for sample in samples2
        )
        post_de_signal_policy = str(config.get("POST_DE_SIGNAL_POLICY", "auto")).strip().lower()
        if post_de_signal_policy == "auto" and int(config.get("MAX_PROJECT_SIZE_BYTES", 0) or 0) <= 0:
            bigwig_inputfolder = master_config["output_folders"][master_config["wig_rule_num"] - 1]
            input_files.extend(
                f"{experiment_dir}/{bigwig_inputfolder}/{sample}.bw"
                for sample in samples2
            )
        if config["THETYPE"] == "CHIP" and str(config.get("BROAD_MODE", "off")).strip().lower() == "genebody":
            input_files.append(
                os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "annotation", "genes.gtf")
            )
    return input_files


rule analyze_peaks_de:
    input:
        analyze_peaks_de_input
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['analyzepeaksde_rule_num'] - 1]}/extra_{master_config['analyzepeaksde_rule_num']}.tmp"
    params:
        thetype=config["THETYPE"],
        broad_mode=str(config.get("BROAD_MODE", "off")).strip().lower(),
        diffuse_merge_gap=int(config.get("CHIP_DIFFUSE_MERGE_GAP", 10000)),
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['analyzepeaksde_rule_num'] - 1]}",
        de_outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['dechrom_rule_num'] - 1]}",
        prede_outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['analyzepeaks_rule_num'] - 1]}",
        bam_inputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num'] - 1]}",
        bigwig_inputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num'] - 1]}",
        post_de_signal_policy=str(config.get("POST_DE_SIGNAL_POLICY", "auto")).strip().lower(),
        post_de_signal_auto_dependencies=int(config.get("MAX_PROJECT_SIZE_BYTES", 0) or 0) <= 0,
        post_de_motif_max_sets=int(master_config.get("post_de_motif_max_sets", 6) or 0),
        post_de_motif_max_peaks=int(master_config.get("post_de_motif_max_peaks", 100) or 100),
        post_de_motif_timeout_seconds=int(master_config.get("post_de_motif_timeout_seconds", 1200) or 1200),
        post_de_motif_threads=int(master_config.get("post_de_motif_threads", 1) or 1),
        post_de_motif_window_bp=int(master_config.get("post_de_motif_window_bp", 200) or 200),
        post_de_motif_database=str(master_config.get("post_de_motif_database", "auto")),
        genome_version=config["THEGENOME"],
        genome_assembly_dir=config["GENOME_ASSEMBLY_DIR"],
        genome_fasta=os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "fasta", "genome.fa"),
        gtf_file=os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "annotation", "genes.gtf"),
    threads:
        Threads_Per_Rule[str(master_config["analyzepeaksde_rule_num"])]
    resources:
        mem_mb=Memory_Per_Rule[str(master_config["analyzepeaksde_rule_num"])],
        partition=master_config["partition"],
        runtime=Runtime_Per_Rule[str(master_config["analyzepeaksde_rule_num"])]
    run:
        tracking = begin_step_sample(master_config["analyzepeaksde_rule_num"], "aggregate", "analyze_peaks_de")
        log_once(
            logfile,
            "step16.header",
            "Analyzing post-DE chromatin peak sets (DE subsets, optional unique sets, signal, motifs)...",
            f"EXECUTING STEP {master_config['analyzepeaksde_rule_num']}",
        )
        log_once(logfile, "step16.inputfolder", f"Input DE folder: {params.de_outputfolder}")
        log_once(logfile, "step16.outputfolder", f"Output folder: {params.outputfolder}")

        def quote(path):
            return shlex.quote(path)

        def ensure_dir(path):
            os.makedirs(path, exist_ok=True)
            return path

        def write_tmp_file():
            tmp_path = os.path.join(
                params.outputfolder,
                f"extra_{master_config['analyzepeaksde_rule_num']}.tmp",
            )
            with open(tmp_path, "w") as handle:
                handle.write("necessity file for analyze_peaks_de. can delete this.\n")

        def tool_available(executable_name):
            return shutil.which(executable_name) is not None

        def parse_peak_coord(gene_id):
            text = str(gene_id)
            match = re.match(r"^([^:]+):([0-9]+)-([0-9]+)$", text)
            if match:
                chrom, start_text, end_text = match.groups()
            else:
                parts = text.rsplit("_", 2)
                if len(parts) != 3:
                    return None
                chrom, start_text, end_text = parts
                if not start_text.isdigit() or not end_text.isdigit():
                    return None
            start = int(start_text)
            end = int(end_text)
            if end <= start:
                return None
            return chrom, start, end

        def load_peak_region_map(peak_metadata_path):
            region_map = {}
            if not os.path.isfile(peak_metadata_path):
                return region_map
            with open(peak_metadata_path, newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    peak_id = str(row.get("underscore", "")).strip()
                    if not peak_id:
                        continue
                    region_map[peak_id] = {
                        "genomic_region": str(row.get("genomic_region", "NA")).strip() or "NA",
                        "assigned_genes": str(row.get("assigned_genes", "NA")).strip() or "NA",
                        "nearest_gene": str(row.get("nearest_gene", "NA")).strip() or "NA",
                        "nearest_promoter_gene": str(
                            row.get("nearest_promoter_gene", row.get("nearest_gene", "NA"))
                        ).strip() or "NA",
                        "distance_to_nearest_promoter_bp": str(
                            row.get(
                                "distance_to_nearest_promoter_bp",
                                row.get("distance_to_nearest_gene_bp", "NA"),
                            )
                        ).strip() or "NA",
                        "chrom": str(row.get("chrom", row.get("#chrom", ""))).strip(),
                        "start": str(row.get("start", "")).strip(),
                        "end": str(row.get("end", "")).strip(),
                    }
            return region_map

        def write_bed_from_rows(rows, bed_path):
            with open(bed_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                for row in rows:
                    writer.writerow(row)

        def load_promoter_map(gtf_file, promoter_upstream=1000, promoter_downstream=1000):
            if not os.path.isfile(gtf_file):
                return {}
            return promoter_intervals_by_gene(
                gtf_file,
                promoter_upstream=promoter_upstream,
                promoter_downstream=promoter_downstream,
            )

        def feature_rows_from_peak_metadata(peak_region_map):
            rows = []
            for peak_id, info in sorted(peak_region_map.items()):
                try:
                    chrom = str(info.get("chrom", "")).strip()
                    start = int(info.get("start", ""))
                    end = int(info.get("end", ""))
                except (TypeError, ValueError):
                    parsed = parse_peak_coord(peak_id)
                    if parsed is None:
                        continue
                    chrom, start, end = parsed
                region = info.get("genomic_region", "NA")
                assigned = info.get("assigned_genes", "NA")
                nearest = info.get("nearest_gene", "NA")
                rows.append([chrom, start, end, peak_id, "0.000000", ".", region, assigned, nearest])
            return rows

        def labels_for_feature(region_info):
            assigned = str(region_info.get("assigned_genes", "NA")).strip()
            nearest = str(region_info.get("nearest_gene", "NA")).strip()
            labels = []
            if assigned and assigned != "NA":
                labels.extend([item.strip() for item in assigned.split(",") if item.strip() and item.strip() != "NA"])
            if nearest and nearest != "NA":
                labels.append(nearest)
            return sorted(set(labels))

        def collect_de_tables():
            file_stem = "peaks"
            if params.thetype == "CHIP" and params.broad_mode == "genebody":
                file_stem = "gene_bodies"
            elif params.thetype == "CHIP" and params.broad_mode == "diffuse":
                file_stem = "bins"
            candidate_paths = sorted(
                glob.glob(os.path.join(params.de_outputfolder, "*", "*", f"*.sig_diff_{file_stem}.DESeq2.txt"))
            )
            if candidate_paths:
                return candidate_paths
            return sorted(
                glob.glob(os.path.join(params.de_outputfolder, "*", "*", f"*.diff_{file_stem}.DESeq2.txt"))
            )

        def classify_region(region_text):
            region = str(region_text).strip().lower()
            if region.startswith("promoter"):
                return "promoter"
            if region in {"intergenic", "distal"}:
                return "distal"
            if region in {"intron", "exon", "intronic", "exonic", "gene_body", "body"}:
                return "gene_body"
            return "other"

        def ranked_rows_from_de_table(table_path, peak_region_map, broad_mode, promoter_map, top_n=500):
            ranked = []
            with open(table_path, newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames is None:
                    return {}
                for row in reader:
                    peak_id = str(row.get("gene_id", "")).strip()
                    if not peak_id:
                        continue
                    parsed = parse_peak_coord(peak_id)
                    if parsed is None:
                        continue
                    chrom, start, end = parsed
                    try:
                        lfc = float(row.get("log2FoldChange", "0"))
                    except (TypeError, ValueError):
                        lfc = 0.0
                    padj_text = row.get("padj", row.get("pvalue", "1"))
                    try:
                        padj = float(padj_text)
                    except (TypeError, ValueError):
                        padj = 1.0
                    if not math.isfinite(padj):
                        padj = 1.0
                    region_info = peak_region_map.get(peak_id, {})
                    region = region_info.get("genomic_region", "NA")
                    assigned = region_info.get("assigned_genes", "NA")
                    nearest = region_info.get("nearest_gene", "NA")
                    bed_row = [chrom, start, end, peak_id, f"{lfc:.6f}", ".", region, assigned, nearest]
                    ranked.append(
                        {
                            "rank_key": (padj, -abs(lfc), chrom, start, end),
                            "bed_row": bed_row,
                            "region_class": classify_region(region),
                            "region_info": region_info,
                            "lfc": lfc,
                        }
                    )
            ranked = sorted(ranked, key=lambda item: item["rank_key"])
            top_all = ranked[:top_n]
            top_promoter = [item for item in ranked if item["region_class"] == "promoter"][:top_n]
            top_distal = [item for item in ranked if item["region_class"] == "distal"][:top_n]
            top_sets = {
                "top_de_ranked": [item["bed_row"] for item in top_all],
                "top_de_ranked_promoter": [item["bed_row"] for item in top_promoter],
                "top_de_ranked_distal": [item["bed_row"] for item in top_distal],
            }
            if broad_mode == "genebody":
                promoter_region_rows = []
                for item in top_all:
                    peak_id = item["bed_row"][3]
                    for gene_label in labels_for_feature(item["region_info"]):
                        promoter_intervals = promoter_map.get(gene_label, [])
                        if not promoter_intervals:
                            continue
                        for prom_chrom, prom_start, prom_end in promoter_intervals:
                            promoter_region_rows.append(
                                [
                                    prom_chrom,
                                    prom_start,
                                    prom_end,
                                    peak_id,
                                    f"{item['lfc']:.6f}",
                                    ".",
                                    item["bed_row"][6],
                                    item["bed_row"][7],
                                    item["bed_row"][8],
                                ]
                            )
                top_sets["top_de_ranked_promoter_regions"] = sorted(
                    {tuple(row) for row in promoter_region_rows},
                    key=lambda x: (x[0], x[1], x[2], x[3]),
                )
            return top_sets

        def build_de_sets(table_path, peak_region_map, sets_dir, broad_mode, promoter_map):
            rel_name = os.path.relpath(table_path, params.de_outputfolder)
            parts = rel_name.split(os.sep)
            if len(parts) < 3:
                return []
            de_subdir = parts[0]
            contrast_dir = parts[1]
            file_stem = "peaks"
            if broad_mode == "genebody":
                file_stem = "gene_bodies"
            elif broad_mode == "diffuse":
                file_stem = "bins"
            contrast_label = (
                parts[2]
                .replace(f".sig_diff_{file_stem}.DESeq2.txt", "")
                .replace(f".diff_{file_stem}.DESeq2.txt", "")
            )

            rows_all = []
            rows_up = []
            rows_down = []
            rows_promoter = []
            rows_distal = []
            rows_promoter_regions = []
            rows_up_promoter_regions = []
            rows_down_promoter_regions = []
            diffuse_rows_all = []
            diffuse_rows_up = []
            diffuse_rows_down = []

            with open(table_path, newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames is None:
                    return []
                for row in reader:
                    peak_id = str(row.get("gene_id", "")).strip()
                    if not peak_id:
                        continue
                    parsed = parse_peak_coord(peak_id)
                    if parsed is None:
                        continue
                    chrom, start, end = parsed
                    lfc_text = row.get("log2FoldChange")
                    try:
                        lfc = float(lfc_text)
                    except (TypeError, ValueError):
                        lfc = 0.0

                    region_info = peak_region_map.get(peak_id, {})
                    region = region_info.get("genomic_region", "NA")
                    assigned = region_info.get("assigned_genes", "NA")
                    nearest = region_info.get("nearest_gene", "NA")
                    row_out = [chrom, start, end, peak_id, f"{lfc:.6f}", ".", region, assigned, nearest]
                    rows_all.append(row_out)
                    diffuse_row = {
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "feature_id": peak_id,
                        "lfc": lfc,
                        "region": region,
                        "assigned": assigned,
                        "nearest": nearest,
                    }
                    diffuse_rows_all.append(diffuse_row)
                    if lfc > 0:
                        rows_up.append(row_out)
                        diffuse_rows_up.append(diffuse_row)
                    if lfc < 0:
                        rows_down.append(row_out)
                        diffuse_rows_down.append(diffuse_row)

                    region_class = classify_region(region)
                    if region_class == "promoter":
                        rows_promoter.append(row_out)
                    elif region_class == "distal":
                        rows_distal.append(row_out)

                    if broad_mode == "genebody":
                        promoter_rows_for_feature = []
                        for gene_label in labels_for_feature(region_info):
                            promoter_intervals = promoter_map.get(gene_label, [])
                            if not promoter_intervals:
                                continue
                            for prom_chrom, prom_start, prom_end in promoter_intervals:
                                promoter_rows_for_feature.append(
                                    [prom_chrom, prom_start, prom_end, peak_id, f"{lfc:.6f}", ".", region, assigned, nearest]
                                )
                        rows_promoter_regions.extend(promoter_rows_for_feature)
                        if lfc > 0:
                            rows_up_promoter_regions.extend(promoter_rows_for_feature)
                        if lfc < 0:
                            rows_down_promoter_regions.extend(promoter_rows_for_feature)

            out_rows = []
            prefix = f"{de_subdir}__{contrast_dir}__{contrast_label}"
            ranked_source_table = table_path.replace(f".sig_diff_{file_stem}.DESeq2.txt", f".diff_{file_stem}.DESeq2.txt")
            if not os.path.exists(ranked_source_table):
                ranked_source_table = table_path
            set_specs = [
                ("de_significant", rows_all),
                ("de_up", rows_up),
                ("de_down", rows_down),
                ("de_significant_promoter", rows_promoter),
                ("de_significant_distal", rows_distal),
            ]
            ranked_sets = ranked_rows_from_de_table(
                ranked_source_table,
                peak_region_map,
                broad_mode,
                promoter_map,
            )
            set_specs.extend(
                [
                    ("top_de_ranked", ranked_sets.get("top_de_ranked", [])),
                    ("top_de_ranked_promoter", ranked_sets.get("top_de_ranked_promoter", [])),
                    ("top_de_ranked_distal", ranked_sets.get("top_de_ranked_distal", [])),
                ]
            )
            if broad_mode == "genebody":
                set_specs.extend(
                    [
                        ("de_significant_promoter_regions", rows_promoter_regions),
                        ("de_up_promoter_regions", rows_up_promoter_regions),
                        ("de_down_promoter_regions", rows_down_promoter_regions),
                        ("top_de_ranked_promoter_regions", ranked_sets.get("top_de_ranked_promoter_regions", [])),
                    ]
                )
            for set_name, set_rows in set_specs:
                out_path = os.path.join(ensure_dir(os.path.join(sets_dir, set_name)), f"{prefix}.{set_name}.bed")
                if broad_mode == "genebody" and set_name.endswith("_promoter_regions"):
                    set_rows = sorted({tuple(row) for row in set_rows}, key=lambda x: (x[0], x[1], x[2], x[3]))
                write_bed_from_rows(set_rows, out_path)
                out_rows.append(
                    {
                        "set_name": f"{prefix}.{set_name}",
                        "set_type": set_name,
                        "contrast_group": de_subdir,
                        "contrast_label": contrast_label,
                        "source_table": table_path,
                        "peak_bed": out_path,
                        "peak_count": len(set_rows),
                    }
                )
            if broad_mode == "diffuse":
                out_rows.extend(
                    write_diffuse_domain_sets(
                        prefix=prefix,
                        contrast_label=contrast_label,
                        de_subdir=de_subdir,
                        source_table=table_path,
                        sets_dir=sets_dir,
                        rows_all=diffuse_rows_all,
                        rows_up=diffuse_rows_up,
                        rows_down=diffuse_rows_down,
                        merge_gap=params.diffuse_merge_gap,
                    )
                )
            return out_rows

        def merge_diffuse_rows_into_domains(rows, merge_gap):
            if not rows:
                return []
            rows_sorted = sorted(rows, key=lambda item: (item["chrom"], item["start"], item["end"]))
            domains = []
            current = None
            for row in rows_sorted:
                if current is None:
                    current = {"chrom": row["chrom"], "start": row["start"], "end": row["end"], "rows": [row]}
                    continue
                if row["chrom"] == current["chrom"] and row["start"] <= current["end"] + int(merge_gap):
                    current["end"] = max(current["end"], row["end"])
                    current["rows"].append(row)
                else:
                    domains.append(current)
                    current = {"chrom": row["chrom"], "start": row["start"], "end": row["end"], "rows": [row]}
            if current is not None:
                domains.append(current)
            return domains

        def write_diffuse_domain_sets(prefix, contrast_label, de_subdir, source_table, sets_dir, rows_all, rows_up, rows_down, merge_gap):
            domain_specs = [
                ("de_significant_domains", rows_all),
                ("de_up_domains", rows_up),
                ("de_down_domains", rows_down),
            ]
            out_rows = []
            domains_summary_dir = ensure_dir(os.path.join(os.path.dirname(sets_dir), "domains"))
            for set_name, source_rows in domain_specs:
                merged_domains = merge_diffuse_rows_into_domains(source_rows, merge_gap)
                bed_path = os.path.join(sets_dir, set_name, f"{prefix}.{set_name}.bed")
                summary_path = os.path.join(domains_summary_dir, f"{prefix}.{set_name}.tsv")
                bed_rows = []
                summary_rows = []
                for idx, domain in enumerate(merged_domains, start=1):
                    domain_rows = domain["rows"]
                    lfcs = [item["lfc"] for item in domain_rows]
                    regions = [str(item["region"]).strip() or "NA" for item in domain_rows]
                    assigned_vals = [str(item["assigned"]).strip() for item in domain_rows if str(item["assigned"]).strip() not in {"", "NA"}]
                    nearest_vals = [str(item["nearest"]).strip() for item in domain_rows if str(item["nearest"]).strip() not in {"", "NA"}]
                    dominant_region = max(sorted(set(regions)), key=regions.count) if regions else "NA"
                    assigned_join = ",".join(sorted(set(assigned_vals))) if assigned_vals else "NA"
                    nearest_join = ",".join(sorted(set(nearest_vals))) if nearest_vals else "NA"
                    domain_id = f"{contrast_label}.domain_{idx}"
                    mean_lfc = sum(lfcs) / len(lfcs)
                    max_abs_lfc = max(abs(value) for value in lfcs)
                    bed_rows.append([domain["chrom"], domain["start"], domain["end"], domain_id, f"{mean_lfc:.6f}", ".", dominant_region, assigned_join, nearest_join])
                    summary_rows.append([domain_id, domain["chrom"], domain["start"], domain["end"], domain["end"] - domain["start"], len(domain_rows), f"{mean_lfc:.6f}", f"{max_abs_lfc:.6f}", dominant_region, assigned_join, nearest_join])
                write_bed_from_rows(bed_rows, bed_path)
                with open(summary_path, "w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["domain_id", "chrom", "start", "end", "width_bp", "bin_count", "mean_log2FoldChange", "max_abs_log2FoldChange", "dominant_region", "assigned_genes", "nearest_genes"])
                    writer.writerows(summary_rows)
                out_rows.append(
                    {
                        "set_name": f"{prefix}.{set_name}",
                        "set_type": set_name,
                        "contrast_group": de_subdir,
                        "contrast_label": contrast_label,
                        "source_table": source_table,
                        "peak_bed": bed_path,
                        "peak_count": len(bed_rows),
                    }
                )
            return out_rows

        def build_all_feature_sets(peak_region_map, sets_dir, broad_mode):
            rows_all = feature_rows_from_peak_metadata(peak_region_map)
            set_name = "all_features"
            if broad_mode == "genebody":
                set_name = "all_gene_bodies"
            elif broad_mode == "diffuse":
                set_name = "all_bins"
            out_path = os.path.join(ensure_dir(os.path.join(sets_dir, set_name)), f"{set_name}.bed")
            write_bed_from_rows(rows_all, out_path)
            out_rows = [
                {
                    "set_name": set_name,
                    "set_type": set_name,
                    "contrast_group": "global",
                    "contrast_label": set_name,
                    "source_table": os.path.join(params.de_outputfolder, "peak_metadata.tsv"),
                    "peak_bed": out_path,
                    "peak_count": len(rows_all),
                }
            ]
            if broad_mode != "genebody":
                region_specs = [
                    ("all_promoter", "promoter"),
                    ("all_distal", "distal"),
                    ("all_gene_body", "gene_body"),
                ]
                for split_name, split_class in region_specs:
                    split_rows = [row for row in rows_all if classify_region(row[6]) == split_class]
                    split_path = os.path.join(ensure_dir(os.path.join(sets_dir, split_name)), f"{split_name}.bed")
                    write_bed_from_rows(split_rows, split_path)
                    out_rows.append(
                        {
                            "set_name": split_name,
                            "set_type": split_name,
                            "contrast_group": "global",
                            "contrast_label": split_name,
                            "source_table": os.path.join(params.de_outputfolder, "peak_metadata.tsv"),
                            "peak_bed": split_path,
                            "peak_count": len(split_rows),
                        }
                    )
            return out_rows

        def collect_unique_sets_from_prede(sets_dir):
            unique_input_dir = os.path.join(params.prede_outputfolder, "analyze_peaks", "intersections", "unique")
            unique_out_dir = ensure_dir(os.path.join(sets_dir, "unique_from_prede"))
            out_rows = []
            if not os.path.isdir(unique_input_dir):
                return out_rows
            for src in sorted(glob.glob(os.path.join(unique_input_dir, "*.unique.bed"))):
                dst = os.path.join(unique_out_dir, os.path.basename(src))
                shutil.copy2(src, dst)
                peak_count = int(
                    subprocess.check_output(
                        f"wc -l {quote(dst)} | awk '{{print $1}}'",
                        shell=True,
                        executable="/bin/bash",
                        text=True,
                    ).strip()
                )
                out_rows.append(
                    {
                        "set_name": os.path.basename(src).replace(".bed", ""),
                        "set_type": "unique_from_prede",
                        "contrast_group": "pre_de",
                        "contrast_label": "pre_de_unique",
                        "source_table": src,
                        "peak_bed": dst,
                        "peak_count": peak_count,
                    }
                )
            return out_rows

        def run_deeptools(set_manifest_rows, signal_dir):
            required = ["computeMatrix", "plotHeatmap"]
            if not all(tool_available(x) for x in required):
                log_it(logfile, "deepTools executables not fully available. Skipping signal plots in analyze_peaks_de.")
                return

            signal_summary_path = os.path.join(signal_dir, "signal_runs.tsv")
            signal_summary_rows = []
            eligible_set_types = {
                "de_significant",
                "de_up",
                "de_down",
                "de_significant_promoter",
                "de_significant_distal",
                "de_significant_domains",
                "de_up_domains",
                "de_down_domains",
                "de_significant_promoter_regions",
                "de_up_promoter_regions",
                "de_down_promoter_regions",
                "top_de_ranked",
                "top_de_ranked_promoter",
                "top_de_ranked_distal",
                "top_de_ranked_promoter_regions",
            }
            min_signal_regions = 10
            max_signal_regions = 250
            if params.thetype == "CHIP" and params.broad_mode in {"genebody", "diffuse"}:
                max_signal_regions = 100
            deeptools_threads = 1

            def write_signal_summary():
                with open(signal_summary_path, "w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["set_name", "status", "source_bed", "signal_bed", "matrix", "heatmap", "reason"])
                    writer.writerows(signal_summary_rows)

            def regions_label_for_set(item):
                set_type = str(item.get("set_type", "features"))
                label_map = {
                    "all_features": "all peaks",
                    "all_promoter": "promoter peaks",
                    "all_distal": "distal peaks",
                    "all_gene_body": "gene-body peaks",
                    "all_gene_bodies": "gene bodies",
                    "all_bins": "bins",
                    "de_significant": "DE features",
                    "de_up": "DE up features",
                    "de_down": "DE down features",
                    "de_significant_promoter": "DE promoter peaks",
                    "de_significant_distal": "DE distal peaks",
                    "de_significant_domains": "DE domains",
                    "de_up_domains": "DE up domains",
                    "de_down_domains": "DE down domains",
                    "de_significant_promoter_regions": "DE promoter regions",
                    "de_up_promoter_regions": "DE up promoter regions",
                    "de_down_promoter_regions": "DE down promoter regions",
                    "top_de_ranked": "top ranked DE features",
                    "top_de_ranked_promoter": "top ranked promoter peaks",
                    "top_de_ranked_distal": "top ranked distal peaks",
                    "top_de_ranked_promoter_regions": "top ranked promoter regions",
                    "unique_from_prede": "unique pre-DE peaks",
                }
                return label_map.get(set_type, set_type.replace("_", " "))

            def copy_bed_limit(source_bed, target_bed, limit):
                copied = 0
                with open(source_bed, "r", newline="") as source, open(target_bed, "w", newline="") as target:
                    for line in source:
                        if copied >= limit:
                            break
                        if not line.strip():
                            continue
                        target.write(line)
                        copied += 1
                return copied

            def file_sha256(path):
                digest = hashlib.sha256()
                with open(path, "rb") as handle:
                    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                        digest.update(chunk)
                return digest.hexdigest()

            def existing_bigwig_for_sample(sample):
                candidates = [
                    os.path.join(params.bigwig_inputfolder, f"{sample}.bw"),
                    os.path.join(params.bigwig_inputfolder, f"{sample}.CPM.bw"),
                ]
                for candidate in candidates:
                    if os.path.isfile(candidate):
                        return candidate
                return None

            def render_profile_from_matrix(matrix_path, profile_path, region_label, sample_labels):
                try:
                    import matplotlib
                    matplotlib.use("Agg")
                    import matplotlib.pyplot as plt
                except Exception as exc:
                    raise RuntimeError(f"matplotlib import failed: {exc}") from exc

                with gzip.open(matrix_path, "rt") as handle:
                    header_line = handle.readline().strip()
                    if not header_line.startswith("@"):
                        raise RuntimeError("deepTools matrix header is missing")
                    matrix_header = json.loads(header_line[1:])

                    sample_boundaries = [int(x) for x in matrix_header.get("sample_boundaries", [])]
                    if len(sample_boundaries) < 2:
                        raise RuntimeError("deepTools matrix header lacks sample boundaries")
                    n_samples = len(sample_boundaries) - 1
                    if len(sample_labels) != n_samples:
                        sample_labels = [str(x) for x in matrix_header.get("sample_labels", sample_labels)]
                    if len(sample_labels) != n_samples:
                        sample_labels = [f"sample_{index + 1}" for index in range(n_samples)]

                    bin_counts = [sample_boundaries[i + 1] - sample_boundaries[i] for i in range(n_samples)]
                    if not bin_counts or any(count <= 0 for count in bin_counts):
                        raise RuntimeError("deepTools matrix sample boundaries are invalid")

                    sums = [[0.0 for _ in range(bin_counts[i])] for i in range(n_samples)]
                    counts = [[0 for _ in range(bin_counts[i])] for i in range(n_samples)]
                    for line in handle:
                        if not line.strip():
                            continue
                        fields = line.rstrip("\n").split("\t")
                        values = fields[6:]
                        for sample_index in range(n_samples):
                            start = sample_boundaries[sample_index]
                            end = sample_boundaries[sample_index + 1]
                            for bin_index, raw_value in enumerate(values[start:end]):
                                try:
                                    value = float(raw_value)
                                except ValueError:
                                    continue
                                if not math.isfinite(value):
                                    continue
                                sums[sample_index][bin_index] += value
                                counts[sample_index][bin_index] += 1

                upstream = int(matrix_header.get("upstream", [1000])[0] or 1000)
                downstream = int(matrix_header.get("downstream", [1000])[0] or 1000)
                bin_size = int(matrix_header.get("bin size", [100])[0] or 100)
                n_bins = bin_counts[0]
                x_values = [(-upstream + (index + 0.5) * bin_size) / 1000.0 for index in range(n_bins)]
                profiles = []
                for sample_index in range(n_samples):
                    sample_profile = []
                    for bin_index in range(bin_counts[sample_index]):
                        if counts[sample_index][bin_index] == 0:
                            sample_profile.append(float("nan"))
                        else:
                            sample_profile.append(sums[sample_index][bin_index] / counts[sample_index][bin_index])
                    profiles.append(sample_profile)

                color_cycle = [
                    "#1f77b4",
                    "#2ca02c",
                    "#ff7f0e",
                    "#d62728",
                    "#9467bd",
                    "#17becf",
                    "#bcbd22",
                    "#8c564b",
                ]
                fig = plt.figure(figsize=(9.5, 6.5))
                grid = fig.add_gridspec(
                    nrows=2,
                    ncols=1,
                    height_ratios=[4.8, 1.0],
                    hspace=0.28,
                    left=0.12,
                    right=0.98,
                    top=0.91,
                    bottom=0.08,
                )
                ax = fig.add_subplot(grid[0])
                legend_ax = fig.add_subplot(grid[1])
                legend_ax.axis("off")

                lines = []
                for sample_index, sample_profile in enumerate(profiles):
                    line, = ax.plot(
                        x_values[: len(sample_profile)],
                        sample_profile,
                        linewidth=2.4,
                        color=color_cycle[sample_index % len(color_cycle)],
                        label=sample_labels[sample_index],
                    )
                    lines.append(line)

                ax.axvline(0, color="#555555", linewidth=1.0, alpha=0.6)
                ax.set_title(region_label, fontsize=16, pad=12)
                ax.set_xlabel("distance from center (kb)", fontsize=12)
                ax.set_ylabel("normalized signal", fontsize=12)
                ax.set_xlim(-upstream / 1000.0, downstream / 1000.0)
                ax.set_xticks([-upstream / 1000.0, 0, downstream / 1000.0])
                ax.set_xticklabels([f"-{upstream / 1000.0:.1f}", "center", f"{downstream / 1000.0:.1f}"])
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.grid(True, axis="y", color="#dddddd", linewidth=0.8)
                legend_ax.legend(
                    handles=lines,
                    labels=sample_labels,
                    loc="center",
                    ncol=min(3, max(1, len(sample_labels))),
                    frameon=False,
                    fontsize=11,
                )
                fig.savefig(profile_path)
                plt.close(fig)

            matrices_dir = ensure_dir(os.path.join(signal_dir, "matrices"))
            heatmaps_dir = ensure_dir(os.path.join(signal_dir, "heatmaps"))
            profiles_dir = ensure_dir(os.path.join(signal_dir, "profiles"))
            regions_dir = ensure_dir(os.path.join(signal_dir, "regions"))

            if params.post_de_signal_policy == "skip":
                reason = "post-DE signal plotting skipped by --post-de-signal-policy skip"
                log_it(logfile, f"Skipping analyze_peaks_de signal plots because {reason}.")
                for item in set_manifest_rows:
                    signal_summary_rows.append([item["set_name"], "SKIP", item["peak_bed"], "", "", "", reason])
                write_signal_summary()
                return

            bigwigs = []
            sample_labels = []
            missing_bigwig_samples = []
            for sample in samples2:
                bw_path = existing_bigwig_for_sample(sample)
                if bw_path:
                    bigwigs.append(bw_path)
                    sample_labels.append(sample)
                else:
                    missing_bigwig_samples.append(sample)
            if missing_bigwig_samples:
                missing_text = ", ".join(missing_bigwig_samples)
                reason = f"missing existing BigWigs for samples: {missing_text}"
                if params.post_de_signal_policy == "require" or (
                    params.post_de_signal_policy == "auto" and params.post_de_signal_auto_dependencies
                ):
                    raise RuntimeError(
                        f"Post-DE signal plots require BigWigs, but {reason}. "
                        "Run step 8 first or use --post-de-signal-policy skip."
                    )
                log_it(logfile, f"Skipping analyze_peaks_de signal plots because {reason}. Run step 8 before step 16 to enable these plots.")
                for item in set_manifest_rows:
                    signal_summary_rows.append([item["set_name"], "SKIP", item["peak_bed"], "", "", "", reason])
                write_signal_summary()
                return
            if not bigwigs:
                reason = "no existing BigWigs found"
                if params.post_de_signal_policy == "require" or (
                    params.post_de_signal_policy == "auto" and params.post_de_signal_auto_dependencies
                ):
                    raise RuntimeError(
                        "Post-DE signal plots require BigWigs, but no existing BigWigs were found. "
                        "Run step 8 first or use --post-de-signal-policy skip."
                    )
                log_it(logfile, "No existing BigWigs found for analyze_peaks_de deepTools. Skipping signal plots.")
                for item in set_manifest_rows:
                    signal_summary_rows.append([item["set_name"], "SKIP", item["peak_bed"], "", "", "", reason])
                write_signal_summary()
                return

            signal_artifact_cache = {}
            for item in set_manifest_rows:
                set_type = str(item.get("set_type", ""))
                if set_type not in eligible_set_types:
                    signal_summary_rows.append([item["set_name"], "SKIP", item["peak_bed"], "", "", "", "set type is not configured for signal plotting"])
                    continue
                bed_path = item["peak_bed"]
                if not os.path.isfile(bed_path):
                    signal_summary_rows.append([item["set_name"], "SKIP", bed_path, "", "", "", "BED file missing"])
                    continue
                if int(item["peak_count"]) < min_signal_regions:
                    signal_summary_rows.append([item["set_name"], "SKIP", bed_path, "", "", "", f"fewer than {min_signal_regions} regions"])
                    continue
                safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", item["set_name"])
                region_label = regions_label_for_set(item)
                signal_bed_path = os.path.join(regions_dir, f"{safe_name}.signal_regions.bed")
                signal_region_count = copy_bed_limit(bed_path, signal_bed_path, max_signal_regions)
                if signal_region_count < min_signal_regions:
                    signal_summary_rows.append([item["set_name"], "SKIP", bed_path, signal_bed_path, "", "", f"fewer than {min_signal_regions} copied regions"])
                    continue
                matrix_path = os.path.join(matrices_dir, f"{safe_name}.matrix.gz")
                heatmap_path = os.path.join(heatmaps_dir, f"{safe_name}.heatmap.pdf")
                profile_path = os.path.join(profiles_dir, f"{safe_name}.profile.pdf")
                cap_reason = ""
                if int(item["peak_count"]) > signal_region_count:
                    cap_reason = f"signal input capped at first {signal_region_count} of {item['peak_count']} regions"
                signal_signature = (
                    file_sha256(signal_bed_path),
                    tuple(bigwigs),
                    tuple(sample_labels),
                    region_label,
                    "reference-point:center:-1000:+1000:bin100:missingDataAsZero",
                )
                cached = signal_artifact_cache.get(signal_signature)
                if cached:
                    shutil.copy2(cached["matrix"], matrix_path)
                    shutil.copy2(cached["heatmap"], heatmap_path)
                    shutil.copy2(cached["profile"], profile_path)
                    reason_parts = [f"reused signal plots from {cached['set_name']}"]
                    if cap_reason:
                        reason_parts.append(cap_reason)
                    signal_summary_rows.append(
                        [item["set_name"], "REUSE", bed_path, signal_bed_path, matrix_path, heatmap_path, "; ".join(reason_parts)]
                    )
                    continue
                shell(
                    f"computeMatrix reference-point --referencePoint center -b 1000 -a 1000 "
                    f"-R {quote(signal_bed_path)} -S {' '.join(quote(x) for x in bigwigs)} "
                    f"--missingDataAsZero --binSize 100 --numberOfProcessors {deeptools_threads} -o {quote(matrix_path)}"
                )
                shell(
                    f"plotHeatmap -m {quote(matrix_path)} -out {quote(heatmap_path)} "
                    f"--whatToShow 'heatmap and colorbar' --sortRegions descend "
                    f"--regionsLabel {quote(region_label)} "
                    f"--samplesLabel {' '.join(quote(label) for label in sample_labels)} "
                    f"--xAxisLabel 'distance from center (bp)' --refPointLabel center"
                )
                try:
                    render_profile_from_matrix(matrix_path, profile_path, region_label, sample_labels)
                except Exception as exc:
                    if not tool_available("plotProfile"):
                        raise RuntimeError(f"Custom profile rendering failed for {item['set_name']}, and plotProfile is not available: {exc}") from exc
                    log_it(logfile, f"Custom profile rendering failed for {item['set_name']}: {exc}. Falling back to deepTools plotProfile.")
                    shell(
                        f"plotProfile -m {quote(matrix_path)} -out {quote(profile_path)} "
                        f"--perGroup --regionsLabel {quote(region_label)} "
                        f"--samplesLabel {' '.join(quote(label) for label in sample_labels)} "
                        f"--refPointLabel center "
                        f"--yAxisLabel 'normalized signal' --legendLocation upper-right "
                        f"--plotWidth 10 --plotHeight 6"
                    )
                signal_artifact_cache[signal_signature] = {
                    "set_name": item["set_name"],
                    "matrix": matrix_path,
                    "heatmap": heatmap_path,
                    "profile": profile_path,
                }
                signal_summary_rows.append([item["set_name"], "OK", bed_path, signal_bed_path, matrix_path, heatmap_path, cap_reason])

            write_signal_summary()

        def run_meme_motifs(set_manifest_rows, motifs_dir):
            if os.path.isdir(motifs_dir):
                shutil.rmtree(motifs_dir)
            motifs_dir = ensure_dir(motifs_dir)
            motif_summary_path = os.path.join(motifs_dir, "motif_runs.tsv")
            eligible_set_types = {
                "de_significant",
                "de_up",
                "de_down",
                "de_significant_promoter",
                "de_significant_distal",
                "unique_from_prede",
                "de_significant_promoter_regions",
                "de_up_promoter_regions",
                "de_down_promoter_regions",
                "top_de_ranked",
                "top_de_ranked_promoter",
                "top_de_ranked_distal",
                "top_de_ranked_promoter_regions",
            }
            min_motif_peaks = 10
            max_motif_peaks = max(10, int(params.post_de_motif_max_peaks))
            max_motif_sets = max(0, int(params.post_de_motif_max_sets))
            motif_timeout_seconds = max(60, int(params.post_de_motif_timeout_seconds))
            motif_threads = max(1, min(int(params.post_de_motif_threads), int(threads)))
            motif_window_bp = max(50, int(params.post_de_motif_window_bp))
            motif_summary = []

            def write_motif_summary():
                with open(motif_summary_path, "w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["set_name", "method", "status", "peak_bed", "input_fasta", "output_dir", "motif_database", "report_pdf", "reason"])
                    writer.writerows(motif_summary)

            required_tools = ["bedtools", "sea"]
            missing_tools = [tool_name for tool_name in required_tools if not tool_available(tool_name)]
            if missing_tools:
                reason = "missing required executable(s): " + ", ".join(missing_tools)
                log_it(logfile, f"Skipping MEME motif analysis in analyze_peaks_de because {reason}.")
                motif_summary.append(["ALL", "MEME", "SKIP", "", "", "", "", "", reason])
                write_motif_summary()
                return motif_summary_path

            def resolve_meme_motif_database():
                try:
                    from omnomnomics.genomes import resolve_meme_motif_database as resolve_reference_motif_database
                except Exception:
                    return None

                try:
                    return resolve_reference_motif_database(
                        params.genome_assembly_dir,
                        requested=params.post_de_motif_database,
                    )
                except Exception as exc:
                    log_it(logfile, f"Could not resolve MEME motif database in analyze_peaks_de: {exc}")
                    return None
                return None

            motif_database = resolve_meme_motif_database()
            if not motif_database:
                reason = (
                    "no MEME motif database found or cached; run 'omnomnomics genomes motifs' "
                    "or set post_de_motif_database to a MEME-format motif file"
                )
                log_it(logfile, f"Skipping MEME motif analysis in analyze_peaks_de because {reason}.")
                motif_summary.append(["ALL", "MEME", "SKIP", "", "", "", "", "", reason])
                write_motif_summary()
                return motif_summary_path

            priority_by_set_type = {
                "de_significant_promoter": 0,
                "de_significant_promoter_regions": 0,
                "de_up_promoter_regions": 1,
                "de_down_promoter_regions": 1,
                "de_significant_distal": 2,
                "top_de_ranked_promoter": 3,
                "top_de_ranked_promoter_regions": 3,
                "top_de_ranked_distal": 4,
                "unique_from_prede": 5,
                "de_significant": 6,
                "de_up": 7,
                "de_down": 7,
                "top_de_ranked": 8,
            }

            candidate_indices = []
            for index, item in enumerate(set_manifest_rows):
                if item["set_type"] not in eligible_set_types:
                    continue
                if int(item["peak_count"]) < min_motif_peaks:
                    continue
                candidate_indices.append(index)
            candidate_indices = sorted(
                candidate_indices,
                key=lambda index: (
                    priority_by_set_type.get(set_manifest_rows[index]["set_type"], 99),
                    set_manifest_rows[index]["contrast_group"],
                    set_manifest_rows[index]["contrast_label"],
                    set_manifest_rows[index]["set_name"],
                ),
            )
            selected_indices = set(candidate_indices[:max_motif_sets]) if max_motif_sets > 0 else set(candidate_indices)
            if max_motif_sets > 0 and len(candidate_indices) > max_motif_sets:
                log_it(
                    logfile,
                    f"Motif analysis capped at {max_motif_sets} highest-priority sets "
                    f"out of {len(candidate_indices)} eligible sets.",
                )
            log_it(
                logfile,
                f"Motif analysis settings: max_sets={max_motif_sets or 'unlimited'}, "
                f"max_peaks={max_motif_peaks}, window_bp={motif_window_bp}, "
                f"timeout_seconds={motif_timeout_seconds}, threads={motif_threads}, motif_database={motif_database}",
            )

            def write_centered_bed_and_scores(source_bed, target_bed, limit, window_bp):
                copied = 0
                scores = {}
                half_window = window_bp // 2
                with open(source_bed, "r", newline="") as source, open(target_bed, "w", newline="") as target:
                    for line in source:
                        if copied >= limit:
                            break
                        if not line.strip() or line.startswith("#"):
                            continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) < 3:
                            continue
                        try:
                            start = int(fields[1])
                            end = int(fields[2])
                        except ValueError:
                            continue
                        if end <= start:
                            continue
                        chrom = fields[0]
                        center = (start + end) // 2
                        window_start = max(0, center - half_window)
                        window_end = window_start + window_bp
                        raw_name = fields[3] if len(fields) > 3 and fields[3] else f"{chrom}:{start}-{end}"
                        safe_name = re.sub(r"[^A-Za-z0-9._:-]+", "_", raw_name)
                        fasta_name = f"{copied + 1}_{safe_name}"
                        score_value = float(copied + 1)
                        scores[fasta_name] = score_value
                        target.write("\t".join([chrom, str(window_start), str(window_end), fasta_name, f"{score_value:.6f}"]) + "\n")
                        copied += 1
                return copied, scores

            def rewrite_fasta_headers_with_scores(fasta_path, scores):
                tmp_path = f"{fasta_path}.tmp"
                with open(fasta_path, "r", newline="") as source, open(tmp_path, "w", newline="") as target:
                    for line in source:
                        if not line.startswith(">"):
                            target.write(line)
                            continue
                        header = line[1:].strip()
                        name = header.split("::", 1)[0].split()[0]
                        score = scores.get(name, 0.0)
                        target.write(f">{name} {score:.6f}\n")
                os.replace(tmp_path, fasta_path)

            def run_meme_command(cmd, timeout_seconds, env):
                completed = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    timeout=timeout_seconds,
                    env=env,
                )
                if completed.stdout:
                    for line in completed.stdout.splitlines():
                        log_it(logfile, line)
                if completed.returncode != 0:
                    raise RuntimeError(f"{cmd[0]} exited with status {completed.returncode}")
                return completed

            def run_meme_text_command(cmd, output_path, timeout_seconds, env):
                completed = subprocess.run(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=timeout_seconds,
                    env=env,
                )
                if completed.stderr:
                    for line in completed.stderr.splitlines():
                        log_it(logfile, line)
                if completed.stdout:
                    with open(output_path, "w", newline="") as output_handle:
                        output_handle.write(completed.stdout)
                if completed.returncode != 0:
                    raise RuntimeError(f"{cmd[0]} exited with status {completed.returncode}")
                return completed

            def parse_meme_text_table(table_path):
                if not os.path.isfile(table_path):
                    return [], []
                header = None
                rows = []
                with open(table_path, "r", newline="") as handle:
                    for line in handle:
                        line = line.rstrip("\n")
                        if not line.strip() or line.startswith("#"):
                            continue
                        fields = line.split("\t")
                        if header is None:
                            upper = {field.upper() for field in fields}
                            if "RANK" in upper and ("ID" in upper or "MOTIF_ID" in upper):
                                header = fields
                                continue
                            return [], []
                        if len(fields) < len(header):
                            fields.extend([""] * (len(header) - len(fields)))
                        rows.append(dict(zip(header, fields)))
                return header or [], rows

            def count_meme_text_rows(table_path):
                _header, rows = parse_meme_text_table(table_path)
                return len(rows)

            def safe_float(value, default=float("nan")):
                try:
                    out = float(value)
                except (TypeError, ValueError):
                    return default
                return out if math.isfinite(out) else default

            def render_motif_report(method, set_label, table_path, report_pdf):
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt

                header, rows = parse_meme_text_table(table_path)
                os.makedirs(os.path.dirname(report_pdf), exist_ok=True)
                if not rows:
                    fig, ax = plt.subplots(figsize=(8.5, 4.5))
                    ax.axis("off")
                    ax.text(
                        0.5,
                        0.62,
                        f"{method} motif enrichment",
                        ha="center",
                        va="center",
                        fontsize=16,
                        fontweight="bold",
                    )
                    ax.text(
                        0.5,
                        0.44,
                        f"{set_label}\nNo enriched motifs reported by {method}.",
                        ha="center",
                        va="center",
                        fontsize=11,
                    )
                    fig.savefig(report_pdf, bbox_inches="tight")
                    plt.close(fig)
                    return report_pdf

                upper_to_col = {col.upper(): col for col in header}
                label_col = upper_to_col.get("ALT_ID") or upper_to_col.get("ID") or upper_to_col.get("MOTIF_ID")
                q_col = (
                    upper_to_col.get("QVALUE")
                    or upper_to_col.get("ADJ_PVALUE")
                    or upper_to_col.get("PVALUE")
                    or upper_to_col.get("P-VALUE")
                )
                effect_col = upper_to_col.get("ENR_RATIO") or upper_to_col.get("SCORE") or upper_to_col.get("TP%")
                if not label_col or not q_col:
                    fig, ax = plt.subplots(figsize=(8.5, 4.5))
                    ax.axis("off")
                    ax.text(
                        0.5,
                        0.5,
                        f"{method} motif enrichment rows were found, but required plotting columns are missing.",
                        ha="center",
                        va="center",
                        fontsize=11,
                    )
                    fig.savefig(report_pdf, bbox_inches="tight")
                    plt.close(fig)
                    return report_pdf

                plot_rows = []
                for row in rows:
                    q_value = safe_float(row.get(q_col))
                    if not math.isfinite(q_value) or q_value <= 0:
                        continue
                    label = str(row.get(label_col, "")).strip() or "motif"
                    effect = safe_float(row.get(effect_col), default=float("nan")) if effect_col else float("nan")
                    plot_rows.append(
                        {
                            "label": label,
                            "q_value": q_value,
                            "score": -math.log10(q_value),
                            "effect": effect,
                        }
                    )
                plot_rows = sorted(plot_rows, key=lambda item: item["q_value"])[:15]
                if not plot_rows:
                    return render_motif_report(method, set_label, "", report_pdf)

                labels = [item["label"] for item in reversed(plot_rows)]
                scores = [item["score"] for item in reversed(plot_rows)]
                colors = ["#4C78A8" if not math.isfinite(item["effect"]) else "#D95F02" for item in reversed(plot_rows)]
                fig_height = max(4.5, 0.35 * len(labels) + 1.8)
                fig, ax = plt.subplots(figsize=(9.5, fig_height))
                ax.barh(labels, scores, color=colors)
                ax.set_xlabel("-log10(q-value)")
                ax.set_title(f"{set_label}: top {method} motif enrichments")
                ax.grid(True, axis="x", color="#dddddd", linewidth=0.8)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                fig.savefig(report_pdf, bbox_inches="tight")
                plt.close(fig)
                return report_pdf

            def write_combined_motif_report(motifs_dir, motif_summary_rows):
                import matplotlib
                matplotlib.use("Agg")
                import matplotlib.pyplot as plt

                report_pdf = os.path.join(motifs_dir, "motif_summary.pdf")
                status_counts = {}
                for row in motif_summary_rows:
                    status = row[2] if len(row) > 2 else "NA"
                    status_counts[status] = status_counts.get(status, 0) + 1
                labels = sorted(status_counts)
                values = [status_counts[label] for label in labels]
                fig, ax = plt.subplots(figsize=(8, 5))
                ax.bar(labels, values, color="#4C78A8")
                ax.set_title("Motif analysis summary")
                ax.set_ylabel("Run count")
                ax.tick_params(axis="x", rotation=30)
                ax.grid(True, axis="y", color="#dddddd", linewidth=0.8)
                fig.savefig(report_pdf, bbox_inches="tight")
                plt.close(fig)
                return report_pdf

            for index, item in enumerate(set_manifest_rows):
                if item["set_type"] not in eligible_set_types:
                    motif_summary.append([item["set_name"], "MEME", "SKIP", item["peak_bed"], "", "", motif_database, "", "set type is not configured for motif analysis"])
                    write_motif_summary()
                    continue
                if int(item["peak_count"]) < min_motif_peaks:
                    motif_summary.append([item["set_name"], "MEME", "SKIP", item["peak_bed"], "", "", motif_database, "", f"fewer than {min_motif_peaks} peaks"])
                    write_motif_summary()
                    continue
                if index not in selected_indices:
                    reason = f"motif set cap reached; retained top {max_motif_sets} highest-priority sets"
                    motif_summary.append([item["set_name"], "MEME", "SKIP", item["peak_bed"], "", "", motif_database, "", reason])
                    write_motif_summary()
                    continue
                safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", item["set_name"])
                out_dir = os.path.join(motifs_dir, safe_name)
                if os.path.isdir(out_dir):
                    shutil.rmtree(out_dir)
                out_dir = ensure_dir(out_dir)
                motif_input_bed = os.path.join(out_dir, "motif_input.centered.bed")
                motif_input_fasta = os.path.join(out_dir, "motif_input.fa")
                motif_peak_count, fasta_scores = write_centered_bed_and_scores(
                    item["peak_bed"],
                    motif_input_bed,
                    max_motif_peaks,
                    motif_window_bp,
                )
                if motif_peak_count < min_motif_peaks:
                    motif_summary.append([item["set_name"], "MEME", "SKIP", item["peak_bed"], "", out_dir, motif_database, "", f"fewer than {min_motif_peaks} valid centered windows"])
                    write_motif_summary()
                    continue
                cap_reason = ""
                if int(item["peak_count"]) > motif_peak_count:
                    cap_reason = f"motif input capped at first {motif_peak_count} of {item['peak_count']} peaks"
                fasta_cmd = [
                    "bedtools",
                    "getfasta",
                    "-fi",
                    params.genome_fasta,
                    "-bed",
                    motif_input_bed,
                    "-name",
                    "-fo",
                    motif_input_fasta,
                ]
                env = os.environ.copy()
                for thread_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
                    env[thread_var] = str(motif_threads)
                try:
                    log_it(logfile, f"Preparing MEME motif FASTA for {item['set_name']} with {motif_peak_count} centered windows.")
                    run_meme_command(fasta_cmd, motif_timeout_seconds, env)
                    rewrite_fasta_headers_with_scores(motif_input_fasta, fasta_scores)

                    sea_dir = ensure_dir(os.path.join(out_dir, "sea"))
                    sea_tsv = os.path.join(sea_dir, "sea.tsv")
                    sea_cmd = [
                        "sea",
                        "--p",
                        motif_input_fasta,
                        "--m",
                        motif_database,
                        "--qvalue",
                        "--align",
                        "center",
                        "--text",
                        "--verbosity",
                        "1",
                    ]
                    log_it(logfile, f"Running MEME SEA for {item['set_name']}.")
                    run_meme_text_command(sea_cmd, sea_tsv, motif_timeout_seconds, env)
                    if os.path.isfile(sea_tsv):
                        sea_report_pdf = os.path.join(sea_dir, "motif_enrichment.pdf")
                        render_motif_report("SEA", regions_label_for_set(item), sea_tsv, sea_report_pdf)
                        sea_status = "OK" if count_meme_text_rows(sea_tsv) > 0 else "NO_MOTIFS"
                        sea_reason = cap_reason
                        if sea_status == "NO_MOTIFS":
                            sea_reason = "; ".join(filter(None, ["SEA produced no enriched motif rows", cap_reason]))
                        motif_summary.append([item["set_name"], "SEA", sea_status, item["peak_bed"], motif_input_fasta, sea_dir, motif_database, sea_report_pdf, sea_reason])
                    else:
                        reason = "SEA completed but did not produce sea.tsv"
                        if cap_reason:
                            reason = f"{reason}; {cap_reason}"
                        motif_summary.append([item["set_name"], "SEA", "NO_MOTIFS", item["peak_bed"], motif_input_fasta, sea_dir, motif_database, "", reason])

                    if item["set_type"].startswith("top_de_ranked") and tool_available("ame"):
                        ame_dir = ensure_dir(os.path.join(out_dir, "ame"))
                        ame_tsv = os.path.join(ame_dir, "ame.tsv")
                        ame_cmd = [
                            "ame",
                            "--method",
                            "ranksum",
                            "--scoring",
                            "max",
                            "--text",
                            motif_input_fasta,
                            motif_database,
                        ]
                        log_it(logfile, f"Running MEME AME ranked enrichment for {item['set_name']}.")
                        run_meme_text_command(ame_cmd, ame_tsv, motif_timeout_seconds, env)
                        if os.path.isfile(ame_tsv):
                            ame_report_pdf = os.path.join(ame_dir, "motif_enrichment.pdf")
                            render_motif_report("AME", regions_label_for_set(item), ame_tsv, ame_report_pdf)
                            ame_status = "OK" if count_meme_text_rows(ame_tsv) > 0 else "NO_MOTIFS"
                            ame_reason = cap_reason
                            if ame_status == "NO_MOTIFS":
                                ame_reason = "; ".join(filter(None, ["AME produced no enriched motif rows", cap_reason]))
                            motif_summary.append([item["set_name"], "AME", ame_status, item["peak_bed"], motif_input_fasta, ame_dir, motif_database, ame_report_pdf, ame_reason])
                        else:
                            reason = "AME completed but did not produce ame.tsv"
                            if cap_reason:
                                reason = f"{reason}; {cap_reason}"
                            motif_summary.append([item["set_name"], "AME", "NO_MOTIFS", item["peak_bed"], motif_input_fasta, ame_dir, motif_database, "", reason])
                    elif item["set_type"].startswith("top_de_ranked"):
                        motif_summary.append([item["set_name"], "AME", "SKIP", item["peak_bed"], motif_input_fasta, "", motif_database, "", "ame executable not found"])
                except subprocess.TimeoutExpired as motif_timeout:
                    if motif_timeout.stdout:
                        output_text = motif_timeout.stdout
                        if isinstance(output_text, bytes):
                            output_text = output_text.decode(errors="replace")
                        for line in str(output_text).splitlines():
                            log_it(logfile, line)
                    reason = f"MEME motif command exceeded {motif_timeout_seconds} seconds"
                    if cap_reason:
                        reason = f"{reason}; {cap_reason}"
                    log_it(logfile, f"MEME motif run timed out for {item['set_name']}: {reason}")
                    motif_summary.append([item["set_name"], "MEME", "TIMEOUT", item["peak_bed"], motif_input_fasta, out_dir, motif_database, "", reason])
                except Exception as motif_error:
                    log_it(logfile, f"MEME motif run failed for {item['set_name']}: {motif_error}")
                    motif_summary.append([item["set_name"], "MEME", "FAIL", item["peak_bed"], motif_input_fasta, out_dir, motif_database, "", str(motif_error)])
                write_motif_summary()

            combined_report = write_combined_motif_report(motifs_dir, motif_summary)
            write_motif_summary()
            log_it(logfile, f"Motif summary PDF: {combined_report}")
            log_it(logfile, f"Motif run summary: {motif_summary_path}")
            return motif_summary_path

        try:
            if params.thetype not in {"ATAC", "CHIP"}:
                log_it(logfile, "Analyze peaks DE step is supported for ATAC/CHIP only. Skipping.")
                write_tmp_file()
                finish_step_sample(
                    master_config["analyzepeaksde_rule_num"],
                    "aggregate",
                    "analyze_peaks_de",
                    tracking["start_time"],
                    "OK",
                )
                return
            sanity_check_dir(
                logfile,
                params.de_outputfolder,
                master_config["input_file_types"][master_config["analyzepeaksde_rule_num"] - 1],
                "step16.sanity",
            )

            analyze_de_root = ensure_dir(os.path.join(params.outputfolder, "analyze_peaks_de"))
            summary_dir = ensure_dir(os.path.join(analyze_de_root, "summary"))
            sets_dir = ensure_dir(os.path.join(analyze_de_root, "sets"))
            ensure_dir(os.path.join(sets_dir, "all_features"))
            ensure_dir(os.path.join(sets_dir, "all_gene_bodies"))
            ensure_dir(os.path.join(sets_dir, "all_bins"))
            ensure_dir(os.path.join(sets_dir, "de_significant"))
            ensure_dir(os.path.join(sets_dir, "de_up"))
            ensure_dir(os.path.join(sets_dir, "de_down"))
            ensure_dir(os.path.join(sets_dir, "de_significant_promoter"))
            ensure_dir(os.path.join(sets_dir, "de_significant_distal"))
            ensure_dir(os.path.join(sets_dir, "de_significant_domains"))
            ensure_dir(os.path.join(sets_dir, "de_up_domains"))
            ensure_dir(os.path.join(sets_dir, "de_down_domains"))
            ensure_dir(os.path.join(sets_dir, "de_significant_promoter_regions"))
            ensure_dir(os.path.join(sets_dir, "de_up_promoter_regions"))
            ensure_dir(os.path.join(sets_dir, "de_down_promoter_regions"))
            signal_dir = ensure_dir(os.path.join(analyze_de_root, "signal"))
            motifs_dir = ensure_dir(os.path.join(analyze_de_root, "motifs"))

            peak_metadata_path = os.path.join(params.de_outputfolder, "peak_metadata.tsv")
            feature_metadata_label = "Peak metadata"
            if params.thetype == "CHIP" and params.broad_mode in {"genebody", "diffuse"}:
                feature_metadata_label = "Feature metadata"
            peak_region_map = load_peak_region_map(peak_metadata_path)
            if not peak_region_map:
                log_it(
                    logfile,
                    f"{feature_metadata_label} map is empty or missing ({peak_metadata_path}); region-aware splits may be incomplete.",
                )
            promoter_map = {}
            if params.thetype == "CHIP" and params.broad_mode == "genebody":
                promoter_map = load_promoter_map(params.gtf_file)
                log_it(
                    logfile,
                    f"ChIP gene-body mode: loaded promoter annotations for {len(promoter_map)} genes from {params.gtf_file}",
                )

            de_tables = collect_de_tables()
            if not de_tables:
                feature_file_stem = "peaks"
                if params.thetype == "CHIP" and params.broad_mode == "genebody":
                    feature_file_stem = "gene_bodies"
                elif params.thetype == "CHIP" and params.broad_mode == "diffuse":
                    feature_file_stem = "bins"
                raise FileNotFoundError(
                    f"No DE {feature_file_stem} result tables found. Expected files like "
                    f"*.sig_diff_{feature_file_stem}.DESeq2.txt or *.diff_{feature_file_stem}.DESeq2.txt "
                    f"under {params.de_outputfolder}/*/*/"
                )

            set_manifest_rows = build_all_feature_sets(peak_region_map, sets_dir, params.broad_mode)
            for de_table in de_tables:
                set_manifest_rows.extend(
                    build_de_sets(
                        de_table,
                        peak_region_map,
                        sets_dir,
                        params.broad_mode,
                        promoter_map,
                    )
                )
            if not (params.thetype == "CHIP" and params.broad_mode in {"genebody", "diffuse"}):
                set_manifest_rows.extend(collect_unique_sets_from_prede(sets_dir))

            set_manifest_path = os.path.join(summary_dir, "set_manifest.tsv")
            with open(set_manifest_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    [
                        "set_name",
                        "set_type",
                        "contrast_group",
                        "contrast_label",
                        "source_table",
                        "peak_bed",
                        "peak_count",
                    ]
                )
                for item in set_manifest_rows:
                    writer.writerow(
                        [
                            item["set_name"],
                            item["set_type"],
                            item["contrast_group"],
                            item["contrast_label"],
                            item["source_table"],
                            item["peak_bed"],
                            item["peak_count"],
                        ]
                    )

            run_deeptools(set_manifest_rows, signal_dir)
            motif_summary_path = run_meme_motifs(set_manifest_rows, motifs_dir)

            report_path = os.path.join(summary_dir, "analyze_peaks_de_report.tsv")
            with open(report_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["metric", "value"])
                writer.writerow(["de_table_count", len(de_tables)])
                writer.writerow(["set_count_total", len(set_manifest_rows)])
                writer.writerow(
                    [
                        "set_count_nonempty",
                        sum(1 for item in set_manifest_rows if int(item["peak_count"]) > 0),
                    ]
                )
                writer.writerow(["set_manifest", set_manifest_path])
                writer.writerow(["signal_dir", signal_dir])
                writer.writerow(["motifs_dir", motifs_dir])
                writer.writerow(["motif_runs", motif_summary_path or "NA"])
                writer.writerow(["peak_metadata_path", peak_metadata_path])

            write_tmp_file()
            finish_step_sample(
                master_config["analyzepeaksde_rule_num"],
                "aggregate",
                "analyze_peaks_de",
                tracking["start_time"],
                "OK",
            )
        except Exception:
            finish_step_sample(
                master_config["analyzepeaksde_rule_num"],
                "aggregate",
                "analyze_peaks_de",
                tracking["start_time"],
                "FAIL",
            )
            raise
