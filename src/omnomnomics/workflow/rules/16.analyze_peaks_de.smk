# Rule 16: Analyze Peaks DE (post-DE, chromatin-focused)

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import csv
import glob
import math
import os
import re
import shlex
import shutil
import subprocess


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

        def parse_gtf_attributes(attr_string):
            attrs = {}
            for item in str(attr_string).strip().split(";"):
                item = item.strip()
                if not item or " " not in item:
                    continue
                key, value = item.split(" ", 1)
                attrs[key] = value.strip().strip('"')
            return attrs

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
            promoter_map = {}
            if not os.path.isfile(gtf_file):
                return promoter_map
            with open(gtf_file, "r", encoding="utf-8") as handle:
                for line in handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 9 or fields[2] != "gene":
                        continue
                    chrom = fields[0]
                    start = int(fields[3]) - 1
                    end = int(fields[4])
                    strand = fields[6]
                    attrs = parse_gtf_attributes(fields[8])
                    gene_name = attrs.get("gene_name") or attrs.get("gene_id") or "NA"
                    gene_id = attrs.get("gene_id") or gene_name
                    gene_label = f"{gene_name}|{gene_id}"
                    if strand == "+":
                        prom_start = max(start - promoter_upstream, 0)
                        prom_end = max(start + promoter_downstream, 0)
                    else:
                        prom_start = max(end - promoter_downstream, 0)
                        prom_end = max(end + promoter_upstream, 0)
                    if prom_end <= prom_start:
                        continue
                    promoter_map[gene_label] = (chrom, prom_start, prom_end)
            return promoter_map

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
                        promoter_coords = promoter_map.get(gene_label)
                        if promoter_coords is None:
                            continue
                        prom_chrom, prom_start, prom_end = promoter_coords
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
                            promoter_coords = promoter_map.get(gene_label)
                            if promoter_coords is None:
                                continue
                            prom_chrom, prom_start, prom_end = promoter_coords
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
            required = ["computeMatrix", "plotHeatmap", "plotProfile"]
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

            def existing_bigwig_for_sample(sample):
                candidates = [
                    os.path.join(params.bigwig_inputfolder, f"{sample}.bw"),
                    os.path.join(params.bigwig_inputfolder, f"{sample}.CPM.bw"),
                ]
                for candidate in candidates:
                    if os.path.isfile(candidate):
                        return candidate
                return None

            matrices_dir = ensure_dir(os.path.join(signal_dir, "matrices"))
            heatmaps_dir = ensure_dir(os.path.join(signal_dir, "heatmaps"))
            profiles_dir = ensure_dir(os.path.join(signal_dir, "profiles"))
            regions_dir = ensure_dir(os.path.join(signal_dir, "regions"))

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
                shell(
                    f"computeMatrix reference-point --referencePoint center -b 1000 -a 1000 "
                    f"-R {quote(signal_bed_path)} -S {' '.join(quote(x) for x in bigwigs)} "
                    f"--missingDataAsZero --binSize 100 --numberOfProcessors {deeptools_threads} -o {quote(matrix_path)}"
                )
                shell(
                    f"plotHeatmap -m {quote(matrix_path)} -out {quote(heatmap_path)} "
                    f"--whatToShow 'heatmap and colorbar' --sortRegions descend "
                    f"--plotTitle {quote(region_label)} --regionsLabel {quote(region_label)} "
                    f"--samplesLabel {' '.join(quote(label) for label in sample_labels)} "
                    f"--xAxisLabel 'distance from center (bp)' --refPointLabel center"
                )
                shell(
                    f"plotProfile -m {quote(matrix_path)} -out {quote(profile_path)} "
                    f"--perGroup --plotTitle {quote(region_label)} "
                    f"--regionsLabel {quote(region_label)} "
                    f"--samplesLabel {' '.join(quote(label) for label in sample_labels)} "
                    f"--refPointLabel center"
                )
                reason = ""
                if int(item["peak_count"]) > signal_region_count:
                    reason = f"signal input capped at first {signal_region_count} of {item['peak_count']} regions"
                signal_summary_rows.append([item["set_name"], "OK", bed_path, signal_bed_path, matrix_path, heatmap_path, reason])

            write_signal_summary()

        def run_meme_motifs(set_manifest_rows, motifs_dir):
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
                    writer.writerow(["set_name", "method", "status", "peak_bed", "input_fasta", "output_dir", "motif_database", "reason"])
                    writer.writerows(motif_summary)

            required_tools = ["bedtools", "sea"]
            missing_tools = [tool_name for tool_name in required_tools if not tool_available(tool_name)]
            if missing_tools:
                reason = "missing required executable(s): " + ", ".join(missing_tools)
                log_it(logfile, f"Skipping MEME motif analysis in analyze_peaks_de because {reason}.")
                motif_summary.append(["ALL", "MEME", "SKIP", "", "", "", "", reason])
                write_motif_summary()
                return motif_summary_path

            def resolve_meme_motif_database():
                requested = str(params.post_de_motif_database).strip()
                if requested and requested.lower() != "auto":
                    if os.path.isfile(requested):
                        return requested
                    return None

                for env_var in ("OMNOMNOMICS_MEME_MOTIF_DATABASE", "MEME_MOTIF_DATABASE", "JASPAR_MOTIF_DATABASE"):
                    env_path = os.environ.get(env_var, "").strip()
                    if env_path and os.path.isfile(env_path):
                        return env_path

                search_roots = []
                conda_prefix = os.environ.get("CONDA_PREFIX")
                if conda_prefix:
                    search_roots.extend(
                        [
                            os.path.join(conda_prefix, "share", "meme", "motif_databases"),
                            os.path.join(conda_prefix, "share", "meme", "db", "motif_databases"),
                            os.path.join(conda_prefix, "share", "meme", "db"),
                            os.path.join(conda_prefix, "share", "meme"),
                            *glob.glob(os.path.join(conda_prefix, "share", "meme-*", "db", "motif_databases")),
                            *glob.glob(os.path.join(conda_prefix, "share", "meme-*", "motif_databases")),
                        ]
                    )
                search_roots.extend(["/usr/local/share", "/usr/share"])
                patterns = [
                    "**/JASPAR*CORE*vertebrates*non-redundant*.meme",
                    "**/JASPAR*CORE*vertebrates*.meme",
                    "**/HOCOMOCO*.meme",
                    "**/*.meme",
                ]
                for root in search_roots:
                    if not root or not os.path.isdir(root):
                        continue
                    for pattern in patterns:
                        matches = sorted(glob.glob(os.path.join(root, pattern), recursive=True))
                        for candidate in matches:
                            candidate_parts = set(os.path.normpath(candidate).split(os.sep))
                            if "doc" in candidate_parts or "examples" in candidate_parts:
                                continue
                            if os.path.isfile(candidate) and os.path.getsize(candidate) > 0:
                                return candidate
                return None

            motif_database = resolve_meme_motif_database()
            if not motif_database:
                reason = (
                    "no MEME motif database found; set post_de_motif_database to a MEME-format motif file "
                    "or install the MEME Suite motif databases"
                )
                log_it(logfile, f"Skipping MEME motif analysis in analyze_peaks_de because {reason}.")
                motif_summary.append(["ALL", "MEME", "SKIP", "", "", "", "", reason])
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

            for index, item in enumerate(set_manifest_rows):
                if item["set_type"] not in eligible_set_types:
                    motif_summary.append([item["set_name"], "MEME", "SKIP", item["peak_bed"], "", "", motif_database, "set type is not configured for motif analysis"])
                    write_motif_summary()
                    continue
                if int(item["peak_count"]) < min_motif_peaks:
                    motif_summary.append([item["set_name"], "MEME", "SKIP", item["peak_bed"], "", "", motif_database, f"fewer than {min_motif_peaks} peaks"])
                    write_motif_summary()
                    continue
                if index not in selected_indices:
                    reason = f"motif set cap reached; retained top {max_motif_sets} highest-priority sets"
                    motif_summary.append([item["set_name"], "MEME", "SKIP", item["peak_bed"], "", "", motif_database, reason])
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
                    motif_summary.append([item["set_name"], "MEME", "SKIP", item["peak_bed"], "", out_dir, motif_database, f"fewer than {min_motif_peaks} valid centered windows"])
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
                    sea_cmd = [
                        "sea",
                        "--p",
                        motif_input_fasta,
                        "--m",
                        motif_database,
                        "--oc",
                        sea_dir,
                        "--qvalue",
                        "--noseqs",
                        "--align",
                        "center",
                        "--verbosity",
                        "1",
                    ]
                    log_it(logfile, f"Running MEME SEA for {item['set_name']}.")
                    run_meme_command(sea_cmd, motif_timeout_seconds, env)
                    sea_tsv = os.path.join(sea_dir, "sea.tsv")
                    if os.path.isfile(sea_tsv):
                        motif_summary.append([item["set_name"], "SEA", "OK", item["peak_bed"], motif_input_fasta, sea_dir, motif_database, cap_reason])
                    else:
                        reason = "SEA completed but did not produce sea.tsv"
                        if cap_reason:
                            reason = f"{reason}; {cap_reason}"
                        motif_summary.append([item["set_name"], "SEA", "NO_MOTIFS", item["peak_bed"], motif_input_fasta, sea_dir, motif_database, reason])

                    if item["set_type"].startswith("top_de_ranked") and tool_available("ame"):
                        ame_dir = ensure_dir(os.path.join(out_dir, "ame"))
                        ame_cmd = [
                            "ame",
                            "--oc",
                            ame_dir,
                            "--method",
                            "ranksum",
                            "--scoring",
                            "max",
                            motif_input_fasta,
                            motif_database,
                        ]
                        log_it(logfile, f"Running MEME AME ranked enrichment for {item['set_name']}.")
                        run_meme_command(ame_cmd, motif_timeout_seconds, env)
                        ame_tsv = os.path.join(ame_dir, "ame.tsv")
                        if os.path.isfile(ame_tsv):
                            motif_summary.append([item["set_name"], "AME", "OK", item["peak_bed"], motif_input_fasta, ame_dir, motif_database, cap_reason])
                        else:
                            reason = "AME completed but did not produce ame.tsv"
                            if cap_reason:
                                reason = f"{reason}; {cap_reason}"
                            motif_summary.append([item["set_name"], "AME", "NO_MOTIFS", item["peak_bed"], motif_input_fasta, ame_dir, motif_database, reason])
                    elif item["set_type"].startswith("top_de_ranked"):
                        motif_summary.append([item["set_name"], "AME", "SKIP", item["peak_bed"], motif_input_fasta, "", motif_database, "ame executable not found"])
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
                    motif_summary.append([item["set_name"], "MEME", "TIMEOUT", item["peak_bed"], motif_input_fasta, out_dir, motif_database, reason])
                except Exception as motif_error:
                    log_it(logfile, f"MEME motif run failed for {item['set_name']}: {motif_error}")
                    motif_summary.append([item["set_name"], "MEME", "FAIL", item["peak_bed"], motif_input_fasta, out_dir, motif_database, str(motif_error)])
                write_motif_summary()

            write_motif_summary()
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
