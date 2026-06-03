# Rule 13: Peak QC

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import csv
import glob
import os
import shutil
import subprocess
import tempfile
from shlex import quote


def peak_qc_input(_wildcards):
    input_files = [
        f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num'] - 1]}/extra_{master_config['callpeaks_rule_num']}.tmp"
    ]
    input_files.append(
        os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "annotation", "genes.gtf")
    )
    if config["THETYPE"] == "ATAC":
        input_files.extend(
            f"{experiment_dir}/{master_config['input_folders'][master_config['peakqc_rule_num'] - 1][0]}/{sample}.sorted.dups_marked.filtered.bam"
            for sample in samples2
        )
    elif config["THETYPE"] == "CHIP":
        input_files.extend(
            f"{experiment_dir}/{master_config['input_folders'][master_config['peakqc_rule_num'] - 1][0]}/{sample}.filtered.bam"
            for sample in samples2
        )
    return input_files


rule peak_qc:
    input:
        peak_qc_input
    output:
        extra=f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/extra_{master_config['peakqc_rule_num']}.tmp",
        union_annotation=f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/peak_qc/peak_annotations/{config['THETYPE'].lower()}.all_groups.merged_peaks.annotated.bed"
    params:
        thetype=config["THETYPE"],
        spp_gate_mode=str(config.get("SPP_GATE", "warn")).strip().lower(),
        bam_inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['peakqc_rule_num'] - 1][0]}",
        peak_outputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['peakqc_rule_num'] - 1][1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}",
        gtf_file=os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "annotation", "genes.gtf"),
    threads:
        Threads_Per_Rule[str(master_config["peakqc_rule_num"])]
    resources:
        mem_mb=Memory_Per_Rule[str(master_config["peakqc_rule_num"])],
        partition=master_config["partition"],
        runtime=Runtime_Per_Rule[str(master_config["peakqc_rule_num"])]
    run:
        tracking = begin_step_sample(master_config["peakqc_rule_num"], "aggregate", "peak_qc")
        log_once(logfile, "step13.header", "Calculating peak-level QC metrics (FRiP and library QC)...", f"EXECUTING STEP {master_config['peakqc_rule_num']}")
        log_once(logfile, "step13.inputfolder", f"Input folders: {params.bam_inputfolder} and {params.peak_outputfolder}")
        log_once(logfile, "step13.outputfolder", f"Output folder: {params.outputfolder}")

        def ensure_peak_qc_dir(outputfolder):
            qc_dir = os.path.join(outputfolder, "peak_qc")
            os.makedirs(qc_dir, exist_ok=True)
            return qc_dir

        def resolve_encode_blacklist(genome_name, outputfolder):
            try:
                import genomepy
            except ImportError as exc:
                raise RuntimeError(
                    "ATAC blacklist filtering requires genomepy in the active environment."
                ) from exc

            cache_root = os.path.join(outputfolder, "peak_qc", ".genomepy_blacklist_cache")
            os.makedirs(cache_root, exist_ok=True)
            local_name = f"{genome_name}.omnomnomics.blacklist"
            genomepy.manage_plugins("enable", ["blacklist"])
            genomepy.install_genome(
                genome_name,
                provider="UCSC",
                genomes_dir=cache_root,
                localname=local_name,
                annotation=False,
                threads=1,
                force=False,
            )
            genome_dir = os.path.join(cache_root, local_name)
            candidates = sorted(
                glob.glob(os.path.join(genome_dir, "**", "*blacklist*.bed*"), recursive=True)
            )
            if not candidates:
                raise FileNotFoundError(
                    f"No blacklist BED found after genomepy install for '{genome_name}' in {genome_dir}."
                )
            return candidates[0]

        def filter_peak_bed_file(in_bed, out_bed, keep_standard=True, drop_chrm=True, blacklist_bed=None):
            work = in_bed
            with tempfile.TemporaryDirectory(prefix="omnomnomics_peak_filter_") as tmpdir:
                if keep_standard:
                    std_bed = os.path.join(tmpdir, "std.bed")
                    shell(
                        f"awk 'BEGIN{{OFS=\"\\t\"}} $1 ~ /^chr([0-9]+|X|Y)$/ {{print $0}}' {quote(work)} > {quote(std_bed)}"
                    )
                    work = std_bed
                if drop_chrm:
                    chrm_bed = os.path.join(tmpdir, "nochrm.bed")
                    shell(
                        f"awk 'BEGIN{{OFS=\"\\t\"}} $1 != \"chrM\" && $1 != \"MT\" && $1 != \"chrMT\" {{print $0}}' {quote(work)} > {quote(chrm_bed)}"
                    )
                    work = chrm_bed
                if blacklist_bed and os.path.exists(blacklist_bed):
                    no_blacklist_bed = os.path.join(tmpdir, "noblacklist.bed")
                    shell(
                        f"bedtools intersect -v -a {quote(work)} -b {quote(blacklist_bed)} > {quote(no_blacklist_bed)}"
                    )
                    work = no_blacklist_bed
                shell(f"sort -k1,1 -k2,2n -k3,3n {quote(work)} > {quote(out_bed)}")

        def build_filtered_peak_sets(peak_outputfolder, outputfolder, genome_name):
            filtered_dir = os.path.join(ensure_peak_qc_dir(outputfolder), "filtered_peaks")
            os.makedirs(filtered_dir, exist_ok=True)
            blacklist_bed = resolve_encode_blacklist(genome_name, outputfolder)

            source_peak_beds = sorted(glob.glob(os.path.join(peak_outputfolder, "*.bed")))
            if not source_peak_beds:
                raise FileNotFoundError(f"No peak BED files found in {peak_outputfolder}")

            produced = []
            for src in source_peak_beds:
                dst = os.path.join(filtered_dir, os.path.basename(src))
                filter_peak_bed_file(
                    in_bed=src,
                    out_bed=dst,
                    keep_standard=True,
                    drop_chrm=True,
                    blacklist_bed=blacklist_bed,
                )
                produced.append(dst)

            summary_path = os.path.join(filtered_dir, "peak_filtering_summary.tsv")
            with open(summary_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["source_bed", "filtered_bed", "source_n", "filtered_n"])
                for src in source_peak_beds:
                    dst = os.path.join(filtered_dir, os.path.basename(src))
                    source_n = int(
                        subprocess.check_output(
                            f"wc -l {quote(src)} | awk '{{print $1}}'",
                            shell=True,
                            executable="/bin/bash",
                            text=True,
                        ).strip()
                    )
                    filtered_n = int(
                        subprocess.check_output(
                            f"wc -l {quote(dst)} | awk '{{print $1}}'",
                            shell=True,
                            executable="/bin/bash",
                            text=True,
                        ).strip()
                    )
                    writer.writerow([src, dst, source_n, filtered_n])

            return filtered_dir, blacklist_bed, produced, summary_path

        def parse_gtf_attributes(attr_string):
            attrs = {}
            for item in attr_string.strip().split(";"):
                item = item.strip()
                if not item or " " not in item:
                    continue
                key, value = item.split(" ", 1)
                attrs[key] = value.strip().strip('"')
            return attrs

        def build_gene_annotation_beds(gtf_file, work_dir, promoter_upstream=1000, promoter_downstream=1000):
            genes_bed = os.path.join(work_dir, "genes.unsorted.bed")
            exons_bed = os.path.join(work_dir, "exons.unsorted.bed")
            promoters_bed = os.path.join(work_dir, "promoters.unsorted.bed")
            genes_sorted = os.path.join(work_dir, "genes.sorted.bed")
            exons_sorted = os.path.join(work_dir, "exons.sorted.bed")
            promoters_sorted = os.path.join(work_dir, "promoters.sorted.bed")
            introns_sorted = os.path.join(work_dir, "introns.sorted.bed")

            with open(gtf_file) as gtf_handle, open(genes_bed, "w") as genes_out, open(exons_bed, "w") as exons_out, open(promoters_bed, "w") as promoter_out:
                for line in gtf_handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 9:
                        continue
                    chrom, _, feature, start_s, end_s, _, strand, _, attrs_raw = fields
                    start = int(start_s)
                    end = int(end_s)
                    if end <= start:
                        continue
                    attrs = parse_gtf_attributes(attrs_raw)
                    gene_name = attrs.get("gene_name") or attrs.get("gene_id") or "NA"
                    gene_id = attrs.get("gene_id") or gene_name
                    gene_label = f"{gene_name}|{gene_id}"
                    bed_start = max(start - 1, 0)
                    bed_end = end

                    if feature == "gene":
                        genes_out.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_label}\t0\t{strand}\n")
                        if strand == "+":
                            prom_start = max(bed_start - promoter_upstream, 0)
                            prom_end = max(bed_start + promoter_downstream, 0)
                        else:
                            prom_start = max(bed_end - promoter_downstream, 0)
                            prom_end = max(bed_end + promoter_upstream, 0)
                        if prom_end > prom_start:
                            promoter_out.write(f"{chrom}\t{prom_start}\t{prom_end}\t{gene_label}\t0\t{strand}\n")
                    elif feature == "exon":
                        exons_out.write(f"{chrom}\t{bed_start}\t{bed_end}\t{gene_label}\t0\t{strand}\n")

            shell(f"sort -k1,1 -k2,2n -k3,3n {quote(genes_bed)} > {quote(genes_sorted)}")
            shell(f"sort -k1,1 -k2,2n -k3,3n {quote(exons_bed)} > {quote(exons_sorted)}")
            shell(f"sort -k1,1 -k2,2n -k3,3n {quote(promoters_bed)} > {quote(promoters_sorted)}")
            shell(
                f"bedtools subtract -a {quote(genes_sorted)} -b {quote(exons_sorted)} "
                f"| sort -k1,1 -k2,2n -k3,3n > {quote(introns_sorted)}"
            )

            return {
                "genes": genes_sorted,
                "promoters": promoters_sorted,
                "exons": exons_sorted,
                "introns": introns_sorted,
            }

        def load_peak_records(peak_bed):
            peaks = []
            with open(peak_bed) as handle:
                for idx, line in enumerate(handle, start=1):
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    chrom, start_s, end_s = fields[0:3]
                    start = int(start_s)
                    end = int(end_s)
                    if end <= start:
                        continue
                    peaks.append({
                        "peak_id": f"peak_{idx}",
                        "chrom": chrom,
                        "start": start,
                        "end": end,
                        "underscore": f"{chrom}_{start}_{end}",
                    })
            return peaks

        def annotate_peak_regions(peak_bed, annotation_beds, outputfolder, assay, group, peak_set):
            qc_dir = ensure_peak_qc_dir(outputfolder)
            annotation_dir = os.path.join(qc_dir, "peak_annotations")
            os.makedirs(annotation_dir, exist_ok=True)
            safe_peak_set = peak_set.replace("/", "_")
            out_bed = os.path.join(annotation_dir, f"{assay.lower()}.{safe_peak_set}.annotated.bed")
            summary_tsv = os.path.join(annotation_dir, f"{assay.lower()}.{safe_peak_set}.genomic_distribution.tsv")

            with tempfile.TemporaryDirectory(prefix="omnomnomics_peak_annot_") as tmpdir:
                peaks_indexed = os.path.join(tmpdir, "peaks.indexed.bed")
                with open(peak_bed) as inp, open(peaks_indexed, "w") as out:
                    for idx, line in enumerate(inp, start=1):
                        if not line.strip():
                            continue
                        f = line.rstrip("\n").split("\t")
                        out.write(f"{f[0]}\t{f[1]}\t{f[2]}\tpeak_{idx}\n")

                overlap_maps = {}
                for region_name, bed_path in (
                    ("promoter", annotation_beds["promoters"]),
                    ("exon", annotation_beds["exons"]),
                    ("intron", annotation_beds["introns"]),
                ):
                    overlap_file = os.path.join(tmpdir, f"overlap.{region_name}.tsv")
                    shell(
                        f"bedtools intersect -a {quote(peaks_indexed)} -b {quote(bed_path)} -wa -wb > {quote(overlap_file)}"
                    )
                    region_map = {}
                    with open(overlap_file) as handle:
                        for line in handle:
                            if not line.strip():
                                continue
                            cols = line.rstrip("\n").split("\t")
                            peak_id = cols[3]
                            gene_label = cols[7]
                            region_map.setdefault(peak_id, set()).add(gene_label)
                    overlap_maps[region_name] = region_map

                closest_file = os.path.join(tmpdir, "closest.tsv")
                shell(
                    f"bedtools closest -a {quote(peaks_indexed)} -b {quote(annotation_beds['genes'])} -d -t first > {quote(closest_file)}"
                )
                closest_map = {}
                with open(closest_file) as handle:
                    for line in handle:
                        if not line.strip():
                            continue
                        cols = line.rstrip("\n").split("\t")
                        peak_id = cols[3]
                        gene_label = cols[7] if len(cols) > 7 else "NA"
                        distance = cols[-1] if cols else "NA"
                        closest_map[peak_id] = (gene_label, distance)

            peaks = load_peak_records(peak_bed)
            region_counts = {"promoter": 0, "exon": 0, "intron": 0, "intergenic": 0}
            with open(out_bed, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "#chrom",
                    "start",
                    "end",
                    "underscore",
                    "peak_id",
                    "score",
                    "strand",
                    "genomic_region",
                    "assigned_genes",
                    "nearest_gene",
                    "distance_to_nearest_gene_bp",
                    "assay",
                    "group",
                    "peak_set",
                    "width_bp",
                ])
                for peak in peaks:
                    peak_id = peak["peak_id"]
                    assigned_region = "intergenic"
                    assigned_gene_set = set()
                    for region in ("promoter", "exon", "intron"):
                        genes_here = overlap_maps.get(region, {}).get(peak_id, set())
                        if genes_here:
                            assigned_region = region
                            assigned_gene_set = genes_here
                            break
                    nearest_gene, nearest_distance = closest_map.get(peak_id, ("NA", "NA"))
                    if assigned_region != "intergenic":
                        nearest_distance = "0"
                        if assigned_gene_set:
                            nearest_gene = sorted(assigned_gene_set)[0]
                    assigned_genes = ",".join(sorted(assigned_gene_set)) if assigned_gene_set else "NA"
                    region_counts[assigned_region] += 1
                    writer.writerow([
                        peak["chrom"],
                        peak["start"],
                        peak["end"],
                        peak["underscore"],
                        peak_id,
                        "0",
                        ".",
                        assigned_region,
                        assigned_genes,
                        nearest_gene,
                        nearest_distance,
                        assay,
                        group,
                        peak_set,
                        peak["end"] - peak["start"],
                    ])

            total = sum(region_counts.values())
            with open(summary_tsv, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["assay", "group", "peak_set", "genomic_region", "peak_count", "fraction"])
                for region in ("promoter", "exon", "intron", "intergenic"):
                    count = region_counts[region]
                    fraction = (count / total) if total else 0.0
                    writer.writerow([assay, group, peak_set, region, count, f"{fraction:.6f}"])

            return out_bed, summary_tsv

        def log_peak_qc_versions():
            bedtools_version = subprocess.check_output(["bedtools", "--version"], stderr=subprocess.STDOUT)
            samtools_version = subprocess.check_output(["samtools", "--version"], stderr=subprocess.STDOUT).decode("utf-8").splitlines()[0]
            log_once(logfile, "step13.bedtools_version", "\n" + bedtools_version.decode("utf-8"), "BEDTOOLS VERSION")
            log_once(logfile, "step13.samtools_version", "\n" + samtools_version + "\n", "SAMTOOLS VERSION")
            if shutil.which("run_spp.R"):
                log_once(logfile, "step13.phantompeakqualtools", f"\nrun_spp.R: {shutil.which('run_spp.R')}\n", "PHANTOMPEAKQUALTOOLS")

        def calculate_library_complexity_metrics(bam_path):
            with tempfile.TemporaryDirectory(prefix="omnomnomics_library_qc_") as tmpdir:
                counts_path = os.path.join(tmpdir, "fragment_counts.tsv")
                if config["PAIRED"]:
                    complexity_command = (
                        f"bedtools bamtobed -bedpe -i {quote(bam_path)} | "
                        "awk 'BEGIN{OFS=\"\\t\"} $1==$4 && $1!=\".\" {"
                        "start=($2<$5?$2:$5); end=($3>$6?$3:$6); print $1,start,end"
                        "}' | "
                        f"sort -T {quote(tmpdir)} -k1,1 -k2,2n -k3,3n | uniq -c > {quote(counts_path)}"
                    )
                else:
                    complexity_command = (
                        f"bedtools bamtobed -i {quote(bam_path)} | "
                        "awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3,$6}' | "
                        f"sort -T {quote(tmpdir)} -k1,1 -k2,2n -k3,3n -k4,4 | uniq -c > {quote(counts_path)}"
                    )
                shell(complexity_command)

                total_reads = 0
                distinct_reads = 0
                one_read = 0
                two_reads = 0
                with open(counts_path) as handle:
                    for line in handle:
                        if not line.strip():
                            continue
                        count = int(line.strip().split()[0])
                        total_reads += count
                        distinct_reads += 1
                        if count == 1:
                            one_read += 1
                        elif count == 2:
                            two_reads += 1

            nrf = (distinct_reads / total_reads) if total_reads else 0.0
            pbc1 = (one_read / distinct_reads) if distinct_reads else 0.0
            pbc2 = (one_read / two_reads) if two_reads else ""
            return {
                "total_reads": total_reads,
                "distinct_reads": distinct_reads,
                "one_read_sites": one_read,
                "two_read_sites": two_reads,
                "nrf": nrf,
                "pbc1": pbc1,
                "pbc2": pbc2,
            }

        def calculate_cross_correlation_metrics(bam_path, output_dir, sample_name):
            crosscorr_prefix = os.path.join(output_dir, sample_name)
            crosscorr_table = f"{crosscorr_prefix}.cross_correlation.tsv"
            crosscorr_pdf = f"{crosscorr_prefix}.cross_correlation.pdf"
            if not shutil.which("run_spp.R"):
                return {
                    "est_frag_len": "",
                    "nsc": "",
                    "rsc": "",
                    "crosscorr_table": "",
                    "crosscorr_pdf": "",
                }

            run_spp_command = (
                f"run_spp.R -c={quote(bam_path)} -savp={quote(crosscorr_pdf)} "
                f"-out={quote(crosscorr_table)} -p={threads}"
            )
            log_it(logfile, run_spp_command, "PHANTOMPEAKQUALTOOLS COMMAND")
            shell(run_spp_command)

            metrics_line = ""
            with open(crosscorr_table) as handle:
                for line in handle:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        metrics_line = line
                        break

            parts = [part.strip() for part in metrics_line.split("\t")] if "\t" in metrics_line else [part.strip() for part in metrics_line.split(",")]
            est_frag_len = parts[2] if len(parts) > 2 else ""
            nsc = parts[8] if len(parts) > 8 else ""
            rsc = parts[9] if len(parts) > 9 else ""
            return {
                "est_frag_len": est_frag_len,
                "nsc": nsc,
                "rsc": rsc,
                "crosscorr_table": crosscorr_table,
                "crosscorr_pdf": crosscorr_pdf,
            }

        def calculate_peak_qc_metrics(peak_bed, bam_files):
            with tempfile.TemporaryDirectory(prefix="omnomnomics_peak_qc_") as tmpdir:
                merged_bed = os.path.join(tmpdir, "merged_peaks.bed")
                multicov_output = os.path.join(tmpdir, "multicov.tsv")
                shell(
                    f"sort -k1,1 -k2,2n -k3,3n {quote(peak_bed)} | "
                    f"bedtools merge -i - > {quote(merged_bed)}"
                )
                shell(
                    f"bedtools multicov -bed {quote(merged_bed)} -bams {' '.join(quote(path) for path in bam_files)} "
                    f"> {quote(multicov_output)}"
                )

                peak_count = 0
                total_peak_bp = 0
                reads_in_peaks = 0
                with open(multicov_output, newline="") as handle:
                    reader = csv.reader(handle, delimiter="\t")
                    for row in reader:
                        peak_count += 1
                        total_peak_bp += int(row[2]) - int(row[1])
                        reads_in_peaks += sum(int(value) for value in row[3:])

                total_aligned_reads = 0
                for bam_path in bam_files:
                    total_aligned_reads += int(
                        subprocess.check_output(
                            ["samtools", "view", "-c", "-F", "260", bam_path],
                            stderr=subprocess.STDOUT,
                        ).decode("utf-8").strip()
                    )

            frip = (reads_in_peaks / total_aligned_reads) if total_aligned_reads else 0.0
            return {
                "peak_count": peak_count,
                "total_peak_bp": total_peak_bp,
                "reads_in_peaks": reads_in_peaks,
                "total_aligned_reads": total_aligned_reads,
                "frip": frip,
            }

        def write_sample_qc_outputs(outputfolder, thetype, qc_rows):
            if not qc_rows:
                return None

            qc_dir = ensure_peak_qc_dir(outputfolder)
            qc_table = os.path.join(qc_dir, f"{thetype.lower()}.sample_qc_metrics.tsv")
            with open(qc_table, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "assay",
                    "sample",
                    "bam_file",
                    "total_reads",
                    "distinct_reads",
                    "one_read_sites",
                    "two_read_sites",
                    "nrf",
                    "pbc1",
                    "pbc2",
                    "est_frag_len",
                    "nsc",
                    "rsc",
                    "crosscorr_table",
                    "crosscorr_pdf",
                ])
                for row in qc_rows:
                    writer.writerow([
                        row["assay"],
                        row["sample"],
                        row["bam_file"],
                        row["total_reads"],
                        row["distinct_reads"],
                        row["one_read_sites"],
                        row["two_read_sites"],
                        f"{row['nrf']:.6f}",
                        f"{row['pbc1']:.6f}",
                        row["pbc2"] if row["pbc2"] == "" else f"{float(row['pbc2']):.6f}",
                        row["est_frag_len"],
                        row["nsc"],
                        row["rsc"],
                        row["crosscorr_table"],
                        row["crosscorr_pdf"],
                    ])
            return qc_table

        def write_spp_qc_summary(outputfolder, thetype, sample_rows):
            if not sample_rows:
                return None, []

            qc_dir = ensure_peak_qc_dir(outputfolder)
            summary_path = os.path.join(qc_dir, f"{thetype.lower()}.spp_qc_summary.tsv")
            flagged_rows = []
            with open(summary_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    [
                        "assay",
                        "sample",
                        "nsc",
                        "rsc",
                        "nsc_flag",
                        "rsc_flag",
                        "overall_flag",
                        "rule",
                    ]
                )
                for row in sample_rows:
                    nsc_raw = str(row.get("nsc", "")).strip()
                    rsc_raw = str(row.get("rsc", "")).strip()
                    nsc_val = None
                    rsc_val = None
                    try:
                        nsc_val = float(nsc_raw) if nsc_raw not in ("", "NA") else None
                    except ValueError:
                        nsc_val = None
                    try:
                        rsc_val = float(rsc_raw) if rsc_raw not in ("", "NA") else None
                    except ValueError:
                        rsc_val = None

                    if nsc_val is None:
                        nsc_flag = "missing"
                    elif nsc_val >= 1.05:
                        nsc_flag = "pass"
                    else:
                        nsc_flag = "warn"

                    if rsc_val is None:
                        rsc_flag = "missing"
                    elif rsc_val >= 0.8:
                        rsc_flag = "pass"
                    elif rsc_val >= 0.5:
                        rsc_flag = "warn"
                    else:
                        rsc_flag = "fail"

                    if nsc_flag == "missing" or rsc_flag == "missing":
                        overall = "missing"
                    elif rsc_flag == "fail":
                        overall = "fail"
                    elif nsc_flag == "warn" or rsc_flag == "warn":
                        overall = "warn"
                    else:
                        overall = "pass"

                    writer.writerow(
                        [
                            thetype,
                            row["sample"],
                            nsc_raw,
                            rsc_raw,
                            nsc_flag,
                            rsc_flag,
                            overall,
                            "ENCODE-style NSC>=1.05, RSC>=0.8 (warn: 0.5<=RSC<0.8)",
                        ]
                    )
                    if overall in {"warn", "fail", "missing"}:
                        flagged_rows.append(
                            {
                                "sample": row["sample"],
                                "nsc": nsc_raw,
                                "rsc": rsc_raw,
                                "overall": overall,
                            }
                        )
            return summary_path, flagged_rows

        def write_spp_dropped_samples(outputfolder, flagged_rows):
            qc_dir = ensure_peak_qc_dir(outputfolder)
            spp_dir = os.path.join(qc_dir, "spp_qc")
            os.makedirs(spp_dir, exist_ok=True)
            dropped_path = os.path.join(spp_dir, "dropped_samples.tsv")
            with open(dropped_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["sample_id", "nsc", "rsc", "status"])
                for row in flagged_rows:
                    writer.writerow([row["sample"], row["nsc"], row["rsc"], row["overall"]])
            return dropped_path

        def write_peak_qc_outputs(outputfolder, thetype, qc_rows):
            if not qc_rows:
                return None

            qc_dir = ensure_peak_qc_dir(outputfolder)
            qc_table = os.path.join(qc_dir, f"{thetype.lower()}.peak_qc_metrics.tsv")
            with open(qc_table, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "assay",
                    "group",
                    "peak_set",
                    "peak_file",
                    "bam_count",
                    "peak_count",
                    "total_peak_bp",
                    "reads_in_peaks",
                    "total_aligned_reads",
                    "frip",
                ])
                for row in qc_rows:
                    writer.writerow([
                        row["assay"],
                        row["group"],
                        row["peak_set"],
                        row["peak_file"],
                        row["bam_count"],
                        row["peak_count"],
                        row["total_peak_bp"],
                        row["reads_in_peaks"],
                        row["total_aligned_reads"],
                        f"{row['frip']:.6f}",
                    ])
            return qc_table

        def write_peak_qc_summary_pdf(outputfolder, thetype, sample_rows, peak_rows):
            if not sample_rows and not peak_rows:
                return

            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages

            log_once(logfile, "step13.matplotlib_version", f"\nmatplotlib {matplotlib.__version__}\n", "MATPLOTLIB VERSION")

            qc_dir = ensure_peak_qc_dir(outputfolder)
            pdf_path = os.path.join(qc_dir, f"{thetype.lower()}.peak_qc_summary.pdf")

            plt.style.use("seaborn-v0_8-whitegrid")
            with PdfPages(pdf_path) as pdf:
                if peak_rows:
                    peak_labels = [row["peak_set"] for row in peak_rows]
                    frip_values = [float(row["frip"]) for row in peak_rows]
                    peak_counts = [float(row["peak_count"]) for row in peak_rows]
                    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
                    axes[0].bar(peak_labels, frip_values, color="#4C78A8")
                    axes[0].set_title(f"{thetype} peak QC: FRiP")
                    axes[0].set_ylabel("FRiP")
                    axes[0].tick_params(axis="x", rotation=45, labelsize=9)
                    axes[1].bar(peak_labels, peak_counts, color="#D95F02")
                    axes[1].set_title(f"{thetype} peak QC: peak count")
                    axes[1].set_ylabel("Peaks")
                    axes[1].tick_params(axis="x", rotation=45, labelsize=9)
                    fig.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

                if sample_rows:
                    sample_labels = [row["sample"] for row in sample_rows]
                    nrf_values = [float(row["nrf"]) for row in sample_rows]
                    pbc1_values = [float(row["pbc1"]) for row in sample_rows]
                    pbc2_values = [float(row["pbc2"]) if row["pbc2"] != "" else 0.0 for row in sample_rows]
                    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
                    axes[0].bar(sample_labels, nrf_values, color="#2A9D8F")
                    axes[0].set_title(f"{thetype} library complexity: NRF")
                    axes[0].set_ylabel("NRF")
                    axes[0].tick_params(axis="x", rotation=45, labelsize=9)
                    axes[1].bar(sample_labels, pbc1_values, color="#E9C46A")
                    axes[1].set_title(f"{thetype} library complexity: PBC1")
                    axes[1].set_ylabel("PBC1")
                    axes[1].tick_params(axis="x", rotation=45, labelsize=9)
                    axes[2].bar(sample_labels, pbc2_values, color="#A44A3F")
                    axes[2].set_title(f"{thetype} library complexity: PBC2")
                    axes[2].set_ylabel("PBC2")
                    axes[2].tick_params(axis="x", rotation=45, labelsize=9)
                    fig.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

                    nsc_values = [float(row["nsc"]) for row in sample_rows if row["nsc"] not in ("", "NA")]
                    rsc_values = [float(row["rsc"]) for row in sample_rows if row["rsc"] not in ("", "NA")]
                    if nsc_values and rsc_values:
                        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
                        axes[0].bar(sample_labels, [float(row["nsc"]) if row["nsc"] not in ("", "NA") else 0.0 for row in sample_rows], color="#6A4C93")
                        axes[0].set_title(f"{thetype} enrichment: NSC")
                        axes[0].set_ylabel("NSC")
                        axes[0].tick_params(axis="x", rotation=45, labelsize=9)
                        axes[1].bar(sample_labels, [float(row["rsc"]) if row["rsc"] not in ("", "NA") else 0.0 for row in sample_rows], color="#F28482")
                        axes[1].set_title(f"{thetype} enrichment: RSC")
                        axes[1].set_ylabel("RSC")
                        axes[1].tick_params(axis="x", rotation=45, labelsize=9)
                        fig.tight_layout()
                        pdf.savefig(fig, bbox_inches="tight")
                        plt.close(fig)

            log_it(logfile, f"Peak QC summary PDF: {pdf_path}")

        def collect_sample_qc_rows(inputfolder, thetype):
            sample_rows = []
            qc_dir = ensure_peak_qc_dir(params.outputfolder)
            bam_suffix = ".sorted.dups_marked.filtered.bam" if thetype == "ATAC" else ".filtered.bam"
            for sample in samples2:
                bam_path = os.path.join(inputfolder, f"{sample}{bam_suffix}")
                if not os.path.exists(bam_path):
                    continue
                log_it(logfile, f"Calculating library complexity metrics for {sample}...")
                complexity_metrics = calculate_library_complexity_metrics(bam_path)
                log_it(logfile, f"Calculating cross-correlation metrics for {sample}...")
                crosscorr_metrics = calculate_cross_correlation_metrics(bam_path, qc_dir, sample)
                sample_rows.append({
                    "assay": thetype,
                    "sample": sample_id_for_sample(sample),
                    "bam_file": os.path.basename(bam_path),
                    **complexity_metrics,
                    **crosscorr_metrics,
                })
            return sample_rows

        def grouped_bams_by_sample_type(inputfolder, bam_suffix):
            grouped_bams = {}
            for sample in samples2:
                bam_path = os.path.join(inputfolder, f"{sample}{bam_suffix}")
                if not os.path.exists(bam_path):
                    continue
                group_name = sample_type_for_sample(sample)
                grouped_bams.setdefault(group_name, []).append((sample, bam_path))
            return grouped_bams

        def write_tmp_file(outputfolder):
            shell(f"""echo "necessity file for peak QC. can delete this." > {outputfolder}/extra_{master_config['peakqc_rule_num']}.tmp""")

        try:
            if params.thetype == "RNA":
                log_it(logfile, "Not a ChIP- or ATAC-seq experiment, skipping this step...")
                write_tmp_file(params.outputfolder)
                finish_step_sample(master_config["peakqc_rule_num"], "aggregate", "peak_qc", tracking["start_time"], "OK")
                return

            if not os.path.isfile(params.gtf_file):
                raise FileNotFoundError(f"Genome annotation GTF not found: {params.gtf_file}")

            log_peak_qc_versions()
            assay = params.thetype
            bam_suffix = ".sorted.dups_marked.filtered.bam" if assay == "ATAC" else ".filtered.bam"
            grouped_bams = grouped_bams_by_sample_type(params.bam_inputfolder, bam_suffix)
            peak_qc_rows = []
            generated_peak_annotation_files = []
            generated_distribution_files = []
            peak_bed_source_folder = params.peak_outputfolder
            if assay == "ATAC":
                filtered_dir, blacklist_bed, filtered_peak_beds, filtering_summary = build_filtered_peak_sets(
                    params.peak_outputfolder,
                    params.outputfolder,
                    config["MYGENOME"],
                )
                peak_bed_source_folder = filtered_dir
                log_it(logfile, f"ATAC blacklist BED: {blacklist_bed}")
                log_it(logfile, f"ATAC filtered peak BEDs: {len(filtered_peak_beds)} files in {filtered_dir}")
                log_it(logfile, f"ATAC peak filtering summary: {filtering_summary}")
            with tempfile.TemporaryDirectory(prefix="omnomnomics_gene_anno_") as anno_tmpdir:
                annotation_beds = build_gene_annotation_beds(params.gtf_file, anno_tmpdir)
                for group in sorted(grouped_bams):
                    bams = [bam_path for _, bam_path in grouped_bams[group]]
                    if assay == "ATAC":
                        preferred = os.path.join(peak_bed_source_folder, f"{group}.MACS3.optimized.bed")
                        fallback = os.path.join(peak_bed_source_folder, f"{group}.MACS3.q-0p01.shiftm100.ext200.group_peaks.bed")
                        candidate_peak_sets = [preferred, fallback]
                    else:
                        preferred = os.path.join(peak_bed_source_folder, f"{group}.MACS3.optimized.bed")
                        fallback_glob = sorted(glob.glob(os.path.join(peak_bed_source_folder, f"{group}.MACS3*.bed")))
                        candidate_peak_sets = [preferred, *fallback_glob]
                    for peak_bed in candidate_peak_sets:
                        if not os.path.exists(peak_bed):
                            continue
                        metrics = calculate_peak_qc_metrics(peak_bed, bams)
                        peak_qc_rows.append({
                            "assay": assay,
                            "group": group,
                            "peak_set": os.path.basename(peak_bed).replace(".bed", ""),
                            "peak_file": os.path.basename(peak_bed),
                            "bam_count": len(bams),
                            **metrics,
                        })
                        peak_set_name = os.path.basename(peak_bed).replace(".bed", "")
                        if assay != "ATAC":
                            ann_file, dist_file = annotate_peak_regions(
                                peak_bed=peak_bed,
                                annotation_beds=annotation_beds,
                                outputfolder=params.outputfolder,
                                assay=assay,
                                group=group,
                                peak_set=peak_set_name,
                            )
                            generated_peak_annotation_files.append(ann_file)
                            generated_distribution_files.append(dist_file)

                if assay in {"ATAC", "CHIP"}:
                    union_peak_bed = os.path.join(peak_bed_source_folder, "all_groups.merged_peaks.bed")
                    if not os.path.exists(union_peak_bed):
                        raise FileNotFoundError(
                            f"Expected merged peak union not found for annotation: {union_peak_bed}"
                        )
                    ann_file, dist_file = annotate_peak_regions(
                        peak_bed=union_peak_bed,
                        annotation_beds=annotation_beds,
                        outputfolder=params.outputfolder,
                        assay=assay,
                        group="all_groups",
                        peak_set="all_groups.merged_peaks",
                    )
                    generated_peak_annotation_files.append(ann_file)
                    generated_distribution_files.append(dist_file)

            sample_qc_rows = collect_sample_qc_rows(params.bam_inputfolder, assay)
            peak_qc_table = write_peak_qc_outputs(params.outputfolder, assay, peak_qc_rows)
            sample_qc_table = write_sample_qc_outputs(params.outputfolder, assay, sample_qc_rows)
            spp_qc_table, spp_flagged = write_spp_qc_summary(params.outputfolder, assay, sample_qc_rows)
            if peak_qc_table:
                log_it(logfile, f"Peak QC metrics: {peak_qc_table}")
            if sample_qc_table:
                log_it(logfile, f"Sample QC metrics: {sample_qc_table}")
            if spp_qc_table:
                log_it(logfile, f"SPP cross-correlation QC summary: {spp_qc_table}")
            spp_drop_file = os.path.join(ensure_peak_qc_dir(params.outputfolder), "spp_qc", "dropped_samples.tsv")
            if os.path.exists(spp_drop_file):
                os.remove(spp_drop_file)
            spp_gate = str(params.spp_gate_mode).strip().lower()
            if spp_gate not in {"none", "warn", "drop", "strict"}:
                spp_gate = "warn"
            log_it(logfile, f"SPP gate mode: {spp_gate}")
            if spp_gate != "none" and spp_flagged:
                flagged_text = "; ".join(
                    f"{item['sample']} (NSC={item['nsc']}, RSC={item['rsc']}, status={item['overall']})"
                    for item in spp_flagged
                )
                if spp_gate == "warn":
                    log_it(logfile, f"SPP gate warning: flagged samples detected: {flagged_text}")
                elif spp_gate == "drop":
                    dropped_path = write_spp_dropped_samples(params.outputfolder, spp_flagged)
                    log_it(logfile, f"SPP gate drop mode: flagged samples excluded from downstream count/DE. Drop list: {dropped_path}")
                else:
                    raise RuntimeError(
                        "SPP strict gate failed. Flagged samples: "
                        + flagged_text
                    )
            if generated_peak_annotation_files:
                log_it(
                    logfile,
                    "Peak genomic annotations:\n" + "\n".join(generated_peak_annotation_files),
                )
            if generated_distribution_files:
                log_it(
                    logfile,
                    "Peak genomic distribution summaries:\n" + "\n".join(generated_distribution_files),
                )
            write_peak_qc_summary_pdf(params.outputfolder, assay, sample_qc_rows, peak_qc_rows)
            write_tmp_file(params.outputfolder)
            finish_step_sample(master_config["peakqc_rule_num"], "aggregate", "peak_qc", tracking["start_time"], "OK")
        except Exception:
            finish_step_sample(master_config["peakqc_rule_num"], "aggregate", "peak_qc", tracking["start_time"], "FAIL")
            raise
