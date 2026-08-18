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

from omnomnomics.peak_annotation import build_gtf_annotation_sources


def peak_qc_input(_wildcards):
    input_files = [
        f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num'] - 1]}/all_groups.merged_peaks.bed"
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
        union_annotation=f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/peak_qc/peak_annotations/{config['THETYPE'].lower()}.all_groups.merged_peaks.annotated.bed",
        metrics=f"{experiment_dir}/{master_config['output_folders'][master_config['peakqc_rule_num'] - 1]}/peak_qc/{config['THETYPE'].lower()}.peak_qc_metrics.tsv"
    params:
        thetype=config["THETYPE"],
        broad_mode=str(config.get("BROAD_MODE", "off")).strip().lower(),
        spp_gate_mode=str(config.get("SPP_GATE", "warn")).strip().lower(),
        bam_inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['peakqc_rule_num'] - 1][0]}",
        raw_bam_inputfolder=f"{experiment_dir}/BAM",
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
                from omnomnomics.genomes import resolve_blacklist_bed
            except ImportError as exc:
                raise RuntimeError(
                    "Chromatin peak blacklist filtering requires the omnomnomics genome helper in the active environment."
                ) from exc

            try:
                return str(resolve_blacklist_bed(genome_name, config["GENOME_ASSEMBLY_DIR"]))
            except FileNotFoundError as exc:
                log_it(logfile, f"Chromatin blacklist BED not available for {genome_name}: {exc}")
                return None

        def filter_peak_bed_file(in_bed, out_bed, keep_standard=True, drop_chrm=True, blacklist_bed=None):
            work = in_bed
            with tempfile.TemporaryDirectory(prefix="omnomnomics_peak_filter_") as tmpdir:
                if keep_standard:
                    std_bed = os.path.join(tmpdir, "std.bed")
                    shell(
                        f"awk 'BEGIN{{{{OFS=\"\\t\"}}}} $1 ~ /^chr([0-9]+|X|Y)$/ {{{{print $0}}}}' {quote(work)} > {quote(std_bed)}"
                    )
                    work = std_bed
                if drop_chrm:
                    chrm_bed = os.path.join(tmpdir, "nochrm.bed")
                    shell(
                        f"awk 'BEGIN{{{{OFS=\"\\t\"}}}} $1 != \"chrM\" && $1 != \"MT\" && $1 != \"chrMT\" {{{{print $0}}}}' {quote(work)} > {quote(chrm_bed)}"
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

        def build_gene_annotation_beds(gtf_file, work_dir, promoter_upstream=1000, promoter_downstream=1000):
            source_beds = build_gtf_annotation_sources(
                gtf_file,
                work_dir,
                promoter_upstream=promoter_upstream,
                promoter_downstream=promoter_downstream,
            )
            genes_sorted = os.path.join(work_dir, "genes.sorted.bed")
            exons_sorted = os.path.join(work_dir, "exons.sorted.bed")
            promoters_sorted = os.path.join(work_dir, "promoters.sorted.bed")
            tss_sorted = os.path.join(work_dir, "tss.sorted.bed")
            introns_sorted = os.path.join(work_dir, "introns.sorted.bed")
            shell(f"sort -k1,1 -k2,2n -k3,3n {quote(source_beds['genes'])} > {quote(genes_sorted)}")
            shell(f"sort -k1,1 -k2,2n -k3,3n {quote(source_beds['exons'])} > {quote(exons_sorted)}")
            shell(f"sort -k1,1 -k2,2n -k3,3n {quote(source_beds['promoters'])} > {quote(promoters_sorted)}")
            shell(f"sort -k1,1 -k2,2n -k3,3n {quote(source_beds['tss'])} > {quote(tss_sorted)}")
            shell(
                f"bedtools subtract -a {quote(genes_sorted)} -b {quote(exons_sorted)} "
                f"| sort -k1,1 -k2,2n -k3,3n > {quote(introns_sorted)}"
            )

            return {
                "genes": genes_sorted,
                "promoters": promoters_sorted,
                "exons": exons_sorted,
                "introns": introns_sorted,
                "tss": tss_sorted,
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

                closest_file = os.path.join(tmpdir, "closest_promoter.tsv")
                shell(
                    f"bedtools closest -a {quote(peaks_indexed)} -b {quote(annotation_beds['tss'])} -d -t first > {quote(closest_file)}"
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
                    "nearest_promoter_gene",
                    "distance_to_nearest_promoter_bp",
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
                    nearest_promoter_gene, nearest_promoter_distance = closest_map.get(peak_id, ("NA", "NA"))
                    if nearest_promoter_gene in {"", ".", "-1"}:
                        nearest_promoter_gene = "NA"
                        nearest_promoter_distance = "NA"
                    nearest_gene = nearest_promoter_gene
                    nearest_distance = nearest_promoter_distance
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
                        nearest_promoter_gene,
                        nearest_promoter_distance,
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

        def make_samtools_sample_arg(fraction, seed=13):
            fraction_digits = f"{fraction:.9f}".split(".", 1)[1].rstrip("0") or "0"
            return f" -s {seed}.{fraction_digits}"

        def count_qc_alignments(bam_path, paired=False, proper_pair=False, filter_flags="2820"):
            count_args = ["samtools", "view", "-@", "1", "-c"]
            if paired and proper_pair:
                count_args.extend(["-f", "2"])
            count_args.extend(["-F", str(filter_flags), bam_path])
            return int(subprocess.check_output(count_args, stderr=subprocess.STDOUT).decode("utf-8").strip())

        def prepare_sampled_bam_for_qc(bam_path, tmpdir, label, max_alignments, paired=False, proper_pair=False, filter_flags="2820", seed=13):
            total_alignments = count_qc_alignments(bam_path, paired=paired, proper_pair=proper_pair, filter_flags=filter_flags)
            if max_alignments <= 0 or total_alignments <= max_alignments:
                return {
                    "bam_path": bam_path,
                    "input_alignments": total_alignments,
                    "counted_alignments": total_alignments,
                    "max_alignments": max_alignments,
                    "sampling_fraction": 1.0,
                    "sampled": False,
                    "samtools_sample_arg": "",
                }
            fraction = max_alignments / total_alignments
            sample_arg = make_samtools_sample_arg(fraction, seed=seed)
            sampled_bam = os.path.join(tmpdir, f"{label}.sampled.bam")
            view_command = f"samtools view -@ {threads} -b"
            if paired and proper_pair:
                view_command += " -f 2"
            view_command += f" -F {filter_flags}{sample_arg} {quote(bam_path)} > {quote(sampled_bam)}"
            shell(view_command)
            shell(f"samtools index -@ {threads} {quote(sampled_bam)}")
            counted_alignments = count_qc_alignments(sampled_bam, paired=paired, proper_pair=proper_pair, filter_flags=filter_flags)
            return {
                "bam_path": sampled_bam,
                "input_alignments": total_alignments,
                "counted_alignments": counted_alignments,
                "max_alignments": max_alignments,
                "sampling_fraction": fraction,
                "sampled": True,
                "samtools_sample_arg": sample_arg,
            }

        def complexity_sampling_parameters(bam_path, paired):
            max_alignments = int(config.get("LIBRARY_COMPLEXITY_MAX_READS", 5000000) or 0)
            if max_alignments < 0:
                max_alignments = 0
            total_alignments = count_qc_alignments(bam_path, paired=paired, proper_pair=paired, filter_flags="2820")
            fraction = 1.0
            sample_arg = ""
            if max_alignments and total_alignments > max_alignments:
                fraction = max_alignments / total_alignments
                sample_arg = make_samtools_sample_arg(fraction)
            return {
                "complexity_input_alignments": total_alignments,
                "complexity_max_alignments": max_alignments,
                "complexity_sampling_fraction": fraction,
                "samtools_sample_arg": sample_arg,
            }

        def calculate_library_complexity_metrics(bam_path):
            with tempfile.TemporaryDirectory(prefix="omnomnomics_library_qc_") as tmpdir:
                counts_path = os.path.join(tmpdir, "fragment_counts.tsv")
                sampling = complexity_sampling_parameters(bam_path, config["PAIRED"])
                if config["PAIRED"]:
                    collate_prefix = os.path.join(tmpdir, "collated")
                    bedpe_stderr = os.path.join(tmpdir, "bamtobed_bedpe.stderr")
                    complexity_command = (
                        f"samtools view -@ {threads} -u -f 2 -F 2820{sampling['samtools_sample_arg']} {quote(bam_path)} | "
                        f"samtools collate -@ {threads} -u -O - {quote(collate_prefix)} | "
                        f"bedtools bamtobed -bedpe -i stdin 2> {quote(bedpe_stderr)} | "
                        "awk 'BEGIN{{OFS=\"\\t\"}} $1==$4 && $1!=\".\" {{"
                        "start=($2<$5?$2:$5); end=($3>$6?$3:$6); print $1,start,end"
                        "}}' | "
                        f"sort -T {quote(tmpdir)} -k1,1 -k2,2n -k3,3n | uniq -c > {quote(counts_path)}"
                    )
                else:
                    complexity_command = (
                        f"samtools view -@ {threads} -u -F 2820{sampling['samtools_sample_arg']} {quote(bam_path)} | "
                        "bedtools bamtobed -i stdin | "
                        "awk 'BEGIN{{OFS=\"\\t\"}} {{print $1,$2,$3,$6}}' | "
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
                "complexity_input_alignments": sampling["complexity_input_alignments"],
                "complexity_max_alignments": sampling["complexity_max_alignments"],
                "complexity_sampling_fraction": sampling["complexity_sampling_fraction"],
            }

        def calculate_cross_correlation_metrics(bam_path, output_dir, sample_name):
            crosscorr_prefix = os.path.join(output_dir, sample_name)
            crosscorr_table = f"{crosscorr_prefix}.cross_correlation.tsv"
            crosscorr_pdf = f"{crosscorr_prefix}.cross_correlation.pdf"
            if not shutil.which("run_spp.R"):
                log_once(logfile, "step13.run_spp_missing", "run_spp.R not found on PATH. NSC/RSC metrics will be reported as missing.", "PHANTOMPEAKQUALTOOLS")
                return {
                    "est_frag_len": "",
                    "nsc": "",
                    "rsc": "",
                    "crosscorr_table": "",
                    "crosscorr_pdf": "",
                    "spp_status": "run_spp.R_not_found",
                    "spp_input_alignments": "",
                    "spp_counted_alignments": "",
                    "spp_max_alignments": "",
                    "spp_sampling_fraction": "",
                }

            max_alignments = int(config.get("SPP_MAX_READS", 10000000) or 0)
            if max_alignments < 0:
                max_alignments = 0
            run_spp_command = (
                f"run_spp.R -c={{spp_bam}} -savp={quote(crosscorr_pdf)} "
                f"-out={quote(crosscorr_table)} -p={threads} -rf"
            )
            with tempfile.TemporaryDirectory(prefix="omnomnomics_spp_qc_") as tmpdir:
                spp_input = prepare_sampled_bam_for_qc(
                    bam_path,
                    tmpdir,
                    sample_name,
                    max_alignments,
                    filter_flags="2820",
                    seed=17,
                )
                log_it(
                    logfile,
                    "SPP sampling for "
                    f"{sample_name}: input_alignments={spp_input['input_alignments']}, "
                    f"counted_alignments={spp_input['counted_alignments']}, "
                    f"max_alignments={spp_input['max_alignments']}, "
                    f"fraction={spp_input['sampling_fraction']:.9f}"
                )
                formatted_command = run_spp_command.format(spp_bam=quote(spp_input["bam_path"]))
                log_it(logfile, formatted_command, "PHANTOMPEAKQUALTOOLS COMMAND")
                shell(formatted_command)

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
            spp_status = "parsed" if nsc and rsc else "unparseable_output"
            if spp_status != "parsed":
                log_it(logfile, f"Could not parse NSC/RSC from run_spp.R output for {sample_name}: {crosscorr_table}")
            return {
                "est_frag_len": est_frag_len,
                "nsc": nsc,
                "rsc": rsc,
                "crosscorr_table": crosscorr_table,
                "crosscorr_pdf": crosscorr_pdf,
                "spp_status": spp_status,
                "spp_input_alignments": spp_input["input_alignments"],
                "spp_counted_alignments": spp_input["counted_alignments"],
                "spp_max_alignments": spp_input["max_alignments"],
                "spp_sampling_fraction": spp_input["sampling_fraction"],
            }

        def calculate_peak_qc_metrics(peak_bed, bam_files):
            with tempfile.TemporaryDirectory(prefix="omnomnomics_peak_qc_") as tmpdir:
                merged_bed = os.path.join(tmpdir, "merged_peaks.bed")
                multicov_output = os.path.join(tmpdir, "multicov.tsv")
                max_alignments = int(config.get("FRIP_MAX_READS", 10000000) or 0)
                if max_alignments < 0:
                    max_alignments = 0
                shell(
                    f"sort -k1,1 -k2,2n -k3,3n {quote(peak_bed)} | "
                    f"bedtools merge -i - > {quote(merged_bed)}"
                )
                frip_inputs = []
                for idx, bam_path in enumerate(bam_files, start=1):
                    frip_inputs.append(
                        prepare_sampled_bam_for_qc(
                            bam_path,
                            tmpdir,
                            f"frip_{idx}",
                            max_alignments,
                            filter_flags="260",
                            seed=23 + idx,
                        )
                    )
                bams_for_frip = [entry["bam_path"] for entry in frip_inputs]
                shell(
                    f"bedtools multicov -bed {quote(merged_bed)} -bams {' '.join(quote(path) for path in bams_for_frip)} "
                    f"> {quote(multicov_output)}"
                )

                peak_count = len(load_peak_records(peak_bed))
                merged_interval_count = 0
                total_peak_bp = 0
                reads_in_peaks = 0
                with open(multicov_output, newline="") as handle:
                    reader = csv.reader(handle, delimiter="\t")
                    for row in reader:
                        merged_interval_count += 1
                        total_peak_bp += int(row[2]) - int(row[1])
                        reads_in_peaks += sum(int(value) for value in row[3:])

                frip_input_alignments = sum(entry["input_alignments"] for entry in frip_inputs)
                total_aligned_reads = sum(entry["counted_alignments"] for entry in frip_inputs)
                min_sampling_fraction = min((entry["sampling_fraction"] for entry in frip_inputs), default=1.0)
                frip_estimated = any(entry["sampled"] for entry in frip_inputs)

            frip = (reads_in_peaks / total_aligned_reads) if total_aligned_reads else 0.0
            return {
                "peak_count": peak_count,
                "merged_interval_count": merged_interval_count,
                "total_peak_bp": total_peak_bp,
                "reads_in_peaks": reads_in_peaks,
                "total_aligned_reads": total_aligned_reads,
                "frip_input_alignments": frip_input_alignments,
                "frip_max_alignments": max_alignments,
                "frip_min_sampling_fraction": min_sampling_fraction,
                "frip_estimated": "yes" if frip_estimated else "no",
                "frip": frip,
            }

        def calculate_sample_frip_metrics(peak_bed, sample_bams):
            with tempfile.TemporaryDirectory(prefix="omnomnomics_sample_frip_") as tmpdir:
                merged_bed = os.path.join(tmpdir, "merged_peaks.bed")
                multicov_output = os.path.join(tmpdir, "multicov.tsv")
                max_alignments = int(config.get("FRIP_MAX_READS", 10000000) or 0)
                if max_alignments < 0:
                    max_alignments = 0
                shell(
                    f"sort -k1,1 -k2,2n -k3,3n {quote(peak_bed)} | "
                    f"bedtools merge -i - > {quote(merged_bed)}"
                )
                frip_inputs = []
                for idx, (sample_name, bam_path) in enumerate(sample_bams, start=1):
                    sampled = prepare_sampled_bam_for_qc(
                        bam_path,
                        tmpdir,
                        f"sample_frip_{idx}",
                        max_alignments,
                        filter_flags="260",
                        seed=101 + idx,
                    )
                    sampled["sample_name"] = sample_id_for_sample(sample_name)
                    sampled["bam_file"] = os.path.basename(bam_path)
                    frip_inputs.append(sampled)
                bams_for_frip = [entry["bam_path"] for entry in frip_inputs]
                if not bams_for_frip:
                    return []
                shell(
                    f"bedtools multicov -bed {quote(merged_bed)} -bams {' '.join(quote(path) for path in bams_for_frip)} "
                    f"> {quote(multicov_output)}"
                )

                peak_count = len(load_peak_records(peak_bed))
                merged_interval_count = 0
                total_peak_bp = 0
                reads_by_sample = [0 for _ in frip_inputs]
                with open(multicov_output, newline="") as handle:
                    reader = csv.reader(handle, delimiter="\t")
                    for row in reader:
                        merged_interval_count += 1
                        total_peak_bp += int(row[2]) - int(row[1])
                        for idx, value in enumerate(row[3:]):
                            reads_by_sample[idx] += int(value)

            sample_rows = []
            for idx, entry in enumerate(frip_inputs):
                total_aligned_reads = entry["counted_alignments"]
                reads_in_peaks = reads_by_sample[idx]
                frip = (reads_in_peaks / total_aligned_reads) if total_aligned_reads else 0.0
                sample_rows.append(
                    {
                        "sample": entry["sample_name"],
                        "bam_file": entry["bam_file"],
                        "peak_count": peak_count,
                        "merged_interval_count": merged_interval_count,
                        "total_peak_bp": total_peak_bp,
                        "reads_in_peaks": reads_in_peaks,
                        "total_aligned_reads": total_aligned_reads,
                        "frip_input_alignments": entry["input_alignments"],
                        "frip_max_alignments": max_alignments,
                        "frip_sampling_fraction": entry["sampling_fraction"],
                        "frip_estimated": "yes" if entry["sampled"] else "no",
                        "frip": frip,
                    }
                )
            return sample_rows

        def write_sample_frip_outputs(outputfolder, thetype, sample_frip_rows):
            if not sample_frip_rows:
                return None

            qc_dir = ensure_peak_qc_dir(outputfolder)
            qc_table = os.path.join(qc_dir, f"{thetype.lower()}.sample_frip_metrics.tsv")
            with open(qc_table, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "assay",
                    "group",
                    "peak_set",
                    "peak_file",
                    "sample",
                    "bam_file",
                    "peak_count",
                    "merged_interval_count",
                    "total_peak_bp",
                    "reads_in_peaks",
                    "total_aligned_reads",
                    "frip_input_alignments",
                    "frip_max_alignments",
                    "frip_sampling_fraction",
                    "frip_estimated",
                    "frip",
                ])
                for row in sample_frip_rows:
                    writer.writerow([
                        row["assay"],
                        row["group"],
                        row["peak_set"],
                        row["peak_file"],
                        row["sample"],
                        row["bam_file"],
                        row["peak_count"],
                        row["merged_interval_count"],
                        row["total_peak_bp"],
                        row["reads_in_peaks"],
                        row["total_aligned_reads"],
                        row["frip_input_alignments"],
                        row["frip_max_alignments"],
                        f"{row['frip_sampling_fraction']:.9f}",
                        row["frip_estimated"],
                        f"{row['frip']:.6f}",
                    ])
            return qc_table

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
                    "complexity_bam_file",
                    "complexity_bam_stage",
                    "complexity_input_alignments",
                    "complexity_max_alignments",
                    "complexity_sampling_fraction",
                    "total_reads",
                    "distinct_reads",
                    "one_read_sites",
                    "two_read_sites",
                    "nrf",
                    "pbc1",
                    "pbc2",
                    "spp_input_alignments",
                    "spp_counted_alignments",
                    "spp_max_alignments",
                    "spp_sampling_fraction",
                    "est_frag_len",
                    "nsc",
                    "rsc",
                    "spp_status",
                    "crosscorr_table",
                    "crosscorr_pdf",
                ])
                for row in qc_rows:
                    writer.writerow([
                        row["assay"],
                        row["sample"],
                        row["bam_file"],
                        row["complexity_bam_file"],
                        row["complexity_bam_stage"],
                        row["complexity_input_alignments"],
                        row["complexity_max_alignments"],
                        f"{row['complexity_sampling_fraction']:.9f}",
                        row["total_reads"],
                        row["distinct_reads"],
                        row["one_read_sites"],
                        row["two_read_sites"],
                        f"{row['nrf']:.6f}",
                        f"{row['pbc1']:.6f}",
                        row["pbc2"] if row["pbc2"] == "" else f"{float(row['pbc2']):.6f}",
                        row["spp_input_alignments"],
                        row["spp_counted_alignments"],
                        row["spp_max_alignments"],
                        row["spp_sampling_fraction"] if row["spp_sampling_fraction"] == "" else f"{float(row['spp_sampling_fraction']):.9f}",
                        row["est_frag_len"],
                        row["nsc"],
                        row["rsc"],
                        row.get("spp_status", ""),
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
            unique_rows = []
            seen_rows = set()
            for row in qc_rows:
                row_key = (row["assay"], row["group"], row["peak_set"], row["peak_file"])
                if row_key in seen_rows:
                    continue
                seen_rows.add(row_key)
                unique_rows.append(row)
            with open(qc_table, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow([
                    "assay",
                    "group",
                    "peak_set",
                    "peak_file",
                    "bam_count",
                    "peak_count",
                    "merged_interval_count",
                    "total_peak_bp",
                    "reads_in_peaks",
                    "total_aligned_reads",
                    "frip_input_alignments",
                    "frip_max_alignments",
                    "frip_min_sampling_fraction",
                    "frip_estimated",
                    "frip",
                ])
                for row in unique_rows:
                    writer.writerow([
                        row["assay"],
                        row["group"],
                        row["peak_set"],
                        row["peak_file"],
                        row["bam_count"],
                        row["peak_count"],
                        row["merged_interval_count"],
                        row["total_peak_bp"],
                        row["reads_in_peaks"],
                        row["total_aligned_reads"],
                        row["frip_input_alignments"],
                        row["frip_max_alignments"],
                        f"{row['frip_min_sampling_fraction']:.9f}",
                        row["frip_estimated"],
                        f"{row['frip']:.6f}",
                    ])
            return qc_table

        def write_peak_qc_summary_pdf(outputfolder, thetype, broad_mode, sample_rows, peak_rows, sample_frip_rows):
            if not sample_rows and not peak_rows and not sample_frip_rows:
                return

            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from matplotlib.backends.backend_pdf import PdfPages

            log_once(logfile, "step13.matplotlib_version", f"\nmatplotlib {matplotlib.__version__}\n", "MATPLOTLIB VERSION")

            qc_dir = ensure_peak_qc_dir(outputfolder)
            pdf_path = os.path.join(qc_dir, f"{thetype.lower()}.peak_qc_summary.pdf")
            feature_terms = {
                "genebody": ("gene-body", "Gene bodies"),
                "domain": ("domain", "Domains"),
                "diffuse": ("bin", "Bins"),
            }
            feature_term, feature_term_plural = feature_terms.get(broad_mode, ("peak", "Peaks"))

            plt.style.use("seaborn-v0_8-whitegrid")
            with PdfPages(pdf_path) as pdf:
                if peak_rows:
                    peak_labels = [row["peak_set"] for row in peak_rows]
                    frip_values = [float(row["frip"]) for row in peak_rows]
                    peak_counts = [float(row["peak_count"]) for row in peak_rows]
                    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
                    axes[0].bar(peak_labels, frip_values, color="#4C78A8")
                    axes[0].set_title(f"{thetype} {feature_term} QC: FRiP")
                    axes[0].set_ylabel("FRiP")
                    axes[0].tick_params(axis="x", rotation=45, labelsize=9)
                    axes[1].bar(peak_labels, peak_counts, color="#D95F02")
                    axes[1].set_title(f"{thetype} {feature_term} QC: {feature_term} count")
                    axes[1].set_ylabel(feature_term_plural)
                    axes[1].tick_params(axis="x", rotation=45, labelsize=9)
                    fig.tight_layout()
                    pdf.savefig(fig, bbox_inches="tight")
                    plt.close(fig)

                if sample_frip_rows:
                    frip_sets = sorted(set((row["group"], row["peak_set"]) for row in sample_frip_rows))
                    for group, peak_set in frip_sets:
                        plot_rows = [
                            row for row in sample_frip_rows
                            if row["group"] == group and row["peak_set"] == peak_set
                        ]
                        if not plot_rows:
                            continue
                        sample_labels = [row["sample"] for row in plot_rows]
                        frip_values = [float(row["frip"]) for row in plot_rows]
                        fig, ax = plt.subplots(1, 1, figsize=(max(8, len(sample_labels) * 0.55), 5))
                        ax.bar(sample_labels, frip_values, color="#4C78A8")
                        ax.set_title(f"{thetype} {feature_term} sample FRiP: {group} / {peak_set}")
                        ax.set_ylabel("FRiP")
                        ax.set_xlabel("Sample")
                        ax.tick_params(axis="x", rotation=45, labelsize=9)
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
                raw_bam_path = os.path.join(params.raw_bam_inputfolder, f"{sample}.bam")
                if os.path.exists(raw_bam_path):
                    complexity_bam_path = raw_bam_path
                    complexity_bam_stage = "merged_pre_dedup"
                else:
                    complexity_bam_path = bam_path
                    complexity_bam_stage = "filtered_fallback"
                    log_it(
                        logfile,
                        f"Merged pre-dedup BAM not found for {sample}; library complexity uses filtered BAM fallback.",
                    )
                log_it(logfile, f"Calculating library complexity metrics for {sample}...")
                complexity_metrics = calculate_library_complexity_metrics(complexity_bam_path)
                log_it(
                    logfile,
                    "Library complexity sampling for "
                    f"{sample}: input_alignments={complexity_metrics['complexity_input_alignments']}, "
                    f"max_alignments={complexity_metrics['complexity_max_alignments']}, "
                    f"fraction={complexity_metrics['complexity_sampling_fraction']:.9f}"
                )
                log_it(logfile, f"Calculating cross-correlation metrics for {sample}...")
                crosscorr_metrics = calculate_cross_correlation_metrics(bam_path, qc_dir, sample)
                sample_rows.append({
                    "assay": thetype,
                    "sample": sample_id_for_sample(sample),
                    "bam_file": os.path.basename(bam_path),
                    "complexity_bam_file": os.path.basename(complexity_bam_path),
                    "complexity_bam_stage": complexity_bam_stage,
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
            sample_frip_rows = []
            generated_peak_annotation_files = []
            generated_distribution_files = []
            peak_bed_source_folder = params.peak_outputfolder
            filter_called_peaks = assay == "ATAC" or (
                assay == "CHIP" and params.broad_mode not in {"genebody", "diffuse"}
            )
            if filter_called_peaks:
                filtered_dir, blacklist_bed, filtered_peak_beds, filtering_summary = build_filtered_peak_sets(
                    params.peak_outputfolder,
                    params.outputfolder,
                    config["THEGENOME"],
                )
                peak_bed_source_folder = filtered_dir
                log_it(logfile, f"{assay} blacklist BED: {blacklist_bed if blacklist_bed else 'NA'}")
                log_it(logfile, f"{assay} filtered peak BEDs: {len(filtered_peak_beds)} files in {filtered_dir}")
                log_it(logfile, f"{assay} peak filtering summary: {filtering_summary}")
            with tempfile.TemporaryDirectory(prefix="omnomnomics_gene_anno_") as anno_tmpdir:
                annotation_beds = build_gene_annotation_beds(params.gtf_file, anno_tmpdir)
                for group in sorted(grouped_bams):
                    sample_bams = grouped_bams[group]
                    bams = [bam_path for _, bam_path in sample_bams]
                    if filter_called_peaks:
                        preferred = os.path.join(peak_bed_source_folder, f"{group}.MACS3.optimized.bed")
                        fallback = os.path.join(peak_bed_source_folder, f"{group}.MACS3.q-0p01.shiftm100.ext200.group_peaks.bed")
                        candidate_peak_sets = [preferred, fallback]
                    else:
                        preferred = os.path.join(peak_bed_source_folder, f"{group}.MACS3.optimized.bed")
                        fallback_glob = sorted(glob.glob(os.path.join(peak_bed_source_folder, f"{group}.MACS3*.bed")))
                        candidate_peak_sets = list(dict.fromkeys([preferred, *fallback_glob]))
                    for peak_bed in candidate_peak_sets:
                        if not os.path.exists(peak_bed):
                            continue
                        metrics = calculate_peak_qc_metrics(peak_bed, bams)
                        peak_set_name = os.path.basename(peak_bed).replace(".bed", "")
                        if params.broad_mode == "genebody":
                            peak_set_name = f"{group}.gene_bodies"
                        elif params.broad_mode == "domain":
                            peak_set_name = f"{group}.domains"
                        elif params.broad_mode == "diffuse":
                            peak_set_name = f"{group}.bins"
                        peak_qc_rows.append({
                            "assay": assay,
                            "group": group,
                            "peak_set": peak_set_name,
                            "peak_file": os.path.basename(peak_bed),
                            "bam_count": len(bams),
                            **metrics,
                        })
                        for sample_frip in calculate_sample_frip_metrics(peak_bed, sample_bams):
                            sample_frip_rows.append({
                                "assay": assay,
                                "group": group,
                                "peak_set": peak_set_name,
                                "peak_file": os.path.basename(peak_bed),
                                **sample_frip,
                            })
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
            sample_frip_table = write_sample_frip_outputs(params.outputfolder, assay, sample_frip_rows)
            sample_qc_table = write_sample_qc_outputs(params.outputfolder, assay, sample_qc_rows)
            spp_qc_table, spp_flagged = write_spp_qc_summary(params.outputfolder, assay, sample_qc_rows)
            if peak_qc_table:
                log_it(logfile, f"Peak QC metrics: {peak_qc_table}")
            if sample_frip_table:
                log_it(logfile, f"Sample FRiP metrics: {sample_frip_table}")
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
            write_peak_qc_summary_pdf(
                params.outputfolder,
                assay,
                params.broad_mode,
                sample_qc_rows,
                peak_qc_rows,
                sample_frip_rows,
            )
            write_tmp_file(params.outputfolder)
            finish_step_sample(master_config["peakqc_rule_num"], "aggregate", "peak_qc", tracking["start_time"], "OK")
        except Exception:
            finish_step_sample(master_config["peakqc_rule_num"], "aggregate", "peak_qc", tracking["start_time"], "FAIL")
            raise
