# Rule 10: Call Peaks

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob
import subprocess
import tempfile
import itertools
import statistics
import math
import csv
import shutil
import hashlib
import concurrent.futures
from shlex import quote

def input_function(wildcards):
    input_folder1 = f"{experiment_dir}/{master_config['input_folders'][master_config['callpeaks_rule_num']-1]}"
    input_files = []
    if config['THETYPE'] == "CHIP":
        for sample in samples2:
            input_files.append(f"{input_folder1}/{sample}.filtered.bam")
    elif config['THETYPE'] == "ATAC":
        for sample in samples2:
            input_files.append(f"{input_folder1}/{sample}.sorted.dups_marked.filtered.bam")
    if config['THETYPE'] == "CHIP" and str(config.get("BROAD_MODE", "off")).lower() == "genebody":
        input_files.append(
            os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "annotation", "genes.gtf")
        )
    return input_files

rule call_peaks:
    input:
        input_function
    output:
        extra=f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/extra_10.tmp",
        merged_features=f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/all_groups.merged_peaks.bed"
    params:
        thetype= lambda wildcards: config['THETYPE'],  
        broad_mode= lambda wildcards: str(config.get('BROAD_MODE', 'off')).lower(),
        input_sample= lambda wildcards: config['INPUT'],
        experiment_dir= lambda wildcards: config['EXPERIMENT_DIR'], 
        inputfolder1=lambda wildcards: f"{experiment_dir}/{master_config['input_folders'][master_config['callpeaks_rule_num']-1]}",
        outputfolder= lambda wildcards: f"{experiment_dir}/{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}",
        gtf_file=lambda wildcards: os.path.join(
            config["GENOME_ASSEMBLY_DIR"],
            config["THEGENOME"],
            "annotation",
            "genes.gtf",
        ),
    threads:
        lambda wildcards: Threads_Per_Rule['10']
    resources:
        mem_mb = lambda wildcards: Memory_Per_Rule['10'],
        partition = lambda wildcards: master_config['partition'],
        runtime = lambda wildcards: Runtime_Per_Rule['10']
    run:
        log_once(logfile, "step10.header", "Calling Peaks...", f"EXECUTING STEP {master_config['callpeaks_rule_num']}")
        log_once(logfile, "step10.inputfolder", f"Input folder: {params.inputfolder1}")
        log_once(logfile, "step10.outputfolder", f"Output folder: {params.outputfolder}")

        def chip_style_label(broad_mode):
            if broad_mode == "domain":
                return "broad domain marks"
            if broad_mode == "genebody":
                return "gene-body marks"
            if broad_mode == "diffuse":
                return "diffuse broad marks"
            return "TF / narrow peaks"

        def macs3_chip_presets(broad_mode):
            presets = [
                ("q-0p05", "0.05"),
                ("q-0p01", "0.01"),
                ("q-0p001", "0.001"),
            ]
            if broad_mode == "domain":
                return [(f"{label}.broad", value) for label, value in presets]
            return presets

        def macs3_peak_to_bed_path(peak_path):
            basename = os.path.basename(peak_path)
            if basename.endswith(".narrowPeak"):
                basename = basename[:-11]
            elif basename.endswith(".broadPeak"):
                basename = basename[:-10]
            return os.path.join(os.path.dirname(peak_path), f"{basename}.bed")

        def grouped_bams_by_sample_type(inputfolder, bam_suffix):
            grouped_bams = {}
            for sample in samples2:
                bam_path = os.path.join(inputfolder, f"{sample}{bam_suffix}")
                if not os.path.exists(bam_path):
                    continue
                group_name = sample_type_for_sample(sample)
                grouped_bams.setdefault(group_name, []).append((sample, bam_path))
            return grouped_bams

        def maybe_resolve_blacklist_bed(genome_name, outdir, cache_label):
            try:
                from omnomnomics.genomes import resolve_blacklist_bed

                return str(resolve_blacklist_bed(genome_name, config["GENOME_ASSEMBLY_DIR"]))
            except Exception:
                return None

        def build_chip_diffuse_bin_feature_sets(
            reference_bam,
            genome_name,
            outputfolder,
            groups,
            bin_size,
            keep_standard,
            drop_chrm,
            blacklist_bed,
        ):
            diffuse_root = os.path.join(outputfolder, "chip_diffuse_features")
            os.makedirs(diffuse_root, exist_ok=True)
            genome_sizes = os.path.join(diffuse_root, "genome.sizes")
            raw_bins = os.path.join(diffuse_root, "all_bins.raw.bed")
            filtered_bins = os.path.join(outputfolder, "all_groups.merged_peaks.bed")
            summary_path = os.path.join(diffuse_root, "diffuse_feature_summary.tsv")

            shell(
                f"""samtools view -H {quote(reference_bam)} | awk 'BEGIN{{{{FS="\\t"; OFS="\\t"}}}} /^@SQ/ {{{{
chrom=""; len="";
for (i=1; i<=NF; i++) {{{{
  if ($i ~ /^SN:/) chrom=substr($i,4);
  else if ($i ~ /^LN:/) len=substr($i,4);
}}}}
if (chrom != "" && len != "") print chrom, len;
}}}}' > {quote(genome_sizes)}"""
            )
            shell(f"bedtools makewindows -g {quote(genome_sizes)} -w {int(bin_size)} > {quote(raw_bins)}")

            work = raw_bins
            with tempfile.TemporaryDirectory(prefix="omnomnomics_chip_diffuse_") as tmpdir:
                if keep_standard:
                    std_bins = os.path.join(tmpdir, "std_bins.bed")
                    shell(
                        f"awk 'BEGIN{{{{OFS=\"\\t\"}}}} $1 ~ /^chr([0-9]+|X|Y)$/ {{{{print $0}}}}' {quote(work)} > {quote(std_bins)}"
                    )
                    work = std_bins
                if drop_chrm:
                    no_chrm_bins = os.path.join(tmpdir, "no_chrm_bins.bed")
                    shell(
                        f"awk 'BEGIN{{{{OFS=\"\\t\"}}}} $1 != \"chrM\" && $1 != \"MT\" && $1 != \"chrMT\" {{{{print $0}}}}' {quote(work)} > {quote(no_chrm_bins)}"
                    )
                    work = no_chrm_bins
                if blacklist_bed and os.path.exists(blacklist_bed):
                    no_blacklist_bins = os.path.join(tmpdir, "noblacklist_bins.bed")
                    shell(
                        f"bedtools intersect -v -a {quote(work)} -b {quote(blacklist_bed)} > {quote(no_blacklist_bins)}"
                    )
                    work = no_blacklist_bins
                shell(f"sort -k1,1 -k2,2n -k3,3n {quote(work)} > {quote(filtered_bins)}")

            feature_count = int(
                subprocess.check_output(
                    f"wc -l {quote(filtered_bins)} | awk '{{print $1}}'",
                    shell=True,
                    executable="/bin/bash",
                    text=True,
                ).strip()
            )
            raw_count = int(
                subprocess.check_output(
                    f"wc -l {quote(raw_bins)} | awk '{{print $1}}'",
                    shell=True,
                    executable="/bin/bash",
                    text=True,
                ).strip()
            )

            with open(summary_path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["group", "feature_bed", "feature_count", "feature_type", "bin_size_bp"])
                for group in sorted(groups):
                    group_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                    shutil.copy2(filtered_bins, group_bed)
                    writer.writerow([group, group_bed, feature_count, "diffuse_bin", int(bin_size)])
                writer.writerow(["all_groups", filtered_bins, feature_count, "diffuse_bin", int(bin_size)])
                writer.writerow([])
                writer.writerow(["metric", "value"])
                writer.writerow(["raw_bin_count", raw_count])
                writer.writerow(["filtered_bin_count", feature_count])
                writer.writerow(["bin_size_bp", int(bin_size)])
                writer.writerow(["standard_chroms_only", str(bool(keep_standard)).lower()])
                writer.writerow(["exclude_chrm", str(bool(drop_chrm)).lower()])
                writer.writerow(["blacklist_bed", blacklist_bed if blacklist_bed else "NA"])

            log_it(logfile, f"Diffuse bin feature BED written to {filtered_bins}")
            log_it(logfile, f"Diffuse bin feature summary: {summary_path}")
            log_it(logfile, f"Total diffuse bins retained: {feature_count}")

        def count_mapped_reads_total(bam_paths):
            total = 0
            for bam in bam_paths:
                cmd = f"samtools view -c -F 260 {quote(bam)}"
                try:
                    total += int(subprocess.check_output(cmd, shell=True, executable="/bin/bash", text=True).strip())
                except Exception:
                    pass
            return total

        def bed_peak_count_and_widths(bed_path):
            widths = []
            count = 0
            with open(bed_path) as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 3:
                        continue
                    start = int(f[1])
                    end = int(f[2])
                    if end <= start:
                        continue
                    widths.append(end - start)
                    count += 1
            median_width = statistics.median(widths) if widths else 0.0
            return count, median_width

        def frip_for_group(bed_path, bam_paths):
            if not bam_paths:
                return 0.0, 0, 0
            with tempfile.NamedTemporaryFile(prefix="omnom_multicov_", suffix=".tsv", delete=False) as tmp:
                multicov_path = tmp.name
            try:
                cmd = (
                    f"bedtools multicov -bed {quote(bed_path)} -bams {' '.join(quote(x) for x in bam_paths)} "
                    f"> {quote(multicov_path)}"
                )
                shell(cmd)
                reads_in_peaks = 0
                with open(multicov_path) as handle:
                    for line in handle:
                        if not line.strip():
                            continue
                        fields = line.rstrip("\n").split("\t")
                        if len(fields) <= 3:
                            continue
                        counts = fields[3:]
                        reads_in_peaks += sum(int(x) for x in counts if x.isdigit())
                mapped_total = count_mapped_reads_total(bam_paths)
                frip = (reads_in_peaks / mapped_total) if mapped_total > 0 else 0.0
                return frip, reads_in_peaks, mapped_total
            finally:
                if os.path.exists(multicov_path):
                    os.remove(multicov_path)

        def mean_pairwise_jaccard(sample_beds):
            if len(sample_beds) < 2:
                return 1.0
            vals = []
            for a, b in itertools.combinations(sample_beds, 2):
                cmd = f"bedtools jaccard -a {quote(a)} -b {quote(b)}"
                out = subprocess.check_output(cmd, shell=True, executable="/bin/bash", text=True).strip().splitlines()
                if len(out) > 1:
                    cols = out[1].split("\t")
                    if len(cols) >= 3:
                        try:
                            vals.append(float(cols[2]))
                        except ValueError:
                            pass
            if not vals:
                return 0.0
            return sum(vals) / len(vals)

        def blacklist_overlap_fraction(bed_path, blacklist_bed):
            if not blacklist_bed or not os.path.exists(blacklist_bed):
                return 0.0
            with tempfile.NamedTemporaryFile(prefix="omnom_black_ov_", suffix=".bed", delete=False) as tmp:
                overlap_path = tmp.name
            try:
                shell(f"bedtools intersect -u -a {quote(bed_path)} -b {quote(blacklist_bed)} > {quote(overlap_path)}")
                total_n, _ = bed_peak_count_and_widths(bed_path)
                ov_n, _ = bed_peak_count_and_widths(overlap_path)
                return (ov_n / total_n) if total_n > 0 else 0.0
            finally:
                if os.path.exists(overlap_path):
                    os.remove(overlap_path)

        def zscore_map(metric_map, higher_is_better=True):
            vals = list(metric_map.values())
            if not vals:
                return {k: 0.0 for k in metric_map}
            mean_v = sum(vals) / len(vals)
            var = sum((x - mean_v) ** 2 for x in vals) / max(1, len(vals) - 1)
            sd = math.sqrt(var) if var > 0 else 0.0
            out = {}
            for key, value in metric_map.items():
                z = ((value - mean_v) / sd) if sd > 0 else 0.0
                out[key] = z if higher_is_better else -z
            return out

        def write_optimization_plots(rows, report_dir):
            try:
                import matplotlib.pyplot as plt
            except Exception:
                log_it(logfile, "matplotlib not available; skipping peak-optimization plots.")
                return

            by_group = {}
            for row in rows:
                by_group.setdefault(row["group"], []).append(row)
            for group, grows in by_group.items():
                grows = sorted(grows, key=lambda x: x["composite_score"], reverse=True)
                labels = [x["candidate"] for x in grows]
                scores = [x["composite_score"] for x in grows]
                frips = [x["frip"] for x in grows]
                jacs = [x["replicate_jaccard"] for x in grows]

                fig = plt.figure(figsize=(10, 4))
                ax = fig.add_subplot(1, 1, 1)
                ax.bar(labels, scores, color="#2b8cbe")
                ax.set_title(f"ATAC peak optimization: composite score ({group})")
                ax.set_ylabel("score")
                ax.tick_params(axis="x", rotation=35)
                fig.tight_layout()
                fig.savefig(os.path.join(report_dir, f"{group}.composite_score_barplot.png"), dpi=150)
                plt.close(fig)

                fig = plt.figure(figsize=(6, 5))
                ax = fig.add_subplot(1, 1, 1)
                ax.scatter(frips, jacs, c=scores, cmap="viridis", s=80)
                for i, label in enumerate(labels):
                    ax.annotate(label, (frips[i], jacs[i]), fontsize=8)
                ax.set_xlabel("FRiP")
                ax.set_ylabel("Replicate Jaccard")
                ax.set_title(f"ATAC peak optimization: FRiP vs Jaccard ({group})")
                fig.tight_layout()
                fig.savefig(os.path.join(report_dir, f"{group}.frip_vs_jaccard.png"), dpi=150)
                plt.close(fig)

        def build_macs3_callpeak_command(
            treatment_bams,
            outdir,
            name,
            qvalue,
            control_bam=None,
            broad_mode=False,
            shift=None,
            extsize=None,
            assay_type=None,
            broad_cutoff=None,
            broad_min_length=None,
            broad_max_gap=None,
        ):
            cmd = [
                "macs3", "callpeak",
                "-t", *treatment_bams,
                "--outdir", outdir,
                "-n", name,
                "-q", str(qvalue),
                "--verbose", "0",
            ]
            assay_label = str(assay_type or "").strip().upper()
            if assay_label == "ATAC":
                cmd.extend(["-f", "BAMPE"])
            elif assay_label == "CHIP":
                cmd.extend(["-f", "BAM"])
            if control_bam and control_bam != "NA":
                cmd.extend(["-c", control_bam])
            if broad_mode:
                cmd.extend(["--broad", "--broad-cutoff", str(float(broad_cutoff) if broad_cutoff is not None else 0.1)])
                if broad_min_length not in (None, "", "NA"):
                    cmd.extend(["--min-length", str(int(broad_min_length))])
                if broad_max_gap not in (None, "", "NA"):
                    cmd.extend(["--max-gap", str(int(broad_max_gap))])
            if shift is not None and extsize is not None:
                cmd.extend(["--nomodel", "--shift", str(int(shift)), "--extsize", str(int(extsize))])
            return " ".join(quote(str(x)) for x in cmd)

        def make_deterministic_pseudorep_bams(input_bam, work_dir, out_prefix, split_seed):
            bam1 = os.path.join(work_dir, f"{out_prefix}.ps1.bam")
            bam2 = os.path.join(work_dir, f"{out_prefix}.ps2.bam")
            seed = int(split_seed) % 1000000
            samtools_threads = max(1, min(int(threads), 8))
            shell(
                f"samtools view -@ {samtools_threads} -b -s {seed}.5 "
                f"-o {quote(bam1)} -U {quote(bam2)} {quote(input_bam)}"
            )
            return bam1, bam2

        def read_bed_triplets(bed_path):
            rows = []
            if not os.path.exists(bed_path):
                return rows
            with open(bed_path) as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 3:
                        continue
                    rows.append((fields[0], int(fields[1]), int(fields[2])))
            return rows

        def write_bed_triplets(rows, bed_path):
            sorted_rows = sorted(set(rows), key=lambda x: (x[0], x[1], x[2]))
            with open(bed_path, "w", encoding="utf-8") as handle:
                for chrom, start, end in sorted_rows:
                    handle.write(f"{chrom}\t{start}\t{end}\n")
            return len(sorted_rows)

        def parse_gtf_attributes(attr_string):
            attrs = {}
            for item in str(attr_string).strip().split(";"):
                item = item.strip()
                if not item or " " not in item:
                    continue
                key, value = item.split(" ", 1)
                attrs[key] = value.strip().strip('"')
            return attrs

        def chrom_is_standard(chrom_name):
            text = str(chrom_name)
            if text in {"chrX", "chrY"}:
                return True
            if text.startswith("chr") and text[3:].isdigit():
                return True
            return False

        def build_chip_genebody_feature_sets(gtf_file, outputfolder, groups):
            if not os.path.isfile(gtf_file):
                raise FileNotFoundError(f"Genome annotation GTF not found for gene-body mode: {gtf_file}")

            feature_root = os.path.join(outputfolder, "chip_genebody_features")
            os.makedirs(feature_root, exist_ok=True)
            union_bed = os.path.join(outputfolder, "all_groups.merged_peaks.bed")
            summary_path = os.path.join(feature_root, "genebody_feature_summary.tsv")
            gene_body_rows = []

            with open(gtf_file, "r", encoding="utf-8") as handle:
                for line in handle:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 9 or fields[2] != "gene":
                        continue
                    chrom = fields[0]
                    if not chrom_is_standard(chrom) or chrom in {"chrM", "MT", "chrMT"}:
                        continue
                    start = int(fields[3]) - 1
                    end = int(fields[4])
                    if end <= start:
                        continue
                    attrs = parse_gtf_attributes(fields[8])
                    gene_id = attrs.get("gene_id", "NA")
                    gene_name = attrs.get("gene_name", gene_id)
                    gene_body_rows.append((chrom, start, end, gene_id, gene_name))

            gene_body_rows = sorted(
                {(chrom, start, end, gene_id, gene_name) for chrom, start, end, gene_id, gene_name in gene_body_rows},
                key=lambda x: (x[0], x[1], x[2], x[3], x[4]),
            )
            feature_triplets = [(chrom, start, end) for chrom, start, end, _, _ in gene_body_rows]
            feature_count = write_bed_triplets(feature_triplets, union_bed)

            with open(summary_path, "w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["group", "feature_bed", "feature_count", "feature_type"])
                for group in sorted(groups):
                    group_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                    shutil.copy2(union_bed, group_bed)
                    writer.writerow([group, group_bed, feature_count, "gene_body"])
                writer.writerow(["all_groups", union_bed, feature_count, "gene_body"])

            log_it(logfile, f"Gene-body feature BED written to {union_bed}")
            log_it(logfile, f"Gene-body feature summary: {summary_path}")
            log_it(logfile, f"Total gene-body features retained: {feature_count}")

        def pooled_overlap_support_bed(pooled_bed, support_beds, min_support, output_bed, overlap_fraction=0.5):
            if not support_beds:
                return write_bed_triplets([], output_bed)
            supported_rows = []
            with tempfile.TemporaryDirectory(prefix="omnom_broad_overlap_") as tmpdir:
                for idx, support_bed in enumerate(support_beds, start=1):
                    filtered = os.path.join(tmpdir, f"support_{idx}.bed")
                    shell(
                        f"bedtools intersect -u -f {float(overlap_fraction)} "
                        f"-a {quote(pooled_bed)} -b {quote(support_bed)} > {quote(filtered)}"
                    )
                    supported_rows.extend(read_bed_triplets(filtered))
            counts = {}
            for row in supported_rows:
                counts[row] = counts.get(row, 0) + 1
            kept = [row for row, count in counts.items() if count >= int(min_support)]
            return write_bed_triplets(kept, output_bed)

        def broad_domain_consensus_from_replicates(
            grouped_bams,
            outputfolder,
            control_bam,
            replicate_fraction,
            overlap_fraction,
            relaxed_qvalue,
            broad_cutoff,
            broad_min_length,
            broad_max_gap,
        ):
            broad_root = os.path.join(outputfolder, "chip_broad_domain_concordance")
            os.makedirs(broad_root, exist_ok=True)
            summary_rows = []

            for group in sorted(grouped_bams):
                group_records = grouped_bams[group]
                group_bams = [bam_path for _, bam_path in group_records]
                group_tmp = tempfile.mkdtemp(prefix=f"omnom_chip_broad_{group}_")
                try:
                    pooled_name = f"{group}.MACS3.relaxed_pooled"
                    pooled_broadpeak = os.path.join(group_tmp, f"{pooled_name}_peaks.broadPeak")
                    pooled_bed = os.path.join(group_tmp, f"{pooled_name}.bed")
                    pooled_cmd = build_macs3_callpeak_command(
                        treatment_bams=group_bams,
                        outdir=group_tmp,
                        name=pooled_name,
                        qvalue=relaxed_qvalue,
                        control_bam=control_bam if control_bam != "NA" else None,
                        broad_mode=True,
                        assay_type="CHIP",
                        broad_cutoff=broad_cutoff,
                        broad_min_length=broad_min_length,
                        broad_max_gap=broad_max_gap,
                    )
                    shell(pooled_cmd)
                    shell(f"cut -f1-3 {quote(pooled_broadpeak)} | sort -k1,1 -k2,2n -k3,3n > {quote(pooled_bed)}")

                    replicate_beds = []
                    for sample_name, sample_bam in group_records:
                        rep_name = f"{group}.{sample_name}.MACS3.relaxed_rep"
                        rep_broadpeak = os.path.join(group_tmp, f"{rep_name}_peaks.broadPeak")
                        rep_bed = os.path.join(group_tmp, f"{rep_name}.bed")
                        rep_cmd = build_macs3_callpeak_command(
                            treatment_bams=[sample_bam],
                            outdir=group_tmp,
                            name=rep_name,
                            qvalue=relaxed_qvalue,
                            control_bam=control_bam if control_bam != "NA" else None,
                            broad_mode=True,
                            assay_type="CHIP",
                            broad_cutoff=broad_cutoff,
                            broad_min_length=broad_min_length,
                            broad_max_gap=broad_max_gap,
                        )
                        shell(rep_cmd)
                        shell(f"cut -f1-3 {quote(rep_broadpeak)} | sort -k1,1 -k2,2n -k3,3n > {quote(rep_bed)}")
                        replicate_beds.append((sample_name, rep_bed))

                    pooled_merge_bam = os.path.join(group_tmp, f"{group}.pooled.bam")
                    shell(f"samtools merge -f {quote(pooled_merge_bam)} {' '.join(quote(x) for x in group_bams)}")
                    pooled_seed = int(hashlib.md5(f"{group}|broad|pooled".encode("utf-8")).hexdigest()[:8], 16)
                    ps1_bam, ps2_bam = make_deterministic_pseudorep_bams(pooled_merge_bam, group_tmp, f"{group}.broad.pooled", pooled_seed)
                    pseudorep_beds = []
                    for idx, ps_bam in enumerate((ps1_bam, ps2_bam), start=1):
                        ps_name = f"{group}.MACS3.relaxed_pooled_ps{idx}"
                        ps_broadpeak = os.path.join(group_tmp, f"{ps_name}_peaks.broadPeak")
                        ps_bed = os.path.join(group_tmp, f"{ps_name}.bed")
                        ps_cmd = build_macs3_callpeak_command(
                            treatment_bams=[ps_bam],
                            outdir=group_tmp,
                            name=ps_name,
                            qvalue=relaxed_qvalue,
                            control_bam=control_bam if control_bam != "NA" else None,
                            broad_mode=True,
                            assay_type="CHIP",
                            broad_cutoff=broad_cutoff,
                            broad_min_length=broad_min_length,
                            broad_max_gap=broad_max_gap,
                        )
                        shell(ps_cmd)
                        shell(f"cut -f1-3 {quote(ps_broadpeak)} | sort -k1,1 -k2,2n -k3,3n > {quote(ps_bed)}")
                        pseudorep_beds.append(ps_bed)

                    required_rep_support = max(1, int(math.ceil(float(replicate_fraction) * float(len(replicate_beds)))))
                    rep_concordant_bed = os.path.join(broad_root, f"{group}.replicate_concordant.bed")
                    rep_concordant_count = pooled_overlap_support_bed(
                        pooled_bed=pooled_bed,
                        support_beds=[bed for _, bed in replicate_beds],
                        min_support=required_rep_support,
                        output_bed=rep_concordant_bed,
                        overlap_fraction=overlap_fraction,
                    )

                    pseudorep_concordant_bed = os.path.join(broad_root, f"{group}.pseudorep_concordant.bed")
                    pseudorep_concordant_count = pooled_overlap_support_bed(
                        pooled_bed=pooled_bed,
                        support_beds=pseudorep_beds,
                        min_support=2,
                        output_bed=pseudorep_concordant_bed,
                        overlap_fraction=overlap_fraction,
                    )

                    final_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                    shell(
                        f"cat {quote(rep_concordant_bed)} {quote(pseudorep_concordant_bed)} | "
                        f"sort -k1,1 -k2,2n -k3,3n | uniq > {quote(final_bed)}"
                    )
                    final_count = int(
                        subprocess.check_output(
                            f"wc -l {quote(final_bed)} | awk '{{print $1}}'",
                            shell=True,
                            executable="/bin/bash",
                            text=True,
                        ).strip()
                    )

                    pooled_count = int(
                        subprocess.check_output(
                            f"wc -l {quote(pooled_bed)} | awk '{{print $1}}'",
                            shell=True,
                            executable="/bin/bash",
                            text=True,
                        ).strip()
                    )
                    summary_rows.append(
                        [
                            group,
                            len(replicate_beds),
                            pooled_count,
                            required_rep_support,
                            rep_concordant_count,
                            pseudorep_concordant_count,
                            final_count,
                            final_bed,
                        ]
                    )
                    log_it(
                        logfile,
                        f"Broad-domain concordance for {group}: pooled={pooled_count}, replicate_support>={required_rep_support} -> {rep_concordant_count}, pooled_pseudorep= {pseudorep_concordant_count}, final={final_count}",
                    )
                finally:
                    shutil.rmtree(group_tmp, ignore_errors=True)

            summary_path = os.path.join(broad_root, "broad_domain_concordance_summary.tsv")
            with open(summary_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    [
                        "group",
                        "num_replicates",
                        "pooled_relaxed_peaks",
                        "required_replicate_support",
                        "replicate_concordant_peaks",
                        "pooled_pseudorep_concordant_peaks",
                        "final_union_peaks",
                        "output_bed",
                    ]
                )
                writer.writerows(summary_rows)
            log_it(logfile, f"Broad-domain concordance summary: {summary_path}")

        def run_shell_jobs_parallel(job_specs, max_parallel, phase_label):
            if not job_specs:
                return
            workers = max(1, min(int(max_parallel), len(job_specs)))
            log_it(logfile, f"{phase_label}: dispatching {len(job_specs)} MACS3 jobs with parallel workers={workers}")

            def _run_one(spec):
                cmd = spec["cmd"]
                label = spec.get("label", "job")
                try:
                    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)
                    return (label, None)
                except subprocess.CalledProcessError as exc:
                    return (label, exc)

            failures = []
            with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [pool.submit(_run_one, spec) for spec in job_specs]
                for fut in concurrent.futures.as_completed(futures):
                    label, err = fut.result()
                    if err is not None:
                        failures.append((label, err))

            if failures:
                first_label, first_err = failures[0]
                raise RuntimeError(
                    f"{phase_label}: {len(failures)} MACS3 jobs failed; first failure={first_label}: {first_err}"
                )

        def sort_narrowpeak_by_pvalue(src_path, dst_path):
            shell(f"sort -k8,8nr {quote(src_path)} > {quote(dst_path)}")

        def idr_consensus_from_replicates(grouped_bams, outputfolder, assay_type, control_bam, idr_mode, pair_fraction, pairing_policy):
            if shutil.which("idr") is None:
                raise RuntimeError(
                    "Narrow peak strategy 'idr' requested, but executable 'idr' was not found on PATH."
                )
            if assay_type == "ATAC":
                log_it(logfile, "Narrow peak strategy: idr (ATAC)")
            else:
                log_it(logfile, "Narrow peak strategy: idr (ChIP narrow)")

            idr_root = os.path.join(outputfolder, f"{assay_type.lower()}_idr_peak_calling")
            os.makedirs(idr_root, exist_ok=True)
            idr_summary_rows = []

            def replicate_pairs(replicates, policy):
                if len(replicates) <= 2:
                    return list(itertools.combinations(replicates, 2))
                if policy == "anchor_vs_all":
                    ordered = sorted(replicates, key=lambda x: str(x[0]))
                    anchor = ordered[0]
                    return [(anchor, rep) for rep in ordered[1:]]
                return list(itertools.combinations(replicates, 2))

            def make_pseudorep_bams(input_bam, out_prefix, split_seed):
                return make_deterministic_pseudorep_bams(input_bam, group_tmp, out_prefix, split_seed)

            def run_idr_pair(sorted_a, sorted_b, peak_list, output_path, log_path):
                for stale_path in (output_path, log_path, f"{output_path}.png"):
                    if os.path.exists(stale_path):
                        os.remove(stale_path)
                try:
                    min_input_peaks = int(config.get("IDR_MIN_INPUT_PEAKS", 20))
                except (TypeError, ValueError):
                    min_input_peaks = 20
                min_input_peaks = max(0, min_input_peaks)
                input_counts = {
                    "sample_a": count_peak_file_rows(sorted_a),
                    "sample_b": count_peak_file_rows(sorted_b),
                    "peak_list": count_peak_file_rows(peak_list),
                }
                too_sparse = {name: count for name, count in input_counts.items() if count < min_input_peaks}
                if too_sparse:
                    reason = ", ".join(f"{name}={count}" for name, count in too_sparse.items())
                    message = (
                        f"Skipping sparse IDR comparison for {os.path.basename(output_path)}: "
                        f"{reason}; required >= {min_input_peaks} peaks."
                    )
                    log_it(logfile, message, "WARNING")
                    open(output_path, "w", encoding="utf-8").close()
                    with open(log_path, "w", encoding="utf-8") as handle:
                        handle.write(message + "\n")
                    return False
                idr_cmd = (
                    f"idr --samples {quote(sorted_a)} {quote(sorted_b)} "
                    f"--peak-list {quote(peak_list)} "
                    f"--input-file-type narrowPeak --rank p.value "
                    f"--output-file {quote(output_path)} --output-file-type narrowPeak "
                    f"--plot --log-output-file {quote(log_path)} --soft-idr-threshold 0.05"
                )
                result = subprocess.run(
                    idr_cmd,
                    shell=True,
                    executable="/bin/bash",
                    text=True,
                    capture_output=True,
                )
                if result.returncode == 0:
                    return True

                combined_output = "\n".join(
                    part for part in [result.stdout, result.stderr] if part
                )
                if os.path.exists(log_path):
                    try:
                        with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
                            combined_output += "\n" + handle.read()
                    except OSError:
                        pass
                sparse_failure = (
                    "Peak files must contain at least" in combined_output
                    or "post-merge" in combined_output
                )
                if sparse_failure:
                    message = (
                        f"Skipping sparse IDR comparison for {os.path.basename(output_path)} "
                        "after IDR reported too few post-merge peaks."
                    )
                    log_it(logfile, message, "WARNING")
                    open(output_path, "w", encoding="utf-8").close()
                    with open(log_path, "a", encoding="utf-8") as handle:
                        handle.write("\n" + message + "\n")
                        if combined_output:
                            handle.write(combined_output + "\n")
                    return False

                raise RuntimeError(
                    f"IDR failed for {os.path.basename(output_path)} with exit code {result.returncode}."
                )

            def count_peak_file_rows(path):
                count = 0
                if not os.path.exists(path):
                    return count
                with open(path, "r", encoding="utf-8", errors="replace") as handle:
                    for line in handle:
                        if line.strip() and not line.startswith("#"):
                            count += 1
                return count

            def idr_row_passes_threshold(row, threshold=0.05):
                candidate_vals = []
                for idx in (10, 11):
                    if idx < len(row):
                        value = row[idx].strip()
                        if value in {"", ".", "NA", "nan", "NaN"}:
                            continue
                        try:
                            numeric = float(value)
                        except ValueError:
                            continue
                        candidate_vals.append(numeric)
                for numeric in candidate_vals:
                    if numeric > 1.0:
                        idr_val = 10.0 ** (-numeric)
                    else:
                        idr_val = numeric
                    if idr_val <= threshold:
                        return True
                if len(row) > 4:
                    try:
                        score = float(row[4])
                        if 0.0 <= score <= 1000.0:
                            approx_idr = 10.0 ** (-(score / 125.0))
                            return approx_idr <= threshold
                    except ValueError:
                        pass
                return False

            def extract_idr_pass_bed(idr_output_path, pass_bed_path, threshold=0.05):
                pass_rows = []
                with open(idr_output_path, "r", encoding="utf-8") as in_handle:
                    for raw_line in in_handle:
                        if not raw_line or raw_line.startswith("#"):
                            continue
                        row = raw_line.rstrip("\n").split("\t")
                        if len(row) < 3:
                            continue
                        if idr_row_passes_threshold(row, threshold=threshold):
                            pass_rows.append((row[0], int(row[1]), int(row[2])))
                pass_rows = sorted(set(pass_rows), key=lambda x: (x[0], x[1], x[2]))
                with open(pass_bed_path, "w", encoding="utf-8") as out_handle:
                    for chrom, start, end in pass_rows:
                        out_handle.write(f"{chrom}\t{start}\t{end}\n")
                return len(pass_rows)

            for group in sorted(grouped_bams):
                group_records = grouped_bams[group]
                group_bams = [bam_path for _, bam_path in group_records]
                group_tmp = tempfile.mkdtemp(prefix=f"omnom_{assay_type.lower()}_idr_{group}_")
                try:
                    pooled_name = f"{group}.MACS3.pooled"
                    pooled_narrow = os.path.join(group_tmp, f"{pooled_name}_peaks.narrowPeak")
                    if assay_type == "ATAC":
                        pooled_cmd = build_macs3_callpeak_command(
                            treatment_bams=group_bams,
                            outdir=group_tmp,
                            name=pooled_name,
                            qvalue="0.01",
                            control_bam=None,
                            broad_mode=False,
                            shift=-100,
                            extsize=200,
                            assay_type=assay_type,
                        )
                    else:
                        pooled_cmd = build_macs3_callpeak_command(
                            treatment_bams=group_bams,
                            outdir=group_tmp,
                            name=pooled_name,
                            qvalue="0.01",
                            control_bam=control_bam if control_bam != "NA" else None,
                            broad_mode=False,
                            shift=None,
                            extsize=None,
                            assay_type=assay_type,
                        )
                    shell(pooled_cmd)
                    pooled_sorted = os.path.join(group_tmp, f"{pooled_name}.sorted.narrowPeak")
                    sort_narrowpeak_by_pvalue(pooled_narrow, pooled_sorted)

                    replicate_sorted = []
                    for sample_name, sample_bam in group_records:
                        rep_name = f"{group}.{sample_name}.MACS3.rep"
                        rep_narrow = os.path.join(group_tmp, f"{rep_name}_peaks.narrowPeak")
                        if assay_type == "ATAC":
                            rep_cmd = build_macs3_callpeak_command(
                                treatment_bams=[sample_bam],
                                outdir=group_tmp,
                                name=rep_name,
                                qvalue="0.01",
                                control_bam=None,
                                broad_mode=False,
                                shift=-100,
                                extsize=200,
                                assay_type=assay_type,
                            )
                        else:
                            rep_cmd = build_macs3_callpeak_command(
                                treatment_bams=[sample_bam],
                                outdir=group_tmp,
                                name=rep_name,
                                qvalue="0.01",
                                control_bam=control_bam if control_bam != "NA" else None,
                                broad_mode=False,
                                shift=None,
                                extsize=None,
                                assay_type=assay_type,
                            )
                        shell(rep_cmd)
                        rep_sorted = os.path.join(group_tmp, f"{rep_name}.sorted.narrowPeak")
                        sort_narrowpeak_by_pvalue(rep_narrow, rep_sorted)
                        replicate_sorted.append((sample_name, rep_sorted))

                    if len(replicate_sorted) < 2:
                        log_it(logfile, f"Group {group} has <2 replicates. Falling back to pooled MACS3 peaks.")
                        final_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                        shell(f"cut -f1-3 {quote(pooled_narrow)} | sort -k1,1 -k2,2n -k3,3n > {quote(final_bed)}")
                        idr_summary_rows.append([group, "pooled_fallback", "NA", "NA", "NA", "NA", final_bed, "fewer_than_two_replicates"])
                        continue

                    pair_pass_beds = []
                    rep_pairs = replicate_pairs(replicate_sorted, pairing_policy)
                    log_it(logfile, f"IDR replicate pairing policy for group {group}: {pairing_policy} ({len(rep_pairs)} pair(s))")
                    for (rep1_name, rep1_sorted), (rep2_name, rep2_sorted) in rep_pairs:
                        pair_tag = f"{group}__{rep1_name}__{rep2_name}"
                        pair_hash = hashlib.md5(pair_tag.encode("utf-8")).hexdigest()[:12]
                        idr_out = os.path.join(idr_root, f"{group}.pair.{pair_hash}.idr.narrowPeak")
                        idr_log = os.path.join(idr_root, f"{group}.pair.{pair_hash}.idr.log")
                        pair_pass_bed = os.path.join(idr_root, f"{group}.pair.{pair_hash}.idr.pass.bed")
                        if run_idr_pair(rep1_sorted, rep2_sorted, pooled_sorted, idr_out, idr_log):
                            pass_count = extract_idr_pass_bed(idr_out, pair_pass_bed, threshold=0.05)
                            pair_pass_beds.append(pair_pass_bed)
                            idr_summary_rows.append([group, "true_pair_idr", rep1_name, rep2_name, "0.05", pass_count, pair_pass_bed, "ok"])
                        else:
                            open(pair_pass_bed, "w", encoding="utf-8").close()
                            idr_summary_rows.append([group, "true_pair_idr_skipped_sparse", rep1_name, rep2_name, "0.05", 0, pair_pass_bed, "too_few_peaks_for_idr"])

                    num_pairs = len(pair_pass_beds)
                    min_pairs = max(1, int(math.ceil(float(pair_fraction) * float(num_pairs))))
                    final_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                    if num_pairs == 0:
                        shell(f"cut -f1-3 {quote(pooled_narrow)} | sort -k1,1 -k2,2n -k3,3n > {quote(final_bed)}")
                        retained_count = int(
                            subprocess.check_output(
                                f"wc -l {quote(final_bed)} | awk '{{print $1}}'",
                                shell=True,
                                executable="/bin/bash",
                                text=True,
                            ).strip()
                        )
                        log_it(
                            logfile,
                            f"Group {group} had no usable true-replicate IDR comparisons. Falling back to pooled MACS3 peaks.",
                            "WARNING",
                        )
                        idr_summary_rows.append([group, "pooled_fallback_no_usable_idr", "all", "all", "NA", retained_count, final_bed, "all_true_pair_idr_comparisons_skipped"])
                    elif num_pairs == 1:
                        shell(f"cp {quote(pair_pass_beds[0])} {quote(final_bed)}")
                    else:
                        pair_join = " ".join(quote(p) for p in pair_pass_beds)
                        shell(
                            f"bedtools multiinter -i {pair_join} | "
                            f"awk 'BEGIN{{{{OFS=\"\\t\"}}}} $4>={min_pairs} {{{{print $1,$2,$3}}}}' | "
                            f"sort -k1,1 -k2,2n -k3,3n > {quote(final_bed)}"
                        )
                        retained_count = int(
                            subprocess.check_output(
                                f"wc -l {quote(final_bed)} | awk '{{print $1}}'",
                                shell=True,
                                executable="/bin/bash",
                                text=True,
                            ).strip()
                        )
                        idr_summary_rows.append([group, "consensus_idr", "all", "all", f"0.05;pair_fraction={pair_fraction};min_pairs={min_pairs}", retained_count, final_bed, "ok"])
                    if num_pairs == 1:
                        retained_count = int(
                            subprocess.check_output(
                                f"wc -l {quote(final_bed)} | awk '{{print $1}}'",
                                shell=True,
                                executable="/bin/bash",
                                text=True,
                            ).strip()
                        )
                        idr_summary_rows.append([group, "consensus_idr", "all", "all", f"0.05;pair_fraction={pair_fraction};min_pairs={min_pairs}", retained_count, final_bed, "ok"])

                    if idr_mode == "encode":
                        pooled_bam = os.path.join(group_tmp, f"{group}.pooled.bam")
                        shell(f"samtools merge -f {quote(pooled_bam)} {' '.join(quote(x) for x in group_bams)}")
                        pooled_seed = int(hashlib.md5(f"{group}|pooled".encode("utf-8")).hexdigest()[:8], 16)
                        ps1_bam, ps2_bam = make_pseudorep_bams(pooled_bam, f"{group}.pooled", pooled_seed)
                        pooled_ps = []
                        for idx, ps_bam in enumerate((ps1_bam, ps2_bam), start=1):
                            ps_name = f"{group}.pooled.ps{idx}"
                            ps_narrow = os.path.join(group_tmp, f"{ps_name}_peaks.narrowPeak")
                            if assay_type == "ATAC":
                                ps_cmd = build_macs3_callpeak_command([ps_bam], group_tmp, ps_name, "0.01", None, False, -100, 200, assay_type=assay_type)
                            else:
                                ps_cmd = build_macs3_callpeak_command([ps_bam], group_tmp, ps_name, "0.01", control_bam if control_bam != "NA" else None, False, None, None, assay_type=assay_type)
                            shell(ps_cmd)
                            ps_sorted = os.path.join(group_tmp, f"{ps_name}.sorted.narrowPeak")
                            sort_narrowpeak_by_pvalue(ps_narrow, ps_sorted)
                            pooled_ps.append(ps_sorted)
                        pooled_idr_out = os.path.join(idr_root, f"{group}.pooled_pseudorep.idr.narrowPeak")
                        pooled_idr_log = os.path.join(idr_root, f"{group}.pooled_pseudorep.idr.log")
                        pooled_pass_bed = os.path.join(idr_root, f"{group}.pooled_pseudorep.idr.pass.bed")
                        if run_idr_pair(pooled_ps[0], pooled_ps[1], pooled_sorted, pooled_idr_out, pooled_idr_log):
                            pooled_pass_count = extract_idr_pass_bed(pooled_idr_out, pooled_pass_bed, threshold=0.05)
                        else:
                            open(pooled_pass_bed, "w", encoding="utf-8").close()
                            pooled_pass_count = 0
                        pooled_note = "ok" if pooled_pass_count > 0 else "zero_or_skipped_sparse"
                        idr_summary_rows.append([group, "pooled_pseudorep_idr", "ps1", "ps2", "0.05", pooled_pass_count, pooled_pass_bed, pooled_note])

                        for sample_name, sample_bam in group_records:
                            self_seed = int(hashlib.md5(f"{group}|{sample_name}|self".encode("utf-8")).hexdigest()[:8], 16)
                            sps1_bam, sps2_bam = make_pseudorep_bams(sample_bam, f"{group}.{sample_name}.self", self_seed)
                            self_sorted = []
                            for idx, ps_bam in enumerate((sps1_bam, sps2_bam), start=1):
                                ps_name = f"{group}.{sample_name}.self.ps{idx}"
                                ps_narrow = os.path.join(group_tmp, f"{ps_name}_peaks.narrowPeak")
                                if assay_type == "ATAC":
                                    ps_cmd = build_macs3_callpeak_command([ps_bam], group_tmp, ps_name, "0.01", None, False, -100, 200, assay_type=assay_type)
                                else:
                                    ps_cmd = build_macs3_callpeak_command([ps_bam], group_tmp, ps_name, "0.01", control_bam if control_bam != "NA" else None, False, None, None, assay_type=assay_type)
                                shell(ps_cmd)
                                ps_sorted = os.path.join(group_tmp, f"{ps_name}.sorted.narrowPeak")
                                sort_narrowpeak_by_pvalue(ps_narrow, ps_sorted)
                                self_sorted.append(ps_sorted)
                            self_idr_out = os.path.join(idr_root, f"{group}.{sample_name}.self_pseudorep.idr.narrowPeak")
                            self_idr_log = os.path.join(idr_root, f"{group}.{sample_name}.self_pseudorep.idr.log")
                            self_pass_bed = os.path.join(idr_root, f"{group}.{sample_name}.self_pseudorep.idr.pass.bed")
                            if run_idr_pair(self_sorted[0], self_sorted[1], pooled_sorted, self_idr_out, self_idr_log):
                                self_pass_count = extract_idr_pass_bed(self_idr_out, self_pass_bed, threshold=0.05)
                            else:
                                open(self_pass_bed, "w", encoding="utf-8").close()
                                self_pass_count = 0
                            self_note = "ok" if self_pass_count > 0 else "zero_or_skipped_sparse"
                            idr_summary_rows.append([group, "self_pseudorep_idr", sample_name, sample_name, "0.05", self_pass_count, self_pass_bed, self_note])
                finally:
                    shutil.rmtree(group_tmp, ignore_errors=True)

            summary_path = os.path.join(idr_root, "idr_selected_peaks.tsv")
            with open(summary_path, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["group", "strategy", "replicate_1", "replicate_2", "idr_threshold", "retained_peaks", "output_bed", "note"])
                writer.writerows(idr_summary_rows)
            log_it(logfile, f"IDR peak summary: {summary_path}")

        def call_peaks(logfile, thetype, inputfolder1, outputfolder, broad_mode, input_sample):
            if thetype == "RNA":
                log_it(logfile, "Not a ChIP- or ATAC-seq experiment, skipping this step...")
                return
            
            log_once(logfile, "step10.calling", "Calling peaks...", f"EXECUTING STEP {master_config['callpeaks_rule_num']}")
            log_once(logfile, "step10.type", f"The type = {thetype}")
            # Report version
            macs3_version = subprocess.check_output(["macs3", "--version"], stderr=subprocess.STDOUT)
            log_once(logfile, "step10.macs3_version", "\n"+macs3_version.decode("utf-8"), "MACS3 VERSION")

            if thetype == "CHIP":
                #If ChIP, call peaks
                log_it(logfile, f"Finding ChIP enriched regions for {chip_style_label(broad_mode)}...")
                if broad_mode == "domain":
                    log_it(logfile, "NSC and RSC are reported for broad histone marks, but they are typically less informative there than for TF / narrow peaks.")
                elif broad_mode == "genebody":
                    log_it(logfile, "Gene-body mode quantifies reads over annotated gene bodies instead of calling peaks.")
                elif broad_mode == "diffuse":
                    log_it(logfile, "Diffuse mode quantifies reads over fixed genomic bins instead of calling peaks.")

                log_it(logfile, "Calling peaks with MACS3...")
                
                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder1, ".filtered.bam")

                grouped_bams = grouped_bams_by_sample_type(inputfolder1, ".filtered.bam")
                chip_groups = sorted(grouped_bams)

                if broad_mode == "genebody":
                    build_chip_genebody_feature_sets(
                        gtf_file=params.gtf_file,
                        outputfolder=outputfolder,
                        groups=chip_groups,
                    )
                    return
                if broad_mode == "diffuse":
                    chip_diffuse_bin_size = int(config.get("CHIP_DIFFUSE_BIN_SIZE", 10000))
                    chip_diffuse_standard_chroms_only = bool(config.get("CHIP_DIFFUSE_STANDARD_CHROMS_ONLY", True))
                    chip_diffuse_exclude_chrm = bool(config.get("CHIP_DIFFUSE_EXCLUDE_CHRM", True))
                    blacklist_bed = maybe_resolve_blacklist_bed(config["THEGENOME"], outputfolder, "step10_diffuse")
                    reference_bam = None
                    for group in chip_groups:
                        if grouped_bams.get(group):
                            reference_bam = grouped_bams[group][0][1]
                            break
                    if reference_bam is None:
                        raise FileNotFoundError("No ChIP BAM files available to derive diffuse bin chromosome sizes.")
                    log_it(
                        logfile,
                        (
                            "Diffuse mode parameters: "
                            f"bin_size={chip_diffuse_bin_size}, "
                            f"standard_chroms_only={chip_diffuse_standard_chroms_only}, "
                            f"exclude_chrm={chip_diffuse_exclude_chrm}, "
                            f"blacklist_bed={blacklist_bed if blacklist_bed else 'NA'}"
                        ),
                    )
                    build_chip_diffuse_bin_feature_sets(
                        reference_bam=reference_bam,
                        genome_name=config["THEGENOME"],
                        outputfolder=outputfolder,
                        groups=chip_groups,
                        bin_size=chip_diffuse_bin_size,
                        keep_standard=chip_diffuse_standard_chroms_only,
                        drop_chrm=chip_diffuse_exclude_chrm,
                        blacklist_bed=blacklist_bed,
                    )
                    return
                if broad_mode == "domain":
                    chip_broad_qvalue = float(config.get("CHIP_BROAD_QVALUE", 0.05))
                    chip_broad_cutoff = float(config.get("CHIP_BROAD_CUTOFF", 0.1))
                    chip_broad_min_length = config.get("CHIP_BROAD_MIN_LENGTH", "NA")
                    chip_broad_max_gap = config.get("CHIP_BROAD_MAX_GAP", "NA")
                    chip_broad_replicate_fraction = float(config.get("CHIP_BROAD_REPLICATE_FRACTION", 1.0))
                    chip_broad_overlap_fraction = float(config.get("CHIP_BROAD_OVERLAP_FRACTION", 0.5))
                    log_it(
                        logfile,
                        (
                            "Broad-domain mode parameters: "
                            f"qvalue={chip_broad_qvalue}, broad_cutoff={chip_broad_cutoff}, "
                            f"min_length={chip_broad_min_length}, max_gap={chip_broad_max_gap}, "
                            f"replicate_fraction={chip_broad_replicate_fraction}, overlap_fraction={chip_broad_overlap_fraction}"
                        ),
                    )
                    broad_domain_consensus_from_replicates(
                        grouped_bams=grouped_bams,
                        outputfolder=outputfolder,
                        control_bam=input_sample,
                        replicate_fraction=chip_broad_replicate_fraction,
                        overlap_fraction=chip_broad_overlap_fraction,
                        relaxed_qvalue=chip_broad_qvalue,
                        broad_cutoff=chip_broad_cutoff,
                        broad_min_length=chip_broad_min_length,
                        broad_max_gap=chip_broad_max_gap,
                    )
                else:
                    narrow_peak_strategy = str(config.get("NARROW_PEAK_STRATEGY", "idr")).strip().lower()
                    if narrow_peak_strategy not in {"idr", "macs3"}:
                        narrow_peak_strategy = "idr"
                    if narrow_peak_strategy == "idr":
                        idr_mode = str(config.get("IDR_MODE", "encode")).strip().lower()
                        if idr_mode not in {"basic", "encode"}:
                            idr_mode = "basic"
                        pair_fraction = config.get("IDR_PAIR_FRACTION", 0.5)
                        pairing_policy = str(config.get("IDR_PAIRING_POLICY", "all_pairs")).strip().lower()
                        if pairing_policy not in {"all_pairs", "anchor_vs_all"}:
                            pairing_policy = "all_pairs"
                        idr_consensus_from_replicates(
                            grouped_bams=grouped_bams,
                            outputfolder=outputfolder,
                            assay_type="CHIP",
                            control_bam=input_sample,
                            idr_mode=idr_mode,
                            pair_fraction=pair_fraction,
                            pairing_policy=pairing_policy,
                        )
                        shell(f"""cat {outputfolder}/*.MACS3.optimized.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {outputfolder}/all_groups.merged_peaks.bed """)
                        log_it(logfile, f"Total merged peaks: {subprocess.getoutput(f'wc -l {outputfolder}/all_groups.merged_peaks.bed').strip().split()[0]}")
                        return
                    opt_mode = str(config.get("CHIP_PEAK_OPT_MODE", "fast")).strip().lower()
                    if opt_mode not in {"none", "fast", "full"}:
                        opt_mode = "fast"
                    log_it(logfile, f"ChIP narrow peak optimization mode: {opt_mode}")
                    log_it(logfile, "ENCODE-style note: IDR is ideal for TF ChIP; this optimizer selects MACS3 parameters for downstream consistency.")

                    if opt_mode == "none":
                        candidate_grid = [("q0p01.model", "0.01", None, None)]
                    elif opt_mode == "fast":
                        candidate_grid = [
                            ("q0p01.model", "0.01", None, None),
                            ("q0p001.model", "0.001", None, None),
                        ]
                    else:
                        candidate_grid = [
                            ("q0p01.model", "0.01", None, None),
                            ("q0p001.model", "0.001", None, None),
                            ("q0p01.shiftm75.ext150", "0.01", -75, 150),
                            ("q0p001.shiftm75.ext150", "0.001", -75, 150),
                            ("q0p01.shiftm100.ext200", "0.01", -100, 200),
                            ("q0p001.shiftm100.ext200", "0.001", -100, 200),
                        ]

                    blacklist_bed = maybe_resolve_blacklist_bed(config["THEGENOME"], outputfolder, "step10_chip_narrow")
                    optimization_dir = os.path.join(outputfolder, "chip_narrow_peak_call_optimization")
                    os.makedirs(optimization_dir, exist_ok=True)
                    score_rows = []
                    selected = {}
                    total_planned_runs = 0
                    total_completed_runs = 0

                    for group in chip_groups:
                        bams = [bam_path for _, bam_path in grouped_bams[group]]
                        group_tmp = tempfile.mkdtemp(prefix=f"omnom_chip_opt_{group}_")
                        group_rows = []
                        try:
                            if opt_mode == "none":
                                cand_name, qval, shift, extsize = candidate_grid[0]
                                name = f"{group}.MACS3.{cand_name}.group"
                                group_narrow = os.path.join(group_tmp, f"{name}_peaks.narrowPeak")
                                group_bed = os.path.join(group_tmp, f"{group}.MACS3.{cand_name}.group.bed")
                                cmd = build_macs3_callpeak_command(
                                    treatment_bams=bams,
                                    outdir=group_tmp,
                                    name=name,
                                    qvalue=qval,
                                    control_bam=input_sample if input_sample != "NA" else None,
                                    broad_mode=False,
                                    shift=shift,
                                    extsize=extsize,
                                )
                                shell(cmd)
                                total_planned_runs += 1
                                total_completed_runs += 1
                                shell(f"cut -f1-3 {quote(group_narrow)} | sort -k1,1 -k2,2n -k3,3n > {quote(group_bed)}")
                                final_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                                shell(f"cp {quote(group_bed)} {quote(final_bed)}")
                                peak_n, median_width = bed_peak_count_and_widths(group_bed)
                                selected[group] = {
                                    "group": group, "candidate": cand_name, "qvalue": float(qval),
                                    "shift": "NA", "extsize": "NA", "composite_score": 0.0,
                                    "frip": float("nan"), "replicate_jaccard": float("nan"),
                                    "blacklist_fraction": float("nan"), "peak_count": int(peak_n),
                                    "median_peak_width_bp": float(median_width),
                                }
                                continue

                            candidate_outputs = {}
                            macs_jobs = []
                            for cand_name, qval, shift, extsize in candidate_grid:
                                group_name = f"{group}.MACS3.{cand_name}.group"
                                group_narrow = os.path.join(group_tmp, f"{group_name}_peaks.narrowPeak")
                                group_bed = os.path.join(group_tmp, f"{group}.MACS3.{cand_name}.group.bed")
                                candidate_outputs[cand_name] = {
                                    "group_narrow": group_narrow,
                                    "group_bed": group_bed,
                                    "sample_beds": [],
                                    "qval": qval,
                                    "shift": shift,
                                    "extsize": extsize,
                                }
                                macs_jobs.append(
                                    {
                                        "label": f"{group}:{cand_name}:group",
                                        "cmd": build_macs3_callpeak_command(
                                            treatment_bams=bams, outdir=group_tmp, name=group_name, qvalue=qval,
                                            control_bam=input_sample if input_sample != "NA" else None,
                                            broad_mode=False, shift=shift, extsize=extsize,
                                        ),
                                    }
                                )
                                for sample, sample_bam in grouped_bams[group]:
                                    sample_name = f"{group}.{sample}.MACS3.{cand_name}.sample"
                                    sample_narrow = os.path.join(group_tmp, f"{sample_name}_peaks.narrowPeak")
                                    sample_bed = os.path.join(group_tmp, f"{group}.{sample}.MACS3.{cand_name}.sample.bed")
                                    macs_jobs.append(
                                        {
                                            "label": f"{group}:{cand_name}:{sample}",
                                            "cmd": build_macs3_callpeak_command(
                                                treatment_bams=[sample_bam], outdir=group_tmp, name=sample_name, qvalue=qval,
                                                control_bam=input_sample if input_sample != "NA" else None,
                                                broad_mode=False, shift=shift, extsize=extsize,
                                            ),
                                        }
                                    )
                                    candidate_outputs[cand_name]["sample_beds"].append((sample_narrow, sample_bed))

                            planned_for_group = len(macs_jobs)
                            total_planned_runs += planned_for_group
                            parallel_workers = min(max(1, int(threads) - 1), planned_for_group)
                            run_shell_jobs_parallel(macs_jobs, max_parallel=parallel_workers, phase_label=f"ChIP narrow optimization {group}")
                            total_completed_runs += planned_for_group

                            for cand_name, out in candidate_outputs.items():
                                shell(f"cut -f1-3 {quote(out['group_narrow'])} | sort -k1,1 -k2,2n -k3,3n > {quote(out['group_bed'])}")
                                for sample_narrow, sample_bed in out["sample_beds"]:
                                    shell(f"cut -f1-3 {quote(sample_narrow)} | sort -k1,1 -k2,2n -k3,3n > {quote(sample_bed)}")

                            for cand_name, out in candidate_outputs.items():
                                qval = out["qval"]
                                group_bed = out["group_bed"]
                                sample_beds = [x[1] for x in out["sample_beds"]]
                                frip, reads_in_peaks, mapped_total = frip_for_group(group_bed, bams)
                                jacc = mean_pairwise_jaccard(sample_beds)
                                peak_n, median_width = bed_peak_count_and_widths(group_bed)
                                black_frac = blacklist_overlap_fraction(group_bed, blacklist_bed)
                                group_rows.append(
                                    {
                                        "group": group, "candidate": cand_name, "qvalue": float(qval),
                                        "shift": int(out["shift"]) if out["shift"] is not None else 0,
                                        "extsize": int(out["extsize"]) if out["extsize"] is not None else 0,
                                        "frip": float(frip), "reads_in_peaks": int(reads_in_peaks),
                                        "mapped_reads": int(mapped_total), "replicate_jaccard": float(jacc),
                                        "peak_count": int(peak_n), "median_peak_width_bp": float(median_width),
                                        "blacklist_fraction": float(black_frac), "group_bed": group_bed,
                                    }
                                )

                            frip_z = zscore_map({r["candidate"]: r["frip"] for r in group_rows}, higher_is_better=True)
                            jacc_z = zscore_map({r["candidate"]: r["replicate_jaccard"] for r in group_rows}, higher_is_better=True)
                            black_z = zscore_map({r["candidate"]: r["blacklist_fraction"] for r in group_rows}, higher_is_better=False)
                            width_penalty_map = {r["candidate"]: abs(math.log10(max(r["median_peak_width_bp"], 1) / 250.0)) for r in group_rows}
                            width_z = zscore_map(width_penalty_map, higher_is_better=False)
                            for row in group_rows:
                                c = row["candidate"]
                                composite = (0.35 * frip_z[c]) + (0.45 * jacc_z[c]) + (0.15 * black_z[c]) + (0.05 * width_z[c])
                                row["score_frip_z"] = frip_z[c]
                                row["score_jaccard_z"] = jacc_z[c]
                                row["score_blacklist_z"] = black_z[c]
                                row["score_width_z"] = width_z[c]
                                row["composite_score"] = float(composite)
                            best_row = sorted(group_rows, key=lambda r: r["composite_score"], reverse=True)[0]
                            selected[group] = best_row
                            final_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                            shell(f"cp {quote(best_row['group_bed'])} {quote(final_bed)}")
                            score_rows.extend(group_rows)
                        finally:
                            shutil.rmtree(group_tmp, ignore_errors=True)

                    score_tsv = os.path.join(optimization_dir, "candidate_scores.tsv")
                    with open(score_tsv, "w", newline="") as handle:
                        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                        writer.writerow(["group", "candidate", "qvalue", "shift", "extsize", "frip", "reads_in_peaks", "mapped_reads", "replicate_jaccard", "peak_count", "median_peak_width_bp", "blacklist_fraction", "score_frip_z", "score_jaccard_z", "score_blacklist_z", "score_width_z", "composite_score"])
                        if opt_mode == "none":
                            for group in sorted(selected):
                                row = selected[group]
                                writer.writerow([row["group"], row["candidate"], row["qvalue"], row["shift"], row["extsize"], "NA", "NA", "NA", "NA", row["peak_count"], row["median_peak_width_bp"], "NA", "NA", "NA", "NA", "NA", row["composite_score"]])
                        else:
                            for row in score_rows:
                                writer.writerow([row["group"], row["candidate"], row["qvalue"], row["shift"], row["extsize"], row["frip"], row["reads_in_peaks"], row["mapped_reads"], row["replicate_jaccard"], row["peak_count"], row["median_peak_width_bp"], row["blacklist_fraction"], row["score_frip_z"], row["score_jaccard_z"], row["score_blacklist_z"], row["score_width_z"], row["composite_score"]])
                    with open(os.path.join(optimization_dir, "selected_parameters.tsv"), "w", newline="") as handle:
                        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                        writer.writerow(["group", "candidate", "qvalue", "shift", "extsize", "composite_score", "frip", "replicate_jaccard", "blacklist_fraction", "peak_count"])
                        for group in sorted(selected):
                            row = selected[group]
                            writer.writerow([group, row["candidate"], row["qvalue"], row["shift"], row["extsize"], row["composite_score"], row["frip"], row["replicate_jaccard"], row["blacklist_fraction"], row["peak_count"]])
                    if opt_mode != "none":
                        write_optimization_plots(score_rows, optimization_dir)
                    log_it(logfile, f"ChIP narrow optimization run summary: planned={total_planned_runs}, completed={total_completed_runs}, slurm_cpus_per_task={threads}, worker_cap={max(1, int(threads) - 1)}")

                shell(f"""cat {outputfolder}/*.MACS3.optimized.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {outputfolder}/all_groups.merged_peaks.bed """)
                log_it(logfile, f"Total merged peaks: {subprocess.getoutput(f'wc -l {outputfolder}/all_groups.merged_peaks.bed').strip().split()[0]}")
            else:
                # If ATAC, call open regions
                log_it(logfile, "Calling ATAC open chromatin peaks...")

                # Sanity check the working dir
                sanity_check_dir(logfile, inputfolder1, ".sorted.dups_marked.filtered.bam")

                grouped_bams = grouped_bams_by_sample_type(inputfolder1, ".sorted.dups_marked.filtered.bam")
                atac_groups = sorted(grouped_bams)
                narrow_peak_strategy = str(config.get("NARROW_PEAK_STRATEGY", "idr")).strip().lower()
                if narrow_peak_strategy not in {"idr", "macs3"}:
                    narrow_peak_strategy = "idr"
                if narrow_peak_strategy == "idr":
                    idr_mode = str(config.get("IDR_MODE", "encode")).strip().lower()
                    if idr_mode not in {"basic", "encode"}:
                        idr_mode = "basic"
                    pair_fraction = config.get("IDR_PAIR_FRACTION", 0.5)
                    pairing_policy = str(config.get("IDR_PAIRING_POLICY", "all_pairs")).strip().lower()
                    if pairing_policy not in {"all_pairs", "anchor_vs_all"}:
                        pairing_policy = "all_pairs"
                    idr_consensus_from_replicates(
                        grouped_bams=grouped_bams,
                        outputfolder=outputfolder,
                        assay_type="ATAC",
                        control_bam="NA",
                        idr_mode=idr_mode,
                        pair_fraction=pair_fraction,
                        pairing_policy=pairing_policy,
                    )
                    shell(f"""cat {outputfolder}/*.MACS3.optimized.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {outputfolder}/all_groups.merged_peaks.bed """)
                    log_it(logfile, f"Total merged peaks: {subprocess.getoutput(f'wc -l {outputfolder}/all_groups.merged_peaks.bed').strip().split()[0]}")
                    return
                opt_mode = str(config.get("ATAC_PEAK_OPT_MODE", "fast")).strip().lower()
                if opt_mode not in {"none", "fast", "full"}:
                    opt_mode = "fast"
                log_it(logfile, f"ATAC peak optimization mode: {opt_mode}")

                if opt_mode == "none":
                    candidate_grid = [("q0p01.shiftm100.ext200", "0.01", -100, 200)]
                elif opt_mode == "fast":
                    candidate_grid = [
                        ("q0p01.shiftm100.ext200", "0.01", -100, 200),
                        ("q0p01.shiftm75.ext150", "0.01", -75, 150),
                    ]
                else:
                    candidate_grid = [
                        ("q0p01.shiftm100.ext200", "0.01", -100, 200),
                        ("q0p05.shiftm100.ext200", "0.05", -100, 200),
                        ("q0p01.shiftm75.ext150", "0.01", -75, 150),
                        ("q0p05.shiftm75.ext150", "0.05", -75, 150),
                        ("q0p01.shiftm50.ext100", "0.01", -50, 100),
                        ("q0p05.shiftm50.ext100", "0.05", -50, 100),
                    ]
                blacklist_bed = maybe_resolve_blacklist_bed(config["THEGENOME"], outputfolder, "step10_atac")
                if blacklist_bed:
                    log_it(logfile, f"ATAC optimization will include blacklist-overlap penalty using: {blacklist_bed}")
                else:
                    log_it(logfile, "ATAC optimization blacklist BED not available; blacklist penalty set to 0.")

                optimization_dir = os.path.join(outputfolder, "atac_peak_call_optimization")
                os.makedirs(optimization_dir, exist_ok=True)
                score_rows = []
                selected = {}
                total_planned_runs = 0
                total_completed_runs = 0

                for group in atac_groups:
                    bams = [bam_path for _, bam_path in grouped_bams[group]]
                    log_it(logfile, f"Optimizing ATAC peak caller settings for group: {group}")
                    log_it(logfile, f"Files in group: {', '.join(bams)}")
                    group_tmp = tempfile.mkdtemp(prefix=f"omnom_atac_opt_{group}_")
                    group_rows = []
                    try:
                        if opt_mode == "none":
                            cand_name, qval, shift, extsize = candidate_grid[0]
                            group_prefix = f"{group}.MACS3.{cand_name}.group"
                            group_narrow = os.path.join(group_tmp, f"{group_prefix}_peaks.narrowPeak")
                            group_bed = os.path.join(group_tmp, f"{group}.MACS3.{cand_name}.group.bed")
                            group_cmd = build_macs3_callpeak_command(
                                treatment_bams=bams,
                                outdir=group_tmp,
                                name=group_prefix,
                                qvalue=qval,
                                control_bam=None,
                                broad_mode=False,
                                shift=shift,
                                extsize=extsize,
                                assay_type="ATAC",
                            )
                            shell(group_cmd)
                            total_planned_runs += 1
                            total_completed_runs += 1
                            shell(
                                f"cut -f1-3 {quote(group_narrow)} | sort -k1,1 -k2,2n -k3,3n > {quote(group_bed)}"
                            )
                            final_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                            shell(f"cp {quote(group_bed)} {quote(final_bed)}")
                            peak_n, median_width = bed_peak_count_and_widths(group_bed)
                            selected[group] = {
                                "group": group,
                                "candidate": cand_name,
                                "qvalue": float(qval),
                                "shift": int(shift),
                                "extsize": int(extsize),
                                "composite_score": 0.0,
                                "frip": float("nan"),
                                "replicate_jaccard": float("nan"),
                                "blacklist_fraction": float("nan"),
                                "peak_count": int(peak_n),
                                "median_peak_width_bp": float(median_width),
                            }
                            log_it(
                                logfile,
                                (
                                    f"ATAC optimization none-mode selected for {group}: {cand_name} "
                                    f"(q={qval}, shift={shift}, extsize={extsize}, peaks={peak_n})"
                                ),
                            )
                            continue

                        candidate_outputs = {}
                        macs_jobs = []
                        for cand_name, qval, shift, extsize in candidate_grid:
                            group_prefix = f"{group}.MACS3.{cand_name}.group"
                            group_narrow = os.path.join(group_tmp, f"{group_prefix}_peaks.narrowPeak")
                            group_bed = os.path.join(group_tmp, f"{group}.MACS3.{cand_name}.group.bed")
                            candidate_outputs[cand_name] = {
                                "group_narrow": group_narrow,
                                "group_bed": group_bed,
                                "sample_beds": [],
                                "qval": qval,
                                "shift": shift,
                                "extsize": extsize,
                            }
                            macs_jobs.append(
                                {
                                    "label": f"{group}:{cand_name}:group",
                                    "cmd": build_macs3_callpeak_command(
                                        treatment_bams=bams,
                                        outdir=group_tmp,
                                        name=group_prefix,
                                        qvalue=qval,
                                        control_bam=None,
                                        broad_mode=False,
                                        shift=shift,
                                        extsize=extsize,
                                        assay_type="ATAC",
                                    ),
                                }
                            )

                            for sample, sample_bam in grouped_bams[group]:
                                sample_prefix = f"{group}.{sample}.MACS3.{cand_name}.sample"
                                sample_narrow = os.path.join(group_tmp, f"{sample_prefix}_peaks.narrowPeak")
                                sample_bed = os.path.join(group_tmp, f"{group}.{sample}.MACS3.{cand_name}.sample.bed")
                                macs_jobs.append(
                                    {
                                        "label": f"{group}:{cand_name}:{sample}",
                                        "cmd": build_macs3_callpeak_command(
                                            treatment_bams=[sample_bam],
                                            outdir=group_tmp,
                                            name=sample_prefix,
                                            qvalue=qval,
                                            control_bam=None,
                                            broad_mode=False,
                                            shift=shift,
                                            extsize=extsize,
                                            assay_type="ATAC",
                                        ),
                                    }
                                )
                                candidate_outputs[cand_name]["sample_beds"].append((sample_narrow, sample_bed))

                        planned_for_group = len(macs_jobs)
                        total_planned_runs += planned_for_group
                        parallel_workers = min(max(1, int(threads) - 1), planned_for_group)
                        log_it(
                            logfile,
                            f"ATAC optimization group {group}: planned MACS3 runs={planned_for_group}, workers={parallel_workers}",
                        )
                        run_shell_jobs_parallel(
                            macs_jobs,
                            max_parallel=parallel_workers,
                            phase_label=f"ATAC optimization {group}",
                        )
                        total_completed_runs += planned_for_group

                        for cand_name, out in candidate_outputs.items():
                            shell(
                                f"cut -f1-3 {quote(out['group_narrow'])} | sort -k1,1 -k2,2n -k3,3n > {quote(out['group_bed'])}"
                            )
                            for sample_narrow, sample_bed in out["sample_beds"]:
                                shell(
                                    f"cut -f1-3 {quote(sample_narrow)} | sort -k1,1 -k2,2n -k3,3n > {quote(sample_bed)}"
                                )

                        for cand_name, out in candidate_outputs.items():
                            qval = out["qval"]
                            shift = out["shift"]
                            extsize = out["extsize"]
                            group_bed = out["group_bed"]
                            sample_beds = [x[1] for x in out["sample_beds"]]
                            frip, reads_in_peaks, mapped_total = frip_for_group(group_bed, bams)
                            jacc = mean_pairwise_jaccard(sample_beds)
                            peak_n, median_width = bed_peak_count_and_widths(group_bed)
                            black_frac = blacklist_overlap_fraction(group_bed, blacklist_bed)

                            group_rows.append(
                                {
                                    "group": group,
                                    "candidate": cand_name,
                                    "qvalue": float(qval),
                                    "shift": int(shift),
                                    "extsize": int(extsize),
                                    "frip": float(frip),
                                    "reads_in_peaks": int(reads_in_peaks),
                                    "mapped_reads": int(mapped_total),
                                    "replicate_jaccard": float(jacc),
                                    "peak_count": int(peak_n),
                                    "median_peak_width_bp": float(median_width),
                                    "blacklist_fraction": float(black_frac),
                                    "group_bed": group_bed,
                                }
                            )

                        frip_z = zscore_map({r["candidate"]: r["frip"] for r in group_rows}, higher_is_better=True)
                        jacc_z = zscore_map({r["candidate"]: r["replicate_jaccard"] for r in group_rows}, higher_is_better=True)
                        black_z = zscore_map({r["candidate"]: r["blacklist_fraction"] for r in group_rows}, higher_is_better=False)
                        width_penalty_map = {
                            r["candidate"]: abs(math.log10(max(r["median_peak_width_bp"], 1) / 250.0))
                            for r in group_rows
                        }
                        width_z = zscore_map(width_penalty_map, higher_is_better=False)

                        for row in group_rows:
                            c = row["candidate"]
                            composite = (
                                0.40 * frip_z[c]
                                + 0.35 * jacc_z[c]
                                + 0.15 * black_z[c]
                                + 0.10 * width_z[c]
                            )
                            row["score_frip_z"] = frip_z[c]
                            row["score_jaccard_z"] = jacc_z[c]
                            row["score_blacklist_z"] = black_z[c]
                            row["score_width_z"] = width_z[c]
                            row["composite_score"] = float(composite)

                        best_row = sorted(group_rows, key=lambda r: r["composite_score"], reverse=True)[0]
                        selected[group] = best_row
                        final_bed = os.path.join(outputfolder, f"{group}.MACS3.optimized.bed")
                        shell(f"cp {quote(best_row['group_bed'])} {quote(final_bed)}")
                        log_it(
                            logfile,
                            (
                                f"ATAC optimizer selected for {group}: {best_row['candidate']} "
                                f"(q={best_row['qvalue']}, shift={best_row['shift']}, extsize={best_row['extsize']}, "
                                f"FRiP={best_row['frip']:.4f}, Jaccard={best_row['replicate_jaccard']:.4f}, "
                                f"blacklist_fraction={best_row['blacklist_fraction']:.4f}, peaks={best_row['peak_count']})"
                            ),
                        )
                        score_rows.extend(group_rows)
                    finally:
                        shutil.rmtree(group_tmp, ignore_errors=True)

                log_it(
                    logfile,
                    f"ATAC optimization run summary: planned MACS3 runs={total_planned_runs}, completed={total_completed_runs}, slurm_cpus_per_task={threads}, worker_cap={max(1, int(threads) - 1)}",
                )

                score_tsv = os.path.join(optimization_dir, "candidate_scores.tsv")
                with open(score_tsv, "w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(
                        [
                            "group",
                            "candidate",
                            "qvalue",
                            "shift",
                            "extsize",
                            "frip",
                            "reads_in_peaks",
                            "mapped_reads",
                            "replicate_jaccard",
                            "peak_count",
                            "median_peak_width_bp",
                            "blacklist_fraction",
                            "score_frip_z",
                            "score_jaccard_z",
                            "score_blacklist_z",
                            "score_width_z",
                            "composite_score",
                        ]
                    )
                    if opt_mode == "none":
                        for group in sorted(selected):
                            row = selected[group]
                            writer.writerow(
                                [
                                    row["group"],
                                    row["candidate"],
                                    row["qvalue"],
                                    row["shift"],
                                    row["extsize"],
                                    "NA",
                                    "NA",
                                    "NA",
                                    "NA",
                                    row["peak_count"],
                                    row["median_peak_width_bp"],
                                    "NA",
                                    "NA",
                                    "NA",
                                    "NA",
                                    "NA",
                                    row["composite_score"],
                                ]
                            )
                    else:
                        for row in score_rows:
                            writer.writerow(
                                [
                                    row["group"],
                                    row["candidate"],
                                    row["qvalue"],
                                    row["shift"],
                                    row["extsize"],
                                    row["frip"],
                                    row["reads_in_peaks"],
                                    row["mapped_reads"],
                                    row["replicate_jaccard"],
                                    row["peak_count"],
                                    row["median_peak_width_bp"],
                                    row["blacklist_fraction"],
                                    row["score_frip_z"],
                                    row["score_jaccard_z"],
                                    row["score_blacklist_z"],
                                    row["score_width_z"],
                                    row["composite_score"],
                                ]
                            )

                selected_tsv = os.path.join(optimization_dir, "selected_parameters.tsv")
                with open(selected_tsv, "w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["group", "candidate", "qvalue", "shift", "extsize", "composite_score", "frip", "replicate_jaccard", "blacklist_fraction", "peak_count"])
                    for group in sorted(selected):
                        row = selected[group]
                        writer.writerow(
                            [
                                group,
                                row["candidate"],
                                row["qvalue"],
                                row["shift"],
                                row["extsize"],
                                row["composite_score"],
                                row["frip"],
                                row["replicate_jaccard"],
                                row["blacklist_fraction"],
                                row["peak_count"],
                            ]
                        )

                if opt_mode != "none":
                    write_optimization_plots(score_rows, optimization_dir)

                # Concatenate the peak files
                shell(f"""cat {outputfolder}/*.MACS3.optimized.bed | sort -k1,1 -k2,2n -k3,3n | bedtools merge -i - > {outputfolder}/all_groups.merged_peaks.bed """)
                log_it(logfile, f"Total merged peaks: {subprocess.getoutput(f'wc -l {outputfolder}/all_groups.merged_peaks.bed').strip().split()[0]}")

        call_peaks(
            logfile,
            thetype=params.thetype,
            inputfolder1=params.inputfolder1,
            outputfolder=params.outputfolder,
            broad_mode=params.broad_mode,
            input_sample=params.input_sample,
        )
        shell(f"""echo "necessity file for callpeaks. can delete this." > {params.outputfolder}/extra_10.tmp""")
