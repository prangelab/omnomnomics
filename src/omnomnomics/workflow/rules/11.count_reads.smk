# Rule 11: Count Reads

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import csv
import hashlib
import os
import shlex
import subprocess
import shutil


def count_reads_input(_wildcards):
    input_files = []
    if config["THETYPE"] == "RNA":
        input_files.extend(
            f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][0]}/{sample}.sorted.dups_marked.filtered.bam"
            for sample in samples2
        )
    elif config["THETYPE"] in {"ATAC", "CHIP"}:
        if master_config.get("analyzepeaks_rule_num") in themode:
            input_files.append(
                f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][1]}/extra_{master_config['analyzepeaks_rule_num']}.tmp"
            )
        elif master_config.get("peakqc_rule_num") in themode:
            input_files.append(
                f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][1]}/extra_{master_config['peakqc_rule_num']}.tmp"
            )
        elif master_config["callpeaks_rule_num"] in themode:
            input_files.append(
                f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][1]}/extra_{master_config['callpeaks_rule_num']}.tmp"
            )
        bam_suffix = ".sorted.dups_marked.filtered.bam" if config["THETYPE"] == "ATAC" else ".filtered.bam"
        input_files.extend(
            f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][0]}/{sample}{bam_suffix}"
            for sample in samples2
        )
        input_files.append(
            f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][1]}/all_groups.merged_peaks.bed"
        )
    return input_files


rule count_reads:
    input:
        count_reads_input
    output:
        (
            f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.raw_read_quant.table.txt",
            f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.featureCounts.summary.txt",
            f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num']-1]}/extra_11.tmp",
        )
        if config["THETYPE"] == "RNA"
        else (
            f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.raw_read_quant.table.txt",
            f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num']-1]}/extra_11.tmp",
        )
    params:
        thetype=config["THETYPE"],
        broad_mode=str(config.get("BROAD_MODE", "off")).strip().lower(),
        genome=config["THEGENOME"],
        experiment_dir=config["EXPERIMENT_DIR"],
        paired=config["PAIRED"],
        bam_input_folder=f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][0]}",
        peak_input_folder=f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num'] - 1]}",
        gtf_file=os.path.join(config["GENOME_ASSEMBLY_DIR"], config["THEGENOME"], "annotation", "genes.gtf")
    threads:
        Threads_Per_Rule[str(master_config["countreads_rule_num"])]
    resources:
        mem_mb=Memory_Per_Rule[str(master_config["countreads_rule_num"])],
        partition=master_config["partition"],
        runtime=Runtime_Per_Rule[str(master_config["countreads_rule_num"])]
    run:
        tracking = begin_step_sample(master_config["countreads_rule_num"], "aggregate", "count_reads")
        log_once(logfile, "step11.header", "Counting Reads...", f"EXECUTING STEP {master_config['countreads_rule_num']}")
        log_once(logfile, "step11.inputfolder", f"Input folder: {params.bam_input_folder} and also {params.peak_input_folder} for ATAC data")
        log_once(logfile, "step11.outputfolder", f"Output folder: {params.outputfolder}")

        def quote(path):
            return shlex.quote(path)

        def samples_after_spp_drop():
            if params.thetype not in {"ATAC", "CHIP"}:
                return list(samples2)
            drop_file = os.path.join(
                experiment_dir,
                master_config["output_folders"][master_config["peakqc_rule_num"] - 1],
                "peak_qc",
                "spp_qc",
                "dropped_samples.tsv",
            )
            if not os.path.exists(drop_file):
                return list(samples2)
            dropped_ids = set()
            with open(drop_file, newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    sample_id = str(row.get("sample_id", "")).strip()
                    if sample_id:
                        dropped_ids.add(sample_id)
            if not dropped_ids:
                return list(samples2)
            filtered = [sample for sample in samples2 if sample_id_for_sample(sample) not in dropped_ids]
            log_it(
                logfile,
                "SPP drop list detected. Excluding samples from count matrix: "
                + ", ".join(sorted(dropped_ids)),
            )
            if not filtered:
                raise RuntimeError("All samples were dropped by SPP gate. Cannot continue with count matrix generation.")
            return filtered

        def write_tmp_file(outputfolder):
            shell(f"""echo "necessity file for count reads. can delete this." > {outputfolder}/extra_11.tmp""")

        def rna_output_path(outputfolder):
            return os.path.join(outputfolder, f"{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt")

        def rna_featurecounts_summary_path(outputfolder):
            return os.path.join(outputfolder, f"{os.path.basename(params.experiment_dir)}.featureCounts.summary.txt")

        def count_cache_dir():
            path = os.path.join(experiment_dir, "run_logs", "count_cache")
            os.makedirs(path, exist_ok=True)
            return path

        def file_state_token(path):
            stat = os.stat(path)
            return f"{path}\t{stat.st_size}\t{stat.st_mtime_ns}"

        def sha256_file(path):
            digest = hashlib.sha256()
            with open(path, "rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            return digest.hexdigest()

        def count_cache_key(feature_bed, bam_files, selected_samples, feature_label, paired, featurecounts_version):
            digest = hashlib.sha256()
            digest.update(f"feature_label={feature_label}\npaired={paired}\n".encode("utf-8"))
            digest.update(f"featurecounts={featurecounts_version}\n".encode("utf-8"))
            digest.update(f"feature_bed_sha256={sha256_file(feature_bed)}\n".encode("utf-8"))
            for sample, bam in zip(selected_samples, bam_files):
                digest.update(f"{sample}\t{sample_id_for_sample(sample)}\t{file_state_token(bam)}\n".encode("utf-8"))
            return digest.hexdigest()

        def prepare_featurecount_inputs(feature_bed, output_folder, feature_label):
            prep_dir = os.path.join(output_folder, "count_features")
            os.makedirs(prep_dir, exist_ok=True)
            cleaned_bed = os.path.join(prep_dir, f"{feature_label.lower()}.cleaned.bed")
            saf_file = os.path.join(prep_dir, f"{feature_label.lower()}.cleaned.saf")

            rows = set()
            with open(feature_bed, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
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
                    rows.add((fields[0], start, end))

            cleaned_rows = sorted(rows, key=lambda item: (item[0], item[1], item[2]))
            with open(cleaned_bed, "w", encoding="utf-8") as bed_handle, open(saf_file, "w", encoding="utf-8") as saf_handle:
                saf_handle.write("GeneID\tChr\tStart\tEnd\tStrand\n")
                for chrom, start, end in cleaned_rows:
                    feature_id = f"{chrom}_{start}_{end}"
                    bed_handle.write(f"{chrom}\t{start}\t{end}\n")
                    saf_handle.write(f"{feature_id}\t{chrom}\t{start + 1}\t{end}\t.\n")

            if not cleaned_rows:
                raise RuntimeError(f"No valid intervals remained after cleaning feature BED: {feature_bed}")
            return cleaned_bed, saf_file, len(cleaned_rows)

        def count_reads_over_features_with_featurecounts(input_folder, peak_bed, output_folder, selected_samples, bam_suffix, feature_label):
            featurecounts_version = subprocess.check_output(["featureCounts", "-v"], stderr=subprocess.STDOUT).decode("utf-8").strip()
            log_once(logfile, "step11.featurecounts_version", "\n" + featurecounts_version, "FEATURECOUNTS VERSION")

            bam_files = [os.path.join(input_folder, f"{sample}{bam_suffix}") for sample in selected_samples]
            if not bam_files:
                raise RuntimeError("No BAM files selected for count matrix generation.")
            for bam in bam_files:
                if not os.path.exists(bam):
                    raise FileNotFoundError(f"Required BAM for count matrix generation not found: {bam}")

            cleaned_bed, saf_file, feature_count = prepare_featurecount_inputs(peak_bed, output_folder, feature_label)
            final_output = rna_output_path(output_folder)
            summary_output = os.path.join(output_folder, f"{os.path.basename(params.experiment_dir)}.featureCounts.summary.txt")
            cache_key = count_cache_key(cleaned_bed, bam_files, selected_samples, feature_label, params.paired, featurecounts_version)
            cache_table = os.path.join(count_cache_dir(), f"{cache_key}.raw_read_quant.table.txt")
            cache_summary = os.path.join(count_cache_dir(), f"{cache_key}.featureCounts.summary.txt")
            cache_manifest = os.path.join(count_cache_dir(), f"{cache_key}.manifest.tsv")

            if os.path.exists(cache_table):
                shutil.copy2(cache_table, final_output)
                if os.path.exists(cache_summary):
                    shutil.copy2(cache_summary, summary_output)
                log_it(logfile, f"Reused cached {feature_label} count matrix: {cache_table}")
                return

            log_it(logfile, f"Cleaned {feature_label} BED for counting: {cleaned_bed} ({feature_count} intervals)")
            log_it(logfile, f"SAF annotation for featureCounts: {saf_file}")
            featurecounts_output = os.path.join(output_folder, f"{os.path.basename(params.experiment_dir)}.featureCounts.tmp.txt")
            paired_flags = "-p --countReadPairs" if params.paired else ""
            featurecounts_command = (
                f"featureCounts -T {threads} -F SAF -a {quote(saf_file)} -o {quote(featurecounts_output)} "
                f"{paired_flags} {' '.join(quote(path) for path in bam_files)}"
            ).strip()
            log_it(logfile, featurecounts_command, "FEATURECOUNTS COMMAND")
            shell(featurecounts_command)

            with open(featurecounts_output, newline="") as source, open(final_output, "w", newline="") as destination:
                reader = csv.reader((line for line in source if not line.startswith("#")), delimiter="\t")
                header = next(reader)
                writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
                writer.writerow([feature_label, *[sample_id_for_sample(sample) for sample in selected_samples]])
                for row in reader:
                    writer.writerow([row[0], *row[6:]])

            featurecounts_summary = f"{featurecounts_output}.summary"
            if os.path.exists(featurecounts_summary):
                os.replace(featurecounts_summary, summary_output)
            os.remove(featurecounts_output)
            shutil.copy2(final_output, cache_table)
            if os.path.exists(summary_output):
                shutil.copy2(summary_output, cache_summary)
            with open(cache_manifest, "w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["key", "value"])
                writer.writerow(["feature_label", feature_label])
                writer.writerow(["feature_bed", peak_bed])
                writer.writerow(["cleaned_bed", cleaned_bed])
                writer.writerow(["feature_count", feature_count])
                writer.writerow(["featurecounts_version", featurecounts_version])
                writer.writerow(["sample_count", len(selected_samples)])
                for sample, bam in zip(selected_samples, bam_files):
                    writer.writerow([f"sample:{sample}", bam])
            log_it(logfile, f"Cached {feature_label} count matrix: {cache_table}")

        def count_reads_rna(input_folder, output_folder, gtf_file, paired):
            log_once(logfile, "step11.rna_mode", "Counting RNA reads from BAMs with featureCounts...")
            sanity_check_dir(logfile, input_folder, master_config["input_file_types"][master_config["countreads_rule_num"] - 1][0], "step11.rna_sanity")

            if not os.path.isfile(gtf_file):
                raise FileNotFoundError(f"Genome annotation GTF not found: {gtf_file}")

            featurecounts_version = subprocess.check_output(["featureCounts", "-v"], stderr=subprocess.STDOUT)
            log_once(logfile, "step11.featurecounts_version", "\n" + featurecounts_version.decode("utf-8"), "FEATURECOUNTS VERSION")

            bam_files = [os.path.join(input_folder, f"{sample}.sorted.dups_marked.filtered.bam") for sample in samples2]
            featurecounts_output = os.path.join(output_folder, f"{os.path.basename(params.experiment_dir)}.featureCounts.tmp.txt")
            paired_flags = "-p --countReadPairs" if paired else ""
            featurecounts_command = (
                f"featureCounts -T {threads} -a {quote(gtf_file)} -o {quote(featurecounts_output)} "
                f"-t exon -g gene_id {paired_flags} {' '.join(quote(path) for path in bam_files)}"
            ).strip()
            log_it(logfile, featurecounts_command, "FEATURECOUNTS COMMAND")
            shell(featurecounts_command)

            final_output = rna_output_path(output_folder)
            with open(featurecounts_output, newline="") as source, open(final_output, "w", newline="") as destination:
                reader = csv.reader((line for line in source if not line.startswith("#")), delimiter="\t")
                header = next(reader)
                sample_headers = [
                    sample_id_for_sample(
                        os.path.basename(path).replace(".sorted.dups_marked.filtered.bam", "")
                    )
                    for path in header[6:]
                ]
                writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
                writer.writerow(["Geneid", *sample_headers])
                for row in reader:
                    writer.writerow([row[0], *row[6:]])

            featurecounts_summary = f"{featurecounts_output}.summary"
            if os.path.exists(featurecounts_summary):
                os.replace(featurecounts_summary, rna_featurecounts_summary_path(output_folder))
            os.remove(featurecounts_output)

        def count_reads_atac(input_folder, peak_folder, output_folder):
            log_once(logfile, "step11.atac_mode", "Counting ATAC reads over peaks with featureCounts...")
            sanity_check_dir(logfile, input_folder, master_config["input_file_types"][master_config["countreads_rule_num"] - 1][0], "step11.atac_bam_sanity")
            sanity_check_dir(logfile, peak_folder, master_config["input_file_types"][master_config["countreads_rule_num"] - 1][1], "step11.atac_peak_sanity")

            filtered_peak_bed = os.path.join(
                experiment_dir,
                master_config["output_folders"][master_config["peakqc_rule_num"] - 1],
                "peak_qc",
                "filtered_peaks",
                "all_groups.merged_peaks.bed",
            )
            peak_bed = filtered_peak_bed if os.path.exists(filtered_peak_bed) else os.path.join(peak_folder, "all_groups.merged_peaks.bed")
            log_it(logfile, f"ATAC peak BED used for counting: {peak_bed}")
            selected_samples = samples_after_spp_drop()
            count_reads_over_features_with_featurecounts(
                input_folder=input_folder,
                peak_bed=peak_bed,
                output_folder=output_folder,
                selected_samples=selected_samples,
                bam_suffix=".sorted.dups_marked.filtered.bam",
                feature_label="Peak",
            )

        def count_reads_chip(input_folder, peak_folder, output_folder, broad_mode):
            if broad_mode == "genebody":
                log_once(logfile, "step11.chip_mode", "Counting ChIP reads over gene-body features with featureCounts...")
            elif broad_mode == "diffuse":
                log_once(logfile, "step11.chip_mode", "Counting ChIP reads over fixed genomic bins with featureCounts...")
            else:
                log_once(logfile, "step11.chip_mode", "Counting ChIP reads over peaks with featureCounts...")
            sanity_check_dir(logfile, input_folder, master_config["input_file_types"][master_config["countreads_rule_num"] - 1][0], "step11.chip_bam_sanity")
            sanity_check_dir(logfile, peak_folder, master_config["input_file_types"][master_config["countreads_rule_num"] - 1][1], "step11.chip_peak_sanity")

            peak_bed = os.path.join(peak_folder, "all_groups.merged_peaks.bed")
            if broad_mode not in {"genebody", "diffuse"}:
                filtered_peak_bed = os.path.join(
                    experiment_dir,
                    master_config["output_folders"][master_config["peakqc_rule_num"] - 1],
                    "peak_qc",
                    "filtered_peaks",
                    "all_groups.merged_peaks.bed",
                )
                if os.path.exists(filtered_peak_bed):
                    peak_bed = filtered_peak_bed
            if broad_mode == "genebody":
                log_it(logfile, f"ChIP gene-body BED used for counting: {peak_bed}")
            elif broad_mode == "diffuse":
                log_it(logfile, f"ChIP diffuse bin BED used for counting: {peak_bed}")
            else:
                log_it(logfile, f"ChIP peak BED used for counting: {peak_bed}")
            selected_samples = samples_after_spp_drop()
            if broad_mode == "genebody":
                feature_label = "Feature"
            elif broad_mode == "diffuse":
                feature_label = "Bin"
            else:
                feature_label = "Peak"
            count_reads_over_features_with_featurecounts(
                input_folder=input_folder,
                peak_bed=peak_bed,
                output_folder=output_folder,
                selected_samples=selected_samples,
                bam_suffix=".filtered.bam",
                feature_label=feature_label,
            )

        try:
            if params.thetype == "RNA":
                count_reads_rna(params.bam_input_folder, params.outputfolder, params.gtf_file, params.paired)
            elif params.thetype == "ATAC":
                count_reads_atac(params.bam_input_folder, params.peak_input_folder, params.outputfolder)
            elif params.thetype == "CHIP":
                count_reads_chip(params.bam_input_folder, params.peak_input_folder, params.outputfolder, params.broad_mode)
            else:
                log_once(logfile, "step11.chip_note", "For ChIP experiments, first determine optimal peak caller settings and quantify peaks with your chosen downstream workflow before continuing.")

            write_tmp_file(params.outputfolder)
            finish_step_sample(master_config["countreads_rule_num"], "aggregate", "count_reads", tracking["start_time"], "OK")
        except Exception:
            finish_step_sample(master_config["countreads_rule_num"], "aggregate", "count_reads", tracking["start_time"], "FAIL")
            raise
