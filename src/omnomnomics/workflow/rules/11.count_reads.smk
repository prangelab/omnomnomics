# Rule 11: Count Reads

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import csv
import os
import shlex
import subprocess


def count_reads_input(_wildcards):
    input_files = []
    if config["THETYPE"] == "RNA":
        input_files.extend(
            f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][0]}/{sample}.sorted.dups_marked.filtered.bam"
            for sample in samples2
        )
    elif config["THETYPE"] == "ATAC":
        if master_config["callpeaks_rule_num"] in themode:
            input_files.append(
                f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][1]}/extra_{master_config['callpeaks_rule_num']}.tmp"
            )
        input_files.extend(
            f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num'] - 1][0]}/{sample}.sorted.dups_marked.filtered.bam"
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
        log_once(logfile, "step11.header", "Counting Reads...", f"EXECUTING STEP {master_config['countreads_rule_num']}")
        log_once(logfile, "step11.inputfolder", f"Input folder: {params.bam_input_folder} and also {params.peak_input_folder} for ATAC data")
        log_once(logfile, "step11.outputfolder", f"Output folder: {params.outputfolder}")

        def quote(path):
            return shlex.quote(path)

        def write_tmp_file(outputfolder):
            shell(f"""echo "necessity file for count reads. can delete this." > {outputfolder}/extra_11.tmp""")

        def rna_output_path(outputfolder):
            return os.path.join(outputfolder, f"{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt")

        def rna_featurecounts_summary_path(outputfolder):
            return os.path.join(outputfolder, f"{os.path.basename(params.experiment_dir)}.featureCounts.summary.txt")

        def count_reads_rna(input_folder, output_folder, gtf_file, paired):
            log_it(logfile, "Counting RNA reads from BAMs with featureCounts...")
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
                sample_headers = [os.path.basename(path).replace(".sorted.dups_marked.filtered.bam", "") for path in header[6:]]
                writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
                writer.writerow(["Geneid", *sample_headers])
                for row in reader:
                    writer.writerow([row[0], *row[6:]])

            featurecounts_summary = f"{featurecounts_output}.summary"
            if os.path.exists(featurecounts_summary):
                os.replace(featurecounts_summary, rna_featurecounts_summary_path(output_folder))
            os.remove(featurecounts_output)

        def count_reads_atac(input_folder, peak_folder, output_folder):
            log_it(logfile, "Counting ATAC reads from BAMs with bedtools multicov...")
            sanity_check_dir(logfile, input_folder, master_config["input_file_types"][master_config["countreads_rule_num"] - 1][0], "step11.atac_bam_sanity")
            sanity_check_dir(logfile, peak_folder, master_config["input_file_types"][master_config["countreads_rule_num"] - 1][1], "step11.atac_peak_sanity")

            bedtools_version = subprocess.check_output(["bedtools", "--version"], stderr=subprocess.STDOUT)
            log_once(logfile, "step11.bedtools_version", "\n" + bedtools_version.decode("utf-8"), "BEDTOOLS VERSION")

            peak_bed = os.path.join(peak_folder, "all_groups.merged_peaks.bed")
            bam_files = [os.path.join(input_folder, f"{sample}.sorted.dups_marked.filtered.bam") for sample in samples2]
            multicov_output = os.path.join(output_folder, f"{os.path.basename(params.experiment_dir)}.multicov.tmp.txt")
            multicov_command = (
                f"bedtools multicov -bed {quote(peak_bed)} -bams {' '.join(quote(path) for path in bam_files)} "
                f"> {quote(multicov_output)}"
            )
            log_it(logfile, multicov_command, "BEDTOOLS MULTICOV COMMAND")
            shell(multicov_command)

            final_output = rna_output_path(output_folder)
            with open(multicov_output, newline="") as source, open(final_output, "w", newline="") as destination:
                reader = csv.reader(source, delimiter="\t")
                writer = csv.writer(destination, delimiter="\t", lineterminator="\n")
                writer.writerow(["Peak", *samples2])
                for row in reader:
                    peak_name = f"{row[0]}_{row[1]}_{row[2]}"
                    writer.writerow([peak_name, *row[3:]])

            os.remove(multicov_output)

        if params.thetype == "RNA":
            count_reads_rna(params.bam_input_folder, params.outputfolder, params.gtf_file, params.paired)
        elif params.thetype == "ATAC":
            count_reads_atac(params.bam_input_folder, params.peak_input_folder, params.outputfolder)
        else:
            log_it(logfile, "For ChIP experiments, first determine optimal peak caller settings and quantify peaks with your chosen downstream workflow before continuing.")

        write_tmp_file(params.outputfolder)
