# Rule 7: Get alignment QC summaries

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import csv
import math
import os
import statistics
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PAIRED = 0x1
PROPER_PAIR = 0x2
UNMAP = 0x4
MATE_UNMAP = 0x8
READ1 = 0x40
DUP = 0x400
SECONDARY = 0x100
SUPPLEMENTARY = 0x800
PRIMARY_FILTER = SECONDARY | SUPPLEMENTARY
PRIMARY_MAPPED_FILTER = UNMAP | SECONDARY | SUPPLEMENTARY
DISCORDANT_FILTER = UNMAP | MATE_UNMAP | PROPER_PAIR | SECONDARY | SUPPLEMENTARY


rule bam_stats:
    input:
        pre_filter_bam=lambda wildcards: (
            f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}/{wildcards.sample}.bam"
            if os.path.exists(f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}/{wildcards.sample}.bam")
            else []
        ),
        filtered_BAM=lambda wildcards: (
            f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}/{wildcards.sample}.sorted.dups_marked.filtered.bam"
            if config['THETYPE'] != "CHIP" and os.path.exists(f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}/{wildcards.sample}.sorted.dups_marked.filtered.bam")
            else (
                f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}/{wildcards.sample}.filtered.bam"
                if config['THETYPE'] == "CHIP" and os.path.exists(f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}/{wildcards.sample}.filtered.bam")
                else []
            )
        )
    output:
        stats_tsv=f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.stats.txt" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}.filtered.bam.stats.txt",
        summary_pdf=f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.qc_summary.pdf" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}.filtered.bam.qc_summary.pdf",
        summary_svg=f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.qc_summary.svg" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}.filtered.bam.qc_summary.svg"
    params:
        inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}",
        prefilterfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}",
        duplicate_handling=config["DUPLICATE_HANDLING"],
        paired=config["PAIRED"]
    threads:
        Threads_Per_Rule['7']
    resources:
        mem_mb=Memory_Per_Rule['7'],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule['7']
    run:
        log_it(logfile, "Generating pre/post-filter alignment QC summaries...", f"EXECUTING STEP {master_config['stats_rule_num']}")
        log_it(logfile, f"Pre-filter BAM folder: {params.prefilterfolder}")
        log_it(logfile, f"Post-filter BAM folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")
        log_it(logfile, f"Duplicate handling: {params.duplicate_handling}")

        samtools_version = subprocess.check_output(["samtools", "--version"], stderr=subprocess.STDOUT).decode("utf-8").splitlines()[:2]
        log_it(logfile, "\n" + "\n".join(samtools_version) + "\n", "SAMTOOLS VERSION")

        pre_filter_bam = input.pre_filter_bam if isinstance(input.pre_filter_bam, str) else None
        filtered_bam = input.filtered_BAM if isinstance(input.filtered_BAM, str) else None

        if pre_filter_bam:
            sanity_check_dir(logfile, params.prefilterfolder, ".bam")
        if filtered_bam:
            sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['stats_rule_num']-1])
        if not pre_filter_bam and not filtered_bam:
            raise FileNotFoundError(
                f"No pre-filter BAM in {params.prefilterfolder} and no post-filter BAM in {params.inputfolder} for sample {wildcards.sample}"
            )

        def samtools_count(bam_path, include_flags=None, exclude_flags=None):
            command = ["samtools", "view", "-c"]
            if include_flags is not None:
                command.extend(["-f", str(include_flags)])
            if exclude_flags is not None:
                command.extend(["-F", str(exclude_flags)])
            command.append(bam_path)
            return int(subprocess.check_output(command, stderr=subprocess.STDOUT).decode("utf-8").strip())

        def collect_alignment_metrics(bam_path, paired_end):
            metrics = {
                "total_primary_reads": samtools_count(bam_path, exclude_flags=PRIMARY_FILTER),
                "mapped_primary_reads": samtools_count(bam_path, exclude_flags=PRIMARY_MAPPED_FILTER),
                "duplicate_primary_reads": samtools_count(bam_path, include_flags=DUP, exclude_flags=PRIMARY_MAPPED_FILTER),
            }

            if paired_end:
                metrics.update({
                    "total_primary_pairs": samtools_count(bam_path, include_flags=PAIRED | READ1, exclude_flags=PRIMARY_FILTER),
                    "mapped_primary_pairs": samtools_count(bam_path, include_flags=PAIRED | READ1, exclude_flags=PRIMARY_MAPPED_FILTER),
                    "properly_paired_templates": samtools_count(bam_path, include_flags=PAIRED | PROPER_PAIR | READ1, exclude_flags=PRIMARY_MAPPED_FILTER),
                    "discordant_templates": samtools_count(bam_path, include_flags=PAIRED | READ1, exclude_flags=DISCORDANT_FILTER),
                    "singleton_mapped_reads": samtools_count(bam_path, include_flags=PAIRED | MATE_UNMAP, exclude_flags=PRIMARY_MAPPED_FILTER),
                    "duplicate_primary_pairs": samtools_count(bam_path, include_flags=PAIRED | READ1 | DUP, exclude_flags=PRIMARY_MAPPED_FILTER),
                })
            return metrics

        def percent(numerator, denominator):
            if denominator == 0:
                return 0.0
            return round((numerator / denominator) * 100, 4)

        pre_metrics = collect_alignment_metrics(pre_filter_bam, params.paired) if pre_filter_bam else {}
        post_metrics = collect_alignment_metrics(filtered_bam, params.paired) if filtered_bam else {}

        derived_metrics = {}
        if pre_metrics and post_metrics:
            derived_metrics = {
                "mapped_primary_reads_removed": pre_metrics["mapped_primary_reads"] - post_metrics["mapped_primary_reads"],
                "mapped_primary_reads_retained_pct": percent(post_metrics["mapped_primary_reads"], pre_metrics["mapped_primary_reads"]),
                "duplicate_primary_reads_removed": pre_metrics["duplicate_primary_reads"] - post_metrics["duplicate_primary_reads"],
            }

            if params.paired:
                derived_metrics.update({
                    "properly_paired_templates_removed": pre_metrics["properly_paired_templates"] - post_metrics["properly_paired_templates"],
                    "properly_paired_templates_retained_pct": percent(post_metrics["properly_paired_templates"], pre_metrics["properly_paired_templates"]),
                    "discordant_templates_removed": pre_metrics["discordant_templates"] - post_metrics["discordant_templates"],
                    "duplicate_primary_pairs_removed": pre_metrics["duplicate_primary_pairs"] - post_metrics["duplicate_primary_pairs"],
                })

        outfile = output.stats_tsv
        log_it(logfile, f"Generating alignment QC summary for {wildcards.sample}")

        with open(outfile, "w", newline="") as handle:
            handle.write(f"# sample\t{wildcards.sample}\n")
            handle.write(f"# paired_end\t{str(bool(params.paired)).lower()}\n")
            handle.write(f"# duplicate_handling\t{params.duplicate_handling}\n")
            handle.write(f"# pre_filter_bam_present\t{str(bool(pre_filter_bam)).lower()}\n")
            handle.write(f"# post_filter_bam_present\t{str(bool(filtered_bam)).lower()}\n")
            if params.paired:
                handle.write("# pair-level metrics count primary read1 records to avoid double counting templates\n")
                handle.write("# singleton_mapped_reads remains a read-level metric because a singleton has no complete pair\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", "stage", "metric", "unit", "value"])

            all_metric_names = set(pre_metrics) | set(post_metrics) | set(derived_metrics)
            if not all_metric_names:
                all_metric_names = {
                    "total_primary_reads",
                    "mapped_primary_reads",
                    "duplicate_primary_reads",
                }
                if params.paired:
                    all_metric_names.update({
                        "total_primary_pairs",
                        "mapped_primary_pairs",
                        "properly_paired_templates",
                        "discordant_templates",
                        "singleton_mapped_reads",
                        "duplicate_primary_pairs",
                        "properly_paired_templates_removed",
                        "properly_paired_templates_retained_pct",
                        "discordant_templates_removed",
                    })
                all_metric_names.update({
                    "mapped_primary_reads_removed",
                    "mapped_primary_reads_retained_pct",
                    "duplicate_primary_reads_removed",
                })

            for metric in sorted(all_metric_names):
                if metric in pre_metrics or metric in post_metrics:
                    unit = "pairs" if metric.endswith("pairs") or metric.endswith("templates") else "reads"
                else:
                    unit = "percent" if metric.endswith("_pct") else "pairs" if metric.endswith("templates_removed") or metric.endswith("pairs_removed") else "reads"
                writer.writerow([wildcards.sample, "pre_filter", metric, unit, pre_metrics.get(metric, "NA")])
                writer.writerow([wildcards.sample, "post_filter", metric, unit, post_metrics.get(metric, "NA")])
                writer.writerow([wildcards.sample, "derived", metric, unit, derived_metrics.get(metric, "NA")])

        def load_metrics_from_tsv(tsv_path):
            metrics = {"pre_filter": {}, "post_filter": {}, "derived": {}}
            with open(tsv_path, newline="") as handle:
                reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
                for row in reader:
                    value = row["value"]
                    if value == "NA":
                        parsed_value = None
                    else:
                        try:
                            parsed_value = float(value)
                        except ValueError:
                            parsed_value = value
                    metrics[row["stage"]][row["metric"]] = parsed_value
            return metrics

        def plot_value_label(ax, x_pos, value):
            if value is None:
                ax.text(x_pos, 0.5, "NA", ha="center", va="bottom", fontsize=10, color="#666666")
            else:
                ax.text(x_pos, value, f"{value:,.0f}" if value >= 10 else f"{value:.2f}", ha="center", va="bottom", fontsize=9)

        def write_qc_summary_plots(tsv_path, pdf_path, svg_path, sample_name, paired_end):
            metrics = load_metrics_from_tsv(tsv_path)
            pre = metrics["pre_filter"]
            post = metrics["post_filter"]
            derived = metrics["derived"]

            fig, axes = plt.subplots(2, 2, figsize=(12, 8))
            fig.patch.set_facecolor("white")
            fig.suptitle(f"{sample_name} alignment QC summary", fontsize=16, fontweight="bold")

            colors = {"pre": "#c9d6df", "post": "#4f6d7a"}

            ax = axes[0, 0]
            read_labels = ["Mapped reads", "Duplicate reads"]
            pre_read_values = [pre.get("mapped_primary_reads"), pre.get("duplicate_primary_reads")]
            post_read_values = [post.get("mapped_primary_reads"), post.get("duplicate_primary_reads")]
            x_positions = range(len(read_labels))
            width = 0.35
            pre_heights = [value or 0 for value in pre_read_values]
            post_heights = [value or 0 for value in post_read_values]
            ax.bar([x - width / 2 for x in x_positions], pre_heights, width=width, color=colors["pre"], label="Pre-filter")
            ax.bar([x + width / 2 for x in x_positions], post_heights, width=width, color=colors["post"], label="Post-filter")
            for idx, value in enumerate(pre_read_values):
                plot_value_label(ax, idx - width / 2, value)
            for idx, value in enumerate(post_read_values):
                plot_value_label(ax, idx + width / 2, value)
            ax.set_xticks(list(x_positions))
            ax.set_xticklabels(read_labels, rotation=10)
            ax.set_ylabel("Reads")
            ax.set_title("Read-level counts")
            ax.legend(frameon=False)

            ax = axes[0, 1]
            retained_metrics = [
                ("Mapped reads retained", derived.get("mapped_primary_reads_retained_pct")),
            ]
            if paired_end:
                retained_metrics.append(("Proper pairs retained", derived.get("properly_paired_templates_retained_pct")))
            retained_labels = [label for label, _ in retained_metrics]
            retained_values = [value for _, value in retained_metrics]
            ax.bar(range(len(retained_labels)), [value or 0 for value in retained_values], color="#2a9d8f")
            for idx, value in enumerate(retained_values):
                if value is None:
                    ax.text(idx, 2, "NA", ha="center", va="bottom", fontsize=10, color="#666666")
                else:
                    ax.text(idx, value, f"{value:.1f}%", ha="center", va="bottom", fontsize=10)
            ax.set_ylim(0, 105)
            ax.set_xticks(range(len(retained_labels)))
            ax.set_xticklabels(retained_labels, rotation=10)
            ax.set_ylabel("Percent")
            ax.set_title("Retention after filtering")

            ax = axes[1, 0]
            removed_metrics = [
                ("Mapped reads removed", derived.get("mapped_primary_reads_removed")),
                ("Duplicate reads removed", derived.get("duplicate_primary_reads_removed")),
            ]
            if paired_end:
                removed_metrics.extend([
                    ("Proper pairs removed", derived.get("properly_paired_templates_removed")),
                    ("Discordant pairs removed", derived.get("discordant_templates_removed")),
                ])
            removed_labels = [label for label, _ in removed_metrics]
            removed_values = [value for _, value in removed_metrics]
            ax.bar(range(len(removed_labels)), [value or 0 for value in removed_values], color="#e76f51")
            for idx, value in enumerate(removed_values):
                plot_value_label(ax, idx, value)
            ax.set_xticks(range(len(removed_labels)))
            ax.set_xticklabels(removed_labels, rotation=15)
            ax.set_ylabel("Reads / pairs")
            ax.set_title("Counts removed by filtering")

            ax = axes[1, 1]
            ax.axis("off")
            summary_lines = [
                f"Duplicate handling: {params.duplicate_handling}",
                f"Pre-filter BAM present: {'yes' if pre_filter_bam else 'no'}",
                f"Post-filter BAM present: {'yes' if filtered_bam else 'no'}",
            ]
            if paired_end:
                summary_lines.append("Pair-level metrics count primary read1 records only")
            summary_lines.append("Singleton counts remain read-level")
            ax.text(
                0.02,
                0.98,
                "\n".join(summary_lines),
                va="top",
                ha="left",
                fontsize=11,
                bbox={"facecolor": "#f7f7f7", "edgecolor": "#d0d0d0", "boxstyle": "round,pad=0.5"},
            )

            for axis in axes.flat[:3]:
                axis.spines["top"].set_visible(False)
                axis.spines["right"].set_visible(False)

            fig.tight_layout(rect=[0, 0.02, 1, 0.95])
            fig.savefig(pdf_path)
            fig.savefig(svg_path)
            plt.close(fig)

        write_qc_summary_plots(outfile, output.summary_pdf, output.summary_svg, wildcards.sample, params.paired)
