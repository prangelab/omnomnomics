from pathlib import Path


RULES_DIR = Path(__file__).parents[1] / "src" / "omnomnomics" / "workflow" / "rules"


def test_peak_qc_enters_through_durable_merged_peaks():
    source = (RULES_DIR / "13.peak_qc.smk").read_text()
    input_block = source.split("def peak_qc_input", 1)[1].split("rule peak_qc", 1)[0]

    assert "all_groups.merged_peaks.bed" in input_block
    assert "extra_{master_config['callpeaks_rule_num']}" not in input_block


def test_pre_de_analysis_enters_through_annotation_not_marker():
    source = (RULES_DIR / "14.analyze_peaks.smk").read_text()
    input_block = source.split("def analyze_peaks_input", 1)[1].split("rule analyze_peaks", 1)[0]

    assert "all_groups.merged_peaks.annotated.bed" in input_block
    assert "extra_{master_config['peakqc_rule_num']}" not in input_block


def test_non_peak_chromatin_de_waits_for_pre_de_completion_marker():
    source = (RULES_DIR / "15.call_DE_chrom.smk").read_text()
    dependency_block = source.split(
        "def _chrom_de_peak_metadata_dependency", 1
    )[1].split("chrom_de_peak_metadata_file", 1)[0]

    assert '{"genebody", "diffuse"}' in dependency_block
    assert "extra_{master_config['analyzepeaks_rule_num']}.tmp" in dependency_block
    assert "return _chrom_de_peak_metadata_file()" not in dependency_block


def test_genebody_feature_builder_uses_transcript_aware_gene_spans():
    source = (RULES_DIR / "10.call_peaks.smk").read_text()
    builder = source.split("def build_chip_genebody_feature_sets", 1)[1].split(
        "def pooled_overlap_support_bed", 1
    )[0]

    assert "build_gtf_annotation_sources" in builder
    assert 'annotation_sources["genes"]' in builder
    assert 'fields[2] != "gene"' not in builder
    assert "if feature_count == 0:" in builder


def test_genebody_features_preserve_source_gene_labels():
    peak_calling = (RULES_DIR / "10.call_peaks.smk").read_text()
    pre_de = (RULES_DIR / "14.analyze_peaks.smk").read_text()

    builder = peak_calling.split("def build_chip_genebody_feature_sets", 1)[1].split(
        "def pooled_overlap_support_bed", 1
    )[0]
    metadata = pre_de.split("def load_genebody_source_labels", 1)[1].split(
        "def write_feature_summary", 1
    )[0]

    assert "feature_labels.setdefault((chrom, start, end), set()).add(label)" in builder
    assert "','.join(sorted(labels))" in builder
    assert '"gene_body"' in metadata
    assert 'row.get("assigned_genes", "")' in metadata


def test_genebody_post_de_omits_peak_style_promoter_distal_sets():
    source = (RULES_DIR / "16.analyze_peaks_de.smk").read_text()
    de_sets = source.split("def build_de_sets", 1)[1].split(
        "def copy_optional_unique_sets", 1
    )[0]

    assert 'if broad_mode != "genebody":' in de_sets
    assert '("de_significant_promoter_regions", rows_promoter_regions)' in de_sets
    assert '("top_de_ranked_promoter_regions", ranked_sets.get(' in de_sets


def test_chromatin_plot_labels_are_readable_for_long_feature_and_sample_names():
    pre_de = (RULES_DIR / "14.analyze_peaks.smk").read_text()
    post_de = (RULES_DIR / "16.analyze_peaks_de.smk").read_text()
    chrom_de = (RULES_DIR.parent / "templates" / "de_core_chrom.R.tmpl").read_text()

    assert "def render_scaled_genebody_profile" in pre_de
    assert 'ax.set_xticklabels(["-1 kb", "TSS", "TES", "+1 kb"])' in pre_de
    assert 're.sub(r"_(rep[0-9]+)$", r"\\n\\1", label)' in post_de
    assert "ylim = c(0, y_lim * 1.25)" in chrom_de
    assert "box.padding = 0.6" in chrom_de
    assert 'compact_volcano_label <- function(label)' in chrom_de
    assert '"-", end_mb, " Mb"' in chrom_de


def test_peak_backed_chip_modes_use_filtered_peak_tree_downstream():
    peak_qc = (RULES_DIR / "13.peak_qc.smk").read_text()
    pre_de = (RULES_DIR / "14.analyze_peaks.smk").read_text()
    counting = (RULES_DIR / "11.count_reads.smk").read_text()

    assert 'assay == "CHIP" and params.broad_mode not in {"genebody", "diffuse"}' in peak_qc
    assert 'chip_broad_mode not in {"genebody", "diffuse"}' in pre_de
    assert 'if broad_mode not in {"genebody", "diffuse"}' in counting
    assert '"peak_qc",\n                    "filtered_peaks"' in counting


def test_post_de_sets_apply_recorded_contrast_thresholds():
    source = (RULES_DIR / "16.analyze_peaks_de.smk").read_text()

    assert 'os.path.join(analysis_dir, "contrast_summary.tsv")' in source
    assert 'return float(row["alpha"]), float(row["lfc_threshold"])' in source
    assert "padj >= alpha or abs(lfc) < lfc_threshold" in source
    assert "sig_diff_" not in source


def test_reruns_clear_pipeline_owned_de_result_directories():
    templates_dir = RULES_DIR.parent / "templates"
    rna_template = (templates_dir / "de_core.R.tmpl").read_text()
    chrom_template = (templates_dir / "de_core_chrom.R.tmpl").read_text()
    post_de = (RULES_DIR / "16.analyze_peaks_de.smk").read_text()

    assert "unlink(contrast_dir, recursive = TRUE, force = TRUE)" in rna_template
    assert "unlink(contrast_dir, recursive = TRUE, force = TRUE)" in chrom_template
    assert "shutil.rmtree(analyze_de_root)" in post_de


def test_aggregate_completion_markers_survive_successful_runs():
    snakefile = (RULES_DIR.parent / "Snakefile.smk").read_text()
    cleanup = snakefile.split("# Clean up tmp files for workflow", 1)[1].split(
        "# Remove old BAM files with lane info", 1
    )[0]

    assert "extra_9.tmp" not in cleanup
    assert "extra_10.tmp" not in cleanup
    assert "extra_11.tmp" not in cleanup
    assert "analyzepeaks_rule_num" not in cleanup
    assert "analyzepeaksde_rule_num" not in cleanup
