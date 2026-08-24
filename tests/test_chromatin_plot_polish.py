from pathlib import Path


ROOT = Path(__file__).parents[1]
TEMPLATES_DIR = ROOT / "src" / "omnomnomics" / "workflow" / "templates"
RULES_DIR = ROOT / "src" / "omnomnomics" / "workflow" / "rules"


def test_de_templates_collapse_equivalent_metadata_partitions():
    for template_name in ("de_core.R.tmpl", "de_core_chrom.R.tmpl"):
        source = (TEMPLATES_DIR / template_name).read_text()
        assert "groups <- match(values, unique(values))" in source


def test_de_distance_plots_use_sample_identifiers_only():
    for template_name in ("de_core.R.tmpl", "de_core_chrom.R.tmpl"):
        source = (TEMPLATES_DIR / template_name).read_text()
        distance_block = source.split("if (distance_plot_enabled)", 1)[1].split("if (variable_heatmap_enabled)", 1)[0]
        assert "sample_labels <- rownames(col_data)" in distance_block
        assert "label_columns" not in distance_block


def test_combined_summaries_keep_full_terms_outside_compact_plot_labels():
    for template_name in ("de_core.R.tmpl", "de_core_chrom.R.tmpl"):
        source = (TEMPLATES_DIR / template_name).read_text()
        assert "path_df$term_plot <- make.unique(compact_pathway_label(path_df$term)" in source
        assert "aes(x = term_plot, y = score, fill = source)" in source


def test_de_enrichment_plots_wrap_terms_and_zero_padj_volcanoes_are_finite():
    for template_name in ("de_core.R.tmpl", "de_core_chrom.R.tmpl"):
        source = (TEMPLATES_DIR / template_name).read_text()
        assert 'vdf$neglog10padj <- -log10(pmax(vdf$padj, 1e-300))' in source
        assert "wrap_enrichment_axis_labels <- function(labels, width = 36)" in source
        assert "dot_plot <- polish_enrichment_axis(dot_plot)" in source
        assert "bar_plot <- polish_enrichment_axis(bar_plot)" in source
        assert "ridge_plot <- polish_enrichment_axis(ridge_plot)" in source


def test_post_de_profile_and_motif_renderers_use_bounded_labels():
    source = (RULES_DIR / "16.analyze_peaks_de.smk").read_text()
    assert "fig.savefig(profile_path)" in source
    assert 'fig.savefig(profile_path, bbox_inches="tight")' not in source
    assert 'ax.set_title(f"{set_label}: top {method} motif enrichments")' in source
    assert 'regions_label_for_set(item), ame_tsv' in source
    assert source.count("def regions_label_for_set(item):") == 1
    assert source.index("def regions_label_for_set(item):") < source.index("def run_deeptools(")
    assert source.index("def regions_label_for_set(item):") < source.index("def run_meme_motifs(")
    assert 'upper_to_col.get("MOTIF_ALT_ID")' in source
    assert 'upper_to_col.get("ADJ_P-VALUE")' in source


def test_centered_signal_heatmaps_label_distance_in_kilobases():
    pre_de = (RULES_DIR / "14.analyze_peaks.smk").read_text()
    post_de = (RULES_DIR / "16.analyze_peaks_de.smk").read_text()

    assert "--xAxisLabel 'distance from center (kb)'" in pre_de
    assert "--xAxisLabel 'distance from center (kb)'" in post_de
    assert "--xAxisLabel 'distance from center (bp)'" not in pre_de
    assert "--xAxisLabel 'distance from center (bp)'" not in post_de


def test_domain_mode_uses_domain_labels_and_readable_chromatin_plots():
    de_rule = (RULES_DIR / "15.call_DE_chrom.smk").read_text()
    chrom_template = (TEMPLATES_DIR / "de_core_chrom.R.tmpl").read_text()
    post_de = (RULES_DIR / "16.analyze_peaks_de.smk").read_text()

    assert 'params.broad_mode == "domain"' in de_rule
    assert 'feature_label_plural = "domains"' in de_rule
    assert 'facet_wrap(~ peak_set, scales = "free", drop = TRUE)' in chrom_template
    assert 'head(label_df, min(max_volcano_labels, 10L))' in chrom_template
    assert '"de_significant": "DE domains"' in post_de


def test_library_complexity_prefers_pre_deduplication_bams():
    peak_qc = (RULES_DIR / "13.peak_qc.smk").read_text()

    assert 'raw_bam_inputfolder=f"{experiment_dir}/BAM"' in peak_qc
    assert 'complexity_bam_stage = "merged_pre_dedup"' in peak_qc
    assert 'complexity_bam_stage = "filtered_fallback"' in peak_qc
    assert 'calculate_library_complexity_metrics(complexity_bam_path)' in peak_qc


def test_peak_qc_reports_input_features_separately_from_merged_intervals():
    peak_qc = (RULES_DIR / "13.peak_qc.smk").read_text()

    assert peak_qc.count("peak_count = len(load_peak_records(peak_bed))") == 2
    assert '"merged_interval_count": merged_interval_count' in peak_qc
    assert peak_qc.count('"merged_interval_count",') >= 2


def test_peak_qc_multipage_pdf_uses_stable_page_bounds():
    peak_qc = (RULES_DIR / "13.peak_qc.smk").read_text()
    summary_block = peak_qc.split("def write_peak_qc_summary_pdf", 1)[1].split(
        "def write_sample_frip_outputs", 1
    )[0]

    assert 'pdf.savefig(fig, bbox_inches="tight")' not in summary_block
    assert summary_block.count("pdf.savefig(fig)") == 4


def test_chromatin_boxplots_use_accessibility_axis_label():
    chrom_template = (TEMPLATES_DIR / "de_core_chrom.R.tmpl").read_text()

    assert 'ggplot2::ylab("VST accessibility")' in chrom_template
    assert 'ggplot2::ylab("VST expression")' not in chrom_template


def test_motif_peak_caps_use_ranked_de_regions():
    post_de = (RULES_DIR / "16.analyze_peaks_de.smk").read_text()

    assert 'de_rank_by_peak[peak_id] = (padj, -abs(lfc), chrom, start, end)' in post_de
    assert "rows_promoter = rank_de_rows(rows_promoter)" in post_de
    assert '"highest-ranked"' in post_de
    assert "item[\"set_type\"] != \"unique_from_prede\"" in post_de
