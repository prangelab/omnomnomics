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
