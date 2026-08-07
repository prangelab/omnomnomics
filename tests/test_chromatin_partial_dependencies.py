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
