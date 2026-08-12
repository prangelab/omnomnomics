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
