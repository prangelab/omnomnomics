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
