from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RULE = ROOT / "src" / "omnomnomics" / "workflow" / "rules" / "11.count_reads.smk"
SITE = ROOT / "src" / "omnomnomics" / "workflow" / "config" / "site.yaml"
WORKFLOW = ROOT / "src" / "omnomnomics" / "workflow" / "config" / "workflow.yaml"


def test_featurecounts_stages_bams_and_temporary_files():
    source = RULE.read_text()

    assert "def stage_bams_for_counting" in source
    assert "shutil.copy2(bam_file, local_bam)" in source
    assert "--tmpDir {quote(local_workdir)}" in source
    assert "shutil.rmtree(local_workdir, ignore_errors=True)" in source


def test_count_reads_supports_a_site_specific_heavy_io_constraint():
    source = RULE.read_text()

    assert 'master_config.get("heavy_io_constraint")' in source
    assert "heavy_io_constraint:" in SITE.read_text()


def test_rna_count_runtime_is_four_hours():
    lines = WORKFLOW.read_text().splitlines()
    start = lines.index("rule_runtime:")
    runtimes = [line.strip()[2:] for line in lines[start + 1 : start + 18]]

    assert runtimes[10] == "240"
