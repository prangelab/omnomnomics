from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RULE = ROOT / "src" / "omnomnomics" / "workflow" / "rules" / "5.touchup_bam.smk"


def test_touchup_uses_layout_specific_samtools_pipeline():
    source = RULE.read_text()

    assert "paired=config['PAIRED']" in source
    assert "stage_count = 5 if paired else 3" in source
    assert "stage_threads = max(1, min(8, samcores // stage_count))" in source
    assert "if paired:" in source
    assert "samtools collate -O" in source
    assert "samtools fixmate -mu" in source
    assert "samtools sort -u -@ {stage_threads} {quote(local_input)}" in source


def test_single_end_touchup_does_not_run_fixmate():
    source = RULE.read_text()
    single_end_branch = source.split("else:\n                    command = f\"\"\"", 1)[1]

    assert "samtools sort" in single_end_branch
    assert "samtools markdup" in single_end_branch
    assert "samtools fixmate" not in single_end_branch.split('"""', 1)[0]
