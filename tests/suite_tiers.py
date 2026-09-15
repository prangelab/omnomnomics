from __future__ import annotations

from pathlib import Path


BIOINFORMATICS_TEST_MODULES = frozenset(
    {
        "test_chip_analysis.py",
        "test_chip_analysis_config.py",
        "test_chip_input_benchmark.py",
        "test_chip_inputs.py",
        "test_chip_review_regressions.py",
    }
)


def test_paths(tests_dir: Path, tier: str) -> list[Path]:
    paths = sorted(tests_dir.glob("test_*.py"))
    if tier == "all":
        return paths
    if tier == "bioinformatics":
        return [path for path in paths if path.name in BIOINFORMATICS_TEST_MODULES]
    if tier == "lightweight":
        return [path for path in paths if path.name not in BIOINFORMATICS_TEST_MODULES]
    raise ValueError(f"Unknown test tier: {tier}")
