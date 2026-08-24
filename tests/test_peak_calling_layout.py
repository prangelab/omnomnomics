from pathlib import Path


RULE = (
    Path(__file__).parents[1]
    / "src"
    / "omnomnomics"
    / "workflow"
    / "rules"
    / "10.call_peaks.smk"
)


def test_atac_macs3_format_tracks_library_layout():
    source = RULE.read_text()
    helper = source.split("def atac_macs3_format_args", 1)[1].split(
        "def idr_bam_suffix", 1
    )[0]

    assert 'return ["-f", "BAMPE"]' in helper
    assert (
        'return ["-f", "BAM", "--nomodel", "--shift", "-100", '
        '"--extsize", "200"]'
    ) in helper


def test_all_atac_peak_calling_paths_use_layout_helper():
    source = RULE.read_text()

    assert source.count("cmd.extend(atac_macs3_format_args())") == 2
    assert 'cmd.extend(["-f", "BAMPE"])' not in source
