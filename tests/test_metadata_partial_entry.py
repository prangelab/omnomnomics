from omnomnomics.cli import validate_metadata_sample_roots


def test_partial_entry_accepts_technical_and_merged_bam_names():
    rows = [
        {
            "filename": "sample_a_T01",
            "filename_key": "sample_a_T01",
            "sample_id": "sample_a",
        },
        {
            "filename": "sample_a_T02",
            "filename_key": "sample_a_T02",
            "sample_id": "sample_a",
        },
        {
            "filename": "sample_b_T01",
            "filename_key": "sample_b_T01",
            "sample_id": "sample_b",
        },
    ]

    validate_metadata_sample_roots(
        rows,
        ["sample_a_T01", "sample_a_T02", "sample_a", "sample_b_T01", "sample_b"],
    )


def test_partial_entry_rejects_unmatched_bam_name():
    rows = [
        {
            "filename": "sample_a_T01",
            "filename_key": "sample_a_T01",
            "sample_id": "sample_a",
        }
    ]

    try:
        validate_metadata_sample_roots(rows, ["sample_a", "unexpected_sample"])
    except Exception as exc:
        assert "unexpected_sample" in str(exc)
    else:
        raise AssertionError("Unmatched BAM names must fail metadata validation.")
