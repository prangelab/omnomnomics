from pathlib import Path

from omnomnomics.peak_annotation import (
    build_gtf_annotation_sources,
    promoter_intervals_by_gene,
)


def _read_bed(path: str) -> list[list[str]]:
    return [line.rstrip("\n").split("\t") for line in Path(path).read_text().splitlines()]


def test_transcript_only_gtf_derives_gene_span_promoters_and_tss(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr1\ttest\ttranscript\t101\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; gene_name \"Gene1\";\n"
        "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; gene_name \"Gene1\";\n"
        "chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; gene_name \"Gene1\";\n"
        "chr2\ttest\ttranscript\t1001\t1600\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\"; gene_name \"Gene2\";\n"
        "chr2\ttest\texon\t1001\t1100\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\"; gene_name \"Gene2\";\n"
    )

    paths = build_gtf_annotation_sources(gtf, tmp_path / "annotation")

    assert _read_bed(paths["genes"]) == [
        ["chr1", "100", "500", "Gene1|g1", "0", "+"],
        ["chr2", "1000", "1600", "Gene2|g2", "0", "-"],
    ]
    assert _read_bed(paths["tss"]) == [
        ["chr1", "100", "101", "Gene1|g1", "0", "+"],
        ["chr2", "1599", "1600", "Gene2|g2", "0", "-"],
    ]
    assert _read_bed(paths["promoters"]) == [
        ["chr1", "0", "1100", "Gene1|g1", "0", "+"],
        ["chr2", "600", "2600", "Gene2|g2", "0", "-"],
    ]


def test_native_gene_span_is_preferred_and_alternative_tss_are_retained(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr3\ttest\tgene\t1001\t3000\t.\t+\t.\tgene_id \"g3\"; gene_name \"Gene3\";\n"
        "chr3\ttest\ttranscript\t1201\t2000\t.\t+\t.\tgene_id \"g3\"; transcript_id \"t3a\"; gene_name \"Gene3\";\n"
        "chr3\ttest\ttranscript\t1501\t2500\t.\t+\t.\tgene_id \"g3\"; transcript_id \"t3b\"; gene_name \"Gene3\";\n"
        "chr3\ttest\texon\t1201\t1300\t.\t+\t.\tgene_id \"g3\"; transcript_id \"t3a\"; gene_name \"Gene3\";\n"
    )

    paths = build_gtf_annotation_sources(gtf, tmp_path / "annotation", 100, 50)

    assert _read_bed(paths["genes"]) == [["chr3", "1000", "3000", "Gene3|g3", "0", "+"]]
    assert _read_bed(paths["tss"]) == [
        ["chr3", "1200", "1201", "Gene3|g3", "0", "+"],
        ["chr3", "1500", "1501", "Gene3|g3", "0", "+"],
    ]
    assert promoter_intervals_by_gene(gtf, 100, 50)["Gene3|g3"] == [
        ("chr3", 1100, 1250),
        ("chr3", 1400, 1550),
    ]


def test_promoter_map_retains_duplicate_labels_on_separate_contigs(tmp_path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr1\ttest\ttranscript\t101\t500\t.\t+\t.\tgene_id \"g1\"; gene_name \"Gene1\";\n"
        "chrUn\ttest\ttranscript\t201\t600\t.\t+\t.\tgene_id \"g1\"; gene_name \"Gene1\";\n"
    )

    assert promoter_intervals_by_gene(gtf, 100, 50)["Gene1|g1"] == [
        ("chr1", 0, 150),
        ("chrUn", 100, 250),
    ]
