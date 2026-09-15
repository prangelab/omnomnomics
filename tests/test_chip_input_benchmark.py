"""Check simulation truth against paired-read unique-assignment geometry."""

import importlib.util
import gzip
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

spec = importlib.util.spec_from_file_location("chip_benchmark", Path(__file__).resolve().parents[1] / "scripts/validation/chip_input_benchmark.py")
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


class BenchmarkTruthTests(unittest.TestCase):
    def test_quality_selection_preserves_original_and_failed_selection(self):
        real_spec = importlib.util.spec_from_file_location("chip_real", Path(__file__).resolve().parents[1] / "scripts/validation/chip_input_real_validation.py")
        real = importlib.util.module_from_spec(real_spec)
        with patch.dict(sys.modules, {"chip_input_benchmark": benchmark}):
            real_spec.loader.exec_module(real)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.fastq.gz"
            bad = b"@bad\n" + b"A" * 20 + b"\n+\n" + b"!" * 20 + b"\n"
            good = b"@good\n" + b"A" * 20 + b"\n+\n" + b"I" * 20 + b"\n"
            source.write_bytes(gzip.compress(bad + good + bad + good))
            destination = root / "FASTQ/input.fastq.gz"
            destination.parent.mkdir()
            original = gzip.compress(bad)
            destination.write_bytes(original)
            with self.assertRaises(RuntimeError):
                real.quality_subset(destination, str(source), limit=2, max_records=2)
            self.assertEqual(destination.read_bytes(), original)
            self.assertFalse(list(root.rglob("*.preparing")))
            result = real.quality_subset(destination, str(source), limit=2, max_records=4)
            self.assertEqual(result["scanned_records"], 4)
            self.assertEqual(gzip.decompress(destination.read_bytes()), good + good)
            self.assertEqual((root / "source_selection_archive/GM12878_input2.original_prefix.fastq.gz").read_bytes(), original)
            self.assertEqual(real.quality_subset(destination, str(source), limit=2, max_records=4), result)

    def test_integer_start_assignment_matches_explicit_mate_overlaps(self):
        regions = [("chr1", 100, 110, "a"), ("chr1", 180, 190, "b"), ("chr1", 350, 360, "c")]
        expected = np.zeros(3)
        for insert, probability in zip(benchmark.INSERTS, benchmark.INSERT_PROBS):
            for start in range(500):
                overlaps = [index for index, (_, left, right, _) in enumerate(regions) if (start < right and start + 50 > left) or (start + insert - 50 < right and start + insert > left)]
                if len(overlaps) == 1:
                    expected[overlaps[0]] += probability / 500
        actual = benchmark.assignment(regions, np.array([0]), np.array([500]), np.array([1.]))
        np.testing.assert_allclose(actual, expected, atol=1e-14)

    def test_streams_are_stable_and_independent(self):
        np.testing.assert_array_equal(benchmark.stream("B1", 0, "WT0").integers(1000, size=10), benchmark.stream("B1", 0, "WT0").integers(1000, size=10))
        self.assertFalse(np.array_equal(benchmark.stream("B1", 0, "WT0").integers(1000, size=10), benchmark.stream("B2", 0, "WT0").integers(1000, size=10)))

    def test_scenarios_keep_binding_and_background_truth_separate(self):
        _, wt, ko, bg, bgko = benchmark.model("B5")
        np.testing.assert_array_equal(wt, ko)
        self.assertEqual(np.count_nonzero(bg != bgko), 200)
        _, wt, ko, _, _ = benchmark.model("B3")
        self.assertEqual(np.count_nonzero(ko > wt), 100)
        self.assertEqual(np.count_nonzero(ko < wt), 100)
        _, wt, ko, _, _ = benchmark.model("B4")
        self.assertEqual(np.count_nonzero(ko == 0), 200)
        _, wt, ko, _, _ = benchmark.model("B6")
        np.testing.assert_array_equal(ko, wt * .5)

    def test_featurecounts_truth_geometry(self):
        from omnomnomics.chip_analysis import count_libraries
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            positions = np.arange(80, 130)
            inserts = np.full(len(positions), 200)
            bam = root / "reads.bam"
            benchmark.write_pairs(bam, positions, inserts, root, 1)
            regions = [("chr1", 100, 110, "a"), ("chr1", 270, 280, "b")]
            counts, depths = count_libraries(regions, {"reads": {"bam": bam, "paired": True}}, root / "counts")
            explicit = [0, 0]
            for start in positions:
                hits = [index for index, (_, left, right, _) in enumerate(regions) if (start < right and start + 50 > left) or (start + 150 < right and start + 200 > left)]
                if len(hits) == 1:
                    explicit[hits[0]] += 1
            self.assertEqual(counts["reads"], explicit)
            self.assertEqual(depths["reads"], len(positions))

    def test_global_loss_has_separate_declared_references(self):
        sites, *_ = benchmark.model("B6")
        regions = [("chr1", int(sites[20] - 500), int(sites[20] + 500), "site")]
        values = benchmark.expectations(regions, "B6")
        self.assertEqual(values[1][0], values[0][0] * .5)
        self.assertLess(values[5][0], -.9)
        self.assertGreater(values[6][0], 0)
        self.assertLess(values[6][0], .03)

    def test_report_preserves_null_false_calls_and_missing_padj(self):
        report_spec = importlib.util.spec_from_file_location("chip_report", Path(__file__).resolve().parents[1] / "scripts/validation/chip_input_benchmark_report.py")
        report = importlib.util.module_from_spec(report_spec)
        with patch.dict(sys.modules, {"chip_input_benchmark": benchmark}):
            report_spec.loader.exec_module(report)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            metric = root / "B1_seed0/metrics/U"
            metric.mkdir(parents=True)
            (metric / "run.json").write_text('{"seconds": 1}\n')
            (metric / "truth.tsv").write_text("region_id\tstart\tend\tspecific_WT\tspecific_KO\texpected_total_log2FC\nchr1_9500_10500\t9500\t10500\t20\t20\t0\nchr1_19500_20500\t19500\t20500\t0\t0\t0\n")
            (metric / "unshrunk.tsv").write_text("region_id\tpadj\tlog2FoldChange\nchr1_9500_10500\t0.01\t2\nchr1_19500_20500\tNA\t0\n")
            (metric / "display_results.tsv").write_text("gene_id\tpadj\tlog2FoldChange\nchr1_9500_10500\t0.01\t0.5\nchr1_19500_20500\tNA\t0\n")
            (metric / "chip_regional_enrichment.tsv").write_text("reliability\nlow_input\n")
            (metric / "size_factors.tsv").write_text("library\tsize_factor\nKO0\t1\nWT0\t1\n")
            (metric / "chip_testing_metadata.tsv").write_text("support_status\nsupported\njoint_only\n")
            result = json.loads(json.dumps(report.summarise(root, "B1", 0, "U")))
            self.assertEqual(result["significant"], 1)
            self.assertEqual(result["dataset_FDP"], 1)
            self.assertEqual(result["padj_unavailable"], 1)
            self.assertEqual(result["displayed_hits"], 0)
            self.assertEqual(result["wrong_total_reference_direction"], 0)


if __name__ == "__main__":
    unittest.main()
