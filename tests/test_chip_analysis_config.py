import contextlib
import io
from unittest.mock import patch

import yaml
import pysam

import test_chip_inputs as fixtures


class AnalysisConfigTests(fixtures.LibraryFixture):
    run_workflow = fixtures.SnakemakePreprocessingTests.run_workflow

    def fixture(self):
        self.bam("a.filtered.bam", folder="filtered_BAM")
        self.bam("i.filtered.bam", paired=True, folder="filtered_BAM")
        bed = self.project / "atlas.bed"
        bed.write_text("chr1\t100\t250\n")
        return [fixtures.row("a.filtered.bam", "chip", "a"), fixtures.row("i.filtered.bam", "input", "i")], bed

    def override(self, values):
        original = fixtures.cli.load_pipeline_config
        def load(*args):
            config = original(*args)
            config.update(values)
            return config
        return patch.object(fixtures.cli, "load_pipeline_config", side_effect=load)

    def test_invalid_site_analysis_settings_fail_before_forced_cleanup(self):
        rows, bed = self.fixture()
        source = self.project / "filtered_BAM/a.filtered.bam"
        before = source.read_bytes()
        for values in ({"chip_input_tracks": "fold"}, {"chip_enrichment_min_input": 0}, {"chip_joint_peak_q": None}, {"chip_se_fragment_length": None}, {"chip_effective_genome_size": float("nan")}):
            with self.subTest(values=values), self.override(values), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as failure:
                self.cli_config(rows, "13", extra=("--chip-regions-bed", str(bed), "--rerun-selected-steps"))
            self.assertEqual(failure.exception.code, 1)
            self.assertEqual(source.read_bytes(), before)

    def test_site_enrichment_threshold_is_forwarded_and_invalidates_summaries(self):
        import csv
        rows, bed = self.fixture()
        raw_before = None
        for minimum, status in ((5, "low_input"), (1, "ok")):
            with self.override({"chip_enrichment_min_input": minimum}):
                config = self.cli_config(rows, "13", extra=("--chip-regions-bed", str(bed)))
            self.assertEqual(yaml.safe_load(config.read_text())["CHIP_ENRICHMENT_MIN_INPUT"], minimum)
            self.run_workflow(config, dry_run=False)
            count_dir = self.project / self.config["output_folders"][self.config["countreads_rule_num"] - 1]
            raw = (count_dir / "project.raw_read_quant.table.txt").read_bytes()
            if raw_before is not None:
                self.assertEqual(raw, raw_before)
            raw_before = raw
            with (count_dir / "chip_regional_enrichment.tsv").open() as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["reliability"], status)

    def test_capped_run_logs_optional_track_omission(self):
        rows, _bed = self.fixture()
        for path in (self.project / "filtered_BAM").glob("*.bam"):
            pysam.index(str(path))
        config = self.cli_config(rows, "8", extra=("--max-project-size", "1M"))
        output = self.run_workflow(config)
        self.assertNotIn("rule chip_input_track:", output)
        records = "\n".join(path.read_text() for path in (self.project / "run_logs").glob("*.log"))
        self.assertIn("Optional input-relative track generation omitted", records)
