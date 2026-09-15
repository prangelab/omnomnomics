import contextlib
import fcntl
import glob
import io
import json
import os
from pathlib import Path
import tempfile
import textwrap
from types import SimpleNamespace
import unittest

import yaml

from omnomnomics.cli import DEFAULT_WORKFLOW_CONFIG, delete_outputs_to_be_updated
from omnomnomics.storage import (
    capture_source_protection,
    remove_generated_path,
    source_protects_path,
)


ROOT = Path(__file__).resolve().parents[1]
SNAKEFILE = ROOT / "src/omnomnomics/workflow/Snakefile.smk"


def load_cleanup_functions(namespace):
    """Exercise the workflow's Python cleanup functions without a scheduler."""
    lines = SNAKEFILE.read_text().splitlines(keepends=True)
    names = (
        "library_cleanup_ready",
        "cleanup_generated_output", "safe_cleanup_for_size_limit",
        "evaluate_post_step_size_cleanup", "terminal_output_dirs",
        "retention_keep_dirs", "apply_retention_policy",
    )
    for name in names:
        start = next(i for i, line in enumerate(lines) if line.startswith(f"def {name}("))
        end = start + 1
        while end < len(lines):
            line = lines[end]
            if line.strip() and not line[0].isspace() and not line.startswith("#"):
                break
            end += 1
        exec(compile("".join(lines[start:end]), str(SNAKEFILE), "exec"), namespace)
    return namespace


class SourceProtectionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.config = yaml.safe_load(DEFAULT_WORKFLOW_CONFIG.read_text())

    def put(self, relative, contents=b"original input\n"):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(contents)
        return path

    def capture(self, additional_sources=()):
        return capture_source_protection(self.root, self.config, additional_sources)

    def cleanup_namespace(self, protected):
        markers = self.root / "markers"
        markers.mkdir(exist_ok=True)
        for step in (4, 5):
            (markers / f"{step}.done").touch()
        self.messages = []
        return load_cleanup_functions({
            "os": os, "fcntl": fcntl, "master_config": self.config,
            "experiment_dir": str(self.root), "protected_source_paths": protected,
            "remove_generated_path": remove_generated_path,
            "logfile": "test.log", "log_marker_dir": str(markers),
            "log_it": lambda _log, message, *args: self.messages.append(message),
            "cache_flow_qc_metrics_from_folder": lambda *args: None,
            "step_tracking_paths": lambda step: {"finished_marker": str(markers / f"{step}.done")},
            "max_project_size_bytes": 1, "project_size_bytes": lambda: 1000,
            "format_bytes": str,
        })

    def test_fastq_sources_keep_both_reads_and_allow_downstream_cleanup(self):
        r1 = self.put("FASTQ/a_R1.fastq.gz")
        r2 = self.put("FASTQ/a_R2.fastq.gz")
        trimmed = self.put("trimmed_FASTQ/a_R1.trimmed.fastq.gz")
        bam = self.put("BAM/a.bam")
        _, protected = self.capture()
        self.assertTrue(source_protects_path(r1, protected))
        self.assertTrue(source_protects_path(r2, protected))
        self.assertTrue(remove_generated_path(trimmed.parent, protected))
        self.assertTrue(remove_generated_path(bam.parent, protected))
        self.assertEqual(r1.read_bytes(), b"original input\n")
        self.assertEqual(r2.read_bytes(), b"original input\n")

    def test_each_sequence_entry_stage_is_preserved_by_retention(self):
        for policy in ("pruned", "minimal"):
            for folder, name in (
                ("trimmed_FASTQ", "a.trimmed.fq.gz"),
                ("BAM", "a.bam"),
                ("filtered_BAM", "a.filtered.bam"),
                ("BigWigs", "a.bw"),
            ):
                with self.subTest(policy=policy, folder=folder), tempfile.TemporaryDirectory() as directory:
                    original_root = self.root
                    self.root = Path(directory).resolve()
                    try:
                        source = self.put(f"{folder}/{name}")
                        _, protected = self.capture()
                        ns = self.cleanup_namespace(protected)
                        # Later stages make the source directory eligible under old policies.
                        ns["apply_retention_policy"]("test.log", policy, [10, 11], "CHIP")
                        self.assertEqual(source.read_bytes(), b"original input\n")
                    finally:
                        self.root = original_root

    def test_bam_indexes_and_lane_named_sources_are_protected(self):
        source = self.put("BAM/a_L001.bam")
        _, protected = self.capture()
        index = self.put("BAM/a_L001.bam.bai")
        for path in (source, index, source.parent):
            self.assertFalse(remove_generated_path(path, protected))
        self.assertEqual(source.read_bytes(), b"original input\n")

    def test_mixed_samples_preserve_bam_without_its_own_fastq(self):
        self.put("FASTQ/a.fastq.gz")
        generated = self.put("BAM/a.bam")
        source = self.put("BAM/b.bam")
        _, protected = self.capture()
        self.assertTrue(remove_generated_path(generated, protected))
        self.assertFalse(remove_generated_path(source, protected))
        self.assertFalse(remove_generated_path(source.parent, protected))

    def test_all_lanes_at_first_available_stage_are_preserved(self):
        sources = [self.put(f"BAM/a_L00{lane}.bam") for lane in (1, 2)]
        _, protected = self.capture()
        for path in sources:
            self.assertFalse(remove_generated_path(path, protected))

    def test_end_of_run_lane_cleanup_preserves_source_and_removes_generated_lane(self):
        source = self.put("BAM/source_L001.bam")
        _, protected = self.capture()
        generated = self.put("BAM/generated_L001.bam")
        ns = self.cleanup_namespace(protected)
        ns.update({"themode": [4], "glob": glob, "library_runtime": None})
        lines = SNAKEFILE.read_text().split("# Remove old BAM files with lane info", 1)[1]
        condition = "        if 4 in themode and not library_runtime:"
        block = condition + lines.split(condition, 1)[1]
        block = block.split("        apply_retention_policy", 1)[0]
        exec(compile(textwrap.dedent(block), str(SNAKEFILE), "exec"), ns)
        self.assertEqual(source.read_bytes(), b"original input\n")
        self.assertFalse(generated.exists())

    def test_recorded_bam_remains_protected_on_later_runs(self):
        source = self.put("BAM/a.bam")
        manifest, _ = self.capture()
        self.put("FASTQ/a.fastq.gz")
        self.put("filtered_BAM/a.filtered.bam")
        self.put("BigWigs/a.bw")
        _, protected = self.capture()
        self.assertTrue(source_protects_path(source, protected))
        self.assertIn("BAM/a.bam", {record["path"] for record in json.loads(Path(manifest).read_text())["sources"]})

    def test_failed_or_unfinished_manifest_consumer_prevents_cleanup(self):
        self.put("FASTQ/a.fastq.gz")
        generated = self.put("BAM/a.bam")
        _, protected = self.capture()
        ns = self.cleanup_namespace(protected)
        states = {}
        for step in (4, 5, 7):
            failed = self.root / "states" / str(step) / "failed"
            failed.mkdir(parents=True)
            finished = failed.parent / "finished.marker"
            finished.touch()
            states[step] = {"failed_dir": str(failed), "finished_marker": str(finished)}
        ns.update({"library_runtime": SimpleNamespace(work={str(step): {} for step in states}),
                   "library_manifest": {"folders": {"bam": str(generated.parent), "trim": str(self.root / "trimmed_FASTQ")}},
                   "ensure_step_tracking_dirs": lambda step: states[step]})
        failure = Path(states[7]["failed_dir"]) / "a"
        failure.touch()
        self.assertFalse(ns["cleanup_generated_output"](str(generated.parent)))
        failure.unlink()
        Path(states[7]["finished_marker"]).unlink()
        self.assertFalse(ns["cleanup_generated_output"](str(generated.parent)))
        Path(states[7]["finished_marker"]).touch()
        self.assertTrue(ns["cleanup_generated_output"](str(generated.parent)))

    def test_unchanged_capture_does_not_rewrite_manifest(self):
        self.put("BAM/a.bam")
        manifest, first = self.capture()
        timestamp = Path(manifest).stat().st_mtime_ns
        _, second = self.capture()
        self.assertEqual(first, second)
        self.assertEqual(timestamp, Path(manifest).stat().st_mtime_ns)

    def test_size_pressure_and_post_step_cleanup_preserve_sources(self):
        source = self.put("BAM/a.bam")
        _, protected = self.capture()
        generated = self.put("trimmed_FASTQ/temporary.trimmed.fastq.gz")
        ns = self.cleanup_namespace(protected)
        ns["safe_cleanup_for_size_limit"]("test.log")
        self.assertFalse(generated.exists())
        self.assertEqual(source.read_bytes(), b"original input\n")
        ns["evaluate_post_step_size_cleanup"]("test.log", 5)
        self.assertEqual(source.read_bytes(), b"original input\n")
        self.assertTrue(any("source input protection" in message for message in self.messages))

    def test_forced_cleanup_preserves_source_but_removes_generated_bam(self):
        source = self.put("BAM/source.bam")
        self.put("FASTQ/generated.fastq.gz")
        generated = self.put("BAM/generated.bam")
        _, protected = self.capture()
        with contextlib.redirect_stdout(io.StringIO()):
            delete_outputs_to_be_updated([3], self.config, str(self.root), protected)
        self.assertFalse(generated.exists())
        self.assertEqual(source.read_bytes(), b"original input\n")

    def test_direct_forced_cleanup_captures_sources_first(self):
        source = self.put("BAM/source.bam")
        with contextlib.redirect_stdout(io.StringIO()):
            delete_outputs_to_be_updated([3, 4], self.config, str(self.root))
        self.assertEqual(source.read_bytes(), b"original input\n")

    def test_table_only_entry_is_preserved(self):
        source = self.put("DE_calling/project.raw_read_quant.table.txt")
        _, protected = self.capture()
        with contextlib.redirect_stdout(io.StringIO()):
            delete_outputs_to_be_updated([11], self.config, str(self.root), protected)
        self.assertFalse(remove_generated_path(source.parent, protected))
        self.assertEqual(source.read_bytes(), b"original input\n")

    def test_legacy_control_and_metadata_paths_are_explicitly_protected(self):
        self.put("FASTQ/a.fastq.gz")
        control = self.put("BAM/a.bam")
        metadata = self.put("DE_calling/metadata.tsv")
        _, protected = self.capture((str(control), str(metadata)))
        self.assertFalse(remove_generated_path(control, protected))
        self.assertFalse(remove_generated_path(metadata.parent, protected))

    def test_symlink_source_and_target_are_preserved(self):
        target = self.put("external/source.bam")
        link = self.root / "BAM/a.bam"
        link.parent.mkdir()
        link.symlink_to(target)
        _, protected = self.capture()
        for path in (link, target, target.parent, link.parent):
            self.assertFalse(remove_generated_path(path, protected))
        self.assertEqual(target.read_bytes(), b"original input\n")

    def test_corrupt_manifest_stops_cleanup(self):
        source = self.put("BAM/a.bam")
        self.put("run_configs/source_protection.json", b"not valid JSON")
        with self.assertRaises(ValueError), contextlib.redirect_stdout(io.StringIO()):
            delete_outputs_to_be_updated([3], self.config, str(self.root))
        self.assertEqual(source.read_bytes(), b"original input\n")

    def test_old_run_without_snapshot_retains_outputs(self):
        source = self.put("BAM/a.bam")
        ns = self.cleanup_namespace(None)
        ns["apply_retention_policy"]("test.log", "minimal", [8, 9], "CHIP")
        self.assertEqual(source.read_bytes(), b"original input\n")
        self.assertTrue(any("missing source protection snapshot" in message for message in self.messages))


if __name__ == "__main__":
    unittest.main()
