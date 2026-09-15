import contextlib
import copy
import csv
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pysam
import yaml

from omnomnomics import cli
from omnomnomics.chip_inputs import resolve_chip_inputs, validate_preparation_scope
from omnomnomics.library_manifest import build_library_manifest, save_library_manifest
from omnomnomics.metadata import MetadataError, derive_metadata_rows, fastq_unit_identity
from omnomnomics.preprocessing import LibraryRuntime


def resolve(rows, selector="cell", name="library"):
    fields = list(rows[0])
    _, derived, _ = derive_metadata_rows(fields, rows, name)
    return resolve_chip_inputs(fields, derived, selector)


def row(filename, role, library=None, cell="A", replicate="1", **extra):
    return {"filename": filename, "role": role, "library": library or filename.split('.')[0], "cell": cell, "replicate": replicate, **extra}


class MatchingTests(unittest.TestCase):
    def test_shared_and_replicate_matching(self):
        rows = [row("a", "chip"), row("b", "chip", replicate="2"), row("c", "input"), row("d", "input", replicate="2")]
        shared = resolve(rows)
        self.assertEqual(len(shared["associations"]), 4)
        matched = resolve(rows, "cell,replicate")
        self.assertEqual([(a["chip_id"], a["input_id"]) for a in matched["associations"]], [("a", "c"), ("b", "d")])
        indexed = resolve(rows, "4,5")
        self.assertEqual(matched["associations"], indexed["associations"])

    def test_explicit_groups_and_unused_inputs(self):
        rows = [row("a", "chip", input_group="shared"), row("b", "input", input_group="shared"), row("c", "input", input_group="unused")]
        result = resolve(rows, "input_group")
        self.assertEqual(result["unassigned_inputs"], ["c"])
        self.assertEqual(result["associations"][0]["input_id"], "b")

    def test_technical_inputs_are_not_duplicated(self):
        rows = [row("a", "chip", technical_replicate="1"), row("b1", "input", "b", technical_replicate="1"), row("b2", "input", "b", technical_replicate="2")]
        self.assertEqual(len(resolve(rows)["associations"]), 1)
        rows[2]["role"] = "chip"
        with self.assertRaisesRegex(MetadataError, "conflicting 'role'"):
            resolve(rows)

    def test_missing_invalid_and_zero_match_values(self):
        for role, cell, pattern in (("", "A", "role"), ("control", "A", "role"), ("input", "", "missing"), ("input", "B", "No input matches")):
            with self.subTest(role=role, cell=cell), self.assertRaisesRegex(MetadataError, pattern):
                resolve([row("a", "chip"), row("b", role, cell=cell)])
        with self.assertRaisesRegex(MetadataError, "--input-match"):
            resolve([row("a", "chip"), row("b", "input")], None)

    def test_matching_uses_original_values(self):
        with self.assertRaisesRegex(MetadataError, "No input matches"):
            resolve([row("a", "chip", cell="cell A"), row("b", "input", cell="cell_A")])

    def test_input_only_and_legacy_conflict(self):
        rows = [row("a", "input")]
        self.assertEqual(resolve(rows, None)["associations"], [])
        fields = list(rows[0])
        _, derived, _ = derive_metadata_rows(fields, rows, "library")
        with self.assertRaisesRegex(MetadataError, "cannot be combined"):
            resolve_chip_inputs(fields, derived, None, "/control.bam")

    def test_downstream_scope_is_explicit(self):
        validate_preparation_scope([1, 3, 4, 5, 6, 7])
        for steps in ([8], [10], [11, 15]):
            validate_preparation_scope(steps)
        with self.assertRaisesRegex(MetadataError, "HOMER"):
            validate_preparation_scope([7], create_homer_tagdirs=True)

    def test_fastq_identity_preserves_lane_and_accepts_common_suffixes(self):
        for filename in ("sample_L001_R1_001.fastq.gz", "sample_L001_R1.fq", "sample_L001_1_001.fastq", "sample_L001_1.fq.gz"):
            self.assertEqual(fastq_unit_identity(filename), ("sample_L001", "R1"))


class LibraryFixture(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.project = self.root / "project"
        self.project.mkdir()
        self.assembly = self.root / "genomes"
        genome = self.assembly / "test"
        (genome / "fasta").mkdir(parents=True)
        (genome / "fasta/genome.fa").write_text(">chr1\n" + "A" * 1000 + "\n")
        (genome / "fasta/genome.fa.fai").write_text("chr1\t1000\t6\t1000\t1001\n")
        (genome / "aux").mkdir()
        (genome / "aux/test_chrom_sizes.2_column").write_text("chr1\t1000\n")
        self.config = yaml.safe_load(cli.DEFAULT_WORKFLOW_CONFIG.read_text())
        self.config["map_tool"] = "hisat2"

    def bam(self, name, paired=False, folder="BAM", length=1000):
        path = self.project / folder / name
        path.parent.mkdir(exist_ok=True)
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": length}]}
        with pysam.AlignmentFile(str(path), "wb", header=header) as handle:
            for index in range(4):
                read = pysam.AlignedSegment()
                read.query_name = f"read{index // 2 if paired else index}"
                read.query_sequence = "A" * 30
                read.flag = (99 if index % 2 == 0 else 147) if paired else 0
                read.reference_id = 0
                read.reference_start = 100 + index * 40
                read.mapping_quality = 60
                read.cigarstring = "30M"
                read.query_qualities = pysam.qualitystring_to_array("I" * 30)
                if paired:
                    read.next_reference_id = 0
                    read.next_reference_start = read.reference_start + (40 if index % 2 == 0 else -40)
                    read.template_length = 70 if index % 2 == 0 else -70
                handle.write(read)
        return path

    def fastq(self, unit, paired=False, folder="FASTQ", suffix=".fastq.gz"):
        directory = self.project / folder
        directory.mkdir(exist_ok=True)
        paths = []
        for read in (["_R1", "_R2"] if paired else [""]):
            path = directory / f"{unit}{read}{suffix}"
            content = "@read\n" + "A" * 80 + "\n+\n" + "I" * 80 + "\n"
            if path.name.endswith(".gz"):
                with gzip.open(path, "wt") as handle:
                    handle.write(content)
            else:
                path.write_text(content)
            paths.append(path)
        return paths

    def manifest(self, rows, steps=range(1, 8), previous=None):
        return build_library_manifest(self.project, resolve(rows), list(steps), "test", self.assembly, self.config, previous)

    def cli_config(self, rows, steps="1-7", extra=(), assay="chip"):
        metadata = self.project / "metadata.tsv"
        with metadata.open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        site = yaml.safe_load(cli.DEFAULT_SITE_CONFIG.read_text())
        site["genome_assembly_dir"] = str(self.assembly)
        site["cellranger_reference_dir"] = str(self.root / "references")
        site["cores_per_node"] = 16
        site["nodes_in_partition"] = 1
        site["max_nodes"] = 1
        site_path = self.root / "site.yaml"
        site_path.write_text("## Omnomnomics pipeline config ##\n" + yaml.safe_dump(site))
        args = ["omnomnomics", assay, "-i", str(self.project), "-g", "test", "-j", steps, "-m", str(metadata), "--sample-name", "library", "--sample-type", "cell", "--site-config", str(site_path), "--dry-run", "--no-multiqc", "--trim-tool", "fastp", *extra]
        if assay == "chip" and "role" in rows[0]:
            args.extend(["--input-match", "cell"])
        with patch("sys.argv", args), patch.object(cli.subprocess, "run") as dispatch, contextlib.redirect_stdout(io.StringIO()):
            cli.main()
        command = dispatch.call_args.args[0]
        self.dispatch_command = command
        return Path(command[command.index("--config") + 1].split("=", 1)[1])


class ManifestTests(LibraryFixture):
    def test_unknown_existing_target_is_protected_before_dispatch(self):
        self.bam("facility.bam")
        original = self.bam("renamed.bam")
        self.bam("b.bam")
        before = original.read_bytes()
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            self.cli_config([row("facility.bam", "chip", "renamed"), row("b.bam", "input")], "4-5")
        self.assertEqual(original.read_bytes(), before)

    def test_missing_pair_and_declared_layout(self):
        reads = self.fastq("a", paired=True)
        reads[1].unlink()
        self.bam("b.bam")
        rows = [row("a", "chip"), row("b.bam", "input")]
        with self.assertRaisesRegex(MetadataError, "Only R1"):
            self.manifest(rows)
        rows[0]["read_layout"] = "PE"
        with self.assertRaisesRegex(MetadataError, "Missing read pair"):
            self.manifest(rows)
        rows[0]["read_layout"] = "SE"
        self.assertFalse(self.manifest(rows)["libraries"][0]["paired"])
        rows[1]["read_layout"] = "PE"
        with self.assertRaisesRegex(MetadataError, "disagrees"):
            self.manifest(rows)

    def test_cross_library_unit_collision(self):
        self.bam("a.bam")
        self.bam("b.bam")
        with self.assertRaisesRegex(MetadataError, "collides"):
            self.manifest([row("a.bam", "chip", "b"), row("b.bam", "input", "c")], [4, 5])

    def test_explicit_bam_overrides_existing_merged_product(self):
        self.bam("renamed.bam")
        supplied = self.bam("facility.bam", paired=True)
        self.bam("b.bam")
        manifest = self.manifest([row("facility.bam", "chip", "renamed"), row("b.bam", "input")], [4, 5])
        lib = next(lib for lib in manifest["libraries"] if lib["sample_id"] == "renamed")
        self.assertEqual(lib["source_paths"], [str(supplied)])
        self.assertTrue(lib["paired"])

    def test_invalid_metadata_force_does_not_delete_outputs(self):
        self.fastq("a")
        source = self.bam("b.bam")
        output = self.bam("a.filtered.bam", folder="filtered_BAM")
        contents = output.read_bytes()
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            self.cli_config([row("a", "chip"), row("b.bam", "input", cell="B")], extra=("--rerun-selected-steps",))
        self.assertEqual(output.read_bytes(), contents)
        self.assertTrue(source.exists())

    def test_generated_provenance_survives_partial_restart(self):
        self.fastq("a")
        self.bam("b.bam")
        rows = [row("a", "chip", "renamed"), row("b.bam", "input")]
        first = self.manifest(rows)
        self.bam("renamed.filtered.bam", folder="filtered_BAM")
        self.bam("b.filtered.bam", folder="filtered_BAM")
        later = self.manifest(rows, [6], previous=first)
        self.assertIn(str(self.project / "BAM/renamed.bam"), later["generated_paths"])
        self.assertEqual(later["libraries"][1]["original_sources"], first["libraries"][1]["original_sources"])

    def test_mixed_fastq_and_bam_and_read_layouts(self):
        self.fastq("a", paired=True)
        bam = self.bam("b.bam")
        manifest = self.manifest([row("a", "chip"), row("b.bam", "input")])
        runtime = LibraryRuntime(manifest)
        self.assertEqual(set(manifest["work"]["3"]), {"a"})
        self.assertEqual(set(manifest["work"]["5"]), {"a", "b"})
        self.assertTrue(runtime.paired("a"))
        self.assertFalse(runtime.paired("b"))
        self.assertEqual(runtime.mapped_units("b"), [])
        self.assertEqual(runtime.bam("b"), str(bam))

    def test_reverse_mix_and_input_only_preparation(self):
        self.bam("a.bam", paired=True)
        self.fastq("b")
        manifest = self.manifest([row("a.bam", "chip"), row("b", "input")])
        self.assertEqual(set(manifest["work"]["1"]), {"b"})
        declarations = resolve([row("b", "input")], None)
        result = build_library_manifest(self.project, declarations, [1, 3], "test", self.assembly, self.config)
        self.assertEqual(result["associations"], [])

    def test_technical_lanes_merge_and_keep_control_independent(self):
        self.fastq("a_L001", paired=True)
        self.fastq("a_L002", paired=True)
        self.bam("b.bam")
        result = self.manifest([row("a", "chip"), row("b.bam", "input")])
        runtime = LibraryRuntime(result)
        self.assertEqual(runtime.input_units("a"), ["a_L001", "a_L002"])
        self.assertTrue(runtime.libraries["a"]["merge_output"])

    def test_ambiguous_source_requires_explicit_type(self):
        self.fastq("a")
        self.bam("a.bam")
        self.bam("b.bam")
        with self.assertRaisesRegex(MetadataError, "Ambiguous"):
            self.manifest([row("a", "chip"), row("b.bam", "input")])
        self.manifest([row("a.fastq.gz", "chip"), row("b.bam", "input")])

    def test_bad_reference_and_missing_intermediate_fail_before_running(self):
        self.bam("a.bam", length=999)
        self.bam("b.bam")
        rows = [row("a.bam", "chip"), row("b.bam", "input")]
        with self.assertRaisesRegex(MetadataError, "reference"):
            self.manifest(rows, [5, 6])
        self.bam("a.bam")
        with self.assertRaisesRegex(MetadataError, "existing filtered BAM"):
            self.manifest(rows, [6, 7])

    def test_filtered_entry_does_not_require_original_bams(self):
        self.bam("a.filtered.bam", folder="filtered_BAM")
        self.bam("b.filtered.bam", folder="filtered_BAM")
        result = self.manifest([row("a", "chip"), row("b", "input")], [6, 7])
        self.assertEqual(result["work"].keys(), {"6", "7"})
        self.assertTrue(all(not lib["units"] for lib in result["libraries"]))

    def test_manifest_and_associations_are_stable(self):
        self.bam("a.bam")
        self.bam("b.bam")
        rows = [row("a.bam", "chip"), row("b.bam", "input")]
        result = self.manifest(rows, [5, 6])
        first, associations = save_library_manifest(self.project, "one", result)
        timestamp = Path(first).stat().st_mtime_ns
        second, _ = save_library_manifest(self.project, "two", self.manifest(rows, [5, 6], previous=result))
        self.assertEqual(first, second)
        self.assertEqual(Path(first).stat().st_mtime_ns, timestamp)
        self.assertIn("a\tb\tA", Path(associations).read_text())

    def test_cli_writes_manifest_before_dispatch(self):
        self.fastq("a", paired=True)
        self.bam("b.bam")
        path = self.cli_config([row("a", "chip"), row("b.bam", "input")])
        config = yaml.safe_load(path.read_text())
        self.assertTrue(Path(config["LIBRARY_MANIFEST_FILE"]).is_file())
        self.assertTrue(Path(config["INPUT_ASSOCIATIONS_FILE"]).is_file())
        self.assertEqual(config["INPUT_FOLDER"], "metadata")
        with Path(config["EXPERIMENT_METADATA_FILE"]).open() as handle:
            experiments = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([entry["role"] for entry in experiments], ["chip"])
        totals = {entry["step_num"]: entry["total"] for entry in cli.load_step_monitor_rows(self.project)}
        self.assertEqual(totals[3], 1)
        self.assertEqual(totals[5], 2)


class SnakemakePreprocessingTests(LibraryFixture):
    def run_workflow(self, config_path, dry_run=True):
        env = dict(os.environ)
        env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
        java_home = Path(sys.prefix) / "lib/jvm"
        if (java_home / "bin/java").is_file():
            env["JAVA_HOME"] = str(java_home)
            env["PATH"] = str(java_home / "bin") + os.pathsep + env["PATH"]
        env["PYTHONPATH"] = str(cli.PACKAGE_ROOT.parent)
        command = [sys.executable, "-m", "snakemake", "--snakefile", str(cli.WORKFLOW_ROOT / "Snakefile.smk"), "--config", f"config_file={config_path}", "--cores", "2", "--scheduler", "greedy", "--rerun-triggers", "mtime", "--latency-wait", "1", "--nocolor"]
        if dry_run:
            command.append("--dry-run")
        if "--forcerun" in self.dispatch_command:
            command.extend(self.dispatch_command[self.dispatch_command.index("--forcerun"):])
        result = subprocess.run(
            command,
            cwd=self.project,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=int(os.environ.get("OMNOM_TEST_TIMEOUT", "120")),
        )
        self.assertEqual(result.returncode, 0, result.stdout[-14000:])
        return result.stdout

    def test_mixed_entry_dag_has_no_alignment_for_supplied_input(self):
        self.fastq("a", paired=True)
        self.bam("b.bam")
        path = self.cli_config([row("a", "chip"), row("b.bam", "input")])
        output = self.run_workflow(path)
        self.assertIn("rule run_hisat2:", output)
        self.assertIn("sample3=a", output)
        self.assertNotIn("sample3=b", output)
        self.assertNotIn("b.extra_3.tmp", output)

    def test_single_end_and_paired_end_fastqs_have_distinct_rules(self):
        self.fastq("a", paired=True)
        self.fastq("b")
        path = self.cli_config([row("a", "chip"), row("b", "input")], "1-3", extra=("--rerun-selected-steps",))
        output = self.run_workflow(path)
        self.assertIn("rule run_fastp_se:", output)
        self.assertIn("rule run_fastp:", output)
        self.assertIn("rule run_fastqc_se:", output)

    def test_hisat2_preserves_reference_contig_names(self):
        if not shutil.which("hisat2-build"):
            self.skipTest("hisat2-build is not available")
        genome = self.assembly / "test"
        fasta = genome / "fasta/genome.fa"
        fasta.write_text(">1\n" + "ACGT" * 250 + "\n")
        (genome / "fasta/genome.fa.fai").write_text("1\t1000\t3\t1000\t1001\n")
        index_dir = genome / "hisat2"
        index_dir.mkdir()
        subprocess.run(
            ["hisat2-build", str(fasta), str(index_dir / "test")],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        self.fastq("a")
        config = self.cli_config([row("a", "chip", read_layout="SE")], "1-3")
        self.run_workflow(config, dry_run=False)
        with pysam.AlignmentFile(str(self.project / "BAM/a.bam"), "rb") as bam:
            self.assertEqual(bam.references, ("1",))

    def test_existing_rna_and_atac_bam_entry_dags(self):
        self.bam("a.bam")
        rows = [{"filename": "a", "library": "a", "cell": "A"}]
        for assay in ("rna", "atac"):
            with self.subTest(assay=assay):
                config = self.cli_config(rows, "5-6", assay=assay)
                output = self.run_workflow(config)
                self.assertIn("rule touchup_bam:", output)
                self.assertNotIn("rule run_hisat2:", output)

    def test_real_mixed_layout_fastq_preparation(self):
        sources = self.fastq("a", paired=True) + self.fastq("b")
        before = [path.read_bytes() for path in sources]
        config = self.cli_config([row("a", "chip"), row("b", "input")], "1-2")
        self.run_workflow(config, dry_run=False)
        for name in ("a_R1", "a_R2", "b"):
            self.assertTrue((self.project / "fastqc_reports" / f"{name}.trimmed_fastqc.html").is_file())
        self.assertEqual([path.read_bytes() for path in sources], before)

    def test_fastqc_from_noncanonical_trimmed_filenames(self):
        self.fastq("a_R1_001", folder="trimmed_FASTQ", suffix=".fq.gz")
        self.bam("b.bam")
        rows = [row("a", "chip", read_layout="SE"), row("b.bam", "input", read_layout="SE")]
        config = self.cli_config(rows, "2")
        self.run_workflow(config, dry_run=False)
        self.assertTrue((self.project / "fastqc_reports/a.trimmed_fastqc.html").exists())

    def test_real_bam_preparation_preserves_sources_and_reuses_outputs(self):
        sources = [self.bam("a.bam", paired=True), self.bam("b.bam")]
        original_bytes = [path.read_bytes() for path in sources]
        rows = [row("a.bam", "chip"), row("b.bam", "input")]
        path = self.cli_config(rows, "4-7")
        self.run_workflow(path, dry_run=False)
        for name in ("a", "b"):
            result = self.project / "filtered_BAM" / f"{name}.filtered.bam"
            self.assertTrue(result.is_file())
            self.assertTrue(Path(f"{result}.bai").is_file())
            pysam.quickcheck(str(result))
            stats = Path(f"{result}.stats.txt").read_text()
            self.assertIn(f"# role\t{'chip' if name == 'a' else 'input'}", stats)
        self.assertEqual([path.read_bytes() for path in sources], original_bytes)
        path = self.cli_config(rows, "4-7")
        output = self.run_workflow(path)
        self.assertIn("Nothing to be done", output)
        self.bam("c.bam")
        path = self.cli_config([*rows, row("c.bam", "chip")], "5-6")
        self.run_workflow(path, dry_run=False)
        progress = cli.load_step_monitor_rows(self.project)
        self.assertTrue(all(entry["state"] == "DONE" and entry["completed"] == 3 for entry in progress if entry["step_num"] in (5, 6)))
        self.assertIn("REUSED", (self.project / "run_logs/steps/step05.summary.tsv").read_text())

    def test_filtered_bam_import_and_supplied_index(self):
        a = self.bam("a.filtered.bam", folder="filtered_BAM")
        b = self.bam("facility.filtered.bam", folder="filtered_BAM")
        pysam.index(str(a))
        sources = [a, Path(f"{a}.bai"), b]
        original = [path.read_bytes() for path in sources]
        config = self.cli_config([row("a.filtered.bam", "chip", "a"), row("facility.filtered.bam", "input", "b")], "6")
        dag = self.run_workflow(config)
        self.assertIn("rule import_filtered_bam:", dag)
        self.assertNotIn("wildcards: sample=a", dag)
        self.run_workflow(config, dry_run=False)
        self.assertEqual([path.read_bytes() for path in sources], original)
        self.assertTrue((self.project / "filtered_BAM/b.filtered.bam.bai").exists())
        config = self.cli_config([row("a.filtered.bam", "chip", "a"), row("facility.filtered.bam", "input", "b")], "6")
        self.assertIn("Nothing to be done", self.run_workflow(config))


if __name__ == "__main__":
    unittest.main()
