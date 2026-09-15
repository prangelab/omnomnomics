import csv
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest

import pysam
import yaml

from omnomnomics.chip_analysis import (
    ControlSets, bam_fragments, enrichment_rows, macs_command,
    read_regions, support_annotations,
)
from omnomnomics.metadata import MetadataError
import test_chip_inputs as fixtures


class RegionTests(unittest.TestCase):
    def test_bed_validation_and_names(self):
        with tempfile.TemporaryDirectory() as temporary:
            bed = Path(temporary) / "regions.bed"
            bed.write_text("chr1\t0\t10\tfirst\nchr1\t10\t20\tsecond\n")
            self.assertEqual(read_regions(bed, {"chr1": 100})[0], ("chr1", 0, 10, "first"))
            for text in ("chr1\t0\t10\nchr1\t5\t15\n", "chr1\t0\t10\nchr1\t0\t10\n", "chr1\t-1\t10\n", "chr2\t0\t10\n", "chr1\t0\t101\n", "chr1 0 10\n", ""):
                bed.write_text(text)
                with self.subTest(text=text), self.assertRaises(MetadataError):
                    read_regions(bed, {"chr1": 100})

    def test_overlap_is_annotation_and_bases_are_not_double_counted(self):
        with tempfile.TemporaryDirectory() as temporary:
            bed = Path(temporary) / "peaks.bed"
            bed.write_text("chr1\t0\t15\nchr1\t10\t20\n")
            regions = [("chr1", 0, 30, ""), ("chr1", 40, 50, "")]
            header, rows = support_annotations(regions, {"WT": bed, "KO": None})
            result = [dict(zip(header, row)) for row in rows]
            self.assertEqual(result[0]["WT.overlap_bp"], 20)
            self.assertEqual(result[0]["support_groups"], "WT")
            self.assertEqual(result[1]["support_status"], "unassessed")
            self.assertEqual(len(rows), len(regions))
            _, assessed = support_annotations(regions, {"WT": bed})
            self.assertEqual(assessed[1][4], "joint_only")

    def test_unique_input_membership_and_depth_weighted_enrichment(self):
        manifest = {"libraries": [{"runtime_id": name, "role": role, "filtered_bam": name + ".bam"} for name, role in (("a", "chip"), ("b", "chip"), ("i1", "input"), ("i2", "input"))], "associations": [{"chip_id": chip, "input_id": inp} for chip, inp in (("a", "i1"), ("a", "i2"), ("b", "i1"))]}
        controls = ControlSets(manifest)
        self.assertEqual(controls.ids(["a", "b", "a"]), ("i1", "i2"))
        counts = {"a": [20, 10], "b": [10, 5], "i1": [2, 0], "i2": [8, 1]}
        depths = {"a": 100, "b": 100, "i1": 100, "i2": 900}
        rows = enrichment_rows([("chr1", 0, 10, ""), ("chr1", 20, 30, "")], ["a", "b"], controls, counts, depths)
        self.assertEqual(rows[0][9], 20)
        self.assertEqual(rows[1][-1], "low_input")
        self.assertEqual(rows[1][9], "NA")
        self.assertEqual(counts["a"], [20, 10])

    def test_macs_settings_distinguish_discovery_and_tracks(self):
        command = macs_command(["a", "b"], ["i1", "i2"], "out", "joint", broad=True, format_args=["-f", "BEDPE"])
        self.assertEqual(command[command.index("-q") + 1], "0.1")
        self.assertEqual(command[command.index("--broad-cutoff") + 1], "0.1")
        self.assertEqual(command[command.index("-c") + 1:command.index("--broad")], ["i1", "i2"])
        self.assertNotIn("--SPMR", macs_command(["a"], ["i"], "out", "q", tracks=True))


class AnalysisWorkflowTests(fixtures.LibraryFixture):
    run_workflow = fixtures.SnakemakePreprocessingTests.run_workflow

    def rows(self, paired_input=False):
        self.bam("a.filtered.bam", folder="filtered_BAM")
        self.bam("b.filtered.bam", paired=paired_input, folder="filtered_BAM")
        return [fixtures.row("a.filtered.bam", "chip", "a"), fixtures.row("b.filtered.bam", "input", "b")]

    def read_table(self, name):
        with (self.project / "counts" / name).open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def counts_dir(self, config):
        data = yaml.safe_load(config.read_text())
        return self.project / self.config["output_folders"][self.config["countreads_rule_num"] - 1]

    def test_custom_bed_counting_mixed_layout_and_reuse(self):
        rows = self.rows(paired_input=True)
        bed = self.project / "custom.bed"
        bed.write_text("chr1\t100\t250\tregion_a\nchr1\t500\t600\tregion_b\n")
        before = {path: path.read_bytes() for path in (bed, self.project / "filtered_BAM/a.filtered.bam", self.project / "filtered_BAM/b.filtered.bam")}
        config = self.cli_config(rows, "13", extra=("--chip-regions-bed", str(bed)))
        dag = self.run_workflow(config)
        self.assertNotIn("idr_pooled_macs3", dag)
        self.assertNotIn("rule call_peaks:", dag)
        self.run_workflow(config, dry_run=False)
        count_dir = self.counts_dir(config)
        with (count_dir / "project.raw_read_quant.table.txt").open() as handle:
            counts = list(csv.reader(handle, delimiter="\t"))
        self.assertEqual(counts, [["Peak", "a"], ["chr1_100_250", "4"], ["chr1_500_600", "0"]])
        with (count_dir / "chip_input_counts.tsv").open() as handle:
            inputs = list(csv.reader(handle, delimiter="\t"))
        self.assertEqual(inputs, [["Peak", "b"], ["chr1_100_250", "2"], ["chr1_500_600", "0"]])
        again = self.run_workflow(config, dry_run=False)
        self.assertIn("Nothing to be done", again)
        for path, content in before.items():
            self.assertEqual(path.read_bytes(), content)
        bed.write_text("chr1\t500\t600\tchanged\n")
        changed = self.cli_config(rows, "13", extra=("--chip-regions-bed", str(bed)))
        dag = self.run_workflow(changed)
        self.assertIn("chip_testing_regions", dag)
        self.assertIn("count_reads", dag)

    def test_fragment_conversion_counts_pairs_once(self):
        rows = self.rows(paired_input=True)
        output = self.root / "fragments.bedpe"
        self.assertEqual(bam_fragments(self.project / "filtered_BAM/b.filtered.bam", output, True), 2)
        self.assertEqual(output.read_text().splitlines(), ["chr1\t100\t170", "chr1\t180\t250"])

    def large_fixture(self):
        length = 100000
        genome = self.assembly / "test"
        (genome / "fasta/genome.fa.fai").write_text(f"chr1\t{length}\t6\t{length}\t{length + 1}\n")
        (genome / "aux/test_chrom_sizes.2_column").write_text(f"chr1\t{length}\n")
        (genome / "annotation").mkdir()
        (genome / "annotation/genes.gtf").write_text('chr1\ttest\tgene\t19000\t23000\t.\t+\t.\tgene_id "ENSG000001"; gene_name "TEST";\n')
        rows = []
        for sample, role, condition, signal in (("w1", "chip", "WT", 1200), ("w2", "chip", "WT", 1200), ("k1", "chip", "KO", 200), ("k2", "chip", "KO", 200), ("i", "input", "input", 0)):
            path = self.project / "filtered_BAM" / f"{sample}.filtered.bam"
            path.parent.mkdir(exist_ok=True)
            rng = random.Random(sample)
            positions = [rng.randrange(100, length - 300) for _ in range(5000)] + [rng.randrange(20000, 20200) for _ in range(signal)]
            with pysam.AlignmentFile(str(path), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": length}]}) as bam:
                for index, position in enumerate(sorted(positions)):
                    read = pysam.AlignedSegment()
                    read.query_name = f"{sample}_{index}"
                    read.query_sequence = "A" * 30
                    read.flag = 0
                    read.reference_id = 0
                    read.reference_start = position
                    read.mapping_quality = 60
                    read.cigarstring = "30M"
                    read.query_qualities = pysam.qualitystring_to_array("I" * 30)
                    bam.write(read)
            rows.append(fixtures.row(path.name, role, sample, condition=condition))
        return rows

    def test_real_joint_discovery_group_calls_counts_and_metadata(self):
        rows = self.large_fixture()
        config = self.cli_config(rows, "10,13", extra=("--sample-type", "condition", "--narrow-peak-strategy", "macs3", "--chip-peak-opt-mode", "none"))
        self.run_workflow(config, dry_run=False)
        peak_dir = self.project / self.config["output_folders"][9]
        regions = read_regions(peak_dir / "chip_testing/regions.bed", {"chr1": 100000})
        self.assertTrue(any(start <= 20100 < end for _, start, end, _ in regions))
        self.assertTrue((peak_dir / "WT.MACS3.optimized.bed").is_file())
        self.assertTrue((peak_dir / "KO.MACS3.optimized.bed").is_file())
        # Build metadata independently of running the complete DE report.
        metadata = self.counts_dir(config) / "chip_testing_metadata.tsv"
        command = [os.sys.executable, "-m", "snakemake", str(metadata), "--snakefile", str(fixtures.cli.WORKFLOW_ROOT / "Snakefile.smk"), "--config", f"config_file={config}", "--cores", "2", "--scheduler", "greedy", "--rerun-triggers", "mtime"]
        env = dict(os.environ, PYTHONPATH=str(fixtures.cli.PACKAGE_ROOT.parent))
        env["PATH"] = str(Path(os.sys.executable).parent) + os.pathsep + env.get("PATH", "")
        result = subprocess.run(command, cwd=self.project, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout[-12000:])
        with metadata.open() as handle:
            annotations = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(annotations), len(regions))
        self.assertIn("support_groups", annotations[0])
        self.assertIn("assigned_genes", annotations[0])
        before = (peak_dir / "chip_testing/regions.bed").read_bytes()
        timestamp = (peak_dir / "chip_testing/regions.bed").stat().st_mtime_ns
        regrouped = self.cli_config(rows, "10,13", extra=("--sample-type", "cell", "--narrow-peak-strategy", "macs3", "--chip-peak-opt-mode", "none"))
        self.run_workflow(regrouped, dry_run=False)
        self.assertEqual((peak_dir / "chip_testing/regions.bed").read_bytes(), before)
        self.assertEqual((peak_dir / "chip_testing/regions.bed").stat().st_mtime_ns, timestamp)
        changed = self.cli_config(rows, "13", extra=("--chip-joint-peak-q", "0.05"))
        self.assertIn("chip_testing_regions", self.run_workflow(changed))

    def test_real_input_fold_enrichment_track(self):
        rows = self.large_fixture()
        config = self.cli_config(rows, "8-9", extra=("--sample-type", "condition"))
        # Existing coverage jobs require indexes at this selective entry point.
        for path in (self.project / "filtered_BAM").glob("*.bam"):
            pysam.index(str(path))
        self.run_workflow(config, dry_run=False)
        track_dir = self.project / self.config["output_folders"][7] / "input_relative"
        self.assertEqual(len(list(track_dir.glob("*.bw"))), 4)
        provenance = json.loads((track_dir / "w1.fold_enrichment.json").read_text())
        self.assertEqual(provenance["inputs"], ["i"])
        hub_dir = self.project / self.config["output_folders"][8]
        containers = list(hub_dir.glob("*.input_relative.fold_enrichment.hub"))
        self.assertEqual(len(containers), 2)
        self.assertEqual(sum(len(list(path.rglob("*.bw"))) for path in containers), 4)

    def test_real_broad_and_split_idr_control_paths(self):
        rows = self.large_fixture()
        for options in (("--broad-mode", "domain"), ("--narrow-peak-strategy", "idr", "--idr-mode", "encode", "--rerun-selected-steps")):
            with self.subTest(options=options):
                config = self.cli_config(rows, "10,13", extra=("--sample-type", "condition", *options))
                self.run_workflow(config, dry_run=False)
        peak_dir = self.project / self.config["output_folders"][9]
        self.assertTrue((peak_dir / "chip_idr_peak_calling/idr_selected_peaks.tsv").is_file())
        command_records = [line for path in peak_dir.rglob("*.xls") for line in path.read_text(errors="replace").splitlines() if line.startswith("# Command line:")]
        self.assertGreater(len(command_records), 10)
        for command in command_records:
            self.assertIn("-c", command)
            self.assertIn("i.bedpe", command)
            self.assertIn("-g 100000", command)

    def test_custom_bed_de_uses_testing_metadata_without_peak_qc(self):
        rows = self.large_fixture()
        bed = self.project / "atlas.bed"
        bed.write_text("".join(f"chr1\t{start}\t{start+900}\tregion_{start}\n" for start in range(1000, 98000, 1500)))
        config = self.cli_config(rows, "13", extra=("--sample-type", "condition", "--chip-regions-bed", str(bed)))
        self.run_workflow(config, dry_run=False)
        de_config = self.project / "de.yaml"
        de_config.write_text(yaml.safe_dump({"version": 1, "design": {"formula": "~ condition"}, "enrichment": {"enabled": False}, "deseq2": {"fit_type": "local", "sf_type": "poscounts", "lfc_shrink": {"enabled": False}}, "plots": {"sig_heatmap": {"enabled": False}}, "qc": {"variable_gene_heatmap": False}}))
        config = self.cli_config(rows, "14", extra=("--sample-type", "condition", "--chip-regions-bed", str(bed), "--de-config", str(de_config)))
        dag = self.run_workflow(config)
        self.assertNotIn("rule peak_qc:", dag)
        self.assertNotIn("rule call_peaks:", dag)
        self.run_workflow(config, dry_run=False)
        de_dir = self.project / self.config["output_folders"][14]
        self.assertTrue((de_dir / "project.chrom.results.zip").is_file())
        with (de_dir / "peak_metadata.tsv").open() as handle:
            metadata = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(row["support_status"] == "unassessed" for row in metadata))
        results = next(de_dir.rglob("*.diff_peaks.DESeq2.txt"))
        def statistics():
            with results.open() as handle:
                return {row["gene_id"]: tuple(row[column] for column in ("log2FoldChange", "pvalue", "padj")) for row in csv.DictReader(handle, delimiter="\t")}
        before = statistics()
        cache = next(de_dir.rglob("chip_model_cache.rds"))
        cache_time = cache.stat().st_mtime_ns
        peak_dir = self.project / self.config["output_folders"][9]
        (peak_dir / "WT.MACS3.optimized.bed").write_text(bed.read_text())
        (peak_dir / "KO.MACS3.optimized.bed").write_text("")
        self.run_workflow(config, dry_run=False)
        self.assertEqual(statistics(), before)
        self.assertEqual(cache.stat().st_mtime_ns, cache_time)
        with (de_dir / "peak_metadata.tsv").open() as handle:
            metadata = list(csv.DictReader(handle, delimiter="\t"))
        self.assertTrue(all(row["support_groups"] == "WT" for row in metadata))
        # A DE-only restart must not need pruned BAMs or the peak tree.
        import shutil
        shutil.rmtree(self.project / "filtered_BAM")
        shutil.rmtree(self.project / self.config["output_folders"][9])
        config = self.cli_config(rows, "14", extra=("--sample-type", "condition", "--chip-regions-bed", str(bed), "--de-config", str(de_config)))
        self.assertIn("Nothing to be done", self.run_workflow(config, dry_run=False))

    def test_optional_input_track_modes(self):
        rows = [row for row in self.large_fixture() if row["library"] in {"w1", "i"}]
        for path in (self.project / "filtered_BAM").glob("*.bam"):
            pysam.index(str(path))
        for mode in ("log2_ratio", "qpois"):
            with self.subTest(mode=mode):
                config = self.cli_config(rows, "8", extra=("--chip-input-tracks", mode))
                self.run_workflow(config, dry_run=False)
                track_dir = self.project / self.config["output_folders"][7] / "input_relative"
                self.assertEqual(len(list(track_dir.glob(f"*.{mode}.bw"))), 1)

    def test_gene_body_and_bins_have_input_enrichment_without_joint_macs(self):
        rows = self.large_fixture()
        annotation = self.assembly / "test/annotation/genes.gtf"
        annotation.write_text(annotation.read_text() + 'chr1\ttest\tgene\t21000\t25000\t.\t+\t.\tgene_id "ENSG000002"; gene_name "OVERLAP";\n')
        for mode in ("genebody", "diffuse"):
            with self.subTest(mode=mode):
                config = self.cli_config(rows, "10,13", extra=("--broad-mode", mode, "--rerun-selected-steps"))
                self.run_workflow(config, dry_run=False)
                peak_dir = self.project / self.config["output_folders"][9]
                self.assertFalse((peak_dir / "chip_testing/joint_peaks.narrowPeak").exists())
                with (self.counts_dir(config) / "chip_regional_enrichment.tsv").open() as handle:
                    enrichment = list(csv.DictReader(handle, delimiter="\t"))
                self.assertTrue(enrichment)
                self.assertTrue(all(row["input_ids"] == "i" for row in enrichment))

    def test_legacy_shared_input_and_explicit_no_input_counting(self):
        rows = self.rows(paired_input=True)
        chip_row = dict(rows[0])
        del chip_row["role"]
        bed = self.project / "atlas.bed"
        bed.write_text("chr1\t100\t250\n")
        control = self.project / "filtered_BAM/b.filtered.bam"
        for options, input_ids, reliability in ((("-I", str(control)), "legacy_input", "low_input"), ((), "", "no_input")):
            with self.subTest(options=options):
                if not options:
                    control.rename(self.root / "external.input.bam")
                config = self.cli_config([chip_row], "13", extra=("--chip-regions-bed", str(bed), *options, "--rerun-selected-steps"))
                self.run_workflow(config, dry_run=False)
                with (self.counts_dir(config) / "chip_regional_enrichment.tsv").open() as handle:
                    enrichment = next(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(enrichment["input_ids"], input_ids)
                self.assertEqual(enrichment["reliability"], reliability)
