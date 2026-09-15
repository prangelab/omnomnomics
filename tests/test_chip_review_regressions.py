import copy
import csv
import json
import shutil
import unittest
from pathlib import Path

import yaml

from omnomnomics.chip_analysis import ControlSets, region_specs_match
import test_chip_analysis as analysis
import test_chip_inputs as fixtures


class RegionSpecTests(unittest.TestCase):
    def spec(self, root):
        return {"controls": {"samples": ["a"], "libraries": [{"id": "a", "bam": f"{root}/filtered_BAM/a.filtered.bam", "layout": True}], "legacy_input": "/external/input.bam", "associations": []}, "bed": f"{root}/atlas.bed", "bed_sha256": "same", "q": 0.1}

    def test_project_relocation_and_old_archive_migration(self):
        old, new = self.spec("/old/project"), self.spec("/new/project")
        self.assertTrue(region_specs_match(old, new, "/new/project", "/old/project"))
        self.assertTrue(region_specs_match(old, new, "/new/project"))
        for mutation in (lambda s: s.update(q=0.05), lambda s: s.update(bed_sha256="different"), lambda s: s["controls"].update(legacy_input="/another/input.bam"), lambda s: s["controls"]["libraries"][0].update(layout=False), lambda s: s["controls"]["libraries"][0].update(bam="/new/project/filtered_BAM/other.bam"), lambda s: s["controls"].update(associations=[{"chip_id": "a", "input_id": "different"}])):
            changed = copy.deepcopy(new)
            mutation(changed)
            with self.subTest(changed=changed):
                self.assertFalse(region_specs_match(old, changed, "/new/project", "/old/project"))
                self.assertFalse(region_specs_match(old, changed, "/new/project"))

    def test_legacy_control_identity_avoids_all_library_names(self):
        names = ["legacy_input", "omnomnomics_legacy_input_1"]
        controls = ControlSets(legacy_input="/input.bam", sample_ids=names)
        self.assertEqual(controls.ids(names), ("omnomnomics_legacy_input_2",))
        self.assertEqual(controls.paths(names), ["/input.bam"])


class ReviewWorkflowTests(fixtures.LibraryFixture):
    run_workflow = fixtures.SnakemakePreprocessingTests.run_workflow
    rows = analysis.AnalysisWorkflowTests.rows
    counts_dir = analysis.AnalysisWorkflowTests.counts_dir
    def review_fixture(self):
        rows = self.rows(paired_input=True)
        self.bam("c.filtered.bam", folder="filtered_BAM")
        rows.append(fixtures.row("c.filtered.bam", "chip", "c"))
        for row in rows:
            row["condition"] = "KO" if row["library"] == "c" else "WT"
        bed = self.project / "atlas.bed"
        bed.write_text("chr1\t100\t250\n")
        annotation = self.assembly / "test/annotation"
        annotation.mkdir()
        (annotation / "genes.gtf").write_text('chr1\ttest\tgene\t100\t300\t.\t+\t.\tgene_id "G1"; gene_name "TEST";\n')
        return rows, bed

    def raw_samples(self, config):
        with (self.counts_dir(config) / "project.raw_read_quant.table.txt").open() as handle:
            return next(csv.reader(handle, delimiter="\t"))[1:]

    def test_chip_named_legacy_input_keeps_its_own_counts(self):
        rows, bed = self.review_fixture()
        (self.project / "filtered_BAM/c.filtered.bam").unlink()
        chip = self.project / "filtered_BAM/a.filtered.bam"
        chip.rename(chip.with_name("legacy_input.filtered.bam"))
        row = dict(rows[0], filename="legacy_input.filtered.bam", library="legacy_input")
        del row["role"]
        config = self.cli_config([row], "13", extra=("--chip-regions-bed", str(bed), "-I", str(self.project / "filtered_BAM/b.filtered.bam")))
        self.run_workflow(config, dry_run=False)
        directory = self.counts_dir(config)
        with (directory / "project.raw_read_quant.table.txt").open() as handle:
            self.assertEqual(list(csv.reader(handle, delimiter="\t")), [["Peak", "legacy_input"], ["chr1_100_250", "4"]])
        with (directory / "chip_input_counts.tsv").open() as handle:
            self.assertEqual(list(csv.reader(handle, delimiter="\t")), [["Peak", "omnomnomics_legacy_input_1"], ["chr1_100_250", "2"]])

    def test_spp_drop_to_none_restores_samples_and_de_rejects_old_exclusions(self):
        rows, bed = self.review_fixture()
        drop = self.project / "filtered_BAM/peak_qc/spp_qc/dropped_samples.tsv"
        drop.parent.mkdir(parents=True)
        drop.write_text("sample_id\na\n")
        options = ("--chip-regions-bed", str(bed), "--de-formula", "~ condition")
        config = self.cli_config(rows, "13", extra=(*options, "--spp-gate", "drop"))
        self.run_workflow(config, dry_run=False)
        self.assertEqual(self.raw_samples(config), ["c"])
        config = self.cli_config(rows, "14", extra=(*options, "--spp-gate", "none"))
        with self.assertRaisesRegex(AssertionError, "Include public step 13"):
            self.run_workflow(config)
        for gate in ("none", "warn", "strict"):
            with self.subTest(gate=gate):
                config = self.cli_config(rows, "13", extra=(*options, "--spp-gate", gate))
                self.run_workflow(config, dry_run=False)
                self.assertEqual(self.raw_samples(config), ["a", "c"])
                provenance = json.loads((self.counts_dir(config) / "chip_count_provenance.json").read_text())
                self.assertEqual(provenance["spp_excluded_samples"], [])
        self.assertTrue(drop.is_file())
        config = self.cli_config(rows, "13", extra=(*options, "--spp-gate", "drop"))
        self.run_workflow(config, dry_run=False)
        self.assertEqual(self.raw_samples(config), ["c"])
        drop.unlink()
        self.run_workflow(config, dry_run=False)
        self.assertEqual(self.raw_samples(config), ["a", "c"])

    def test_de_only_after_relocation_reuses_counts_and_checks_changed_bed(self):
        rows, bed = self.review_fixture()
        options = ("--chip-regions-bed", str(bed), "--de-formula", "~ condition")
        config = self.cli_config(rows, "13", extra=options)
        self.run_workflow(config, dry_run=False)
        original = self.project
        for old_archive in (False, True):
            with self.subTest(old_archive=old_archive):
                self.project = self.root / f"relocated_{old_archive}/project"
                shutil.copytree(original, self.project)
                shutil.rmtree(self.project / "filtered_BAM")
                provenance = self.counts_dir(config) / "chip_count_provenance.json"
                if old_archive:
                    data = json.loads(provenance.read_text())
                    del data["project_root"]
                    provenance.write_text(json.dumps(data))
                bed = self.project / "atlas.bed"
                config = self.cli_config(rows, "14", extra=("--chip-regions-bed", str(bed), "--de-formula", "~ condition"))
                dag = self.run_workflow(config)
                self.assertNotIn("rule count_reads:", dag)
                self.assertNotIn("rule chip_testing_regions:", dag)
                self.assertIn("rule call_DE_chrom:", dag)
                bed.write_text("chr1\t110\t250\n")
                config = self.cli_config(rows, "14", extra=("--chip-regions-bed", str(bed), "--de-formula", "~ condition"))
                with self.assertRaisesRegex(AssertionError, "Include public step 13"):
                    self.run_workflow(config)


class DefaultDesignTests(fixtures.LibraryFixture):
    def test_default_de_design_ignores_inputs_and_role_column(self):
        rows = []
        for sample, role, condition in (("w1", "chip", "WT"), ("w2", "chip", "WT"), ("k1", "chip", "KO"), ("k2", "chip", "KO"), ("i", "input", "WT")):
            self.bam(f"{sample}.filtered.bam", folder="filtered_BAM")
            row = fixtures.row(f"{sample}.filtered.bam", role, sample, condition=condition)
            row["role"] = row.pop("role")
            rows.append(row)
        config = self.cli_config(rows, "14")
        data = yaml.safe_load(config.read_text())
        self.assertEqual(data["RESOLVED_DE_FORMULA"], "~ condition")
        with Path(data["DERIVED_METADATA_FILE"]).open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 5)
