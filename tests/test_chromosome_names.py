import json
from pathlib import Path
import tempfile
import unittest

from omnomnomics.chromosomes import ChromosomeNameError, filter_primary_bed, normalize_external_bed, read_reference_fai, resolve_chromosome
from omnomnomics.chip_analysis import read_regions, write_regions
from omnomnomics.genomes import resolve_blacklist_bed, validate_installed_intervals
from omnomnomics.metadata import MetadataError


class ChromosomeNameTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_primary_peak_filter_preserves_reference_names(self):
        for name in ("1", "chr1"):
            with self.subTest(name=name):
                source, filtered = self.root / "peaks.bed", self.root / "filtered.bed"
                source.write_text(f"{name}\t10\t20\nalt1\t30\t40\n")
                self.assertEqual(filter_primary_bed(source, filtered, {name: 100, "alt1": 100}), (2, 1))
                self.assertEqual(filtered.read_text(), f"{name}\t10\t20\n")

    def test_primary_peak_filter_rejects_name_mismatch(self):
        source = self.root / "peaks.bed"
        source.write_text("1\t10\t20\n")
        with self.assertRaisesRegex(ChromosomeNameError, "do not match"):
            filter_primary_bed(source, self.root / "filtered.bed", {"chr1": 100})

    def test_blacklist_aliases_preserve_source_and_refresh_output(self):
        for reference_name, source_name in (("1", "chr1"), ("chr1", "1")):
            with self.subTest(reference_name=reference_name):
                assembly = self.root / reference_name / "test"
                (assembly / "fasta").mkdir(parents=True)
                (assembly / "aux").mkdir()
                (assembly / "fasta/genome.fa.fai").write_text(f"{reference_name}\t100\t0\t100\t101\n")
                source = assembly / "aux/test.blacklist.bed"
                source.write_text(f"{source_name}\t10\t20\tblacklist\nother\t5\t15\n")
                result = resolve_blacklist_bed("test", assembly.parent, download_if_missing=False)
                self.assertEqual(result.read_text(), f"{reference_name}\t10\t20\tblacklist\n")
                self.assertEqual(source.read_text(), f"{source_name}\t10\t20\tblacklist\nother\t5\t15\n")
                provenance = json.loads(Path(f"{result}.source.json").read_text())
                self.assertEqual(provenance["assembly"], "test")
                self.assertEqual(provenance["chromosome_mapping"], {source_name: reference_name})
                self.assertEqual(provenance["skipped_records"], 1)
                before = result.stat().st_mtime_ns
                resolve_blacklist_bed("test", assembly.parent, download_if_missing=False)
                self.assertEqual(result.stat().st_mtime_ns, before)
                source.write_text(f"{source_name}\t30\t40\n")
                resolve_blacklist_bed("test", assembly.parent, download_if_missing=False)
                self.assertEqual(result.read_text(), f"{reference_name}\t30\t40\n")
                prior_reference_sha = json.loads(Path(f"{result}.source.json").read_text())["reference_fai_sha256"]
                (assembly / "fasta/genome.fa.fai").write_text(f"{reference_name}\t101\t0\t101\t102\n")
                resolve_blacklist_bed("test", assembly.parent, download_if_missing=False)
                current_reference_sha = json.loads(Path(f"{result}.source.json").read_text())["reference_fai_sha256"]
                self.assertNotEqual(current_reference_sha, prior_reference_sha)

    def test_blacklist_rejects_collisions_and_zero_usable_intervals(self):
        fai, source = self.root / "genome.fa.fai", self.root / "source.bed"
        fai.write_text("1\t100\t0\t100\t101\n")
        source.write_text("1\t1\t10\nchr1\t20\t30\n")
        with self.assertRaisesRegex(ChromosomeNameError, "both map"):
            normalize_external_bed(source, self.root / "normalized.bed", fai, unknown_policy="skip")
        source.write_text("unknown\t1\t10\n")
        with self.assertRaisesRegex(ChromosomeNameError, "No usable"):
            normalize_external_bed(source, self.root / "normalized.bed", fai, unknown_policy="skip")

    def test_aliases_do_not_guess_alt_contigs_or_ambiguous_mitochondria(self):
        self.assertEqual(resolve_chromosome("chr1", {"1": 100}), "1")
        self.assertIsNone(resolve_chromosome("chr1_alt", {"1_alt": 100}))
        with self.assertRaisesRegex(ChromosomeNameError, "Ambiguous"):
            resolve_chromosome("chrM", {"M": 100, "MT": 100})

    def test_custom_bed_uses_reference_names_and_records_mapping(self):
        source = self.root / "source.bed"
        source.write_text("chr1\t10\t20\tmy_region\n")
        regions, mapping = read_regions(source, {"1": 100}, allow_aliases=True, return_mapping=True)
        bed, manifest = self.root / "regions.bed", self.root / "regions.tsv"
        write_regions(regions, bed, manifest, source, "GRCh38.p14", mapping, "fai-checksum")
        self.assertEqual(source.read_text(), "chr1\t10\t20\tmy_region\n")
        self.assertEqual(bed.read_text(), "1\t10\t20\t1_10_20\n")
        self.assertIn("my_region", manifest.read_text())
        provenance = json.loads(Path(f"{manifest}.source.json").read_text())
        self.assertEqual(provenance["chromosome_mapping"], {"chr1": "1"})
        self.assertEqual(provenance["reference_fai_sha256"], "fai-checksum")

    def test_custom_bed_rejects_collision_and_out_of_bounds(self):
        source = self.root / "source.bed"
        source.write_text("1\t1\t10\nchr1\t20\t30\n")
        with self.assertRaisesRegex(MetadataError, "both map"):
            read_regions(source, {"1": 100}, allow_aliases=True)
        source.write_text("chr1\t90\t110\n")
        with self.assertRaisesRegex(MetadataError, "outside the selected reference"):
            read_regions(source, {"1": 100}, allow_aliases=True)

    def test_install_checks_sizes_and_annotation(self):
        fai, sizes, gtf = (self.root / name for name in ("genome.fa.fai", "chrom.sizes", "genes.gtf"))
        fai.write_text("1\t100\t0\t100\t101\n")
        sizes.write_text("1\t100\n")
        gtf.write_text('1\tsource\tgene\t10\t20\t.\t+\t.\tgene_id "a";\n')
        self.assertEqual(read_reference_fai(fai), {"1": 100})
        validate_installed_intervals(fai, sizes, gtf)
        sizes.write_text("chr1\t100\n")
        with self.assertRaisesRegex(ValueError, "do not match"):
            validate_installed_intervals(fai, sizes, gtf)
        sizes.write_text("1\t100\n")
        gtf.write_text('chr1\tsource\tgene\t10\t20\t.\t+\t.\tgene_id "a";\n')
        with self.assertRaisesRegex(ValueError, "Primary gene-annotation contigs"):
            validate_installed_intervals(fai, sizes, gtf)
