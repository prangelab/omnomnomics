import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from omnomnomics.genomes import (
    annotation_gene_intervals,
    build_hisat2_index,
    build_star_index,
    configure_genomepy_reference_plugins,
    create_gene_bed_files,
    has_complete_hisat2_index,
    star_sa_index_nbases,
    validate_installed_intervals,
)


class GenomeIndexBuilderTests(unittest.TestCase):
    def test_transcript_exon_gtf_produces_gene_beds_without_changing_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            gtf = root / "genes.gtf"
            original = ('chr1\tgenomepy\ttranscript\t10\t40\t.\t-\t.\tgene_id "a"; gene_name "A"; transcript_id "t1";\n'
                        'chr1\tgenomepy\texon\t20\t30\t.\t-\t.\tgene_id "a"; transcript_id "t1";\n'
                        'chr1\tgenomepy\ttranscript\t15\t60\t.\t-\t.\tgene_id "a"; gene_name "A"; transcript_id "t2";\n')
            gtf.write_text(original)
            fai, sizes = root / "genome.fa.fai", root / "chrom.sizes"
            fai.write_text("chr1\t100\t0\t100\t101\n")
            sizes.write_text("chr1\t100\n")
            validate_installed_intervals(fai, sizes, gtf)
            create_gene_bed_files(gtf, "mm39", root)
            self.assertEqual((root / "mm39_genes.sorted.bed").read_text(), "chr1\t9\t60\ta\t0\t-\tA\n")
            self.assertEqual((root / "mm39_refseq_genes.sorted.bed").read_text(), "chr1\t9\t60\ta\t0\t-\n")
            self.assertEqual(gtf.read_text(), original)

    def test_explicit_genes_take_precedence_and_inference_keeps_contigs_and_strands(self):
        with tempfile.TemporaryDirectory() as temporary:
            gtf = Path(temporary) / "genes.gtf"
            gtf.write_text('chr1\tsource\tgene\t10\t50\t.\t+\t.\tgene_id "a";\n'
                           'chr1\tsource\texon\t20\t30\t.\t+\t.\tgene_id "a";\n'
                           'chr1\tsource\texon\t70\t80\t.\t-\t.\tgene_id "a";\n'
                           'chr2\tsource\texon\t5\t15\t.\t+\t.\tgene_id "a";\n')
            self.assertEqual(annotation_gene_intervals(gtf), [
                ("chr1", 9, 50, "a", "0", "+", "a"),
                ("chr1", 69, 80, "a", "0", "-", "a"),
                ("chr2", 4, 15, "a", "0", "+", "a"),
            ])

    def test_exon_only_annotation_still_rejects_mismatches_and_invalid_intervals(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fai, sizes, gtf = (root / name for name in ("genome.fa.fai", "chrom.sizes", "genes.gtf"))
            fai.write_text("chr1\t100\t0\t100\t101\n")
            sizes.write_text("chr1\t100\n")
            for record, error in [
                ('1\tsource\texon\t10\t20\t.\t+\t.\tgene_id "a";\n', "Primary gene-annotation"),
                ('chr1\tsource\texon\t10\t110\t.\t+\t.\tgene_id "a";\n', "exceeds reference"),
                ('chr1\tsource\texon\t0\t20\t.\t+\t.\tgene_id "a";\n', "Invalid gene annotation"),
                ('chr1\tsource\texon\t10\t20\t.\t+\t.\ttranscript_id "t";\n', "no gene-associated"),
            ]:
                with self.subTest(record=record):
                    gtf.write_text(record)
                    with self.assertRaisesRegex(ValueError, error):
                        validate_installed_intervals(fai, sizes, gtf)

    def test_genomepy_reference_plugins_disable_aligner_indexes(self):
        class FakeGenomepy:
            def __init__(self):
                self.calls = []

            def manage_plugins(self, action, plugins):
                self.calls.append((action, plugins))

        fake_genomepy = FakeGenomepy()

        configure_genomepy_reference_plugins(fake_genomepy, cache_blacklist=True)

        self.assertEqual(
            fake_genomepy.calls,
            [("disable", ["hisat2", "star"]), ("enable", ["blacklist"])],
        )

    def test_hisat2_completion_accepts_large_indexes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            index_dir = Path(tmpdir)
            for i in range(1, 9):
                (index_dir / f"hg38.{i}.ht2l").write_text("")

            self.assertTrue(has_complete_hisat2_index(index_dir, "hg38"))

    def test_star_sa_index_nbases_is_capped_for_large_genomes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fai_path = Path(tmpdir) / "genome.fa.fai"
            fai_path.write_text("chr1\t3200000000\t0\t80\t81\n")

            self.assertEqual(star_sa_index_nbases(fai_path), 14)

    @patch("omnomnomics.genomes.shutil.which", return_value="/usr/bin/STAR")
    @patch("omnomnomics.genomes.subprocess.run")
    def test_star_builder_uses_normalized_fasta_and_gtf(self, run_mock, _which_mock):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "genome.fa"
            gtf = root / "genes.gtf"
            fai = root / "genome.fa.fai"
            fasta.write_text(">chr1\nACGT\n")
            gtf.write_text("")
            fai.write_text("chr1\t4\t0\t4\t5\n")

            build_star_index(root / "star", fasta, gtf, fai, 16)

        command = run_mock.call_args.args[0]
        self.assertIn("--runMode", command)
        self.assertIn("genomeGenerate", command)
        self.assertIn("--genomeFastaFiles", command)
        self.assertIn(str(fasta), command)
        self.assertIn("--sjdbGTFfile", command)
        self.assertIn(str(gtf), command)

    @patch("omnomnomics.genomes.shutil.which", return_value="/usr/bin/hisat2-build")
    @patch("omnomnomics.genomes.subprocess.run")
    def test_hisat2_builder_uses_normalized_fasta(self, run_mock, _which_mock):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "genome.fa"
            fasta.write_text(">chr1\nACGT\n")

            build_hisat2_index(root / "hisat2", fasta, "hg38", 16)

        self.assertEqual(
            run_mock.call_args.args[0],
            ["hisat2-build", "-p", "16", str(fasta), str(root / "hisat2" / "hg38")],
        )


if __name__ == "__main__":
    unittest.main()
