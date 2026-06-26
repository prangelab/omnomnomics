import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from omnomnomics.genomes import (
    build_hisat2_index,
    build_star_index,
    has_complete_hisat2_index,
    star_sa_index_nbases,
)


class GenomeIndexBuilderTests(unittest.TestCase):
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
