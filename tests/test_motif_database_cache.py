import tempfile
import unittest
from pathlib import Path

from omnomnomics.genomes import (
    DEFAULT_MEME_MOTIF_DATABASE_NAME,
    find_cached_meme_motif_database,
    looks_like_meme_database,
    meme_motif_database_path,
    resolve_meme_motif_database,
)


MEME_HEADER = """MEME version 5

ALPHABET= ACGT

strands: + -
"""


class MotifDatabaseCacheTests(unittest.TestCase):
    def test_default_cache_path_is_under_reference_root(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = meme_motif_database_path(tmpdir)

        self.assertEqual(path.name, DEFAULT_MEME_MOTIF_DATABASE_NAME)
        self.assertEqual(path.parent.name, "motif_databases")

    def test_cached_database_requires_meme_format(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_path = meme_motif_database_path(tmpdir)
            cache_path.parent.mkdir(parents=True)
            cache_path.write_text("not a motif db\n")

            self.assertFalse(looks_like_meme_database(cache_path))
            self.assertIsNone(find_cached_meme_motif_database(tmpdir))

    def test_resolve_uses_cached_database_without_download(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cache_path = meme_motif_database_path(tmpdir)
            cache_path.parent.mkdir(parents=True)
            cache_path.write_text(MEME_HEADER)

            resolved = resolve_meme_motif_database(tmpdir, download=False)

        self.assertEqual(Path(resolved), cache_path)


if __name__ == "__main__":
    unittest.main()
