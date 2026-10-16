import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "scripts/release"))
from build_environment_locks import checked_records
from inventory_environments import inventory
from verify_environment import verify


class ReleaseLockTests(unittest.TestCase):
    def record(self):
        return {"name": "example", "version": "1.0", "build": "0", "build_string": "0",
                "subdir": "noarch", "url": "https://conda.anaconda.org/conda-forge/noarch/example-1.0-0.conda",
                "md5": "0" * 32, "sha256": "1" * 64, "depends": None}

    def test_lock_rejects_unapproved_sources_and_incomplete_plans(self):
        record = self.record()
        plan = {"success": True, "actions": {"LINK": [record]}}
        self.assertEqual(checked_records(plan)[0]["depends"], [])
        for updates in [
            {"url": "https://example.org/example.conda"},
            {"url": "https://token@conda.anaconda.org/conda-forge/noarch/example.conda"},
            {"url": record["url"] + "?token=secret"},
            {"sha256": ""}, {"depends": ["missing >=1"]}, {"subdir": "osx-arm64"},
        ]:
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                checked_records({"success": True, "actions": {"LINK": [record | updates]}})
        with self.assertRaises(ValueError):
            checked_records({"success": False})

    def test_post_install_check_rejects_changed_archives_and_extra_packages(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            (prefix / "conda-meta").mkdir()
            record = self.record()
            metadata = prefix / "conda-meta/example.json"
            metadata.write_text(json.dumps(record))
            verify(prefix, [record])
            metadata.write_text(json.dumps(record | {"md5": "2" * 32}))
            with self.assertRaises(ValueError):
                verify(prefix, [record])
            metadata.write_text(json.dumps(record))
            (prefix / "conda-meta/extra.json").write_text(json.dumps(record | {"name": "extra"}))
            with self.assertRaises(ValueError):
                verify(prefix, [record])

    def test_inventory_resolves_metadata_directory_symlinks(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix = root / "environment"
            directory = prefix / "lib/python3.11/site-packages/example-1.0.dist-info"
            directory.mkdir(parents=True)
            (directory / "METADATA").write_text("Metadata-Version: 2.1\nName: example\nVersion: 1.0\n")
            (prefix / "lib/python3.1").symlink_to("python3.11")
            (prefix / "conda-meta").mkdir()
            record = self.record() | {"depends": [], "files": [str((directory / "METADATA").relative_to(prefix))]}
            (prefix / "conda-meta/example.json").write_text(json.dumps(record))
            inventory(prefix, root / "inventory")
            self.assertEqual(json.loads((root / "inventory/pip-only.json").read_text()), [])


if __name__ == "__main__":
    unittest.main()
