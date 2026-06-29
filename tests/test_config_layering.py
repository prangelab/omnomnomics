import tempfile
import unittest
from pathlib import Path

from omnomnomics.cli import DEFAULT_WORKFLOW_CONFIG, load_pipeline_config


class ConfigLayeringTests(unittest.TestCase):
    def test_sparse_site_config_inherits_packaged_site_defaults(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            site_config = Path(tmpdir) / "site.yaml"
            site_config.write_text(
                "\n".join(
                    [
                        "## Omnomnomics pipeline config ##",
                        'genome_assembly_dir: "/scratch/test/genomes/assemblies"',
                        'cellranger_reference_dir: "/scratch/test/genomes/cellranger"',
                    ]
                )
                + "\n"
            )

            config = load_pipeline_config(DEFAULT_WORKFLOW_CONFIG, site_config)

        self.assertEqual(config["partition"], "rome")
        self.assertEqual(config["default_runtime"], 120)
        self.assertEqual(config["genome_assembly_dir"], "/scratch/test/genomes/assemblies")
        self.assertEqual(config["cellranger_reference_dir"], "/scratch/test/genomes/cellranger")


if __name__ == "__main__":
    unittest.main()
