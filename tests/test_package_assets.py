import unittest
from importlib.resources import files


class PackageAssetTests(unittest.TestCase):
    def test_required_workflow_assets_are_packaged(self):
        package_root = files("omnomnomics")
        required = (
            "workflow/Snakefile.smk",
            "workflow/config/workflow.yaml",
            "workflow/rules/16.analyze_peaks_de.smk",
            "workflow/templates/de_core.R.tmpl",
            "workflow/R/shiny_app/app.R",
            "workflow/R/shiny_app/www/explorer.css",
            "workflow/R/shiny_app/www/explorer.js",
            "workflow/R/shiny_app/modules/mod_settings.R",
            "workflow/R/de_shared/explorer_colors.R",
            "workflow/slurm_profile/config.yaml",
        )

        missing = [relative for relative in required if not package_root.joinpath(relative).is_file()]
        self.assertEqual(missing, [])
