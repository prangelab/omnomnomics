import ast
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


RULE = Path(__file__).parents[1] / "src/omnomnomics/workflow/rules/16.analyze_peaks_de.smk"


class PostDERerunTests(unittest.TestCase):
    def workflow(self, directory, assay, selected, force):
        source = RULE.read_text().split("\nrule analyze_peaks_de:", 1)[0]
        function = next(node for node in ast.parse(source).body
                        if isinstance(node, ast.FunctionDef) and node.name == "analyze_peaks_de_input")
        inputs = ast.get_source_segment(source, function)
        project = directory / "project"
        project.mkdir()
        config = {
            "THETYPE": assay, "EXPERIMENT_DIR": str(project),
            "THEMODE": selected, "RERUN_SELECTED_STEPS": bool(force),
        }
        (project / "config.json").write_text(json.dumps(config))
        workflow = """import os
from pathlib import Path
experiment_dir = str(Path.cwd())
master_config = dict(dechrom_rule_num=15, touchup_rule_num=5, wig_rule_num=8,
                     output_folders=['filtered_BAM'] * 7 + ['BigWigs'] * 7 + ['DE_calling'] * 2)
samples2 = ['sample']
""" + inputs + """

rule all:
    input: 'DE_calling/extra_16.tmp'

rule call_DE_chrom:
    input: 'DE_calling/raw_counts.tsv'
    output: 'DE_calling/project.chrom.results.zip', 'DE_calling/peak_metadata.tsv'
    shell: 'touch {output}; echo DE >> executed.txt'

rule analyze_peaks_de:
    input: analyze_peaks_de_input
    output: 'DE_calling/extra_16.tmp'
    shell: 'touch {output}; echo post-DE >> executed.txt'
"""
        (project / "Snakefile").write_text(workflow)
        bam = "sample.filtered.bam" if assay == "CHIP" else "sample.sorted.dups_marked.filtered.bam"
        files = ["DE_calling/raw_counts.tsv", f"filtered_BAM/{bam}", "BigWigs/sample.bw",
                 "DE_calling/project.chrom.results.zip", "DE_calling/peak_metadata.tsv", "DE_calling/extra_16.tmp"]
        for index, name in enumerate(files):
            path = project / name
            path.parent.mkdir(exist_ok=True)
            path.write_text("retained output\n")
            os.utime(path, (1600000000 + index, 1600000000 + index))
        command = [sys.executable, "-m", "snakemake", "--configfile", "config.json",
                   "--cores", "1", "--scheduler", "greedy", "--rerun-triggers", "mtime", "--nocolor"]
        if force:
            command.extend(["--forcerun", *force])
        return project, command

    def run_workflow(self, project, command, dry_run=True):
        result = subprocess.run(command + (["--dry-run"] if dry_run else []), cwd=project,
                                text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout)
        return result.stdout

    def test_forced_de_and_post_de_both_run_with_existing_outputs(self):
        for assay in ("CHIP", "ATAC"):
            with self.subTest(assay=assay), tempfile.TemporaryDirectory() as temporary:
                project, command = self.workflow(Path(temporary), assay, [15, 16],
                                                 ["call_DE_chrom", "analyze_peaks_de"])
                retained = {name: (project / name).read_bytes()
                            for name in ("DE_calling/raw_counts.tsv", "BigWigs/sample.bw")}
                plan = self.run_workflow(project, command)
                self.assertIn("rule call_DE_chrom:", plan)
                self.assertIn("rule analyze_peaks_de:", plan)
                self.assertFalse((project / "executed.txt").exists())
                self.run_workflow(project, command, dry_run=False)
                self.assertEqual((project / "executed.txt").read_text().splitlines(), ["DE", "post-DE"])
                for name, content in retained.items():
                    self.assertEqual((project / name).read_bytes(), content)

    def test_completed_outputs_are_reused(self):
        for assay in ("CHIP", "ATAC"):
            with self.subTest(assay=assay), tempfile.TemporaryDirectory() as temporary:
                project, command = self.workflow(Path(temporary), assay, [15, 16], [])
                self.assertIn("Nothing to be done", self.run_workflow(project, command))

    def test_forcing_only_post_de_reuses_existing_de(self):
        for assay in ("CHIP", "ATAC"):
            with self.subTest(assay=assay), tempfile.TemporaryDirectory() as temporary:
                project, command = self.workflow(Path(temporary), assay, [16], ["analyze_peaks_de"])
                plan = self.run_workflow(project, command)
                self.assertIn("rule analyze_peaks_de:", plan)
                self.assertNotIn("rule call_DE_chrom:", plan)


if __name__ == "__main__":
    unittest.main()
