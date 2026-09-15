"""Read-only access to a resolved library manifest inside Snakemake."""

from __future__ import annotations

from pathlib import Path
import re


class LibraryRuntime:
    def __init__(self, manifest):
        self.manifest = manifest
        self.libraries = {library["runtime_id"]: library for library in manifest["libraries"]}
        self.units = {unit["unit_id"]: unit for library in self.libraries.values() for unit in library["units"]}
        self.unit_libraries = {unit["unit_id"]: library for library in self.libraries.values() for unit in library["units"]}
        self.work = manifest["work"]

    def library(self, sample):
        return self.libraries.get(sample) or self.unit_libraries[sample]

    def paired(self, sample):
        return self.library(sample)["paired"]

    def role(self, sample):
        return self.library(sample)["role"]

    def pattern(self, step, paired=None):
        samples = sorted(self.work.get(str(step), {}))
        if step == 6:
            samples = [sample for sample in samples if not self.library(sample).get("protected_index")]
        if paired is not None:
            samples = [sample for sample in samples if self.paired(sample) == paired]
        return "|".join(re.escape(sample) for sample in samples) or r"$.^"

    def targets(self):
        return sorted({path for jobs in self.work.values() for outputs in jobs.values() for path in outputs})

    def trimmed(self, sample, read):
        paths = self.units[sample]["trimmed"]
        return paths[0] if read == 1 else (paths[1] if len(paths) == 2 else [])

    def raw(self, sample, read):
        paths = self.units[sample]["fastqs"]
        return paths[0] if read == 1 else (paths[1] if len(paths) == 2 else [])

    def input_units(self, sample):
        return [unit["unit_id"] for unit in self.library(sample)["units"]]

    def mapped_units(self, sample):
        return [unit for unit in self.input_units(sample) if unit in self.work.get("3", {})]

    def bam(self, sample):
        return self.library(sample)["bam"]

    def merge_marker(self, sample):
        if sample not in self.work.get("4", {}):
            return []
        return str(Path(self.manifest["folders"]["bam"]) / f"{sample}.extra_4.tmp")

    def filtered_imports(self):
        return {
            sample: library["filtered_source"]
            for sample, library in self.libraries.items()
            if library.get("filtered_source") and library["filtered_source"] != library["filtered_bam"]
        }
