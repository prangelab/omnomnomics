"""Resolve source files and preprocessing targets before a ChIP run starts."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile

from omnomnomics.metadata import (
    FASTQ_EXTENSIONS, MetadataError, fastq_unit_identity, normalize_metadata_filename,
)
from omnomnomics.storage import source_protects_path


def write_stable_text(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text() == content:
        return
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def load_manifest(path):
    if not path or str(path) == "NA" or not Path(path).exists():
        return None
    data = json.loads(Path(path).read_text())
    if data.get("version") != 1 or not isinstance(data.get("libraries"), list):
        raise MetadataError(f"Invalid library manifest: {path}")
    return data


def reference_lengths(assembly_dir, genome):
    root = Path(assembly_dir) / genome
    candidates = [root / "fasta/genome.fa.fai", root / "aux" / f"{genome}_chrom_sizes.2_column"]
    for path in candidates:
        if path.is_file():
            with path.open() as handle:
                return {fields[0]: int(fields[1]) for line in handle if len(fields := line.split()) >= 2}
    raise MetadataError(f"BAM validation needs genome.fa.fai or chromosome sizes for {genome} under {root}.")


def inspect_bam(path, reference, max_records=10000):
    import pysam

    try:
        pysam.quickcheck(str(path))
        with pysam.AlignmentFile(str(path), "rb") as bam:
            mismatches = [name for name, length in zip(bam.references, bam.lengths) if reference.get(name) != length]
            if not bam.references or mismatches:
                raise MetadataError(f"BAM reference does not match the selected genome: {path}; contigs: {', '.join(mismatches[:5])}")
            layouts = set()
            count = 0
            for read in bam.fetch(until_eof=True):
                count += 1
                if not read.is_secondary and not read.is_supplementary:
                    layouts.add(bool(read.is_paired))
                if count >= max_records or len(layouts) > 1:
                    break
            if len(layouts) != 1:
                raise MetadataError(f"Cannot determine one SE/PE layout from primary alignments in {path}.")
            return layouts.pop()
    except (OSError, ValueError, pysam.SamtoolsError) as exc:
        raise MetadataError(f"Invalid BAM {path}: {exc}") from exc


def _fastq_inventory(folder):
    inventory = {}
    if not folder.is_dir():
        return inventory
    for path in sorted(folder.iterdir()):
        if not path.is_file() or not path.name.endswith(FASTQ_EXTENSIONS):
            continue
        unit, read = fastq_unit_identity(path.name)
        record = inventory.setdefault(unit, {})
        if read in record:
            raise MetadataError(f"Multiple FASTQs resolve to {unit}/{read} in {folder}.")
        record[read] = str(path)
    return inventory


def _read_pair(unit, reads, declared_layout=""):
    if set(reads) == {"R1"} and not declared_layout:
        raise MetadataError(f"Only R1 is present for '{unit}'; supply R2 or declare read_layout=SE.")
    if declared_layout == "PE" and set(reads) != {"R1", "R2"}:
        raise MetadataError(f"Missing read pair for PE library unit '{unit}'.")
    if declared_layout == "SE" and set(reads) == {"R1", "R2"}:
        raise MetadataError(f"SE read_layout conflicts with paired FASTQs for '{unit}'.")
    if set(reads) == {"R1", "R2"}:
        return [reads["R1"], reads["R2"]], True
    if set(reads) in ({"SE"}, {"R1"}):
        return [next(iter(reads.values()))], False
    raise MetadataError(f"Incomplete or conflicting FASTQ reads for '{unit}': {', '.join(reads)}")


def build_library_manifest(root, resolved, steps, genome, assembly_dir, config, previous=None):
    root = Path(root).resolve()
    steps = sorted(steps)
    folders = {name: root / config["output_folders"][config[rule] - 1] for name, rule in (
        ("trim", "trim_rule_num"), ("bam", "merge_rule_num"), ("filtered", "touchup_rule_num"), ("qc", "qc_rule_num"),
    )}
    raw_folder = root / config["input_folders"][config["trim_rule_num"] - 1]
    raw = _fastq_inventory(raw_folder) if min(steps) < 4 else {}
    trimmed = _fastq_inventory(folders["trim"]) if min(steps) < 4 else {}
    previous_libs = {library["sample_id"]: library for library in (previous or {}).get("libraries", [])}
    known_outputs = {path for library in previous_libs.values() for paths in library.get("targets", {}).values() for path in paths}
    known_outputs.update(unit["bam"] for library in previous_libs.values() for unit in library["units"])
    known_outputs.update((previous or {}).get("generated_paths", []))
    bam_reference = None

    def bam_layout(path):
        nonlocal bam_reference
        if bam_reference is None:
            bam_reference = reference_lengths(assembly_dir, genome)
        return inspect_bam(path, bam_reference)

    libraries = []
    claimed_units = {}
    for declaration in resolved["libraries"]:
        sample = declaration["sample_id"]
        rows = declaration["rows"]
        declared_layouts = {row.get("read_layout", "").strip().upper() for row in rows}
        if len(declared_layouts) != 1 or not declared_layouts <= {"", "SE", "PE"}:
            raise MetadataError(f"Library '{sample}' has invalid or conflicting read_layout; use SE or PE.")
        declared_layout = declared_layouts.pop()
        canonical_bam = str(folders["bam"] / f"{sample}.bam")
        canonical_filtered = str(folders["filtered"] / f"{sample}.filtered.bam")
        library = {key: declaration[key] for key in ("sample_id", "role", "metadata")}
        library.update({"runtime_id": sample, "units": [], "bam": canonical_bam, "filtered_bam": canonical_filtered, "source_paths": [], "targets": {}})
        prior = previous_libs.get(sample, {})
        if steps == [15]:
            library.update({"paired": prior.get("paired", declared_layout == "PE"), "entry_stage": 15,
                            "original_sources": prior.get("original_sources", []), "source_paths": []})
            library["units"] = []
            libraries.append(library)
            continue
        original_sources = list(prior.get("original_sources", []))
        filtered_candidates = [Path(canonical_filtered)]
        if len(rows) == 1:
            filtered_candidates.append(folders["filtered"] / f"{rows[0]['filename_key']}.filtered.bam")
            if str(rows[0]["filename"]).endswith(".filtered.bam"):
                filtered_candidates.insert(0, folders["filtered"] / Path(rows[0]["filename"]).name)
        filtered_source = next((path for path in filtered_candidates if path.is_file()), None)
        if min(steps) >= 6:
            if filtered_source is None:
                raise MetadataError(f"Library '{sample}' needs an existing filtered BAM for entry at step {min(steps)}; include step 5 to prepare it.")
            library.update({"paired": bam_layout(filtered_source), "entry_stage": 6, "filtered_source": str(filtered_source), "source_paths": [str(filtered_source)]})
        elif min(steps) >= 4 and Path(canonical_bam).is_file() and all(not str(row["filename"]).endswith(".bam") or Path(row["filename"]).name == Path(canonical_bam).name for row in rows):
            library["units"] = [{"unit_id": sample, "filename_key": rows[0]["filename_key"], "source_type": "bam", "source_paths": [canonical_bam], "bam": canonical_bam, "paired": bam_layout(canonical_bam), "trimmed": [], "fastqs": []}]
            library["entry_stage"] = 4
        else:
            for row in rows:
                key = row["filename_key"]
                raw_units = {unit: reads for unit, reads in raw.items() if normalize_metadata_filename(unit) == key}
                trim_units = {unit: reads for unit, reads in trimmed.items() if normalize_metadata_filename(unit) == key}
                bam_files = sorted(path for path in folders["bam"].glob("*.bam") if normalize_metadata_filename(path.name) == key)
                explicit_bam = str(row["filename"]).endswith(".bam")
                explicit_fastq = str(row["filename"]).endswith(FASTQ_EXTENSIONS)
                if explicit_bam:
                    exact = folders["bam"] / Path(row["filename"]).name
                    bam_files = [exact] if exact.is_file() else []
                if min(steps) >= 4 or explicit_bam or (not raw_units and not trim_units):
                    if not bam_files:
                        raise MetadataError(f"No BAM found for '{row['filename']}' in {folders['bam']} at this entry stage.")
                    for path in bam_files:
                        library["units"].append({"unit_id": path.stem, "filename_key": key, "source_type": "bam", "source_paths": [str(path)], "bam": str(path), "paired": bam_layout(path), "trimmed": [], "fastqs": []})
                    continue
                if raw_units and bam_files and not explicit_fastq and any(str(path) not in known_outputs for path in bam_files):
                    raise MetadataError(f"Ambiguous FASTQ/BAM sources for '{row['filename']}'. Specify a FASTQ or BAM filename in metadata.")
                source_units = raw_units if min(steps) == 1 and raw_units else trim_units
                if not source_units:
                    raise MetadataError(f"Library '{sample}' needs trimmed FASTQs at step {min(steps)}; include step 1 to trim raw reads.")
                for unit, reads in sorted(source_units.items()):
                    paths, paired = _read_pair(unit, reads, declared_layout)
                    source_type = "fastq" if source_units is raw_units else "trimmed_fastq"
                    trimmed_paths = [str(folders["trim"] / f"{unit}{suffix}.trimmed.fastq.gz") for suffix in (["_R1", "_R2"] if paired else [""])]
                    if 1 not in steps or source_type == "trimmed_fastq":
                        if unit not in trim_units:
                            raise MetadataError(f"No trimmed FASTQs for '{unit}'; include step 1.")
                        trimmed_paths, trimmed_paired = _read_pair(unit, trim_units[unit], declared_layout)
                        if paired != trimmed_paired:
                            raise MetadataError(f"Raw and trimmed read layouts differ for '{unit}'.")
                    library["units"].append({"unit_id": unit, "filename_key": key, "source_type": source_type, "source_paths": paths, "fastqs": paths if source_type == "fastq" else [], "trimmed": trimmed_paths, "bam": str(folders["bam"] / f"{unit}.bam"), "paired": paired})
        units = library["units"]
        if units:
            layouts = {unit["paired"] for unit in units}
            if len(layouts) != 1:
                raise MetadataError(f"Technical units of library '{sample}' have different SE/PE layouts.")
            library["paired"] = layouts.pop()
            library["entry_stage"] = min({"fastq": 1, "trimmed_fastq": 3, "bam": 4}[unit["source_type"]] for unit in units)
            library["source_paths"] = [path for unit in units for path in unit["source_paths"]]
            if 4 not in steps and len(units) == 1:
                library["bam"] = units[0]["bam"]
            if 4 not in steps and len(units) > 1 and any(step >= 5 for step in steps) and not Path(canonical_bam).is_file():
                raise MetadataError(f"Library '{sample}' needs step 4 to merge its technical units.")
        if declared_layout and library["paired"] != (declared_layout == "PE"):
            raise MetadataError(f"Library '{sample}' read_layout disagrees with its BAM/FASTQ layout.")
        for source in library["source_paths"]:
            if not os.access(source, os.R_OK):
                raise MetadataError(f"Unreadable library source: {source}")
            if source not in known_outputs and source not in original_sources:
                original_sources.append(source)
        library["original_sources"] = sorted(set(original_sources))
        index_path = f"{library['filtered_bam']}.bai"
        library["protected_index"] = (
            library["filtered_bam"] in original_sources and Path(index_path).is_file()
        )
        for unit in units:
            name = unit["unit_id"]
            if name in claimed_units:
                raise MetadataError(f"Input unit '{name}' is assigned more than once ({claimed_units[name]}, {sample}).")
            claimed_units[name] = sample
        libraries.append(library)

    library_names = {lib["runtime_id"] for lib in libraries}
    for unit, owner in claimed_units.items():
        if unit in library_names and unit != owner:
            raise MetadataError(f"Input unit '{unit}' collides with another library's sample name.")
    manifest = {"version": 1, "genome": genome, "steps": steps, "folders": {key: str(path) for key, path in folders.items()}, "libraries": libraries,
                "match_columns": resolved["match_columns"], "associations": resolved["associations"], "unassigned_inputs": resolved["unassigned_inputs"]}
    plan_preprocessing(manifest, config)
    source_paths = {str(Path(path).resolve()) for lib in libraries for path in lib["original_sources"]}
    claimed_outputs = {}
    for lib in libraries:
        for stage, outputs in lib["targets"].items():
            for path in outputs:
                if str(Path(path).resolve()) in source_paths:
                    raise MetadataError(f"Step {stage} would overwrite a source input: {path}. Choose distinct sample names or start after that stage.")
                identity = str(Path(path).resolve())
                if identity in claimed_outputs:
                    raise MetadataError(f"Conflicting output paths for {claimed_outputs[identity]} and {lib['sample_id']}: {path}")
                claimed_outputs[identity] = lib["sample_id"]
        imported = lib.get("filtered_source")
        if imported and imported != lib["filtered_bam"] and str(Path(lib["filtered_bam"]).resolve()) in source_paths:
            raise MetadataError(f"Filtered BAM import would overwrite a source input: {lib['filtered_bam']}")
    generated_paths = set((previous or {}).get("generated_paths", []))
    generated_paths.update(path for lib in libraries for outputs in lib["targets"].values() for path in outputs)
    generated_paths.update(lib["filtered_bam"] for lib in libraries if lib.get("filtered_source") and lib["filtered_source"] != lib["filtered_bam"])
    manifest["generated_paths"] = sorted(path for path in generated_paths if str(Path(path).resolve()) not in source_paths and not any(lib.get("protected_index") and path == f"{lib['filtered_bam']}.bai" for lib in libraries))
    return manifest


def validate_protected_targets(manifest, protected_sources):
    for library in manifest["libraries"]:
        imported = library.get("filtered_source")
        destination = library["filtered_bam"]
        if imported and imported != destination and Path(destination).exists() and source_protects_path(destination, protected_sources):
            raise MetadataError(f"Filtered BAM import would overwrite a protected source: {destination}")
        for outputs in library["targets"].values():
            for path in outputs:
                if library.get("protected_index") and path == f"{library['filtered_bam']}.bai":
                    continue
                if Path(path).exists() and source_protects_path(path, protected_sources):
                    raise MetadataError(f"Planned output is a protected source: {path}. Use a distinct library name or enter after the stage that would overwrite it.")


def plan_preprocessing(manifest, config):
    steps = manifest["steps"]
    folders = manifest["folders"]
    work = {str(step): {} for step in steps if step <= 7}
    for lib in manifest["libraries"]:
        sample = lib["runtime_id"]
        targets = lib["targets"]

        def add(step, key, outputs):
            work[str(step)][key] = outputs
            targets.setdefault(str(step), []).extend(outputs)

        for unit in lib["units"]:
            name = unit["unit_id"]
            if unit["source_type"] != "bam":
                if 1 in steps and unit["source_type"] == "fastq":
                    add(1, name, [*unit["trimmed"], str(Path(folders["trim"]) / f"{name}.trim_metrics.tsv")])
                if 2 in steps:
                    stems = [f"{name}{read}.trimmed_fastqc" for read in (["_R1", "_R2"] if unit["paired"] else [""])]
                    add(2, name, [str(Path(folders["qc"]) / f"{stem}.{ext}") for stem in stems for ext in ("html", "zip")])
                if 3 in steps:
                    label = {"hisat2": "HISAT2", "star": "STAR", "star_te": "STAR_TE"}[config["map_tool"]]
                    add(3, name, [unit["bam"], str(Path(folders["bam"]) / f"{name}.{label}_stats.txt"), str(Path(folders["bam"]) / f"{name}.extra_3.tmp")])
            if any(step in steps for step in (4, 5)) and not Path(unit["bam"]).is_file() and name not in work.get("3", {}):
                raise MetadataError(f"Missing BAM for '{name}'; include alignment step 3.")
        if lib["entry_stage"] < 6:
            if 4 in steps:
                merged = str(Path(folders["bam"]) / f"{sample}.bam")
                lib["merge_output"] = len(lib["units"]) != 1 or lib["units"][0]["bam"] != merged
                add(4, sample, ([merged] if lib["merge_output"] else []) + [str(Path(folders["bam"]) / f"{sample}.extra_4.tmp")])
            if 5 in steps:
                add(5, sample, [lib["filtered_bam"], str(Path(folders["filtered"]) / f"{sample}.extra_5.tmp")])
        if any(step in steps for step in (6, 7)) and 5 not in steps and not Path(lib.get("filtered_source", lib["filtered_bam"])).is_file():
            raise MetadataError(f"Missing filtered BAM for '{sample}'; include step 5.")
        if 6 in steps:
            add(6, sample, [f"{lib['filtered_bam']}.bai"])
        if 7 in steps:
            add(7, sample, [f"{lib['filtered_bam']}.{suffix}" for suffix in ("stats.txt", "qc_summary.pdf", "qc_summary.svg")])
    manifest["work"] = work


def save_library_manifest(root, run_date, manifest):
    root = Path(root) / "run_configs"
    content = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    fingerprint = hashlib.sha256(content.encode()).hexdigest()[:16]
    snapshot = root / f"libraries.{fingerprint}.json"
    write_stable_text(snapshot, content)
    write_stable_text(root / "library_manifest.json", content)
    write_stable_text(root / f"omnomnomics.run.{run_date}.libraries.json", content)
    buffer = io.StringIO()
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(["chip_id", "input_id", *manifest["match_columns"]])
    for association in manifest["associations"]:
        writer.writerow([association["chip_id"], association["input_id"], *[association["match_values"][column] for column in manifest["match_columns"]]])
    association_path = root / f"inputs.{fingerprint}.tsv"
    write_stable_text(association_path, buffer.getvalue())
    write_stable_text(root / "input_associations.tsv", buffer.getvalue())
    write_stable_text(root / f"omnomnomics.run.{run_date}.input_associations.tsv", buffer.getvalue())
    return str(snapshot), str(association_path)
