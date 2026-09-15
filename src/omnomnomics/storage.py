"""Preserve the earliest available project inputs during output cleanup."""

from __future__ import annotations

import fcntl
import json
import os
from pathlib import Path
import shutil
import tempfile

from omnomnomics.metadata import normalize_metadata_filename


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
SOURCE_MANIFEST = "source_protection.json"


def _folder(config, key, rule_key):
    value = config[key][config[rule_key] - 1]
    return value[0] if isinstance(value, list) else value


def discover_source_files(experiment_dir, config):
    """Find the first available sequence stage for each filename identity.

    Unrelated or renamed identities are retained conservatively. Discovery does
    not assume that a FASTQ belonging to one sample can replace another's BAM.
    """
    root = Path(experiment_dir)
    stages = [
        (_folder(config, "input_folders", "trim_rule_num"), FASTQ_SUFFIXES),
        (_folder(config, "output_folders", "trim_rule_num"), FASTQ_SUFFIXES),
        (_folder(config, "output_folders", "merge_rule_num"), (".bam",)),
        (_folder(config, "output_folders", "touchup_rule_num"), (".bam",)),
        (_folder(config, "output_folders", "wig_rule_num"), (".bw",)),
    ]
    seen = set()
    sources = set()
    for folder, suffixes in stages:
        folder_path = root / folder
        if not folder_path.is_dir():
            continue
        stage_files = [
            path for path in folder_path.iterdir()
            if path.is_file() and path.name.endswith(suffixes)
        ]
        for path in stage_files:
            key = normalize_metadata_filename(path.name)
            if key not in seen:
                sources.add(path)
        seen.update(normalize_metadata_filename(path.name) for path in stage_files)

    # Table/feature-only entry has no per-library sequence source to preserve.
    if not seen:
        for folders, suffixes in zip(config["input_folders"], config["input_file_types"]):
            folders = folders if isinstance(folders, list) else [folders]
            suffixes = suffixes if isinstance(suffixes, list) else [suffixes]
            candidates = set()
            for folder, suffix in zip(folders, suffixes):
                folder_path = root / folder
                if folder_path.is_dir():
                    candidates.update(
                        path for path in folder_path.iterdir()
                        if path.is_file() and path.name.endswith(suffix)
                    )
            if candidates:
                sources.update(candidates)
                break
    return sources


def capture_source_protection(experiment_dir, config, additional_sources=(), generated_paths=()):
    """Persist source identities before any outputs can be removed.

    Existing protections are never downgraded by a later run. Records use project
    relative names for portability and retain external/symlink targets as well.
    """
    root = Path(experiment_dir).resolve()
    manifest = root / "run_configs" / SOURCE_MANIFEST
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        records = {}
        if manifest.exists():
            data = json.loads(manifest.read_text())
            if data.get("version") != 1 or not isinstance(data.get("sources"), list):
                raise ValueError(f"Invalid source protection manifest: {manifest}")
            for record in data["sources"]:
                if not isinstance(record, dict) or not isinstance(record.get("path"), str):
                    raise ValueError(f"Invalid source protection record in {manifest}")
                records[record["path"]] = record

        sources = discover_source_files(root, config)
        generated = {Path(path).absolute() for path in generated_paths}
        sources.difference_update(generated)
        sources.update(Path(path).absolute() for path in additional_sources if path and path != "NA")
        for source in sorted(sources):
            paths = [source]
            if source.suffix == ".bam":
                paths.extend((Path(f"{source}.bai"), source.with_suffix(".bai"), Path(f"{source}.csi")))
            for path in paths:
                path = path.absolute()
                try:
                    name = str(path.relative_to(root))
                except ValueError:
                    name = str(path)
                record = records.setdefault(name, {"path": name, "targets": []})
                resolved = str(path.resolve())
                if resolved != str(path) and resolved not in record["targets"]:
                    record["targets"].append(resolved)

        data = {"version": 1, "sources": [records[name] for name in sorted(records)]}
        content = json.dumps(data, indent=2, sort_keys=True) + "\n"
        if not manifest.exists() or manifest.read_text() != content:
            fd, temporary_path = tempfile.mkstemp(prefix=".source_protection.", dir=manifest.parent)
            try:
                with os.fdopen(fd, "w") as handle:
                    handle.write(content)
                    handle.flush()
                    os.fsync(handle.fileno())
                os.replace(temporary_path, manifest)
            finally:
                if os.path.exists(temporary_path):
                    os.unlink(temporary_path)

    protected = set()
    for record in data["sources"]:
        path = Path(record["path"])
        protected.add(str(path if path.is_absolute() else root / path))
        protected.update(record["targets"])
    return str(manifest), sorted(protected)


def source_protects_path(path, protected_sources):
    """Protect sources and any directory that contains them, including aliases."""
    if protected_sources is None:
        # Old run configurations have no trustworthy pre-run source snapshot.
        return True
    candidate = Path(path).absolute()
    candidates = {candidate, candidate.resolve()}
    for source in protected_sources:
        source = Path(source).absolute()
        for identity in {source, source.resolve()}:
            for candidate in candidates:
                if candidate == identity or candidate in identity.parents or identity in candidate.parents:
                    return True
    return False


def remove_generated_path(path, protected_sources):
    """Remove an unprotected output; retain mixed source/output directories."""
    path = Path(path)
    if source_protects_path(path, protected_sources):
        return False
    if path.is_symlink() or path.is_file():
        path.unlink(missing_ok=True)
    elif path.is_dir():
        shutil.rmtree(path)
    return True
