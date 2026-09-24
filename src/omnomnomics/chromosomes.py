"""Reference-native chromosome names for external interval files."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path


ALIAS_VERSION = 1
_PRIMARY = re.compile(r"^(?:chr)?(?:[0-9]+|X|Y)$")
_MITOCHONDRIAL = ("M", "MT", "chrM", "chrMT")


class ChromosomeNameError(ValueError):
    pass


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_reference_fai(path):
    reference = {}
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ChromosomeNameError(f"Invalid FASTA index row at {path}:{number}")
            name, length = fields[:2]
            if name in reference or not name:
                raise ChromosomeNameError(f"Duplicate or empty reference contig at {path}:{number}: {name!r}")
            try:
                reference[name] = int(length)
            except ValueError as exc:
                raise ChromosomeNameError(f"Invalid reference length at {path}:{number}") from exc
            if reference[name] <= 0:
                raise ChromosomeNameError(f"Invalid reference length at {path}:{number}")
    if not reference:
        raise ChromosomeNameError(f"Empty FASTA index: {path}")
    return reference


def is_primary_chromosome(name, reference):
    return name in reference and bool(_PRIMARY.fullmatch(name))


def is_primary_name(name):
    return bool(_PRIMARY.fullmatch(name))


def resolve_chromosome(name, reference):
    if name in reference:
        return name
    if name in _MITOCHONDRIAL:
        candidates = [candidate for candidate in _MITOCHONDRIAL if candidate in reference]
    elif _PRIMARY.fullmatch(name):
        alternate = name[3:] if name.startswith("chr") else f"chr{name}"
        candidates = [alternate] if alternate in reference else []
    else:
        candidates = []
    if len(candidates) > 1:
        raise ChromosomeNameError(f"Ambiguous reference contig alias for {name}: {', '.join(candidates)}")
    return candidates[0] if candidates else None


def filter_primary_bed(source, destination, reference):
    source_n = 0
    primary_n = 0
    primary_like_n = 0
    with open(source, encoding="utf-8") as source_handle, open(destination, "w", encoding="utf-8") as target_handle:
        for line in source_handle:
            if not line.strip():
                continue
            source_n += 1
            chrom = line.split("\t", 1)[0]
            primary_like_n += is_primary_name(chrom)
            if is_primary_chromosome(chrom, reference):
                target_handle.write(line)
                primary_n += 1
    if primary_like_n and primary_n == 0:
        raise ChromosomeNameError(
            f"Peak chromosome names do not match the selected reference: {source}; "
            f"source={source_n}, primary_like={primary_like_n}, retained=0"
        )
    return source_n, primary_n


def _open_text(path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def normalize_external_bed(source, destination, reference_fai, *, assembly=None, unknown_policy="error"):
    """Write an external BED in reference names and preserve its source."""
    if unknown_policy not in {"error", "skip"}:
        raise ValueError(f"Invalid unknown-contig policy: {unknown_policy}")
    source, destination, reference_fai = map(Path, (source, destination, reference_fai))
    if source.resolve() == destination.resolve():
        raise ChromosomeNameError("The reference-native BED must not overwrite its source")
    reference = read_reference_fai(reference_fai)
    source_sha256 = file_sha256(source)
    reference_sha256 = file_sha256(reference_fai)
    provenance_path = Path(f"{destination}.source.json")
    if destination.is_file() and provenance_path.is_file():
        try:
            prior = json.loads(provenance_path.read_text(encoding="utf-8"))
            if (
                prior.get("source_sha256") == source_sha256
                and prior.get("reference_fai_sha256") == reference_sha256
                and prior.get("alias_version") == ALIAS_VERSION
                and prior.get("assembly") == assembly
                and prior.get("unknown_policy") == unknown_policy
                and prior.get("output_sha256") == file_sha256(destination)
            ):
                return destination, prior
        except (OSError, ValueError):
            pass

    destination.parent.mkdir(parents=True, exist_ok=True)
    handle_id, temporary_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    mappings = {}
    reverse = {}
    source_records = 0
    written_records = 0
    skipped_records = 0
    try:
        with os.fdopen(handle_id, "w", encoding="utf-8") as output, _open_text(source) as input_handle:
            for number, line in enumerate(input_handle, 1):
                if not line.strip() or line.startswith(("#", "track ", "browser ")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ChromosomeNameError(f"Invalid BED record at {source}:{number}")
                source_records += 1
                chrom = fields[0]
                target = resolve_chromosome(chrom, reference)
                if target is None:
                    if unknown_policy == "skip":
                        skipped_records += 1
                        continue
                    raise ChromosomeNameError(f"Unknown BED contig at {source}:{number}: {chrom}")
                previous = reverse.setdefault(target, chrom)
                if previous != chrom:
                    raise ChromosomeNameError(
                        f"BED contigs {previous} and {chrom} both map to reference contig {target}"
                    )
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError as exc:
                    raise ChromosomeNameError(f"Invalid BED coordinates at {source}:{number}") from exc
                if not 0 <= start < end <= reference[target]:
                    raise ChromosomeNameError(
                        f"BED coordinate outside the selected reference at {source}:{number}: {chrom}:{start}-{end}"
                    )
                mappings[chrom] = target
                fields[0] = target
                output.write("\t".join(fields) + "\n")
                written_records += 1
        if not written_records:
            raise ChromosomeNameError(f"No usable BED records in {source}")
        os.replace(temporary_name, destination)
    finally:
        if os.path.exists(temporary_name):
            os.unlink(temporary_name)

    provenance = {
        "alias_version": ALIAS_VERSION,
        "assembly": assembly,
        "source": str(source.resolve()),
        "source_sha256": source_sha256,
        "reference_fai": str(reference_fai.resolve()),
        "reference_fai_sha256": reference_sha256,
        "unknown_policy": unknown_policy,
        "chromosome_mapping": dict(sorted(mappings.items())),
        "source_records": source_records,
        "written_records": written_records,
        "skipped_records": skipped_records,
        "output_sha256": file_sha256(destination),
    }
    provenance_text = json.dumps(provenance, sort_keys=True, indent=2) + "\n"
    provenance_id, temporary_provenance = tempfile.mkstemp(prefix=f".{provenance_path.name}.", dir=destination.parent)
    try:
        with os.fdopen(provenance_id, "w", encoding="utf-8") as handle:
            handle.write(provenance_text)
        os.replace(temporary_provenance, provenance_path)
    finally:
        if os.path.exists(temporary_provenance):
            os.unlink(temporary_provenance)
    return destination, provenance
