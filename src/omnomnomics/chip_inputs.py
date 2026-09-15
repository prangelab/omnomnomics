"""Resolve declared ChIP libraries and their input-control associations."""

from __future__ import annotations

from collections import defaultdict

from omnomnomics.metadata import MetadataError, parse_column_selector


def resolve_chip_inputs(fieldnames, rows, match_selector=None, legacy_input="NA"):
    if "role" not in fieldnames:
        raise MetadataError("ChIP input matching requires a metadata 'role' column (chip or input).")
    if legacy_input and legacy_input != "NA":
        raise MetadataError("Metadata roles cannot be combined with -I/--input. Use one control declaration method.")
    columns = parse_column_selector(match_selector, fieldnames)
    forbidden = {"role", "filename", "sample_id", "filename_key"}.intersection(columns)
    if forbidden:
        raise MetadataError("Input matching must use shared metadata fields, not " + ", ".join(sorted(forbidden)))
    grouped = defaultdict(list)
    for row in rows:
        if row.get("role") not in {"chip", "input"}:
            raise MetadataError(f"Metadata row '{row['filename']}' must have role chip or input.")
        grouped[row["sample_id"]].append(row)

    libraries = []
    for sample_id, units in sorted(grouped.items()):
        for column in ["role", *columns]:
            values = {row.get(column, "").strip() for row in units}
            if len(values) != 1 or (column in columns and "" in values):
                raise MetadataError(f"Library '{sample_id}' has missing or conflicting '{column}' metadata across its technical units.")
        libraries.append({"sample_id": sample_id, "role": units[0]["role"], "metadata": units[0], "rows": units})

    chips = [library for library in libraries if library["role"] == "chip"]
    inputs = [library for library in libraries if library["role"] == "input"]
    if chips and inputs and not columns:
        raise MetadataError("Select matching metadata columns with --input-match, for example --input-match cell_type,condition.")
    if chips and not inputs and columns:
        raise MetadataError("--input-match was supplied but no input libraries were declared.")

    candidates = defaultdict(list)
    for library in inputs:
        key = tuple(library["metadata"].get(column, "").strip() for column in columns)
        candidates[key].append(library["sample_id"])
    associations = []
    used = set()
    for library in chips:
        if not inputs:
            continue
        key = tuple(library["metadata"][column].strip() for column in columns)
        matches = candidates.get(key, [])
        if not matches:
            values = ", ".join(f"{column}={value}" for column, value in zip(columns, key))
            raise MetadataError(f"No input matches ChIP library '{library['sample_id']}' ({values}).")
        for input_id in sorted(matches):
            associations.append({"chip_id": library["sample_id"], "input_id": input_id, "match_values": dict(zip(columns, key))})
            used.add(input_id)
    return {
        "libraries": libraries,
        "match_columns": columns,
        "associations": associations,
        "unassigned_inputs": sorted(library["sample_id"] for library in inputs if library["sample_id"] not in used),
    }


def validate_preparation_scope(steps, create_homer_tagdirs=False):
    if create_homer_tagdirs:
        raise MetadataError(
            "Optional HOMER tag directories for metadata-declared ChIP inputs are not enabled yet."
        )
