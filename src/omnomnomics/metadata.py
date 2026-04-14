from __future__ import annotations

import csv
import re
from pathlib import Path


class MetadataError(ValueError):
    pass


DERIVED_METADATA_COLUMNS = ("sample_id", "sample_type", "sample_color")
SANITIZE_IDENTIFIER_RE = re.compile(r"[^A-Za-z0-9._-]+")


def read_metadata_table(metadata_path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    metadata_path = Path(metadata_path)
    sample_text = metadata_path.read_text()
    delimiter = "\t"
    try:
        dialect = csv.Sniffer().sniff(sample_text[:4096], delimiters="\t,")
        delimiter = dialect.delimiter
    except csv.Error:
        pass

    with metadata_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = list(reader.fieldnames or [])
        rows = [{key: (value or "").strip() for key, value in row.items()} for row in reader]

    if not fieldnames:
        raise MetadataError(f"Metadata file '{metadata_path}' is empty or malformed.")
    if fieldnames[0] != "filename":
        raise MetadataError(
            f"Metadata file '{metadata_path}' must have 'filename' as its first column."
        )

    reserved = set(DERIVED_METADATA_COLUMNS)
    present_reserved = reserved.intersection(fieldnames)
    if present_reserved:
        raise MetadataError(
            "Metadata file uses reserved derived column names: "
            + ", ".join(sorted(present_reserved))
        )

    if not rows:
        raise MetadataError(f"Metadata file '{metadata_path}' does not contain any sample rows.")

    return fieldnames, rows


def parse_column_selector(selector: str | None, fieldnames: list[str]) -> list[str]:
    if selector is None:
        return []

    tokens = [token.strip() for token in str(selector).split(",") if token.strip()]
    resolved: list[str] = []
    for token in tokens:
        if token.isdigit():
            index = int(token) - 1
            if index < 0 or index >= len(fieldnames):
                raise MetadataError(
                    f"Metadata column index '{token}' is out of range for columns: {', '.join(fieldnames)}"
                )
            column_name = fieldnames[index]
        else:
            if token not in fieldnames:
                raise MetadataError(
                    f"Metadata column '{token}' was not found. Available columns: {', '.join(fieldnames)}"
                )
            column_name = token

        if column_name not in resolved:
            resolved.append(column_name)
    return resolved


def unique_terms(terms: list[str]) -> list[str]:
    ordered_terms: list[str] = []
    for term in terms:
        if term and term not in ordered_terms:
            ordered_terms.append(term)
    return ordered_terms


def sanitize_identifier(value: str) -> str:
    sanitized = SANITIZE_IDENTIFIER_RE.sub("_", value.strip())
    sanitized = re.sub(r"_+", "_", sanitized).strip("._-")
    return sanitized


def join_metadata_values(row: dict[str, str], column_names: list[str]) -> str:
    return "_".join(row[column_name].strip() for column_name in column_names if row[column_name].strip())


def derive_metadata_rows(
    fieldnames: list[str],
    rows: list[dict[str, str]],
    sample_name_selector: str,
    sample_type_selector: str | None = None,
    sample_color_selector: str | None = None,
) -> tuple[list[str], list[dict[str, str]], dict[str, list[str]]]:
    sample_name_columns = parse_column_selector(sample_name_selector, fieldnames)
    if not sample_name_columns:
        raise MetadataError("At least one metadata column must be selected for --sample-name.")

    sample_type_columns = parse_column_selector(sample_type_selector, fieldnames)
    sample_color_columns = parse_column_selector(sample_color_selector, fieldnames)

    derived_rows: list[dict[str, str]] = []
    seen_filenames: set[str] = set()
    seen_sample_ids: set[str] = set()

    for row in rows:
        filename = row["filename"].strip()
        if not filename:
            raise MetadataError("Metadata column 'filename' contains an empty value.")
        if filename in seen_filenames:
            raise MetadataError(f"Metadata column 'filename' contains a duplicate value: {filename}")
        seen_filenames.add(filename)

        sample_id = sanitize_identifier(join_metadata_values(row, sample_name_columns))
        if not sample_id:
            raise MetadataError(
                f"Derived sample_id is empty for metadata row '{filename}'. Check --sample-name columns."
            )
        if sample_id in seen_sample_ids:
            raise MetadataError(
                f"Derived sample_id values are not unique. Duplicate sample_id: {sample_id}"
            )
        seen_sample_ids.add(sample_id)

        sample_type = sanitize_identifier(join_metadata_values(row, sample_type_columns)) if sample_type_columns else ""
        if not sample_type:
            sample_type = "all_samples"

        sample_color = sanitize_identifier(join_metadata_values(row, sample_color_columns)) if sample_color_columns else ""
        if not sample_color:
            sample_color = sample_type if sample_type else "all_samples"

        derived_row = dict(row)
        derived_row["sample_id"] = sample_id
        derived_row["sample_type"] = sample_type
        derived_row["sample_color"] = sample_color
        derived_rows.append(derived_row)

    derived_fieldnames = list(fieldnames) + list(DERIVED_METADATA_COLUMNS)
    selector_map = {
        "sample_name_columns": sample_name_columns,
        "sample_type_columns": sample_type_columns,
        "sample_color_columns": sample_color_columns,
    }
    return derived_fieldnames, derived_rows, selector_map


def choose_default_de_column(
    original_fieldnames: list[str],
    derived_rows: list[dict[str, str]],
) -> str:
    candidate_columns = [column for column in original_fieldnames[1:] if column not in DERIVED_METADATA_COLUMNS]
    for column_name in reversed(candidate_columns):
        values = [row[column_name].strip() for row in derived_rows]
        if any(value == "" for value in values):
            continue
        distinct_values = set(values)
        if len(distinct_values) < 2:
            continue
        if len(distinct_values) == len(values):
            continue
        return column_name
    raise MetadataError(
        "Could not determine a default DE column from metadata. Provide --de-formula or --de-columns."
    )


def resolve_de_metadata(
    original_fieldnames: list[str],
    derived_fieldnames: list[str],
    derived_rows: list[dict[str, str]],
    de_formula: str | None,
    de_columns_selector: str | None,
    de_block_selector: str | None,
    de_interactions: bool,
) -> tuple[list[str], list[dict[str, str]], str, dict[str, str | list[str]]]:
    resolved_formula = (de_formula or "").strip()
    de_columns = parse_column_selector(de_columns_selector, derived_fieldnames)
    de_block = parse_column_selector(de_block_selector, derived_fieldnames)

    if "filename" in de_columns or "filename" in de_block:
        raise MetadataError("The metadata column 'filename' cannot be used in DE design terms.")

    if resolved_formula:
        return derived_fieldnames, derived_rows, resolved_formula, {
            "mode": "explicit_formula",
            "de_columns": de_columns,
            "de_block": unique_terms(de_block),
        }

    if de_columns:
        if len(de_columns) == 1:
            resolved_formula = "~ " + " + ".join(unique_terms([*de_block, de_columns[0]]))
            return derived_fieldnames, derived_rows, resolved_formula, {
                "mode": "assisted_single_column",
                "de_columns": de_columns,
                "de_block": unique_terms(de_block),
            }

        if de_interactions:
            if len(de_columns) != 2:
                raise MetadataError(
                    "--de-interactions currently supports exactly two --de-columns values."
                )
            resolved_formula = "~ " + " + ".join(
                unique_terms([*de_block, *de_columns, f"{de_columns[0]}:{de_columns[1]}"])
            )
            return derived_fieldnames, derived_rows, resolved_formula, {
                "mode": "assisted_interaction",
                "de_columns": de_columns,
                "de_block": unique_terms(de_block),
            }

        updated_rows = []
        for row in derived_rows:
            updated_row = dict(row)
            updated_row["de_group"] = sanitize_identifier(join_metadata_values(row, de_columns))
            if not updated_row["de_group"]:
                raise MetadataError(
                    "Derived DE group is empty for one or more metadata rows. Check --de-columns."
                )
            updated_rows.append(updated_row)
        updated_fieldnames = list(derived_fieldnames)
        if "de_group" not in updated_fieldnames:
            updated_fieldnames.append("de_group")
        resolved_formula = "~ " + " + ".join(unique_terms([*de_block, "de_group"]))
        return updated_fieldnames, updated_rows, resolved_formula, {
            "mode": "assisted_group",
            "de_columns": de_columns,
            "de_block": unique_terms(de_block),
        }

    default_column = choose_default_de_column(original_fieldnames, derived_rows)
    resolved_formula = f"~ {default_column}"
    return derived_fieldnames, derived_rows, resolved_formula, {
        "mode": "default_last_column",
        "de_columns": [default_column],
        "de_block": unique_terms(de_block),
    }


def write_metadata_table(path: str | Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
