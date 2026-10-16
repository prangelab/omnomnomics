#!/usr/bin/env python3
"""Export installed Conda archive provenance without importing installed packages."""

import argparse
import csv
import importlib.metadata
import json
from pathlib import Path
from urllib.parse import urlsplit, urlunsplit


def public_archive_url(record):
    url = record.get("url", "")
    parsed = urlsplit(url)
    allowed = (
        parsed.hostname == "conda.anaconda.org"
        and parsed.path.split("/")[1] in {"conda-forge", "bioconda"}
    ) or (
        parsed.hostname == "repo.anaconda.com"
        and parsed.path.startswith(("/pkgs/main/", "/pkgs/r/"))
    )
    if (not allowed or parsed.scheme != "https" or parsed.username or
            parsed.password or parsed.query or parsed.fragment):
        raise ValueError(f"Unapproved archive source for {record['name']}")
    if not record.get("md5"):
        raise ValueError(f"Missing archive MD5 for {record['name']}")
    return urlunsplit(parsed)


def inventory(prefix, output):
    records = []
    owned_metadata = set()
    for path in sorted((prefix / "conda-meta").glob("*.json")):
        record = json.loads(path.read_text())
        url = public_archive_url(record)
        records.append({key: record.get(key, "") for key in
                        ["name", "version", "build", "subdir", "md5", "sha256", "depends"]} | {"url": url})
        owned_metadata.update(record.get("files", []))
    if not records:
        raise ValueError(f"No Conda records in {prefix}")
    output.mkdir(parents=True, exist_ok=False)
    (output / "conda-records.json").write_text(json.dumps(records, indent=2) + "\n")
    (output / "explicit.txt").write_text(
        "# platform: linux-64\n@EXPLICIT\n" +
        "".join(f"{r['url']}#{r['md5']}\n" for r in records))
    with (output / "packages.tsv").open("w") as handle:
        fields = ["name", "version", "build", "subdir", "url", "md5", "sha256"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    pip_only = []
    seen_metadata = set()
    for site in prefix.glob("lib/python*/site-packages"):
        for distribution in importlib.metadata.distributions(path=[str(site)]):
            # Conda RECORD entries can be relative to the prefix instead of site-packages.
            directory = Path(distribution._path)
            if directory.resolve() in seen_metadata:
                continue
            seen_metadata.add(directory.resolve())
            metadata = [directory / "METADATA", directory / "PKG-INFO", directory]
            if any(str(path.resolve().relative_to(prefix.resolve())) in owned_metadata for path in metadata):
                continue
            direct = distribution.read_text("direct_url.json")
            source = json.loads(direct) if direct else None
            if source and source.get("url", "").startswith("file:"):
                source = {"kind": "local", "editable": source.get("dir_info", {}).get("editable", False)}
            elif source:
                parsed = urlsplit(source.get("url", ""))
                source = {"kind": "remote", "host": parsed.hostname,
                          "has_credentials": bool(parsed.username or parsed.password or parsed.query)}
            pip_only.append({"name": distribution.metadata["Name"], "version": distribution.version,
                             "source": source, "installer": distribution.read_text("INSTALLER")})
    (output / "pip-only.json").write_text(json.dumps(pip_only, indent=2) + "\n")
    print(json.dumps({"environment": prefix.name, "conda_packages": len(records), "pip_only": pip_only}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--environments-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    for name in ["omnomnomics", "omnomnomics-idr", "omnomnomics-spp"]:
        inventory(args.environments_root / name, args.output / name)


if __name__ == "__main__":
    main()
