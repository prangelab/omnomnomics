#!/usr/bin/env python3
"""Write single-platform release locks from reviewed Micromamba dry-run plans."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import urlsplit

from inventory_environments import public_archive_url


NAMES = {"main": "omnomnomics", "idr": "omnomnomics-idr",
         "spp": "omnomnomics-spp", "explorer": "omnomnomics-explorer"}


def checked_records(plan):
    if not plan.get("success"):
        raise ValueError("The environment solve did not succeed")
    records = []
    for item in plan["actions"]["LINK"]:
        record = {key: item.get(key) for key in ["name", "version", "subdir", "url", "md5", "sha256", "depends"]}
        record["depends"] = record["depends"] or []
        record["build"] = item.get("build_string", item.get("build"))
        public_archive_url(record)
        if record["subdir"] not in {"linux-64", "noarch"}:
            raise ValueError(f"Unexpected platform: {record['name']}")
        for key, length in [("md5", 32), ("sha256", 64)]:
            if not re.fullmatch(rf"[0-9a-f]{{{length}}}", record.get(key) or ""):
                raise ValueError(f"Missing or malformed {key}: {record['name']}")
        records.append(record)
    if not records or len({r["name"] for r in records}) != len(records):
        raise ValueError("Empty plan or duplicate package names")
    names = {r["name"] for r in records}
    for record in records:
        for dependency in record["depends"]:
            name = dependency.split()[0]
            if not name.startswith("__") and name not in names:
                raise ValueError(f"Missing dependency {name} of {record['name']}")
    # Link hooks require their dependencies to be installed first.
    return records


def yaml_lock(records, source, source_hash):
    lines = ["# Generated from the checked solver plan; includes every transitive package.",
             "version: 1", "metadata:", "  content_hash:", f"    linux-64: {source_hash}",
             "  channels:", "    - url: https://conda.anaconda.org/conda-forge",
             "      used_env_vars: []", "    - url: https://conda.anaconda.org/bioconda",
             "      used_env_vars: []", "  platforms: [linux-64]", f"  sources: [{json.dumps(source)}]", "package:"]
    for record in records:
        lines.extend([f"  - name: {json.dumps(record['name'])}",
                      f"    version: {json.dumps(record['version'])}", "    manager: conda",
                      "    platform: linux-64", "    dependencies:"])
        if not record["depends"]:
            lines[-1] = "    dependencies: {}"
        for dependency in record["depends"]:
            name, _, specification = dependency.partition(" ")
            lines.append(f"      {json.dumps(name)}: {json.dumps(specification or '*')}")
        lines.extend([f"    url: {json.dumps(record['url'])}", "    hash:",
                      f"      md5: {json.dumps(record['md5'])}", f"      sha256: {json.dumps(record['sha256'])}",
                      "    category: main", "    optional: false"])
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plans", type=Path, required=True)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise ValueError("Output directory must be empty")
    provenance = []
    for component, name in NAMES.items():
        plan = json.loads((args.plans / f"{component}.json").read_text())
        records = checked_records(plan)
        source = "environment.yml" if component == "main" else f"environment.{component}.yml"
        source_hash = hashlib.sha256((args.repository / source).read_bytes()).hexdigest()
        (args.output / f"{name}-lock.yml").write_text(yaml_lock(records, source, source_hash))
        (args.output / f"{name}.packages.json").write_text(json.dumps(records, indent=2) + "\n")
        (args.output / f"{name}.explicit.txt").write_text(
            "# platform: linux-64\n@EXPLICIT\n" + "".join(f"{r['url']}#{r['md5']}\n" for r in records))
        for record in records:
            channel = urlsplit(record["url"]).path.split("/")[1]
            provenance.append({"environment": name, "publisher": f"{channel} community", **record})
        print(f"{name}: {len(records)} exact packages")
    with (args.output / "package_provenance.tsv").open("w") as handle:
        fields = ["environment", "name", "version", "build", "publisher", "subdir", "url", "md5", "sha256"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(provenance)
    files = sorted(path for path in args.output.iterdir() if path.is_file())
    (args.output / "SHA256SUMS").write_text("".join(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n" for path in files))


if __name__ == "__main__":
    main()
