#!/usr/bin/env python3
"""Constrain release environments to the installed dependency closure for solving."""

import argparse
import json
from pathlib import Path
import re


ENVIRONMENTS = {
    "main": ("environment.yml", "omnomnomics"),
    "idr": ("environment.idr.yml", "omnomnomics-idr"),
    "spp": ("environment.spp.yml", "omnomnomics-spp"),
    "explorer": ("environment.explorer.yml", "omnomnomics"),
}


def dependency_name(spec):
    return re.split(r"[ <>=!~\[]", spec, maxsplit=1)[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replace-package", action="append", default=[],
                        help="Re-solve a reviewed package while retaining pins on its dependency closure")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    summary = {}
    for name, (manifest, original) in ENVIRONMENTS.items():
        text = (args.repository / manifest).read_text().split("dependencies:", 1)[1]
        roots = [line.strip()[2:] for line in text.splitlines() if line.strip().startswith("- ")]
        if name == "main":
            roots.extend(["setuptools>=61", "wheel"])
        records = json.loads((args.inventory / original / "conda-records.json").read_text())
        by_name = {record["name"]: record for record in records}
        selected = set()
        missing = set()
        pending = [dependency_name(spec) for spec in roots]
        while pending:
            package = pending.pop()
            if package in selected or package.startswith("__"):
                continue
            if package not in by_name:
                missing.add(package)
                continue
            selected.add(package)
            pending.extend(dependency_name(spec) for spec in by_name[package]["depends"])
        pins = [f"{package}=={by_name[package]['version']}={by_name[package]['build']}"
                for package in sorted(selected) if package not in args.replace_package]
        (args.output / f"{name}.txt").write_text("\n".join([*roots, *pins]) + "\n")
        summary[name] = {"manifest": manifest, "direct_requirements": roots,
                         "pinned_closure": len(selected), "missing_from_original": sorted(missing)}
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
