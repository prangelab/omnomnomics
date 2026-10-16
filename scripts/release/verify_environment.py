#!/usr/bin/env python3
"""Check an installed prefix against the exact Conda archive inventory."""

import argparse
import json
from pathlib import Path


def verify(prefix, expected):
    actual = {}
    for path in (prefix / "conda-meta").glob("*.json"):
        record = json.loads(path.read_text())
        actual[record["name"]] = record
    desired = {record["name"]: record for record in expected}
    if actual.keys() != desired.keys():
        raise ValueError(f"Package-set mismatch; missing={sorted(desired.keys() - actual.keys())}; "
                         f"extra={sorted(actual.keys() - desired.keys())}")
    for name, record in desired.items():
        for key in ["version", "build", "md5"]:
            if actual[name].get(key) != record[key]:
                raise ValueError(f"Installed {name} does not match locked {key}")
        if actual[name].get("sha256") and actual[name]["sha256"] != record["sha256"]:
            raise ValueError(f"Installed {name} does not match locked sha256")
    print(f"Verified {len(actual)} exact Conda packages in {prefix}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--packages", type=Path, required=True)
    args = parser.parse_args()
    verify(args.prefix, json.loads(args.packages.read_text()))


if __name__ == "__main__":
    main()
