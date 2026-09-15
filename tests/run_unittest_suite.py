from __future__ import annotations

import argparse
from pathlib import Path
import sys
import unittest

from suite_tiers import test_paths


def main() -> int:
    parser = argparse.ArgumentParser(description="Run an omnomnomics unittest tier.")
    parser.add_argument(
        "--tier",
        choices=("lightweight", "bioinformatics", "all"),
        default="all",
    )
    args = parser.parse_args()

    tests_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(tests_dir))
    loader = unittest.defaultTestLoader
    suite = unittest.TestSuite()
    paths = test_paths(tests_dir, args.tier)
    for path in paths:
        suite.addTests(loader.loadTestsFromName(path.stem))

    print(f"Test tier: {args.tier}; modules: {len(paths)}")
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
