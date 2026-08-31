from __future__ import annotations

import importlib.util
import inspect
from pathlib import Path
import traceback


def main() -> int:
    tests_dir = Path(__file__).resolve().parent
    failures: list[tuple[str, str, str]] = []
    total = 0

    for path in sorted(tests_dir.glob("test_*.py")):
        spec = importlib.util.spec_from_file_location(path.stem, path)
        if spec is None or spec.loader is None:
            failures.append((path.name, "import", "Could not create an import specification."))
            continue
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
        except Exception:
            failures.append((path.name, "import", traceback.format_exc()))
            continue

        for name, function in inspect.getmembers(module, inspect.isfunction):
            if not name.startswith("test_") or inspect.signature(function).parameters:
                continue
            total += 1
            try:
                function()
            except Exception:
                failures.append((path.name, name, traceback.format_exc()))

    print(f"Standalone regression functions: {total}; failures: {len(failures)}")
    for path, name, detail in failures:
        print(f"FAIL {path}::{name}\n{detail}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
