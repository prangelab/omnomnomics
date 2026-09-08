#!/usr/bin/env python3

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


PACKAGE_ROOT = Path(__file__).resolve().parent
APP_DIR = PACKAGE_ROOT / "workflow" / "R" / "shiny_app"


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Launch the omnomnomics Differential Explorer.",
        allow_abbrev=False,
    )
    parser.add_argument(
        "-i",
        "--project-dir",
        default=".",
        help="Project directory or DE_calling directory to pre-load. Default: current directory",
    )
    parser.add_argument("--host", default="127.0.0.1", help="Host interface. Default: 127.0.0.1")
    parser.add_argument("--port", type=int, default=3838, help="Port. Default: 3838")
    parser.add_argument("--no-browser", action="store_true", help="Do not open a browser window.")
    return parser.parse_args(argv)


def build_launch(args):
    project_dir = str(Path(args.project_dir).expanduser().resolve())
    env = os.environ.copy()
    env["OMNOMNOMICS_DE_APP_PROJECT"] = project_dir

    launch_browser = "FALSE" if args.no_browser else "TRUE"
    r_expr = (
        f"shiny::runApp({json.dumps(str(APP_DIR))}, "
        f"host={json.dumps(str(args.host))}, "
        f"port={int(args.port)}, "
        f"launch.browser={launch_browser})"
    )
    return ["Rscript", "-e", r_expr], env


def main(argv=None):
    args = parse_arguments(argv)
    app_entry = APP_DIR / "app.R"
    if not app_entry.is_file():
        print(f"Differential Explorer entry file '{app_entry}' does not exist. Aborting...", file=sys.stderr)
        return 1

    cmd, env = build_launch(args)
    try:
        completed = subprocess.run(cmd, env=env)
    except FileNotFoundError:
        print("Rscript was not found in PATH. Activate an environment with R installed and retry.", file=sys.stderr)
        return 1
    return completed.returncode


if __name__ == "__main__":
    sys.exit(main())
