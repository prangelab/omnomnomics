from __future__ import annotations

import subprocess
import sys
from importlib import resources


def _run_packaged_script(script_name: str) -> None:
    script_resource = resources.files("omnomnomics") / "workflow" / "bin" / "scripts" / script_name
    with resources.as_file(script_resource) as script_path:
        completed = subprocess.run(["bash", str(script_path), *sys.argv[1:]], check=False)
    raise SystemExit(completed.returncode)


def create_track_color_table_main() -> None:
    _run_packaged_script("createTrackColorTable.sh")


def display_track_color_table_main() -> None:
    _run_packaged_script("displayTrackColorTable.sh")
