from argparse import Namespace
from pathlib import Path

from omnomnomics.de_app import APP_DIR, build_launch, parse_arguments


def test_differential_explorer_assets_are_packaged():
    assert (APP_DIR / "app.R").is_file()
    assert (APP_DIR / "global.R").is_file()


def test_differential_explorer_defaults():
    args = parse_arguments([])
    assert args.project_dir == "."
    assert args.host == "127.0.0.1"
    assert args.port == 3838
    assert args.no_browser is False


def test_differential_explorer_launch_uses_resolved_project(tmp_path):
    args = Namespace(
        project_dir=str(tmp_path),
        host="127.0.0.1",
        port=4848,
        no_browser=True,
    )
    cmd, env = build_launch(args)

    assert cmd[:2] == ["Rscript", "-e"]
    assert "port=4848" in cmd[2]
    assert "launch.browser=FALSE" in cmd[2]
    assert env["OMNOMNOMICS_DE_APP_PROJECT"] == str(Path(tmp_path).resolve())
