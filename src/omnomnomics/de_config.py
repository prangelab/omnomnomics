from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import yaml


class DEConfigError(ValueError):
    pass


DEFAULT_DE_CONFIG: dict = {
    "version": 1,
    "io": {
        "out_dir": "results",
        "write_rendered_r_script": True,
        "rendered_r_script_name": "DE_analysis.rendered.R",
        "save_r_session_image": False,
    },
    "design": {
        "formula": None,
        "reference_levels": {},
    },
    "contrasts": {
        "mode": "auto",
        "auto": {
            "pairwise": True,
            "primary_variable": None,
        },
        "explicit": {
            "items": [],
        },
    },
    "filtering": {
        "enabled": True,
        "method": "min_count_samples",
        "min_count": 10,
        "min_samples": 2,
    },
    "deseq2": {
        "fit_type": "parametric",
        "sf_type": "ratio",
        "beta_prior": False,
        "parallel": True,
        "lfc_shrink": {
            "enabled": True,
            "type": "apeglm",
            "use_for_tables": True,
        },
    },
    "thresholds": {
        "alpha": 0.05,
        "lfc_for_sig": 1.0,
        "lfc_for_heatmap": 1.0,
        "max_volcano_labels": 20,
    },
    "qc": {
        "enabled": True,
        "vst_blind": True,
        "top_variable_genes": 1000,
        "distance_plot": True,
        "variable_gene_heatmap": True,
        "pca": {
            "enabled": True,
            "color_by": ["sample_type", "sample_color", "replicate"],
            "shape_by": ["donor", "replicate"],
            "extra_pairs": [[1, 2], [2, 3]],
        },
    },
    "plots": {
        "ma_plot": True,
        "volcano": {
            "enabled": True,
            "labeled": True,
            "auto_scale_axes": True,
            "max_xlim": 8,
            "max_ylim": 50,
        },
        "sig_heatmap": {
            "enabled": True,
            "cluster_rows": True,
            "cluster_cols": True,
        },
    },
    "tables": {
        "write_full_results": True,
        "write_sig_only_table": True,
        "sig_only_name_suffix": ".sig_only.tsv",
    },
    "runtime": {
        "seed": 1337,
        "stop_on_error": True,
        "verbose": True,
    },
}


def _deep_merge(base: dict, override: dict) -> dict:
    merged = deepcopy(base)
    for key, value in override.items():
        if (
            key in merged
            and isinstance(merged[key], dict)
            and isinstance(value, dict)
        ):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def load_de_config_file(config_path: str | None) -> tuple[str, dict]:
    if not config_path or config_path == "NA":
        return "NA", {}
    resolved = Path(config_path).expanduser().resolve()
    if not resolved.is_file():
        raise DEConfigError(f"DE config file '{resolved}' does not exist.")
    try:
        loaded = yaml.safe_load(resolved.read_text()) or {}
    except yaml.YAMLError as exc:
        raise DEConfigError(f"Could not parse DE config file '{resolved}': {exc}") from exc
    if not isinstance(loaded, dict):
        raise DEConfigError("DE config root must be a YAML mapping.")
    return str(resolved), loaded


def resolve_de_config(
    config_path: str | None,
    resolved_formula: str,
    de_out_dir_override: str | None = None,
) -> tuple[str, dict]:
    resolved_path, user_config = load_de_config_file(config_path)
    resolved = _deep_merge(DEFAULT_DE_CONFIG, user_config)

    configured_formula = resolved.get("design", {}).get("formula")
    if configured_formula in ("", "NA"):
        configured_formula = None
    if configured_formula is None:
        resolved["design"]["formula"] = resolved_formula
    else:
        resolved["design"]["formula"] = str(configured_formula).strip()

    if de_out_dir_override:
        resolved["io"]["out_dir"] = str(de_out_dir_override).strip()

    _validate_resolved_de_config(resolved)
    return resolved_path, resolved


def _validate_resolved_de_config(config: dict) -> None:
    if str(config.get("version", "")) != "1":
        raise DEConfigError("DE config version must be 1.")

    contrast_mode = str(config["contrasts"]["mode"]).strip().lower()
    if contrast_mode not in {"auto", "explicit"}:
        raise DEConfigError("DE config contrasts.mode must be 'auto' or 'explicit'.")
    config["contrasts"]["mode"] = contrast_mode

    formula = str(config["design"]["formula"]).strip()
    if not formula.startswith("~"):
        raise DEConfigError("Resolved DE formula must start with '~'.")

    explicit_items = config["contrasts"]["explicit"]["items"]
    if not isinstance(explicit_items, list):
        raise DEConfigError("DE config contrasts.explicit.items must be a list.")
    for item in explicit_items:
        if not isinstance(item, list) or len(item) != 3:
            raise DEConfigError(
                "Each explicit contrast must be a 3-item list: [factor, numerator, denominator]."
            )

    out_dir = str(config["io"]["out_dir"]).strip()
    if not out_dir:
        raise DEConfigError("DE config io.out_dir must not be empty.")

    top_var = int(config["qc"]["top_variable_genes"])
    if top_var < 10:
        raise DEConfigError("DE config qc.top_variable_genes must be >= 10.")

    alpha = float(config["thresholds"]["alpha"])
    if alpha <= 0 or alpha >= 1:
        raise DEConfigError("DE config thresholds.alpha must be between 0 and 1.")

    min_count = int(config["filtering"]["min_count"])
    min_samples = int(config["filtering"]["min_samples"])
    if min_count < 0:
        raise DEConfigError("DE config filtering.min_count must be >= 0.")
    if min_samples < 1:
        raise DEConfigError("DE config filtering.min_samples must be >= 1.")
