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
        "write_customization_guide": True,
        "customization_guide_script_name": "DE_analysis.customization_guide.R",
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
        "parallel": False,
        "latent_factors": {
            "enabled": False,
            "method": "sva",
            "n_sv": None,
        },
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
            "color_by": ["sample_type", "condition", "sample_color"],
            "shape_by": [],
            "extra_pairs": [],
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
    "enrichment": {
        "enabled": True,
        "clusterprofiler": {
            "enabled": True,
            "run_ora": True,
            "run_gsea": True,
            "msigdb_sets": [
                {"name": "msig_hallmark", "category": "H"},
                {"name": "msig_reactome", "category": "C2", "subcategory": "CP:REACTOME"},
            ],
            "pvalue_cutoff": 0.05,
            "qvalue_cutoff": 0.2,
            "min_gs_size": 10,
            "max_gs_size": 500,
            "top_terms": 20,
            "gsea_permutations": 1000,
        },
        "decoupler": {
            "enabled": True,
            "run_progeny": True,
            "run_tf_network": True,
            "progeny_top": 500,
            "tf_split_complexes": False,
            "minsize": 5,
            "top_features_heatmap": 25,
            "top_regulators_barplot": 25,
            "top_regulators_detail_each_side": 2,
        },
        "custom_modules": {
            "enabled": False,
            "gmt_file": None,
            "name": "custom_modules",
            "run_ora": True,
            "run_gsea": True,
        },
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
    custom_modules_gmt: str | None = None,
    enable_custom_modules: bool = False,
) -> tuple[str, dict]:
    resolved_path, user_config = load_de_config_file(config_path)
    resolved = _deep_merge(DEFAULT_DE_CONFIG, user_config)

    # Backward compatibility for older enrichment.gene_sets boolean config.
    cp_cfg = resolved.get("enrichment", {}).get("clusterprofiler", {})
    if "msigdb_sets" not in cp_cfg:
        legacy_gene_sets = cp_cfg.get("gene_sets", {}) if isinstance(cp_cfg.get("gene_sets", {}), dict) else {}
        msigdb_sets = []
        if legacy_gene_sets.get("msig_hallmark", False):
            msigdb_sets.append({"name": "msig_hallmark", "category": "H"})
        if legacy_gene_sets.get("msig_c2_cp", False):
            msigdb_sets.append({"name": "msig_c2_cp", "category": "C2", "subcategory": "CP"})
        if msigdb_sets:
            cp_cfg["msigdb_sets"] = msigdb_sets

    configured_formula = resolved.get("design", {}).get("formula")
    if configured_formula in ("", "NA"):
        configured_formula = None
    if configured_formula is None:
        resolved["design"]["formula"] = resolved_formula
    else:
        resolved["design"]["formula"] = str(configured_formula).strip()

    if de_out_dir_override:
        resolved["io"]["out_dir"] = str(de_out_dir_override).strip()

    cm_cfg = resolved.setdefault("enrichment", {}).setdefault("custom_modules", {})
    if custom_modules_gmt:
        cm_cfg["gmt_file"] = str(Path(custom_modules_gmt).expanduser().resolve())
        cm_cfg["enabled"] = True
    if enable_custom_modules:
        cm_cfg["enabled"] = True

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
    normalized_explicit_items = []
    for item in explicit_items:
        if isinstance(item, list):
            if len(item) != 3:
                raise DEConfigError(
                    "List-style explicit contrasts must be [factor, numerator, denominator]."
                )
            factor = str(item[0]).strip()
            numerator = str(item[1]).strip()
            denominator = str(item[2]).strip()
            if not factor or not numerator or not denominator:
                raise DEConfigError(
                    "List-style explicit contrasts must have non-empty factor, numerator, denominator."
                )
            normalized_explicit_items.append(
                {
                    "contrast_type": "factor",
                    "factor": factor,
                    "numerator": numerator,
                    "denominator": denominator,
                    "label": f"{numerator}_vs_{denominator}",
                }
            )
            continue

        if isinstance(item, dict):
            contrast_type = str(
                item.get("contrast_type", item.get("type", "factor"))
            ).strip().lower()
            if contrast_type == "factor":
                factor = str(item.get("factor", "")).strip()
                numerator = str(item.get("numerator", "")).strip()
                denominator = str(item.get("denominator", "")).strip()
                if not factor or not numerator or not denominator:
                    raise DEConfigError(
                        "Factor explicit contrast requires non-empty 'factor', 'numerator', and 'denominator'."
                    )
                label = str(item.get("label", f"{numerator}_vs_{denominator}")).strip()
                normalized_explicit_items.append(
                    {
                        "contrast_type": "factor",
                        "factor": factor,
                        "numerator": numerator,
                        "denominator": denominator,
                        "label": label,
                    }
                )
                continue

            if contrast_type == "coefficient":
                coefficient_name = str(
                    item.get("coefficient_name", item.get("name", ""))
                ).strip()
                if not coefficient_name:
                    raise DEConfigError(
                        "Coefficient explicit contrast requires non-empty 'coefficient_name' (or 'name')."
                    )
                label = str(item.get("label", coefficient_name)).strip()
                normalized_explicit_items.append(
                    {
                        "contrast_type": "coefficient",
                        "coefficient_name": coefficient_name,
                        "label": label,
                    }
                )
                continue

            raise DEConfigError(
                "Explicit contrast dict contrast_type must be 'factor' or 'coefficient'."
            )

        raise DEConfigError(
            "Each explicit contrast must be either a 3-item list [factor, numerator, denominator] "
            "or a mapping with contrast_type."
        )
    config["contrasts"]["explicit"]["items"] = normalized_explicit_items

    out_dir = str(config["io"]["out_dir"]).strip()
    if not out_dir:
        raise DEConfigError("DE config io.out_dir must not be empty.")
    rendered_name = str(config["io"].get("rendered_r_script_name", "")).strip()
    if not rendered_name:
        raise DEConfigError("DE config io.rendered_r_script_name must not be empty.")
    guide_name = str(config["io"].get("customization_guide_script_name", "")).strip()
    if not guide_name:
        raise DEConfigError("DE config io.customization_guide_script_name must not be empty.")

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

    latent_cfg = config.get("deseq2", {}).get("latent_factors", {})
    if latent_cfg:
        latent_method = str(latent_cfg.get("method", "sva")).strip().lower()
        if latent_method not in {"sva"}:
            raise DEConfigError("DE config deseq2.latent_factors.method must be 'sva'.")
        latent_cfg["method"] = latent_method
        latent_n_sv = latent_cfg.get("n_sv")
        if latent_n_sv is not None:
            try:
                latent_n_sv = int(latent_n_sv)
            except (TypeError, ValueError) as exc:
                raise DEConfigError("DE config deseq2.latent_factors.n_sv must be an integer or null.") from exc
            if latent_n_sv < 1:
                raise DEConfigError("DE config deseq2.latent_factors.n_sv must be >= 1 when set.")
            latent_cfg["n_sv"] = latent_n_sv

    enrichment_cfg = config.get("enrichment", {})
    cp_cfg = enrichment_cfg.get("clusterprofiler", {})
    if cp_cfg:
        min_gs_size = int(cp_cfg.get("min_gs_size", 10))
        max_gs_size = int(cp_cfg.get("max_gs_size", 500))
        top_terms = int(cp_cfg.get("top_terms", 20))
        pvalue_cutoff = float(cp_cfg.get("pvalue_cutoff", 0.05))
        qvalue_cutoff = float(cp_cfg.get("qvalue_cutoff", 0.2))
        gsea_permutations = int(cp_cfg.get("gsea_permutations", 1000))
        if min_gs_size < 1:
            raise DEConfigError("DE config enrichment.clusterprofiler.min_gs_size must be >= 1.")
        if max_gs_size < min_gs_size:
            raise DEConfigError(
                "DE config enrichment.clusterprofiler.max_gs_size must be >= min_gs_size."
            )
        if top_terms < 1:
            raise DEConfigError("DE config enrichment.clusterprofiler.top_terms must be >= 1.")
        if pvalue_cutoff <= 0 or pvalue_cutoff > 1:
            raise DEConfigError(
                "DE config enrichment.clusterprofiler.pvalue_cutoff must be in (0, 1]."
            )
        if qvalue_cutoff <= 0 or qvalue_cutoff > 1:
            raise DEConfigError(
                "DE config enrichment.clusterprofiler.qvalue_cutoff must be in (0, 1]."
            )
        if gsea_permutations < 100:
            raise DEConfigError(
                "DE config enrichment.clusterprofiler.gsea_permutations must be >= 100."
            )
        msigdb_sets = cp_cfg.get("msigdb_sets", [])
        if not isinstance(msigdb_sets, list) or len(msigdb_sets) == 0:
            raise DEConfigError(
                "DE config enrichment.clusterprofiler.msigdb_sets must be a non-empty list."
            )
        for idx, entry in enumerate(msigdb_sets):
            if not isinstance(entry, dict):
                raise DEConfigError(
                    f"DE config enrichment.clusterprofiler.msigdb_sets[{idx}] must be a mapping."
                )
            name = str(entry.get("name", "")).strip()
            category = str(entry.get("category", "")).strip().upper()
            if not name:
                raise DEConfigError(
                    f"DE config enrichment.clusterprofiler.msigdb_sets[{idx}].name must not be empty."
                )
            if not category:
                raise DEConfigError(
                    f"DE config enrichment.clusterprofiler.msigdb_sets[{idx}].category must not be empty."
                )
            subcategory = entry.get("subcategory")
            if subcategory is not None and not str(subcategory).strip():
                raise DEConfigError(
                    f"DE config enrichment.clusterprofiler.msigdb_sets[{idx}].subcategory must be non-empty when provided."
                )

    dc_cfg = enrichment_cfg.get("decoupler", {})
    if dc_cfg:
        progeny_top = int(dc_cfg.get("progeny_top", 500))
        minsize = int(dc_cfg.get("minsize", 5))
        top_features_heatmap = int(dc_cfg.get("top_features_heatmap", 25))
        top_regulators_barplot = int(dc_cfg.get("top_regulators_barplot", 25))
        top_regulators_detail_each_side = int(dc_cfg.get("top_regulators_detail_each_side", 2))
        if progeny_top < 50:
            raise DEConfigError("DE config enrichment.decoupler.progeny_top must be >= 50.")
        if minsize < 1:
            raise DEConfigError("DE config enrichment.decoupler.minsize must be >= 1.")
        if top_features_heatmap < 2:
            raise DEConfigError(
                "DE config enrichment.decoupler.top_features_heatmap must be >= 2."
            )
        if top_regulators_barplot < 2:
            raise DEConfigError(
                "DE config enrichment.decoupler.top_regulators_barplot must be >= 2."
            )
        if top_regulators_detail_each_side < 1:
            raise DEConfigError(
                "DE config enrichment.decoupler.top_regulators_detail_each_side must be >= 1."
            )

    cm_cfg = enrichment_cfg.get("custom_modules", {})
    if cm_cfg:
        cm_name = str(cm_cfg.get("name", "custom_modules")).strip()
        if not cm_name:
            raise DEConfigError("DE config enrichment.custom_modules.name must not be empty.")
        cm_enabled = bool(cm_cfg.get("enabled", False))
        cm_gmt = cm_cfg.get("gmt_file")
        if cm_enabled:
            if not cm_gmt or not str(cm_gmt).strip():
                raise DEConfigError(
                    "DE config enrichment.custom_modules.gmt_file is required when custom_modules.enabled is true."
                )
            cm_path = Path(str(cm_gmt)).expanduser().resolve()
            if not cm_path.is_file():
                raise DEConfigError(
                    f"DE config enrichment.custom_modules.gmt_file does not exist: {cm_path}"
                )
            if cm_path.suffix.lower() != ".gmt":
                raise DEConfigError(
                    "DE config enrichment.custom_modules.gmt_file must have a .gmt extension."
                )
            cm_cfg["gmt_file"] = str(cm_path)
