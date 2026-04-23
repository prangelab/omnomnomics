# Rule 12 call DE

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import json
import os
import shutil
import zipfile

import yaml


def _r_string(value):
    text = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{text}"'


def _r_bool(value):
    return "TRUE" if bool(value) else "FALSE"


def _r_char_vector(values):
    cleaned = [str(item) for item in values if str(item).strip()]
    if not cleaned:
        return "character(0)"
    return "c(" + ", ".join(_r_string(item) for item in cleaned) + ")"


def _r_contrast_list(items):
    if not items:
        return "list()"
    rows = []
    for item in items:
        if isinstance(item, list) and len(item) == 3:
            rows.append(
                f"list(contrast_type={_r_string('factor')}, factor={_r_string(item[0])}, numerator={_r_string(item[1])}, denominator={_r_string(item[2])}, coefficient_name=NA_character_, label={_r_string(str(item[1]) + '_vs_' + str(item[2]))})"
            )
            continue
        if isinstance(item, dict):
            contrast_type = str(item.get("contrast_type", item.get("type", "factor"))).strip().lower()
            if contrast_type == "factor":
                factor = str(item.get("factor", "")).strip()
                numerator = str(item.get("numerator", "")).strip()
                denominator = str(item.get("denominator", "")).strip()
                if not factor or not numerator or not denominator:
                    continue
                label = str(item.get("label", f"{numerator}_vs_{denominator}")).strip()
                rows.append(
                    f"list(contrast_type={_r_string('factor')}, factor={_r_string(factor)}, numerator={_r_string(numerator)}, denominator={_r_string(denominator)}, coefficient_name=NA_character_, label={_r_string(label)})"
                )
            elif contrast_type == "coefficient":
                coefficient_name = str(item.get("coefficient_name", item.get("name", ""))).strip()
                if not coefficient_name:
                    continue
                label = str(item.get("label", coefficient_name)).strip()
                rows.append(
                    f"list(contrast_type={_r_string('coefficient')}, factor=NA_character_, numerator=NA_character_, denominator=NA_character_, coefficient_name={_r_string(coefficient_name)}, label={_r_string(label)})"
                )
    if not rows:
        return "list()"
    return "list(\n  " + ",\n  ".join(rows) + "\n)"


def _r_numeric_pair_list(items):
    if not items:
        return "list(c(1, 2), c(2, 3))"
    pairs = []
    for item in items:
        if isinstance(item, list) and len(item) == 2:
            try:
                first = int(item[0])
                second = int(item[1])
            except (TypeError, ValueError):
                continue
            pairs.append(f"c({first}, {second})")
    if not pairs:
        return "list(c(1, 2), c(2, 3))"
    return "list(" + ", ".join(pairs) + ")"


def _r_msigdb_set_list(items):
    if not isinstance(items, list) or not items:
        return "list()"
    rows = []
    for item in items:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", "")).strip()
        category = str(item.get("category", "")).strip()
        if not name or not category:
            continue
        subcategory = item.get("subcategory")
        if subcategory is None or not str(subcategory).strip():
            rows.append(
                f"list(name={_r_string(name)}, category={_r_string(category)}, subcategory=NA_character_)"
            )
        else:
            rows.append(
                f"list(name={_r_string(name)}, category={_r_string(category)}, subcategory={_r_string(str(subcategory).strip())})"
            )
    if not rows:
        return "list()"
    return "list(\n  " + ",\n  ".join(rows) + "\n)"


def _add_tree_to_zip(archive, path_on_disk, archive_prefix):
    for root, _, files in os.walk(path_on_disk):
        for file_name in files:
            abs_path = os.path.join(root, file_name)
            rel_path = os.path.relpath(abs_path, path_on_disk)
            archive.write(abs_path, arcname=os.path.join(archive_prefix, rel_path))


def _render_de_template(template_path, replacements):
    with open(template_path) as handle:
        rendered = handle.read()
    for key, value in replacements.items():
        rendered = rendered.replace(key, value)
    return rendered


def _build_de_customization_guide(
    counts_table_path,
    metadata_path,
    result_root,
    resolved_formula,
    de_columns_resolved,
    de_block_resolved,
):
    de_columns_text = ", ".join(str(x) for x in de_columns_resolved) if de_columns_resolved else "<none>"
    de_block_text = ", ".join(str(x) for x in de_block_resolved) if de_block_resolved else "<none>"
    return f"""

# -------------------------------------------------------------------------------------------------
# Customization guide for DE analyses
# -------------------------------------------------------------------------------------------------
# This section is intentionally non-executing. It documents common modifications with project paths.
# Keep the defaults above for exact run reproduction, then adapt options below for custom analyses.
#
# Active run inputs
#   counts table: {counts_table_path}
#   metadata:     {metadata_path}
#   output root:  {result_root}
#   formula:      {resolved_formula}
#   de columns:   {de_columns_text}
#   blocks:       {de_block_text}
#
# 1) Factor contrasts (main effects / grouped effects)
# Example structure:
# contrast_plan <- list(
#   list(contrast_type = "factor", factor = "treatment", numerator = "L", denominator = "C",
#        coefficient_name = NA_character_, label = "treatment_L_vs_C"),
#   list(contrast_type = "factor", factor = "type", numerator = "KD24", denominator = "NT",
#        coefficient_name = NA_character_, label = "type_KD24_vs_NT"),
#   list(contrast_type = "factor", factor = "de_group", numerator = "KD24_L", denominator = "KD24_C",
#        coefficient_name = NA_character_, label = "de_group_KD24_L_vs_KD24_C")
# )
#
# 2) Coefficient contrasts (interaction coefficients)
# Inspect available coefficients written during this run:
#   file.path(qc_dir, "deseq2_results_names.tsv")
#
# Example structure:
# contrast_plan <- list(
#   list(contrast_type = "coefficient", factor = NA_character_, numerator = NA_character_,
#        denominator = NA_character_, coefficient_name = "typeKD24.treatmentL",
#        label = "interaction_typeKD24_treatmentL")
# )
# You can also provide shorthand in de_config coefficient_name (for example "type:treatment").
# The pipeline resolves shorthand to an exact DESeq2 coefficient when unique.
#
# 3) Mixed contrast plan (factor + coefficient)
# contrast_plan <- list(
#   list(contrast_type = "factor", factor = "treatment", numerator = "L", denominator = "C",
#        coefficient_name = NA_character_, label = "treatment_L_vs_C"),
#   list(contrast_type = "coefficient", factor = NA_character_, numerator = NA_character_,
#        denominator = NA_character_, coefficient_name = "typeKD24.treatmentL",
#        label = "interaction_typeKD24_treatmentL")
# )
#
# 4) Custom modules from GMT (already supported in this script)
# custom_modules_enabled <- TRUE
# custom_modules_gmt <- "/absolute/path/to/modules.gmt"
# custom_modules_name <- "custom_modules"
#
# Notes:
# - Factor contrast plots are subset to the numerator/denominator samples.
# - Coefficient contrasts are model-coefficient tests; they do not subset samples by factor levels.
# - For reproducibility, keep this rendered script together with metadata_derived.tsv and de_run_manifest.tsv.
# -------------------------------------------------------------------------------------------------
"""


rule call_DE:
    input:
        counts_table=(
            f"{experiment_dir}/{master_config['input_folders'][master_config['de_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.raw_read_quant.table.txt"
            if config['THETYPE'] == "RNA"
            else []
        )
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['de_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.results.zip"
    params:
        thetype=config['THETYPE'],
        inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['de_rule_num']-1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['de_rule_num']-1]}",
        metadata_file=config.get("DERIVED_METADATA_FILE", "NA"),
        metadata_source=config.get("MYMETADATA", "NA"),
        resolved_formula=config.get("RESOLVED_DE_FORMULA", "NA"),
        design_mode=config.get("DE_DESIGN_MODE", "NA"),
        resolved_de_config_file=config.get("DE_CONFIG_RESOLVED_FILE", "NA"),
        resolved_de_config_json=config.get("DE_CONFIG_RESOLVED_JSON", "{}"),
        de_columns_resolved=config.get("DE_COLUMNS_RESOLVED", []),
        de_block_resolved=config.get("DE_BLOCK_RESOLVED", []),
    threads:
        Threads_Per_Rule[str(master_config['de_rule_num'])]
    resources:
        mem_mb=Memory_Per_Rule[str(master_config['de_rule_num'])],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule[str(master_config['de_rule_num'])]
    run:
        tracking = begin_step_sample(master_config['de_rule_num'], "aggregate", "call_DE")
        log_once(logfile, "step12.header", "Running core DESeq2 analysis...", f"EXECUTING STEP {master_config['de_rule_num']}")
        log_once(logfile, "step12.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step12.outputfolder", f"Output folder: {params.outputfolder}")
        log_once(logfile, "step12.metadata_source", f"Metadata source: {params.metadata_source}")
        log_once(logfile, "step12.derived_metadata", f"Derived metadata: {params.metadata_file}")
        log_once(logfile, "step12.design_formula", f"Resolved DE formula: {params.resolved_formula}")
        log_once(logfile, "step12.de_config", f"Resolved DE config: {params.resolved_de_config_file}")

        try:
            if params.thetype != "RNA":
                raise ValueError("Step 12 is currently supported for RNA only.")

            sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['de_rule_num']-1], "step12.sanity")
            if params.metadata_file == "NA" or not os.path.isfile(params.metadata_file):
                raise FileNotFoundError("Derived metadata file for step 12 was not found.")

            if params.resolved_de_config_file != "NA" and os.path.isfile(params.resolved_de_config_file):
                with open(params.resolved_de_config_file) as handle:
                    resolved_de_config = yaml.safe_load(handle) or {}
            else:
                resolved_de_config = json.loads(params.resolved_de_config_json or "{}")

            if not resolved_de_config:
                raise ValueError("Resolved DE config is missing.")

            io_cfg = resolved_de_config.get("io", {})
            de_subdir = str(io_cfg.get("out_dir", "results")).strip() or "results"
            rendered_script_name = str(io_cfg.get("rendered_r_script_name", "DE_analysis.rendered.R"))
            write_customization_guide = bool(io_cfg.get("write_customization_guide", True))
            customization_guide_script_name = str(
                io_cfg.get("customization_guide_script_name", "DE_analysis.customization_guide.R")
            ).strip() or "DE_analysis.customization_guide.R"
            if os.path.normpath(de_subdir).lower() == "qc":
                raise ValueError("DE output subdirectory cannot be 'qc' because that path is reserved for shared QC outputs.")
            qc_dir = os.path.join(params.outputfolder, "qc")
            de_dir = os.path.join(params.outputfolder, de_subdir)

            os.makedirs(params.outputfolder, exist_ok=True)
            os.makedirs(qc_dir, exist_ok=True)
            os.makedirs(de_dir, exist_ok=True)

            metadata_copy = os.path.join(params.outputfolder, "metadata_derived.tsv")
            shutil.copy2(params.metadata_file, metadata_copy)

            contrasts_cfg = resolved_de_config.get("contrasts", {})
            filtering_cfg = resolved_de_config.get("filtering", {})
            deseq2_cfg = resolved_de_config.get("deseq2", {})
            latent_cfg = deseq2_cfg.get("latent_factors", {})
            lfc_cfg = deseq2_cfg.get("lfc_shrink", {})
            thresholds_cfg = resolved_de_config.get("thresholds", {})
            qc_cfg = resolved_de_config.get("qc", {})
            pca_cfg = qc_cfg.get("pca", {})
            plots_cfg = resolved_de_config.get("plots", {})
            volcano_cfg = plots_cfg.get("volcano", {})
            sig_heatmap_cfg = plots_cfg.get("sig_heatmap", {})
            tables_cfg = resolved_de_config.get("tables", {})
            enrichment_cfg = resolved_de_config.get("enrichment", {})
            enrichment_cp_cfg = enrichment_cfg.get("clusterprofiler", {})
            enrichment_msigdb_sets_cfg = enrichment_cp_cfg.get("msigdb_sets", [])
            enrichment_dc_cfg = enrichment_cfg.get("decoupler", {})
            enrichment_cm_cfg = enrichment_cfg.get("custom_modules", {})
            runtime_cfg = resolved_de_config.get("runtime", {})

            pca_shape_cfg = pca_cfg.get("shape_by", [])
            if isinstance(pca_shape_cfg, list):
                pca_shape_values = [str(item) for item in pca_shape_cfg if str(item).strip()]
            elif pca_shape_cfg:
                pca_shape_values = [str(pca_shape_cfg).strip()]
            else:
                pca_shape_values = []

            template_path = os.path.join(workflow_root, "templates", "de_core.R.tmpl")
            if not os.path.isfile(template_path):
                raise FileNotFoundError(f"DE template file not found: {template_path}")

            template_values = {
                "__COUNTS_PATH__": _r_string(input.counts_table),
                "__METADATA_PATH__": _r_string(metadata_copy),
                "__OUTPUT_ROOT__": _r_string(params.outputfolder),
                "__QC_DIR__": _r_string(qc_dir),
                "__DE_DIR__": _r_string(de_dir),
                "__DESIGN_FORMULA_TEXT__": _r_string(str(resolved_de_config.get("design", {}).get("formula", params.resolved_formula))),
                "__CONTRAST_MODE__": _r_string(str(contrasts_cfg.get("mode", "auto"))),
                "__AUTO_PAIRWISE__": _r_bool(contrasts_cfg.get("auto", {}).get("pairwise", True)),
                "__AUTO_PRIMARY_VARIABLE__": _r_string(str(contrasts_cfg.get("auto", {}).get("primary_variable") or "")),
                "__EXPLICIT_CONTRASTS__": _r_contrast_list(contrasts_cfg.get("explicit", {}).get("items", [])),
                "__FILTERING_ENABLED__": _r_bool(filtering_cfg.get("enabled", True)),
                "__FILTERING_METHOD__": _r_string(str(filtering_cfg.get("method", "min_count_samples"))),
                "__FILTERING_MIN_COUNT__": str(int(filtering_cfg.get("min_count", 10))),
                "__FILTERING_MIN_SAMPLES__": str(int(filtering_cfg.get("min_samples", 2))),
                "__DESEQ_FIT_TYPE__": _r_string(str(deseq2_cfg.get("fit_type", "parametric"))),
                "__DESEQ_SF_TYPE__": _r_string(str(deseq2_cfg.get("sf_type", "ratio"))),
                "__DESEQ_PARALLEL__": _r_bool(deseq2_cfg.get("parallel", True)),
                "__LATENT_ENABLED__": _r_bool(latent_cfg.get("enabled", False)),
                "__LATENT_METHOD__": _r_string(str(latent_cfg.get("method", "sva"))),
                "__LATENT_N_SV__": (
                    str(int(latent_cfg.get("n_sv")))
                    if latent_cfg.get("n_sv") not in (None, "", "NA")
                    else "NA_integer_"
                ),
                "__SHRINK_ENABLED__": _r_bool(lfc_cfg.get("enabled", True)),
                "__SHRINK_TYPE__": _r_string(str(lfc_cfg.get("type", "apeglm"))),
                "__SHRINK_FOR_TABLES__": _r_bool(lfc_cfg.get("use_for_tables", True)),
                "__ALPHA__": str(float(thresholds_cfg.get("alpha", 0.05))),
                "__LFC_FOR_SIG__": str(float(thresholds_cfg.get("lfc_for_sig", 1.0))),
                "__LFC_FOR_HEATMAP__": str(float(thresholds_cfg.get("lfc_for_heatmap", 1.0))),
                "__MAX_VOLCANO_LABELS__": str(int(thresholds_cfg.get("max_volcano_labels", 50))),
                "__QC_ENABLED__": _r_bool(qc_cfg.get("enabled", True)),
                "__VST_BLIND__": _r_bool(qc_cfg.get("vst_blind", True)),
                "__TOP_VARIABLE_GENES__": str(int(qc_cfg.get("top_variable_genes", 1000))),
                "__DISTANCE_PLOT_ENABLED__": _r_bool(qc_cfg.get("distance_plot", True)),
                "__VARIABLE_HEATMAP_ENABLED__": _r_bool(qc_cfg.get("variable_gene_heatmap", True)),
                "__PCA_ENABLED__": _r_bool(pca_cfg.get("enabled", True)),
                "__PCA_COLOR_BY__": _r_char_vector(pca_cfg.get("color_by", [])),
                "__PCA_SHAPE_BY_VALUES__": _r_char_vector(pca_shape_values),
                "__PCA_EXTRA_PAIRS__": _r_numeric_pair_list(pca_cfg.get("extra_pairs", [[1, 2], [2, 3]])),
                "__DE_COLUMNS_RESOLVED__": _r_char_vector(params.de_columns_resolved),
                "__DE_BLOCK_RESOLVED__": _r_char_vector(params.de_block_resolved),
                "__MA_PLOT_ENABLED__": _r_bool(plots_cfg.get("ma_plot", True)),
                "__VOLCANO_ENABLED__": _r_bool(volcano_cfg.get("enabled", True)),
                "__VOLCANO_LABELS_ENABLED__": _r_bool(volcano_cfg.get("labeled", True)),
                "__VOLCANO_AUTO_AXES__": _r_bool(volcano_cfg.get("auto_scale_axes", True)),
                "__VOLCANO_MAX_XLIM__": str(float(volcano_cfg.get("max_xlim", 8))),
                "__VOLCANO_MAX_YLIM__": str(float(volcano_cfg.get("max_ylim", 50))),
                "__SIG_HEATMAP_ENABLED__": _r_bool(sig_heatmap_cfg.get("enabled", True)),
                "__SIG_HEATMAP_CLUSTER_ROWS__": _r_bool(sig_heatmap_cfg.get("cluster_rows", True)),
                "__SIG_HEATMAP_CLUSTER_COLS__": _r_bool(sig_heatmap_cfg.get("cluster_cols", True)),
                "__WRITE_FULL_RESULTS__": _r_bool(tables_cfg.get("write_full_results", True)),
                "__WRITE_SIG_ONLY__": _r_bool(tables_cfg.get("write_sig_only_table", True)),
                "__SIG_SUFFIX__": _r_string(str(tables_cfg.get("sig_only_name_suffix", ".sig_only.tsv"))),
                "__ENRICHMENT_ENABLED__": _r_bool(enrichment_cfg.get("enabled", True)),
                "__ENRICHMENT_CP_ENABLED__": _r_bool(enrichment_cp_cfg.get("enabled", True)),
                "__ENRICHMENT_CP_RUN_ORA__": _r_bool(enrichment_cp_cfg.get("run_ora", True)),
                "__ENRICHMENT_CP_RUN_GSEA__": _r_bool(enrichment_cp_cfg.get("run_gsea", True)),
                "__ENRICHMENT_MSIGDB_SETS__": _r_msigdb_set_list(enrichment_msigdb_sets_cfg),
                "__ENRICHMENT_PVALUE_CUTOFF__": str(float(enrichment_cp_cfg.get("pvalue_cutoff", 0.05))),
                "__ENRICHMENT_QVALUE_CUTOFF__": str(float(enrichment_cp_cfg.get("qvalue_cutoff", 0.2))),
                "__ENRICHMENT_MIN_GS_SIZE__": str(int(enrichment_cp_cfg.get("min_gs_size", 10))),
                "__ENRICHMENT_MAX_GS_SIZE__": str(int(enrichment_cp_cfg.get("max_gs_size", 500))),
                "__ENRICHMENT_TOP_TERMS__": str(int(enrichment_cp_cfg.get("top_terms", 20))),
                "__ENRICHMENT_GSEA_PERMUTATIONS__": str(int(enrichment_cp_cfg.get("gsea_permutations", 1000))),
                "__DECOUPLER_ENABLED__": _r_bool(enrichment_dc_cfg.get("enabled", True)),
                "__DECOUPLER_RUN_PROGENY__": _r_bool(enrichment_dc_cfg.get("run_progeny", True)),
                "__DECOUPLER_RUN_TF_NETWORK__": _r_bool(enrichment_dc_cfg.get("run_tf_network", True)),
                "__DECOUPLER_PROGENY_TOP__": str(int(enrichment_dc_cfg.get("progeny_top", 500))),
                "__DECOUPLER_TF_SPLIT_COMPLEXES__": _r_bool(enrichment_dc_cfg.get("tf_split_complexes", False)),
                "__DECOUPLER_MINSIZE__": str(int(enrichment_dc_cfg.get("minsize", 5))),
                "__DECOUPLER_TOP_FEATURES_HEATMAP__": str(int(enrichment_dc_cfg.get("top_features_heatmap", 25))),
                "__DECOUPLER_TOP_REGULATORS_BARPLOT__": str(int(enrichment_dc_cfg.get("top_regulators_barplot", 25))),
                "__DECOUPLER_TOP_REGULATORS_DETAIL_EACH_SIDE__": str(int(enrichment_dc_cfg.get("top_regulators_detail_each_side", 2))),
                "__CUSTOM_MODULES_ENABLED__": _r_bool(enrichment_cm_cfg.get("enabled", False)),
                "__CUSTOM_MODULES_GMT__": _r_string(str(enrichment_cm_cfg.get("gmt_file") or "")),
                "__CUSTOM_MODULES_NAME__": _r_string(str(enrichment_cm_cfg.get("name", "custom_modules"))),
                "__CUSTOM_MODULES_RUN_ORA__": _r_bool(enrichment_cm_cfg.get("run_ora", True)),
                "__CUSTOM_MODULES_RUN_GSEA__": _r_bool(enrichment_cm_cfg.get("run_gsea", True)),
                "__RUNTIME_SEED__": str(int(runtime_cfg.get("seed", 1337))),
            }

            analysis_script = os.path.join(params.outputfolder, rendered_script_name)
            rendered_r_script = _render_de_template(template_path, template_values)
            with open(analysis_script, "w") as handle:
                handle.write(rendered_r_script)

            customization_guide_script = os.path.join(params.outputfolder, customization_guide_script_name)
            if write_customization_guide:
                guide_script_text = rendered_r_script + _build_de_customization_guide(
                    counts_table_path=input.counts_table,
                    metadata_path=metadata_copy,
                    result_root=params.outputfolder,
                    resolved_formula=str(resolved_de_config.get("design", {}).get("formula", params.resolved_formula)),
                    de_columns_resolved=params.de_columns_resolved,
                    de_block_resolved=params.de_block_resolved,
                )
                with open(customization_guide_script, "w") as handle:
                    handle.write(guide_script_text)
                os.chmod(customization_guide_script, 0o755)

            os.chmod(analysis_script, 0o755)
            de_r_command = f"Rscript {analysis_script}"
            record_step_command(master_config['de_rule_num'], "aggregate", de_r_command)
            shell(de_r_command)

            with zipfile.ZipFile(output[0], "w", compression=zipfile.ZIP_DEFLATED) as archive:
                archive.write(metadata_copy, arcname=os.path.basename(metadata_copy))
                archive.write(analysis_script, arcname=os.path.basename(analysis_script))
                if write_customization_guide and os.path.isfile(customization_guide_script):
                    archive.write(
                        customization_guide_script,
                        arcname=os.path.basename(customization_guide_script),
                    )
                _add_tree_to_zip(archive, qc_dir, os.path.basename(qc_dir))
                _add_tree_to_zip(archive, de_dir, os.path.basename(de_dir))

            log_it(logfile, f"Step 12 DE results folder: {de_dir}")
            log_it(logfile, f"Step 12 shared QC folder: {qc_dir}")
            if write_customization_guide and os.path.isfile(customization_guide_script):
                log_it(logfile, f"Step 12 customization guide: {customization_guide_script}")
            log_it(logfile, f"Step 12 archive: {output[0]}")
            finish_step_sample(master_config['de_rule_num'], "aggregate", "call_DE", tracking["start_time"], "OK")
        except Exception:
            finish_step_sample(master_config['de_rule_num'], "aggregate", "call_DE", tracking["start_time"], "FAIL")
            raise
