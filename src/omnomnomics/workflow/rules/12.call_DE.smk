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
                f"list(factor={_r_string(item[0])}, numerator={_r_string(item[1])}, denominator={_r_string(item[2])})"
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
            de_result_root = os.path.join(params.outputfolder, de_subdir)
            qc_dir = os.path.join(de_result_root, "qc")
            de_dir = os.path.join(de_result_root, "differential_expression")

            os.makedirs(params.outputfolder, exist_ok=True)
            os.makedirs(qc_dir, exist_ok=True)
            os.makedirs(de_dir, exist_ok=True)

            metadata_copy = os.path.join(params.outputfolder, "metadata_derived.tsv")
            shutil.copy2(params.metadata_file, metadata_copy)

            contrasts_cfg = resolved_de_config.get("contrasts", {})
            filtering_cfg = resolved_de_config.get("filtering", {})
            deseq2_cfg = resolved_de_config.get("deseq2", {})
            lfc_cfg = deseq2_cfg.get("lfc_shrink", {})
            thresholds_cfg = resolved_de_config.get("thresholds", {})
            qc_cfg = resolved_de_config.get("qc", {})
            pca_cfg = qc_cfg.get("pca", {})
            plots_cfg = resolved_de_config.get("plots", {})
            volcano_cfg = plots_cfg.get("volcano", {})
            sig_heatmap_cfg = plots_cfg.get("sig_heatmap", {})
            tables_cfg = resolved_de_config.get("tables", {})
            runtime_cfg = resolved_de_config.get("runtime", {})

            template_path = os.path.join(workflow_root, "templates", "de_core.R.tmpl")
            if not os.path.isfile(template_path):
                raise FileNotFoundError(f"DE template file not found: {template_path}")

            template_values = {
                "__COUNTS_PATH__": _r_string(input.counts_table),
                "__METADATA_PATH__": _r_string(metadata_copy),
                "__OUTPUT_ROOT__": _r_string(de_result_root),
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
                "__PCA_SHAPE_BY__": _r_string(str(pca_cfg.get("shape_by") or "")),
                "__PCA_EXTRA_PAIRS__": _r_numeric_pair_list(pca_cfg.get("extra_pairs", [[1, 2], [2, 3]])),
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
                "__RUNTIME_SEED__": str(int(runtime_cfg.get("seed", 1337))),
            }

            analysis_script = os.path.join(params.outputfolder, rendered_script_name)
            rendered_r_script = _render_de_template(template_path, template_values)
            with open(analysis_script, "w") as handle:
                handle.write(rendered_r_script)

            os.chmod(analysis_script, 0o755)
            de_r_command = f"Rscript {analysis_script}"
            record_step_command(master_config['de_rule_num'], "aggregate", de_r_command)
            shell(de_r_command)

            with zipfile.ZipFile(output[0], "w", compression=zipfile.ZIP_DEFLATED) as archive:
                archive.write(metadata_copy, arcname=os.path.basename(metadata_copy))
                archive.write(analysis_script, arcname=os.path.basename(analysis_script))
                _add_tree_to_zip(archive, de_result_root, os.path.basename(de_result_root))

            log_it(logfile, f"Step 12 core DE results: {de_result_root}")
            log_it(logfile, f"Step 12 archive: {output[0]}")
            finish_step_sample(master_config['de_rule_num'], "aggregate", "call_DE", tracking["start_time"], "OK")
        except Exception:
            finish_step_sample(master_config['de_rule_num'], "aggregate", "call_DE", tracking["start_time"], "FAIL")
            raise
