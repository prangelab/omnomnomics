# Rule 12 call DE

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import shutil
import textwrap
import zipfile


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
    threads:
        Threads_Per_Rule[str(master_config['de_rule_num'])]
    resources:
        mem_mb=Memory_Per_Rule[str(master_config['de_rule_num'])],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule[str(master_config['de_rule_num'])]
    run:
        tracking = begin_step_sample(master_config['de_rule_num'], "aggregate", "call_DE")
        log_once(logfile, "step12.header", "Preparing DESeq2 analysis scaffold...", f"EXECUTING STEP {master_config['de_rule_num']}")
        log_once(logfile, "step12.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step12.outputfolder", f"Output folder: {params.outputfolder}")
        log_once(logfile, "step12.metadata_source", f"Metadata source: {params.metadata_source}")
        log_once(logfile, "step12.derived_metadata", f"Derived metadata: {params.metadata_file}")
        log_once(logfile, "step12.design_formula", f"Resolved DE formula: {params.resolved_formula}")

        def render_stub_script(script_path, counts_path, metadata_path, formula, design_mode):
            script_text = textwrap.dedent(
                f"""\
                #!/usr/bin/env Rscript

                counts_path <- {counts_path!r}
                metadata_path <- {metadata_path!r}
                design_formula_text <- {formula!r}
                design_mode <- {design_mode!r}

                message("Omnomnomics DESeq2 scaffold")
                message("Counts table: ", counts_path)
                message("Derived metadata: ", metadata_path)
                message("Resolved design formula: ", design_formula_text)
                message("Design mode: ", design_mode)

                stop(
                  "Step 12 scaffolding is in place, but the project-specific DESeq2 analysis has not been wired in yet."
                )
                """
            )
            with open(script_path, "w") as handle:
                handle.write(script_text)

        try:
            if params.thetype != "RNA":
                raise ValueError("Step 12 is currently supported for RNA only.")

            sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['de_rule_num']-1], "step12.sanity")
            if params.metadata_file == "NA" or not os.path.isfile(params.metadata_file):
                raise FileNotFoundError("Derived metadata file for step 12 was not found.")

            os.makedirs(params.outputfolder, exist_ok=True)
            metadata_copy = os.path.join(params.outputfolder, "metadata_derived.tsv")
            shutil.copy2(params.metadata_file, metadata_copy)

            analysis_script = os.path.join(
                params.outputfolder,
                f"{os.path.basename(config['EXPERIMENT_DIR'])}.deseq2_analysis.R",
            )
            render_stub_script(
                analysis_script,
                input.counts_table,
                metadata_copy,
                params.resolved_formula,
                params.design_mode,
            )

            with zipfile.ZipFile(output[0], "w", compression=zipfile.ZIP_DEFLATED) as archive:
                archive.write(metadata_copy, arcname=os.path.basename(metadata_copy))
                archive.write(analysis_script, arcname=os.path.basename(analysis_script))

            log_it(logfile, f"Step 12 scaffold archive: {output[0]}")
            finish_step_sample(master_config['de_rule_num'], "aggregate", "call_DE", tracking["start_time"], "OK")
        except Exception:
            finish_step_sample(master_config['de_rule_num'], "aggregate", "call_DE", tracking["start_time"], "FAIL")
            raise
