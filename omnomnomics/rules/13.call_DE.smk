# Rule 13 call DE

## Omnomnomics Snake Rule  ##
rule call_DE:
    input:
        f"{master_config['input_folders'][master_config['de_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.raw_read_quant.table.txt" if config['THETYPE'] != "CHIP" else []
    output:
        f"{master_config['output_folders'][master_config['de_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.results.zip"
    params:
        thetype = config['THETYPE'],
        inputfolder = master_config['input_folders'][master_config['de_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['de_rule_num']-1]
    threads:
        Threads_Per_Rule['13']
    resources:
        mem_mb = Memory_Per_Rule['13']
    benchmark:
        f"{master_config['output_folders'][master_config['de_rule_num']-1]}/call_DE_benchmark.tsv"
    run:
        def calling_DE(logfile, thetype, inputfolder, outputfolder):
            if thetype == "RNA":
                log_it(logfile, "Calling DE genes...", f"EXECUTING STEP {master_config['de_rule_num']}")
                log_it(logfile, f"Input folder: {inputfolder}")
                log_it(logfile, f"Output folder: {outputfolder}")
                sanity_check_dir(logfile, inputfolder,  master_config['input_file_types'][master_config['de_rule_num']-1])
                log_it(logfile, "stay tuned for this feature!")
                #....
            elif thetype == "ATAC":
                log_it(logfile, "Calling DE peaks...", f"EXECUTING STEP {master_config['de_rule_num']}")
                log_it(logfile, f"Input folder: {inputfolder}")
                log_it(logfile, f"Output folder: {outputfolder}")
                sanity_check_dir(logfile, inputfolder,  master_config['input_file_types'][master_config['de_rule_num']-1])
                log_it(logfile, "stay tuned for this feature!")
                #....
            else: #type = CHIP
                log_it(logfile, "Calling DE peaks...", f"EXECUTING STEP {master_config['de_rule_num']}")
                log_it(logfile, f"Input folder: {inputfolder}")
                log_it(logfile, f"Output folder: {outputfolder}")
                log_it(logfile, "To call DE peaks for ChIP  data, please first manually determine the best peak calling settings for your experiment and use run_quant_peaks.sh.")
                log_it(logfile, "Then execute 'run_call_DE_peaks.sh' on your optimal peak set.")
        calling_DE(logfile, params.thetype, inputfolder, params.outputfolder)      