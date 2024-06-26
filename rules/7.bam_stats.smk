# Rule 7: Get Bam Stats

## Omnomnomics Snake Rule  ##
import os
import sys

rule bam_stats:
    input:
        filtered_BAM=f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}/{{sample}}.filtered.bam",
        bai_BAM = f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.bai" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}/{{sample}}.filtered.bam.bai"
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.stats.txt" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}.filtered.bam.stats.txt"
    params:
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['stats_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}"
    threads:
        Threads_Per_Rule['7']
    resources:
        mem_mb = Memory_Per_Rule['7']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['stats_rule_num']-1]}/benchmarks/{{sample}}_bamstats_benchmark.tsv"
    run:
        log_it(logfile, f"Getting BAM file statistics...",f"EXECUTING STEP {master_config['stats_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        bamtools_version = subprocess.check_output(["bamtools", "--version"])
        log_it(logfile, "\n"+bamtools_version.decode("utf-8"), "BAMTOOLS VERSION")
        
        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['stats_rule_num']-1])

        outfile = os.path.join(params.outputfolder, f"{wildcards.sample}.sorted.dups_marked.filtered.bam.stats.txt") if config['THETYPE'] != "CHIP" else os.path.join(params.outputfolder, f"{wildcards.sample}.filtered.bam.stats.txt")
        log_it(logfile, f"Generating stats for {wildcards.sample}")
        shell(f"""bamtools stats -in {input.filtered_BAM} > {outfile}""")