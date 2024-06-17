# Rule 7: Get Bam Stats

## Omnomnomics Snake Rule  ##
import os

rule bam_stats:
    input:
        filtered_BAM=f"{master_config['input_folders'][master_config['stats_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{master_config['input_folders'][master_config['stats_rule_num']-1]}/{{sample}}.filtered.bam"
    output:
        f"{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.stats.txt" if config['THETYPE'] != "CHIP" else f"{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.filtered.bam.stats.txt"
    params:
        inputfolder = master_config['input_folders'][master_config['stats_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['stats_rule_num']-1]
    threads:
        Threads_Per_Rule['7']
    resources:
        mem_mb = Memory_Per_Rule['7']
    benchmark:
        f"{master_config['output_folders'][master_config['stats_rule_num']-1]}/{{sample}}_benchmark.tsv"
    run:
        import sys
        print(sys.executable)
        log_it(logfile, f"Getting BAM file statistics...",f"EXECUTING STEP {master_config['stats_rule_num']}")
        log_it(logfile, f"Input folder: filtere_BAM")
        log_it(logfile, f"Output folder: filtered_BAM")

        bamtools_version = subprocess.check_output(["bamtools", "--version"])

        log_it(logfile, "\n"+bamtools_version.decode("utf-8"), "BAMTOOLS VERSION")
        print(bamtools_version.decode("utf-8"))
        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['stats_rule_num']-1])

        outfile = os.path.join(params.outputfolder, f"{wildcards.sample}.sorted.dups_marked.filtered.bam.stats.txt") if config['THETYPE'] != "CHIP" else os.path.join(params.outputfolder, f"{wildcards.sample}.filtered.bam.stats.txt")
        log_it(logfile, f"Generating stats for {wildcards.sample}")
        shell(f"""bamtools stats -in {input.filtered_BAM} > {outfile}""")