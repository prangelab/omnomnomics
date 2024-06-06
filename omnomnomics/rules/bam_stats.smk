# Rule 7: Get BAM file statistics

## Omnomnomics Snake Rule ##
import os

rule bam_stats:
    input:
        filtered_BAM=f"{master_config['input_folders'][master_config['stats_rule_num']-1]}/{{sample}}.filtered.bam"
    output:
        f"{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.filtered.bam.stats.txt"
    params:
        inputfolder = master_config['input_folders'][master_config['stats_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['stats_rule_num']-1]
    threads:
        Threads_Per_Rule['7']
    resources:
        mem_mb = Memory_Per_Rule['7']
    run:
        log_it(logfile, f"Getting BAM file statistics...",f"EXECUTING STEP {master_config['stats_rule_num']}")
        log_it(logfile, f"Input folder: filtere_BAM")
        log_it(logfile, f"Output folder: filtered_BAM")

        outfile = os.path.join(params.outputfolder, f"{wildcards.sample}.filtered.bam.stats.txt")
        log_it(logfile, f"Generating stats for {wildcards.sample}")
        shell(f"""
            module load bamtools && \
            bamtools stats -in {input.filtered_BAM} > {outfile}""")