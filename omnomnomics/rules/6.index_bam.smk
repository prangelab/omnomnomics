# Rule 6 Index BAM

## Omnomnomics Snake Rule  ##
import os
rule index_bam:
    input:
        filtered_BAM= f"{master_config['input_folders'][master_config['index_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.filtered.bam",
    output:
        f"{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.bai" if config['THETYPE'] != "CHIP" else f"{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}.filtered.bam.bai"
    threads:
        Threads_Per_Rule['6']
    params:
        inputfolder = master_config['input_folders'][master_config['index_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['index_rule_num']-1]       
    resources:
        mem_mb = Memory_Per_Rule['6']
    benchmark:
        f"{master_config['output_folders'][master_config['index_rule_num']-1]}/{{sample}}_bamindex_benchmark.tsv"
    run:
        log_it(logfile, f"Index BAM files...",f"EXECUTING STEP {master_config['index_rule_num']}")
        log_it(logfile, f"Input folder: BAM")
        log_it(logfile, f"Output folder: filtered_BAM")
        log_it(logfile, f"Indexing {wildcards.sample}")
        samtools_version = subprocess.check_output("module load samtools && samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+samtools_version.decode("utf-8"), "SAMTOOLS VERSION")
        print(samtools_version.decode("utf-8"))
        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['index_rule_num']-1])
        shell(f"""
            module load samtools && \
            samtools index -@ {threads} {input.filtered_BAM}""")