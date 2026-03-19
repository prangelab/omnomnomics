# Rule 5 touchup BAM

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os

rule touchup_bam:
    input:
        bamfile=f"{experiment_dir}/{master_config['input_folders'][master_config['touchup_rule_num']-1]}/{{sample5}}.bam",
        extrafile = f"{experiment_dir}/{master_config['input_folders'][master_config['touchup_rule_num']-1]}/{{sample5}}.extra_4.tmp" if 4 in themode else []
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num']-1]}/{{sample5}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num']-1]}/{{sample5}}.filtered.bam",
        f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num']-1]}/{{sample5}}.extra_5.tmp"
    params:
        thetype=config['THETYPE'],
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['touchup_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num']-1]}"
    threads:
        Threads_Per_Rule['5']
    resources:
        mem_mb = Memory_Per_Rule['5'],
        partition = master_config['partition']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num']-1]}/benchmarks/{{sample5}}_touchupbam_benchmark.tsv"
    run:
        log_it(logfile, f"Touching up BAM files with samtools (collate | fixmate | sort | markdup | filter)...",f"EXECUTING STEP {master_config['touchup_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")
        samtools_version = subprocess.check_output("module load samtools && samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+samtools_version.decode("utf-8"), "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['touchup_rule_num']-1])

        def touchup_bam_file(input_file, samcores, thetype, sample, outputfolder):

            if thetype == "RNA":
                log_it(logfile, f"Touching up RNA-seq sample {input_file} (sort, mark dups, and keep MAPQ > 15)")
                shell(f"""
                    module load samtools && \
                    samtools collate -O -@ {samcores} {input_file} collate.{sample}.tmp | \
                    samtools fixmate -mu -@ {samcores} - - | \
                    samtools sort -u -@ {samcores} - | \
                    samtools markdup -u -@ {samcores} - - | \
                    samtools view -@ {samcores} -q 15 -b -o {outputfolder}/{sample}.sorted.dups_marked.filtered.bam - """
                )
            elif thetype == "ATAC":
                log_it(logfile, f"Touching up ATAC-seq sample {input_file} (sort, mark dups, and keep MAPQ > 15)")
                shell(f"""
                    module load samtools && \
                    samtools collate -O -@ {samcores} {input_file} collate.{sample}.tmp | \
                    samtools fixmate -mu -@ {samcores} - - | \
                    samtools sort -u -@ {samcores} - | \
                    samtools markdup -u -@ {samcores} - - | \
                    samtools view -@ {samcores} -q 15 -b -o {outputfolder}/{sample}.sorted.dups_marked.filtered.bam - """
                )
                log_it(logfile, f"Getting ATAC statistics for sample {input_file}")
                shell(f"""module load samtools && samtools view {input_file} | awk '{{{{if($3~/chrM/)chrm=chrm+1}}}}END{{{{printf \"%13s\\t%11d\\n\", \"Total Reads:\",NR;printf \"%13s\\t%11d %1s%2.2f%2s\\n\", \"ChrM Reads:\",chrm, \"(\", (chrm/NR)*100,\"%)\"}}}}' > filtered_BAM/{sample}.ATAC_stats.txt""")       

            else:
                log_it(logfile, f"Touching up ChIP-seq sample {input_file} (sort, remove dups, and keep MAPQ > 15)")
                shell(f"""
                    module load samtools && \
                    samtools collate -O -@ {samcores} {input_file} collate.{sample}.tmp | \
                    samtools fixmate -mu -@ {samcores} - - | \
                    samtools sort -u -@ {samcores} - | \
                    samtools markdup -ru -@ {samcores} - - | \
                    samtools view -@ {samcores} -q 15 -b -o {outputfolder}/{sample}.filtered.bam -"""
                )

            shell(f"""echo "necessity file for touchup_bam. can delete this." > {outputfolder}/{sample}.extra_5.tmp""")
        samcores = max(1, threads // 5)
        touchup_bam_file(input.bamfile, samcores, params.thetype, wildcards.sample5, params.outputfolder)
