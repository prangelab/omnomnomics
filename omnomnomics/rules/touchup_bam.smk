# Rule 3 touchup BAM

## Omnomnomics Snake Rule  ##
import os

rule touchup_bam:
    input:
        bamfile=f"{master_config['input_folders'][master_config['touchup_rule_num']-1]}/{{sample2}}.bam"
    output:
        filtered_BAM1=f"{master_config['output_folders'][master_config['touchup_rule_num']-1]}/{{sample2}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{master_config['output_folders'][master_config['touchup_rule_num']-1]}/{{sample}}.filtered.bam"
    params:
        thetype=config['THETYPE'],
        inputfolder = master_config['input_folders'][master_config['touchup_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['touchup_rule_num']-1]
    threads:
        Threads_Per_Rule['5']
    resources:
        mem_mb = Memory_Per_Rule['5']
    run:
        print("HEREEEE!!!!!!!")
        log_it(logfile, f"Touching up BAM files with samtools (collate | fixmate | sort | markdup | filter)...",f"EXECUTING STEP {master_config['touchup_rule_num']}")
        log_it(logfile, f"Input folder: BAM")
        log_it(logfile, f"Output folder: filtered_BAM")
        samtools_version = subprocess.check_output("module load samtools && samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+samtools_version.decode("utf-8"), "SAMTOOLS VERSION")
        print(samtools_version.decode("utf-8"))
        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['touchup_rule_num']-1])

        def touchup_bam_file(input_file, samcores, thetype, sample, outputfolder):
            # input_file_new = input_file.replace("_merged","")
            # shell(f"mv {input_file} {input_file_new}")

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
                shell(
                    f"""samtools view {input_file} | awk '{{if($3~/chrM/)chrm=chrm+1}}END{{printf \"%13s\\t%11d\\n\", \"Total Reads:\",NR;printf \"%13s\\t%11d %1s%2.2f%2s\\n\", \"ChrM Reads:\",chrm, \"(\", (chrm/NR)*100,\"%)\"}}' > filtered_BAM/{sample}.ATAC_stats.txt"""
                )
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
        samcores = max(1, threads // 5)
        touchup_bam_file(input.bamfile, samcores, params.thetype, wildcards.sample2, params.outputfolder)