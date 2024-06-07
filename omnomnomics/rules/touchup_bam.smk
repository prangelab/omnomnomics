# Rule 3 touchup BAM

## Omnomnomics Snake Rule  ##
import os

rule touchup_bam:
    input:
        bamfile=f"{master_config['input_folders'][master_config['touchup_rule_num']-1]}/{{sample}}.bam"
    output:
        filtered_BAM1=f"{master_config['output_folders'][master_config['touchup_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{master_config['output_folders'][master_config['touchup_rule_num']-1]}/{{sample}}.filtered.bam"
    params:
        thetype=config['THETYPE'],
        outputfolder = master_config['output_folders'][master_config['touchup_rule_num']-1]
    threads:
        Threads_Per_Rule['5']
    resources:
        mem_mb = Memory_Per_Rule['5']
    run:
        log_it(f"Touching up BAM files with samtools (collate | fixmate | sort | markdup | filter)...",f"EXECUTING STEP {master_config['touchup_rule_num']}")
        log_it(f"Input folder: BAM")
        log_it(f"Output folder: filtered_BAM")
        # def sanity_check_dir(directory, filetype):
        #     if not os.path.isdir(directory):
        #         raise Exception(f"Directory {directory} does not exist")
        #     if not any(f.endswith(filetype) for f in os.listdir(directory)):
        #         raise Exception(f"No files with type {filetype} in directory {directory}")


        def touchup_bam_file(bamfile, samcores, thetype, sample, outputfolder):
            if thetype == "RNA":
                log_it(logfile, f"Touching up RNA-seq sample {input} (sort, mark dups, and keep MAPQ > 15)")
                shell(f"""
                    module load samtools && \
                    samtools collate -O -@ {samcores} {input} collate.{sample}.tmp | \
                    samtools fixmate -mu -@ {samcores} - - | \
                    samtools sort -u -@ {samcores} - | \
                    samtools markdup -u -@ {samcores} - - | \
                    samtools view -@ {samcores} -q 15 -b -o {outputfolder}/{sample}.sorted.dups_marked.filtered.bam - """
                )
            elif thetype == "ATAC":
                log_it(logfile, f"Touching up ATAC-seq sample {input} (sort, mark dups, and keep MAPQ > 15)")
                shell(f"""
                    module load samtools && \
                    samtools collate -O -@ {samcores} {input} collate.{sample}.tmp | \
                    samtools fixmate -mu -@ {samcores} - - | \
                    samtools sort -u -@ {samcores} - | \
                    samtools markdup -u -@ {samcores} - - | \
                    samtools view -@ {samcores} -q 15 -b -o {outputfolder}/{sample}.sorted.dups_marked.filtered.bam - """
                )
                log_it(logfile, f"Getting ATAC statistics for sample {input}")
                shell(
                    f"""samtools view {input} | awk '{{if($3~/chrM/)chrm=chrm+1}}END{{printf \"%13s\\t%11d\\n\", \"Total Reads:\",NR;printf \"%13s\\t%11d %1s%2.2f%2s\\n\", \"ChrM Reads:\",chrm, \"(\", (chrm/NR)*100,\"%)\"}}' > filtered_BAM/{sample}.ATAC_stats.txt"""
                )
            else:
                log_it(logfile, f"Touching up ChIP-seq sample {input} (sort, remove dups, and keep MAPQ > 15)")
                shell(f"""
                    module load samtools && \
                    samtools collate -O -@ {samcores} {input} collate.{sample}.tmp | \
                    samtools fixmate -mu -@ {samcores} - - | \
                    samtools sort -u -@ {samcores} - | \
                    samtools markdup -ru -@ {samcores} - - | \
                    samtools view -@ {samcores} -q 15 -b -o {outputfolder}/{sample}.filtered.bam -"""
                )
        # sanity_check_dir(BAM, .bam)

        samcores = max(1, threads // 5)
        touchup_bam_file(input.bamfile, samcores, thetype, wildcards.sample, params.outputfolder)