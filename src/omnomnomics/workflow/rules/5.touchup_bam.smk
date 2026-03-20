# Rule 5 touchup BAM

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import shlex
import shutil
import subprocess
import tempfile

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
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['5']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num']-1]}/benchmarks/{{sample5}}_touchupbam_benchmark.tsv"
    run:
        log_it(logfile, f"Touching up BAM files with samtools (collate | fixmate | sort | markdup | filter)...",f"EXECUTING STEP {master_config['touchup_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")
        samtools_version = subprocess.check_output("samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+samtools_version.decode("utf-8"), "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['touchup_rule_num']-1])

        def touchup_bam_file(input_file, samcores, thetype, sample, outputfolder):
            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample}.", dir=tmpdir_root if tmpdir_root else None)
            log_it(logfile, f"Scratch directory: {local_workdir}")

            def quote(path):
                return shlex.quote(path)

            local_input = os.path.join(local_workdir, os.path.basename(input_file))
            stage_input_command = f'cp {quote(input_file)} {quote(local_input)}'
            log_it(logfile, stage_input_command, "STAGE INPUT COMMAND")
            shell(stage_input_command)

            try:
                local_output = os.path.join(
                    local_workdir,
                    f"{sample}.sorted.dups_marked.filtered.bam" if thetype != "CHIP" else f"{sample}.filtered.bam",
                )
                local_extra = os.path.join(local_workdir, f"{sample}.extra_5.tmp")

                if thetype == "RNA":
                    log_it(logfile, f"Touching up RNA-seq sample {input_file} (sort, mark dups, and keep MAPQ > 15)")
                    command = f"""
                        samtools collate -O -@ {samcores} {quote(local_input)} {quote(os.path.join(local_workdir, f'collate.{sample}.tmp'))} | \
                        samtools fixmate -mu -@ {samcores} - - | \
                        samtools sort -u -@ {samcores} - | \
                        samtools markdup -u -@ {samcores} - - | \
                        samtools view -@ {samcores} -q 15 -b -o {quote(local_output)} -
                    """
                elif thetype == "ATAC":
                    log_it(logfile, f"Touching up ATAC-seq sample {input_file} (sort, mark dups, and keep MAPQ > 15)")
                    command = f"""
                        samtools collate -O -@ {samcores} {quote(local_input)} {quote(os.path.join(local_workdir, f'collate.{sample}.tmp'))} | \
                        samtools fixmate -mu -@ {samcores} - - | \
                        samtools sort -u -@ {samcores} - | \
                        samtools markdup -u -@ {samcores} - - | \
                        samtools view -@ {samcores} -q 15 -b -o {quote(local_output)} -
                    """
                else:
                    log_it(logfile, f"Touching up ChIP-seq sample {input_file} (sort, remove dups, and keep MAPQ > 15)")
                    command = f"""
                        samtools collate -O -@ {samcores} {quote(local_input)} {quote(os.path.join(local_workdir, f'collate.{sample}.tmp'))} | \
                        samtools fixmate -mu -@ {samcores} - - | \
                        samtools sort -u -@ {samcores} - | \
                        samtools markdup -ru -@ {samcores} - - | \
                        samtools view -@ {samcores} -q 15 -b -o {quote(local_output)} -
                    """

                command = " ".join(command.split())
                log_it(logfile, command, "SAMTOOLS TOUCHUP COMMAND")
                shell(command)

                if thetype == "ATAC":
                    log_it(logfile, f"Getting ATAC statistics for sample {input_file}")
                    atac_stats_path = os.path.join(outputfolder, f"{sample}.ATAC_stats.txt")
                    shell(f"""samtools view {quote(local_input)} | awk '{{{{if($3~/chrM/)chrm=chrm+1}}}}END{{{{printf \"%13s\\t%11d\\n\", \"Total Reads:\",NR;printf \"%13s\\t%11d %1s%2.2f%2s\\n\", \"ChrM Reads:\",chrm, \"(\", (chrm/NR)*100,\"%)\"}}}}' > {quote(atac_stats_path)}""")

                shell(f"""echo "necessity file for touchup_bam. can delete this." > {quote(local_extra)}""")
                shell(f'cp {quote(local_output)} {quote(os.path.join(outputfolder, os.path.basename(local_output)))}')
                shell(f'cp {quote(local_extra)} {quote(os.path.join(outputfolder, os.path.basename(local_extra)))}')
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)
        samcores = max(1, threads // 5)
        touchup_bam_file(input.bamfile, samcores, params.thetype, wildcards.sample5, params.outputfolder)
