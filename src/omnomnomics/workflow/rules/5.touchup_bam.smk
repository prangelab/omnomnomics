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
        duplicate_handling=config['DUPLICATE_HANDLING'],
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['touchup_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['touchup_rule_num']-1]}"
    threads:
        Threads_Per_Rule['5']
    resources:
        mem_mb = Memory_Per_Rule['5'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['5']
    run:
        log_once(logfile, "step5.header", "Touching up BAM files with samtools (collate | fixmate | sort | markdup | filter)...", f"EXECUTING STEP {master_config['touchup_rule_num']}")
        log_once(logfile, "step5.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step5.outputfolder", f"Output folder: {params.outputfolder}")
        log_once(logfile, "step5.duplicate_handling", f"Duplicate handling: {params.duplicate_handling}")
        samtools_version = subprocess.check_output("samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_once(logfile, "step5.samtools_version", "\n"+samtools_version.decode("utf-8"), "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['touchup_rule_num']-1])

        def touchup_bam_file(input_file, samcores, thetype, duplicate_handling, sample, outputfolder):
            tracking = begin_step_sample(master_config['touchup_rule_num'], sample, "touchup_bam")
            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample}.", dir=tmpdir_root if tmpdir_root else None)
            record_step_note(master_config['touchup_rule_num'], sample, f"scratch_dir={local_workdir}")

            def quote(path):
                return shlex.quote(path)

            local_input = os.path.join(local_workdir, os.path.basename(input_file))
            stage_input_command = f'cp {quote(input_file)} {quote(local_input)}'
            record_step_note(master_config['touchup_rule_num'], sample, "staging_bam_to_scratch")
            shell(stage_input_command)

            try:
                local_output = os.path.join(
                    local_workdir,
                    f"{sample}.sorted.dups_marked.filtered.bam" if thetype != "CHIP" else f"{sample}.filtered.bam",
                )
                local_extra = os.path.join(local_workdir, f"{sample}.extra_5.tmp")
                min_mapq = 15 if thetype == "RNA" else 30
                filter_flags = 2820
                markdup_args = "-ru" if duplicate_handling == "remove" else "-u"

                if thetype == "RNA":
                    record_step_note(master_config['touchup_rule_num'], sample, f"running_rna_touchup duplicate_handling={duplicate_handling} mapq={min_mapq}")
                    command = f"""
                        samtools collate -O -@ {samcores} {quote(local_input)} {quote(os.path.join(local_workdir, f'collate.{sample}.tmp'))} | \
                        samtools fixmate -mu -@ {samcores} - - | \
                        samtools sort -u -@ {samcores} - | \
                        samtools markdup {markdup_args} -@ {samcores} - - | \
                        samtools view -@ {samcores} -q {min_mapq} -F {filter_flags} -b -o {quote(local_output)} -
                    """
                elif thetype == "ATAC":
                    record_step_note(master_config['touchup_rule_num'], sample, f"running_atac_touchup duplicate_handling={duplicate_handling} mapq={min_mapq}")
                    command = f"""
                        samtools collate -O -@ {samcores} {quote(local_input)} {quote(os.path.join(local_workdir, f'collate.{sample}.tmp'))} | \
                        samtools fixmate -mu -@ {samcores} - - | \
                        samtools sort -u -@ {samcores} - | \
                        samtools markdup {markdup_args} -@ {samcores} - - | \
                        samtools view -@ {samcores} -q {min_mapq} -F {filter_flags} -b -o {quote(local_output)} -
                    """
                else:
                    record_step_note(master_config['touchup_rule_num'], sample, f"running_chip_touchup duplicate_handling={duplicate_handling} mapq={min_mapq}")
                    command = f"""
                        samtools collate -O -@ {samcores} {quote(local_input)} {quote(os.path.join(local_workdir, f'collate.{sample}.tmp'))} | \
                        samtools fixmate -mu -@ {samcores} - - | \
                        samtools sort -u -@ {samcores} - | \
                        samtools markdup {markdup_args} -@ {samcores} - - | \
                        samtools view -@ {samcores} -q {min_mapq} -F {filter_flags} -b -o {quote(local_output)} -
                    """

                command = " ".join(command.split())
                record_step_command(master_config['touchup_rule_num'], sample, command)
                shell(command)

                if thetype == "ATAC":
                    record_step_note(master_config['touchup_rule_num'], sample, "writing_atac_mitochondrial_summary")
                    atac_stats_path = os.path.join(outputfolder, f"{sample}.ATAC_stats.txt")
                    shell(f"""samtools view {quote(local_input)} | awk '{{{{if($3~/chrM/)chrm=chrm+1}}}}END{{{{printf \"%34s\\t%11d\\n\", \"Total aligned reads before filtering:\",NR;printf \"%34s\\t%11d %1s%2.2f%2s\\n\", \"chrM aligned reads before filtering:\",chrm, \"(\", (chrm/NR)*100,\"%)\"}}}}' > {quote(atac_stats_path)}""")

                shell(f"""echo "necessity file for touchup_bam. can delete this." > {quote(local_extra)}""")
                record_step_note(master_config['touchup_rule_num'], sample, "copying_filtered_bam_back")
                shell(f'cp {quote(local_output)} {quote(os.path.join(outputfolder, os.path.basename(local_output)))}')
                shell(f'cp {quote(local_extra)} {quote(os.path.join(outputfolder, os.path.basename(local_extra)))}')
                finish_step_sample(master_config['touchup_rule_num'], sample, "touchup_bam", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['touchup_rule_num'], sample, "touchup_bam", tracking["start_time"], "FAIL")
                raise
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)
        samcores = threads
        touchup_bam_file(input.bamfile, samcores, params.thetype, params.duplicate_handling, wildcards.sample5, params.outputfolder)
