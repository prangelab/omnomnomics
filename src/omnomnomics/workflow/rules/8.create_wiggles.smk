# Rule 9 Create Wiggles

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

rule create_wiggles:
    input:
        input_bam=f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.filtered.bam",
        input_bai=f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.sorted.dups_marked.filtered.bam.bai" if config['THETYPE'] != "CHIP" else f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}/{{sample}}.filtered.bam.bai"
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num']-1]}/{{sample}}.extra_8.tmp"
    params:
        thetype=config['THETYPE'],
        inputfolder=f"{experiment_dir}/{master_config['input_folders'][master_config['wig_rule_num']-1]}",
        outputfolder=f"{experiment_dir}/{master_config['output_folders'][master_config['wig_rule_num']-1]}"
    threads:
        Threads_Per_Rule['8']
    resources:
        mem_mb=Memory_Per_Rule['8'],
        partition=master_config['partition'],
        runtime=Runtime_Per_Rule['8']
    run:
        log_once(logfile, "step8.header", "Creating BigWigs...", f"EXECUTING STEP {master_config['wig_rule_num']}")
        log_once(logfile, "step8.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step8.outputfolder", f"Output folder: {params.outputfolder}")
        if params.thetype == "RNA":
            log_once(
                logfile,
                "step8.rna_stranded_mode",
                "RNA BigWigs are split into plus/minus strand tracks with bamCoverage --filterRNAstrand forward/reverse.",
            )

        bamcoverage_version = subprocess.check_output(["bamCoverage", "--version"], stderr=subprocess.STDOUT)
        log_once(logfile, "step8.bamcoverage_version", "\n" + bamcoverage_version.decode("utf-8"), "DEEPTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder, master_config['input_file_types'][master_config['wig_rule_num'] - 1], "step8.sanity")

        def create_wig(input_bam, input_bai, outputfolder, sample):
            tracking = begin_step_sample(master_config['wig_rule_num'], sample, "create_wiggles")
            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample}.", dir=tmpdir_root if tmpdir_root else None)
            record_step_note(master_config['wig_rule_num'], sample, f"scratch_dir={local_workdir}")

            def quote(path):
                return shlex.quote(path)

            def stage_input(path):
                local_path = os.path.join(local_workdir, os.path.basename(path))
                command = f"cp {quote(path)} {quote(local_path)}"
                record_step_note(master_config['wig_rule_num'], sample, f"staging {os.path.basename(path)}")
                shell(command)
                return local_path

            try:
                local_bam = stage_input(input_bam)
                local_bai = stage_input(input_bai)
                if params.thetype == "RNA":
                    local_plus_bw = os.path.join(local_workdir, f"{sample}.plus.bw")
                    local_minus_bw = os.path.join(local_workdir, f"{sample}.minus.bw")
                    output_plus_bw = os.path.join(outputfolder, f"{sample}.plus.bw")
                    output_minus_bw = os.path.join(outputfolder, f"{sample}.minus.bw")

                    plus_command = f"""
                        bamCoverage -b {quote(local_bam)} -o {quote(local_plus_bw)} \
                        --numberOfProcessors {threads} --binSize 10 --normalizeUsing CPM --filterRNAstrand forward
                    """
                    minus_command = f"""
                        bamCoverage -b {quote(local_bam)} -o {quote(local_minus_bw)} \
                        --numberOfProcessors {threads} --binSize 10 --normalizeUsing CPM --filterRNAstrand reverse
                    """
                    plus_command = " ".join(plus_command.split())
                    minus_command = " ".join(minus_command.split())
                    record_step_command(master_config['wig_rule_num'], sample, plus_command)
                    record_step_command(master_config['wig_rule_num'], sample, minus_command)
                    shell(plus_command)
                    shell(minus_command)

                    record_step_note(master_config['wig_rule_num'], sample, "copying_stranded_bigwigs_back")
                    shell(f"cp {quote(local_plus_bw)} {quote(output_plus_bw)}")
                    shell(f"cp {quote(local_minus_bw)} {quote(output_minus_bw)}")
                else:
                    local_bw = os.path.join(local_workdir, f"{sample}.bw")
                    output_bw = os.path.join(outputfolder, f"{sample}.bw")

                    bamcoverage_command = f"""
                        bamCoverage -b {quote(local_bam)} -o {quote(local_bw)} \
                        --numberOfProcessors {threads} --binSize 10 --normalizeUsing CPM
                    """
                    bamcoverage_command = " ".join(bamcoverage_command.split())
                    record_step_command(master_config['wig_rule_num'], sample, bamcoverage_command)
                    shell(bamcoverage_command)

                    copy_bw_command = f"cp {quote(local_bw)} {quote(output_bw)}"
                    record_step_note(master_config['wig_rule_num'], sample, "copying_bigwig_back")
                    shell(copy_bw_command)
                finish_step_sample(master_config['wig_rule_num'], sample, "create_wiggles", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['wig_rule_num'], sample, "create_wiggles", tracking["start_time"], "FAIL")
                raise
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

            shell(f"""echo "necessity file for create wiggles. can delete this." > {outputfolder}/{sample}.extra_8.tmp""")

        create_wig(input.input_bam, input.input_bai, params.outputfolder, wildcards.sample)
