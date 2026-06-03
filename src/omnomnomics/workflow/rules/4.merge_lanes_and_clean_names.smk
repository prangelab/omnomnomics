#Rule 4 merge lanes and clean names

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import glob
import shlex
import shutil
import subprocess
import tempfile

rule merge_bam:
    wildcard_constraints:
        sample4=merged_sample_wildcard_pattern
    input:
        extra_files=lambda wildcards: [
            f"{experiment_dir}/{master_config['input_folders'][master_config['merge_rule_num']-1]}/{lane_sample}.extra_3.tmp"
            for lane_sample in lane_samples_for_merged_sample(wildcards.sample4)
        ] if 3 in themode else []
    output:
        f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}/{{sample4}}.bam",
        f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}/{{sample4}}.extra_4.tmp"
    params:
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['merge_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}"
    threads:
        Threads_Per_Rule['4']
    resources:
        mem_mb = Memory_Per_Rule['4'],
        partition = master_config['partition'],
        runtime = Runtime_Per_Rule['4']
    run:
        log_once(logfile, "step4.header", "Merging BAM inputs...", f"EXECUTING STEP {master_config['merge_rule_num']}")
        log_once(logfile, "step4.inputfolder", f"Input folder: {params.inputfolder}")
        log_once(logfile, "step4.outputfolder", f"Output folder: {params.outputfolder}")
            
        merging_version = subprocess.check_output("samtools --version | head -n2", shell=True, executable='/bin/bash')
        log_once(logfile, "step4.samtools_version", "\n"+merging_version.decode("utf-8"), "SAMTOOLS VERSION")

        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['merge_rule_num']-1])

        def merge_bam_files(threads, inputfolder, sample4):
            tracking = begin_step_sample(master_config['merge_rule_num'], sample4, "merge_bam")
            tmpdir_root = os.environ.get("TMPDIR")
            local_workdir = tempfile.mkdtemp(prefix=f"{sample4}.", dir=tmpdir_root if tmpdir_root else None)
            record_step_note(master_config['merge_rule_num'], sample4, f"scratch_dir={local_workdir}")

            def quote(path):
                return shlex.quote(path)

            grouped_input_units = input_units_for_merged_sample(sample4)
            input_bamfiles = [
                os.path.join(inputfolder, f"{input_unit}.bam")
                for input_unit in grouped_input_units
                if os.path.exists(os.path.join(inputfolder, f"{input_unit}.bam"))
            ]
            output_bam = os.path.join(inputfolder, f"{sample4}.bam")
            technical_replicates = sorted(
                {
                    technical_replicate
                    for technical_replicate in (
                        metadata_value_for_sample(input_unit, "technical_replicate", "")
                        for input_unit in grouped_input_units
                    )
                    if technical_replicate
                }
            )
            try:
                if input_bamfiles:
                    record_step_note(
                        master_config['merge_rule_num'],
                        sample4,
                        "merge_inputs=" + ",".join(os.path.basename(path) for path in input_bamfiles),
                    )
                    if technical_replicates:
                        record_step_note(
                            master_config['merge_rule_num'],
                            sample4,
                            "technical_replicates=" + ",".join(technical_replicates),
                        )
                    if len(input_bamfiles) > 1:
                        local_bam_list = []
                        for bamfile in input_bamfiles:
                            local_bam = os.path.join(local_workdir, os.path.basename(bamfile))
                            shell(f"cp {quote(bamfile)} {quote(local_bam)}")
                            local_bam_list.append(local_bam)
                        local_output = os.path.join(local_workdir, f"{sample4}.bam")
                        merge_command = (
                            f"samtools merge -@ {threads} -o {quote(local_output)} "
                            + " ".join(quote(bam) for bam in local_bam_list)
                        )
                        record_step_note(master_config['merge_rule_num'], sample4, "merging_input_units")
                        record_step_command(master_config['merge_rule_num'], sample4, merge_command)
                        shell(merge_command)
                        shell(f"cp {quote(local_output)} {quote(output_bam)}")
                    else:
                        source_bam = input_bamfiles[0]
                        if os.path.abspath(source_bam) != os.path.abspath(output_bam):
                            record_step_note(master_config['merge_rule_num'], sample4, "renaming_single_input_bam")
                            shell(f"cp {quote(source_bam)} {quote(output_bam)}")
                        else:
                            record_step_note(master_config['merge_rule_num'], sample4, "single_input_bam_already_canonical")
                    shell(f"""echo "necessity file for merge bams. can delete this." > {inputfolder}/{wildcards.sample4}.extra_4.tmp""")
                elif os.path.exists(output_bam):
                    record_step_note(master_config['merge_rule_num'], sample4, "canonical_bam_already_present")
                    shell(f"""echo "necessity file for merge bams. can delete this." > {inputfolder}/{wildcards.sample4}.extra_4.tmp""")
                else:
                    raise FileNotFoundError(
                        f"No input BAMs resolved for merged sample '{sample4}' in {inputfolder}."
                    )
                finish_step_sample(master_config['merge_rule_num'], sample4, "merge_bam", tracking["start_time"], "OK")
            except Exception:
                finish_step_sample(master_config['merge_rule_num'], sample4, "merge_bam", tracking["start_time"], "FAIL")
                raise
            finally:
                shutil.rmtree(local_workdir, ignore_errors=True)

        merge_bam_files(threads, params.inputfolder, wildcards.sample4)
