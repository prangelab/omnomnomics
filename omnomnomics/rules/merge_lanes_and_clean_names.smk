#Rule 4 merge lanes and clean names

## Omnomnomics Snake Rule  ##
import os
import glob

rule merge_bam:
    input:
        bam_files = glob.glob(f"BAM/{{sample}}_L00*.bam")
    output:
        bam="merged_BAM/{sample}.bam"
    params:
        infolder="trimmed_FASTQ",
        outfolder="merged_BAM",
        maxcores=config["MAXCORES"],
        step_merge=4,
        inputfiletype=".bam"
    threads:
        Threads_Per_Rule['4']
    resources:
        mem_mb = Memory_Per_Rule['4']
    run:
        log_it(logfile, "Merging lanes...", f"EXECUTING STEP {master_config['merge_rule_num']}")
        log_it(logfile, "Input folder: BAM")
        log_it(logfile, "Output folder: BAM")

        def merge_bam_files(infolder, outfolder, inputfiletype, maxcores, sample):
            # Check if we have lane info
            if len(glob.glob("BAM/{sample}_L00*_HISAT2.bam")) > 0 or len(glob.glob("BAM/{sample}_L00*_STAR.bam")) > 0 or len(glob.glob("BAM/{sample}_L00*_STAR_TE.bam")) > 0:
                
                log_it(logfile, f"Sample: {sample} with multiple lane IDs (_L00x) detected:")
                # Determine number of cores to use from now based on number of samples
                MYCORES = int(maxcores) // len(sample_list)

                # Iterate over the samples
                for sample in sample_list:
                    # Get number of lanes (L00n)
                    num_lanes = len(set(os.path.basename(f).split('_L00')[1][0:3] for f in glob.glob(f"{infolder}/{sample}*{inputfiletype}")))

                    # Merge only if more than one lane was run
                    if num_lanes > 1:
                        # Collect the set of files per sample
                        bam_list = glob.glob(f"BAM/{sample}*.bam")

                        # Set clean name
                        myname = os.path.basename(bam_list[0]).split('_L00')[0]

                        # Run samtools merge
                        log_it(logfile, f"Merging {sample} lanes...")
                        shell(f"samtools merge -@ {MYCORES} -o BAM/{myname}.bam {' '.join(bam_list)}")
                    else:
                        # If there is only one name we don't need to bother with merging and can simply clean the name
                        from_file = glob.glob(f"BAM/{sample}*.bam")[0]
                        to_file = os.path.join(BAM, os.path.basename(glob.glob(os.path.join(BAM, f"{sample}*.bam"))[0]).split("_L00")[0] + ".bam")
                        log_it(logfile, f"Renaming {from_file} to {to_file}...")
                        shell(f"mv {from_file} {to_file}")

                for bam_file in glob.glob(f"BAM/*_L00*.bam"):
                    shell(f"rm {bam_file}")

            else:
                log_it(logfile, "No split lanes found!")

        merge_bam_files(THEINFOLDER, THEOUTFOLDER, INPUTFILETYPE, params.maxcores, wildcards.sample)