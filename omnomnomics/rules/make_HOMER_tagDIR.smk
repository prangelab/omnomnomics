# Rule 8: Create HOMER tag directories

## Omnomnomics Snake Rule  ##
import os
import glob
import pysam

rule create_homer_tagDir:
    input:
        filtered_BAM=f"{master_config['input_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.filtered.bam"
    output:
        f"{master_config['output_folders'][master_config['tagdir_rule_num']-1]}/{{sample}}.HOMER_tagDir.tar.gz"
    threads:
        Threads_Per_Rule['8']
    resources:
        mem_mb = Memory_Per_Rule['8']
    params:
        genome = config['THEGENOME'],
        inputfolder = master_config['input_folders'][master_config['tagdir_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['tagdir_rule_num']-1]
    run:
        log_it(logfile, f"Creating HOMER tag directories...", f"EXECUTING STEP {master_config['tagdir']}")
        log_it(logfile, f"Input folder: filtered_BAM")
        log_it(logfile, f"Output folder: HOMER_tagDirs")

        version = subprocess.check_output(["perl", os.path.join(OMNOM_HOME, "bin", "homer", "configureHomer.pl"), "-list", "2>", "/dev/null", "|", "grep", "homer"])
        log_it(logfile, "\n"+version.decode("utf-8"), "VERSION")
        print(version.decode("utf-8"))
        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['tagdir_rule_num']-1])

        # Function to create HOMER tag directories
        def create_homer_tagDir(filepath, outfolder, genome):
            basename = os.path.basename(filepath)
            sample_name = basename.replace('.filtered.bam', '')
            tag_dir = os.path.join(outfolder, f"{sample_name}.HOMER_tagDir")
            tar_gz_path = f"{tag_dir}.tar.gz"

            # Remove existing output files if they exist
            if os.path.exists(tar_gz_path):
                log_it(logfile, f"Removing existing compressed tag directory {tar_gz_path}")
                os.remove(tar_gz_path)
            if os.path.exists(tag_dir):
                log_it(logfile, f"Removing existing uncompressed tag directory {tag_dir}")
                shell(f"rm -r {tag_dir}")

            if THETYPE == "RNA":
                with pysam.AlignmentFile(filepath, "rb") as bamfile:
                    paired_end_count = 0
                    for read in bamfile.fetch():
                        if read.is_paired:
                            paired_end_count += 1
                            break

                if paired_end_count > 0:
                    log_it(logfile, f"Running makeTagDirectory on paired-end RNA sample {filepath}")
                    shell(f"""
                        module load samtools && \
                        makeTagDirectory {tag_dir} {filepath} -genome {genome} -sspe -single """)
                else:
                    log_it(logfile, f"Running makeTagDirectory on {filepath}")
                    shell(f"""
                        module load samtools && \
                        makeTagDirectory {tag_dir} {filepath} -genome {genome} -single """)
            else:
                log_it(logfile, f"Running makeTagDirectory on {filepath}")
                shell(f"""
                    module load samtools && \
                    makeTagDirectory {tag_dir} {filepath} -genome {genome} -single""")

            # Compress the tag directory
            log_it(logfile, f"Compressing {tag_dir} into {tag_dir}.tar.gz")
            shell(f"tar czf {tag_dir}.tar.gz  {tag_dir}")
            
            # Remove the uncompressed tag directory
            log_it(logfile, f"Removing uncompressed tag directory {tag_dir}")
            shell(f"rm -r {tag_dir}")

        create_homer_tagDir(input.filtered_BAM, params.outputfolder, params.genome)