# Rule 12: Count Reads

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os

def input_function(wilcards):
    input_files = []
    if config['THETYPE'] != "CHIP":
        if config['THETYPE'] == "ATAC":
            if 11 in themode:
                input_files.append(f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num']-1][1]}/extra_11.tmp")
            for sample in samples2:
                input_files.append(f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num']-1][0]}/{sample}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz")
        else:
            for sample in samples2:
                input_files.append(f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num']-1][0]}/{sample}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz")
    else:
        pass
    return input_files

rule count_reads:
    input:
        input_function
    output:
        #f"{master_config['output_folders'][master_config['countreads_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.raw_read_quant.table.txt" will get produced if type is RNA or ATAC
        f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num']-1]}/extra_12.tmp" #must specify this so that it always performs this rule if asked for, also vor CHIP.
    params:
        thetype=config['THETYPE'],  
        genome=config['THEGENOME'],  
        experiment_dir=config['EXPERIMENT_DIR'], 
        namefields=config['NAMEFIELDS'],  
        separator=config['THESEPARATOR'],
        inputfolder1 = f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num']-1][0]}",
        inputfolder2 = f"{experiment_dir}/{master_config['input_folders'][master_config['countreads_rule_num']-1][1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num']-1]}"
    threads:
        Threads_Per_Rule['12']
    resources:
        mem_mb = Memory_Per_Rule['12'],
        partition = master_config['partition']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['countreads_rule_num']-1]}/benchmarks/counts_reads_benchmark.tsv"
    run:
        log_it(logfile, "Counting Reads...", f"EXECUTING STEP {master_config['countreads_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder1} and also {params.inputfolder2} for ATAC data")
        log_it(logfile, f"Output folder: {params.outputfolder}")
        
        def count_reads_rna(input_folder, output_folder, genome, namefields, separator): 
            sanity_check_dir(logfile, input_folder,  master_config['input_file_types'][master_config['countreads_rule_num']-1][0]) 

            tagdir_files = glob.glob(os.path.join(input_folder, '*tagDir.tar.gz'))
            if len(tagdir_files) != 0:
                log_it(logfile, "Unpacking HOMER tagDir tar balls...")

                for file in tagdir_files:
                    basename = os.path.basename(file)
                    tag_dir = os.path.join(input_folder, basename.replace(".tar.gz", ""))
                    tag_dir_short = basename.replace(".tar.gz", "")

                    # unpack the tar.gz file
                    log_it(logfile, f"Unpacking {basename} in {input_folder}")
                    shell(f"""
                        mkdir -p {tag_dir} && \
                        cd {input_folder} && \
                        tar --strip-components=1 -xzf {basename} -C {tag_dir_short}
                        """)

            shell(f"ls -d {input_folder}/*tagDir/ > TAGDIRlist.txt") # Make a list of tag dirs
            log_it(logfile, "Counting reads with analyzeRepeats.pl...")
            shell(f"""analyzeRepeats.pl rna {genome} -dfile TAGDIRlist.txt -count exons -noadj > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt""")

            shell(f"""sed '1d' {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt | cut -f1,9- | sort -k1,1 | sed 's/ \\+/\\t/g' > clean.tmp""")
            
            #Command in script works better: cat TAGDIRlist.txt | xargs -l basename | cut -f {namefields} -d '{separator}' | awk 'BEGIN{{ORS="\t"}}{{print $0}}' > clean.header.tmp
            shell(f"""{config['SCRIPT_DIR']}/generate_header.sh {params.namefields} '{params.separator}' """)

            shell(f"""sed "1iRefSeq_ID\\t$(cat clean.header.tmp)" clean.tmp > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt""")

            shell("rm TAGDIRlist.txt clean.header.tmp clean.tmp")

            for tag_dir in tagdir_files:
                tag_dir_folder = tag_dir.replace(".tar.gz", "")
                # Remove the uncompressed tag directory
                log_it(logfile, f"Removing uncompressed tag directory {tag_dir_folder}")
                shell(f"rm -r {tag_dir_folder}")
            if os.path.isfile(f"{input_folder1}/genome.tags.tsv"):
                shell(f"rm {input_folder1}/genome.tags.tsv")
            if os.path.isfile(f"{input_folder1}/tagAutocorrelation.txt"):
                shell(f"rm {input_folder1}/tagAutocorrelation.txt")
            if os.path.isfile(f"{input_folder1}/tagCountDistribution.txt"):
                shell(f"rm {input_folder1}/tagCountDistribution.txt")
            if os.path.isfile(f"{input_folder1}/tagInfo.txt"):
                shell(f"rm {input_folder1}/tagInfo.txt")
            if os.path.isfile(f"{input_folder1}/tagLengthDistribution.txt"):
                shell(f"rm {input_folder1}/tagLengthDistribution.txt")


        def count_reads_atac(input_folder1, input_folder2, output_folder, genome, namefields, separator):
            sanity_check_dir(logfile, input_folder1,  master_config['input_file_types'][master_config['countreads_rule_num']-1][0])
            sanity_check_dir(logfile, input_folder2,  master_config['input_file_types'][master_config['countreads_rule_num']-1][1])

            tagdir_files = glob.glob(os.path.join(input_folder, '*tagDir.tar.gz'))
            if len(tagdir_files) != 0:
                log_it(logfile, "Unpacking HOMER tagDir tar balls...")

                for file in tagdir_files:
                    basename = os.path.basename(file)
                    tag_dir = os.path.join(input_folder1, basename.replace(".tar.gz", ""))
                    tag_dir_short = basename.replace(".tar.gz", "")

                    # unpack the tar.gz file
                    log_it(logfile, f"Unpacking {basename} in {input_folder1}")
                    shell(f"""
                        mkdir -p {tag_dir} && \
                        cd {input_folder1} && \
                        tar --strip-components=1 -xzf {basename} -C {tag_dir_short}
                        """)

            shell(f"ls -d {input_folder1}/*tagDir/ > TAGDIRlist.txt") # Make a list of tag dirs

            # Convert peaks
            log_it(logfile, "Converting BED peak file to HOMER POS peak file...")
            shell("bed2pos.pl peak_calling/all_groups.merged_peaks.bed > all_groups.merged_peaks.pos")

            ## Run HOMER's analyze_repeats.pl to count tags
            log_it(logfile, "Counting reads with analyzeRepeats.pl...")
            shell(f"""analyzeRepeats.pl all_groups.merged_peaks.pos {genome} -dfile TAGDIRlist.txt -noadj > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt""")

            # Clean the table
            shell(f"""sed '1d' {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt | cut -f2-4,9- | sort -k1,1 -k2,2n -k3,3n | sed 's/ \\t/_/;s/\\t/_/' | sed 's/ \+/\t/g' > clean.tmp""")

            #cat TAGDIRlist.txt | xargs -l basename | cut -f {namefields} -d '{separator}' | awk 'BEGIN{{ORS="\t"}}{{print $0}}' > clean.header.tmp
            shell(f"""{config['SCRIPT_DIR']}/generate_header.sh {params.namefields} '{params.separator}'""")


            shell(f"""sed "1iPeak\\t$(cat clean.header.tmp)" clean.tmp > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt""")


            shell("rm TAGDIRlist.txt clean.header.tmp clean.tmp all_groups.merged_peaks.pos")

            for tag_dir in tagdir_files:
                tag_dir_folder = tag_dir.replace(".tar.gz", "")
                # Remove the uncompressed tag directory
                log_it(logfile, f"Removing uncompressed tag directory {tag_dir_folder}")
                shell(f"rm -r {tag_dir_folder}")
            if os.path.isfile(f"{input_folder1}/genome.tags.tsv"):
                shell(f"rm {input_folder1}/genome.tags.tsv")
            if os.path.isfile(f"{input_folder1}/tagAutocorrelation.txt"):
                shell(f"rm {input_folder1}/tagAutocorrelation.txt")
            if os.path.isfile(f"{input_folder1}/tagCountDistribution.txt"):
                shell(f"rm {input_folder1}/tagCountDistribution.txt")
            if os.path.isfile(f"{input_folder1}/tagInfo.txt"):
                shell(f"rm {input_folder1}/tagInfo.txt")
            if os.path.isfile(f"{input_folder1}/tagLengthDistribution.txt"):
                shell(f"rm {input_folder1}/tagLengthDistribution.txt")

        if params.thetype == "RNA":
            count_reads_rna(params.inputfolder1, params.outputfolder, params.genome, params.namefields, params.separator)
        elif params.thetype == "ATAC":
            count_reads_atac(params.inputfolder1, params.inputfolder2, params.outputfolder, params.genome, params.namefields, params.separator)
        else:
            log_it(logfile, "For ChIP experiments, first determine optimal peak caller settings, then manually run run_quant_peaks.sh and then continue with the next step!")

        shell(f"""echo "necessity file for count reads. can delete this." > {params.outputfolder}/extra_12.tmp""")
