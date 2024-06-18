# Rule 12: Count Reads

## Omnomnomics Snake Rule  ##
import os

def input_function(wilcards):
    input_files = []
    if config['THETYPE'] != "CHIP":
        if config['THETYPE'] == "ATAC":
            input_files.append(f"{master_config['input_folders'][master_config['countreads_rule_num']-1][1]}/extra.tmp")
            for sample in samples2:
                input_files.append(f"{master_config['input_folders'][master_config['countreads_rule_num']-1][0]}/{sample}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz")
    else:
        for sample in samples2:
            input_files.append(f"{master_config['input_folders'][master_config['countreads_rule_num']-1][0]}/{sample}.filtered.HOMER_tagDir.tar.gz")

rule count_reads:
    input:
        #input_function, #possibly leave this out? or keep it for readability? and makes sure that there is input, although there is also a sanitycheck
        #expand(f"{master_config['input_folders'][master_config['callpeaks_rule_num']-1][0]}/{{sample}}.extra_8.tmp", sample = samples2) if 8 in themode else None,
        #f"{master_config['input_folders'][master_config['countreads_rule_num']-1][1]}/extra_11.tmp" if config['THETYPE'] != "RNA" else None
        #expand(f"{master_config['input_folders'][master_config['countreads_rule_num']-1]}/{{sample}}.HOMER_tagDir.tar.gz")
        extra_input_rule_12_1 = expand(f"{master_config['input_folders'][master_config['countreads_rule_num']-1][0]}/{{sample}}.extra_8.tmp", sample = samples2) if 8 in themode else [],
        extra_input_rule_12_2 = expand(f"{master_config['input_folders'][master_config['countreads_rule_num']-1][1]}/{{sample}}.extra_11.tmp", sample = samples2) if 11 in themode and config['THETYPE'] == "ATAC" else []
    output:
        f"{master_config['output_folders'][master_config['countreads_rule_num']-1]}/{os.path.basename(config['EXPERIMENT_DIR'])}.raw_read_quant.table.txt"
    params:
        thetype=config['THETYPE'],  
        genome=config['THEGENOME'],  
        experiment_dir=config['EXPERIMENT_DIR'], 
        namefields=config['NAMEFIELDS'],  
        separator=config['THESEPARATOR'],
        inputfolder = master_config['input_folders'][master_config['countreads_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['countreads_rule_num']-1]
    threads:
        Threads_Per_Rule['12']
    resources:
        mem_mb = Memory_Per_Rule['12']
    benchmark:
        f"{master_config['output_folders'][master_config['countreads_rule_num']-1]}/counts_reads_benchmark.tsv"
    run:
        log_it(logfile, "Counting Reads...", f"EXECUTING STEP {master_config['countreads_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")
        
        path = os.path.join(OMNOM_HOME, "bin", "homer", "configureHomer.pl")
        version = subprocess.check_output("perl {path} -list 2> /dev/null | grep homer",  shell=True, executable='/bin/bash')
        log_it(logfile, "\n"+version.decode("utf-8"), "VERSION")
        print(version.decode("utf-8"))
        sanity_check_dir(logfile, input_folder,  master_config['input_file_types'][master_config['countreads_rule_num']-1])

        def count_reads_rna(input_folder, output_folder, genome, namefields, separator):  
            tagdir_files = glob.glob(os.path.join(input_folder, '*tagDir.tar.gz'))
            if len(tagdir_files) != 0:
                log_it(logfile, "Unpacking HOMER tagDir tar balls...")
                shell(f"""
                    cd {input_folder} &&
                    for TAGDIR in *tagDir.tar.gz; do
                        tar --strip-components=1 -xzf $TAGDIR &
                    done
                """)

            shell(f"ls -d {input_folder}/*tagDir/ > TAGDIRlist.txt") # Make a list of tag dirs
            log_it(logfile, "Counting reads with analyzeRepeats.pl...")
            shell(f"""analyzeRepeats.pl rna {genome} -dfile TAGDIRlist.txt -count exons -noadj > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt""")

            shell(f"""sed '1d' {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt | cut -f1,9- | sort -k1,1 | sed 's/ \\+/\\t/g' > clean.tmp""")
            
            #cat TAGDIRlist.txt | xargs -l basename | cut -f {namefields} -d '{separator}' | awk 'BEGIN{{ORS="\t"}}{{print $0}}' > clean.header.tmp
            shell(f"""./generate_header.sh {params.namefields} '{params.separator}' """)


            shell(f"""sed "1iRefSeq_ID\\t$(cat clean.header.tmp)" clean.tmp > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt""")


            shell("rm TAGDIRlist.txt clean.header.tmp clean.tmp")


        def count_reads_atac(input_folder, output_folder, genome, namefields, separator):
            tagdir_files = glob.glob(os.path.join(input_folder, '*tagDir.tar.gz'))
            if len(tagdir_files) != 0:
                log_it(logfile, "Unpacking HOMER tagDir tar balls...")
                shell(f"""
                    cd {input_folder} &&
                    for TAGDIR in *tagDir.tar.gz; do
                        tar --strip-components=1 -xzf $TAGDIR &
                    done
                """)

            sanity_check_dir(logfile, "peak_calling", ".merged_peaks.bed")

            shell(f"ls -d {input_folder}/*tagDir/ > TAGDIRlist.txt") # Make a list of tag dirs

            # Convert peaks
            log_it(logfile, "Converting BED peak file to HOMER POS peak file...")
            shell("bed2pos.pl peak_calling/all_groups.merged_peaks.bed > all_groups.merged_peaks.pos")

            ## Run HOMER's analyze_repeats.pl to count tags
            log_it(logfile, "Counting reads with analyzeRepeats.pl...")
            shell(f"""analyzeRepeats.pl all_groups.merged_peaks.pos {genome} -dfile TAGDIRlist.txt -noadj > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt""")

            # Clean the table
            shell(f"""sed '1d' {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt | cut -f2-4,9- | sort -k1,1 -k2,2n -k3,3n | sed 's/ \\t/_/;s/\\t/_/' | sed 's/ \+/\t/g' > clean.tmp""")

            #cat TAGDIRlist.txt | xargs -l basename | cut -f {namefields} -d '{separator}' | awk 'BEGIN{{ORS="\t"}}{{print $0}}' > clean.header.tmp
            shell(f"""./generate_header.sh {params.namefields} '{params.separator}'""")


            shell(f"""sed "1iPeak\\t$(cat clean.header.tmp)" clean.tmp > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt""")


            shell("rm TAGDIRlist.txt clean.header.tmp clean.tmp all_groups.merged_peaks.pos")

        if params.thetype == "RNA":
            count_reads_rna(params.inputfolder, params.outputfolder, params.genome, params.namefields, params.separator)
        elif params.thetype == "ATAC":
            count_reads_atac(params.inputfolder, params.outputfolder, params.genome, params.namefields, params.separator)
        else:
            log_it(logfile, "For ChIP experiments, first determine optimal peak caller settings, then manually run run_quant_peaks.sh to continue!")