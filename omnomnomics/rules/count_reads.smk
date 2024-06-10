# Rule 12: Count Reads

## Omnomnomics Snake Rule  ##
import os

rule count_reads:
    input:
        tagDirs_tar_gz = glob.glob(f"{master_config['input_folders'][master_config['countreads_rule_num']-1]}/{{sample}}.HOMER_tagDir.tar.gz")
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
    run:
        def count_reads_rna(input_folder, output_folder, genome, namefields, separator):
            log_it(logfile, f"Input folder: {input_folder}")
            log_it(logfile, f"Output folder: {output_folder}")

            version = subprocess.check_output(["perl", os.path.join(OMNOM_HOME, "bin", "homer", "configureHomer.pl"), "-list", "2>", "/dev/null", "|", "grep", "homer"])
            log_it(logfile, "\n"+version.decode("utf-8"), "VERSION")
            print(version.decode("utf-8"))
            sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['countreads_rule_num']-1])

            if shell(f"ls {input_folder}/*tagDir.tar.gz | wc -l").strip() != "0":
                log_it(logfile, "Unpacking HOMER tagDir tar balls...")
                shell(f"""
                    cd {input_folder} &&
                    for TAGDIR in *tagDir.tar.gz; do
                        tar --strip-components=1 -xzf $TAGDIR &
                    done
                    wait &&
                    cd ..
                """)


            shell(f"ls -d {input_folder}/*tagDir/ > TAGDIRlist.txt")
            log_it(logfile, "Counting reads with analyzeRepeats.pl...")
            shell(f"""
                analyzeRepeats.pl rna {genome} -dfile TAGDIRlist.txt -count exons -noadj > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt
            """)

            shell(f"""
                sed '1d' {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt | cut -f1,9- | sort -k1,1 | sed 's/ \\+/\\t/g' > clean.tmp
            """)

            shell(f"""
                cat TAGDIRlist.txt | xargs -l basename | cut -f {namefields} -d '{separator}' | awk 'BEGIN{{ORS="\\t"}}{{print $0}}' > clean.header.tmp
            """)


            shell(f"""
                sed "1iRefSeq_ID\t$(cat clean.header.tmp)" clean.tmp > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt
            """)


            shell("rm TAGDIRlist.txt clean.header.tmp clean.tmp")


        def count_reads_atac(input_folder, output_folder, genome, namefields, separator):
            log_it(logfile, f"Input folder: {input_folder}")
            log_it(logfile, f"Output folder: {output_folder}")
            
            version = subprocess.check_output(["perl", os.path.join(OMNOM_HOME, "bin", "homer", "configureHomer.pl"), "-list", "2>", "/dev/null", "|", "grep", "homer"])
            log_it(logfile, "\n"+version.decode("utf-8"), "VERSION")
            print(version.decode("utf-8"))
            sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['countreads_rule_num']-1])

            if shell(f"ls {input_folder}/*tagDir.tar.gz | wc -l").strip() != "0":
                log_it(logfile, "Unpacking HOMER tagDir tar balls...")
                shell(f"""
                    cd {input_folder} &&
                    for TAGDIR in *tagDir.tar.gz; do
                        tar --strip-components=1 -xzf $TAGDIR &
                    done
                    wait &&
                    cd ..
                """)
                shell("cd ..")


            sanity_check_dir("peak_calling", ".merged_peaks.bed")
            shell(f"ls -d {input_folder}/*tagDir/ > TAGDIRlist.txt")
            log_it(logfile, "Converting BED peak file to HOMER POS peak file...")
            shell("bed2pos.pl peak_calling/all_groups.merged_peaks.bed > all_groups.merged_peaks.pos")
            log_it(logfile, f"Peaks: {shell('wc -l all_groups.merged_peaks.pos').strip()}")
            log_it(logfile, "Counting reads with analyzeRepeats.pl...")
            shell(f"""
                analyzeRepeats.pl all_groups.merged_peaks.pos {genome} -dfile TAGDIRlist.txt -noadj > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt
            """)


            shell(f"""
                sed '1d' {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt | cut -f2-4,9- | sort -k1,1 -k2,2n -k3,3n | sed 's/\\t/_/;s/\\t/_/' | sed 's/ \\+/\\t/g' > clean.tmp
            """)


            shell(f"""
                cat TAGDIRlist.txt | xargs -l basename | cut -f {namefields} -d '{separator}' | awk 'BEGIN{{ORS="\\t"}}{{print $0}}' > clean.header.tmp
            """)


            shell(f"""
                sed "1iPeak\t$(cat clean.header.tmp)" clean.tmp > {output_folder}/{os.path.basename(params.experiment_dir)}.raw_read_quant.table.txt
            """)


            shell("rm TAGDIRlist.txt clean.header.tmp clean.tmp all_groups.merged_peaks.pos")


        log_it(logfile, "Counting reads...", f"EXECUTING STEP {master_config['countreads']}")


        if params.thetype == "RNA":
            count_reads_rna(params.inputfolder, params.outputfolder, params.genome, params.namefields, params.separator)
        elif params.thetype == "ATAC":
            count_reads_atac(params.inputfolder, params.outputfolder, params.genome, params.namefields, params.separator)
        else:
            log_it(logfile, "For ChIP experiments, first determine optimal peak caller settings, then manually run run_quant_peaks.sh to continue!")