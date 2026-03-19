#Rule 3 option HISAT2

## Omnomnomics Snake Rule ##
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os

rule run_hisat2:
    input:
        trimmed_fastq1= f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample3}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{experiment_dir}/{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample3}}.trimmed.fastq.gz",
        trimmed_fastq2= f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample3}}_R2.trimmed.fastq.gz" if config['PAIRED'] else []
    output:
        bam=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.bam",
        stats=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.HISAT2_stats.txt",
        extra=f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample3}}.extra_3.tmp"
        #two extra  output files are create if keepunpaired = 1, but not necesarry to specify because will be made then automatically
    params:
        genome_path=os.path.join(config['HISAT2_GENOME_DIR'], f"{config['THEGENOME']}"),
        keepunpaired=config.get("KEEPUNPAIRED", "0"),
        seq_type=config["THETYPE"],
        inputfolder = f"{experiment_dir}/{master_config['input_folders'][master_config['map_rule_num']-1]}",
        outputfolder = f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}"
    threads:
        Threads_Per_Rule['3']
    resources:
        mem_mb = Memory_Per_Rule['3']
    benchmark:
        f"{experiment_dir}/{master_config['output_folders'][master_config['map_rule_num']-1]}/benchmarks/{{sample3}}_hisat2_benchmark.tsv"
    run:
        log_it(logfile, "Mapping reads...", f"EXECUTING STEP {master_config['map_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        hisat2_version = subprocess.check_output(["hisat2", "--version"])
        log_it(logfile, "\n"+hisat2_version.decode("utf-8"), "HISAT2 VERSION")
        
        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['map_rule_num']-1])

        def run_hisat2(seq_type, threads, genome_path, fastq1, fastq2, keepunpaired, inputfolder, outputfolder, sample):
            if seq_type == "RNA":
                if config["PAIRED"]:
                    log_it(logfile, f"Running HISAT2 in Paired End mode on RNA data.")
                    MYNAME = os.path.basename(fastq1)
                    log_it(logfile, f"Launching: {MYNAME}...")
                    if keepunpaired:
                        shell(f"""
                        module load samtools && \
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --dta-cufflinks \
                        --summary-file "{outputfolder}/{sample}.HISAT2_stats.txt" \
                        --un-gz "{outputfolder}/{sample}.unpaired.unaligned.bam" --al-gz "{outputfolder}/{sample}.unpaired.aligned.bam" \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 
                        """)
                    else:
                        shell(f"""
                        module load samtools && \
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --dta-cufflinks \
                        --summary-file "{outputfolder}/{sample}.HISAT2_stats.txt" \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 
                        """) 
                else:
                    log_it(logfile, f"Running HISAT2 in Single End mode on RNA data.")
                    shell(f"""
                    module load samtools && \
                    hisat2 -p {threads} -x {genome_path} \
                    -U {fastq1} --mm --add-chrname --new-summary --dta-cufflinks \
                    --summary-file "{outputfolder}/{sample}.HISAT2_stats.txt" \
                    | samtools view -b - 1> "{outputfolder}/{sample}.bam" 
                    """)
            else: #if ChIP- or ATAC data
                if config["PAIRED"]:
                    log_it(logfile, f"Running HISAT2 in Paired End mode  on ChIP or ATAC data.")
                    MYNAME = os.path.basename(fastq1)
                    log_it(logfile, f"Launching: {MYNAME}...")
                    if keepunpaired:
                        shell(f"""
                        module load samtools && \
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --no-spliced-alignment \
                        --un-gz "BAM/{sample}.unpaired.unaligned.bam" --al-gz "BAM/{sample}.unpaired.aligned.bam"\
                        --summary-file "{outputfolder}/{sample}.HISAT2_stats.txt" \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 
                        """) 
                    else:
                        shell(f"""
                        module load samtools && \
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --no-spliced-alignment \
                        --summary-file "{outputfolder}/{sample}.HISAT2_stats.txt" \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 
                        """)
                else:
                    log_it(logfile, f"Running HISAT2 in Single End mode  on ChIP or ATAC data.")
                    shell(f"""
                    module load samtools && \
                    hisat2 -p {threads} -x {genome_path} \
                    -U {fastq1} --mm --add-chrname --new-summary --no-spliced-alignment \
                    --summary-file "{outputfolder}/{sample}.HISAT2_stats.txt" \
                    | samtools view -b - 1> "{outputfolder}/{sample}.bam" 
                    """) 
            shell(f"""echo "necessity file for aligners. can delete this." > {outputfolder}/{sample}.extra_3.tmp""")
        # Call the function with parameters
        run_hisat2(params.seq_type, threads, params.genome_path, input.trimmed_fastq1, input.trimmed_fastq2, params.keepunpaired, params.inputfolder, params.outputfolder, wildcards.sample3)
