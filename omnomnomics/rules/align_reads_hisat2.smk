#Rule 3 option HISAT2

## Omnomnomics Snake Rule  ##
import os

rule run_hisat2:
    input:
        fastq1=f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R1.fastq.gz",
        fastq2=f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R2.fastq.gz" if config["PAIRED"] else None
    output:
        bam=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample}}.bam",
        stats=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample}}.HISAT2_stats.txt"
    params:
        genome_path=os.path.join(f"{config['OMNOM_HOME']}", "genomes", "HISAT2", f"{config['THEGENOME']}"),
        keepunpaired=config.get("KEEPUNPAIRED", "0"),
        seq_type=config["THETYPE"],
        inputfolder = master_config['input_folders'][master_config['map_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['map_rule_num']-1]
    threads:
        10
    resources:
        mem_mb = (10*4000)
    run:
        log_it(logfile, "Mapping reads...", f"EXECUTING STEP {master_config['map_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        # Load the samtools module
        os.system("module load samtools")

        def run_hisat2(seq_type, paired, threads, genome_path, fastq1, fastq2, keepunpaired, inputfolder, outputfolder, sample):
            if seq_type == "RNA":
                if config["PAIRED"]:
                    log_it(logfile, f"Running HISAT2 in Paired End mode on RNA data.")
                    MYNAME = os.path.basename(input.fastq1)
                    logIt(logfile, f"Launching: {MYNAME}...")
                    if keepunpaired == "1":
                        shell(f"""
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --dta-cufflinks \
                        --un-gz "{outputfolder}/{sample}.unpaired.unaligned.bam" --al-gz "{outputfolder}/{sample}.unpaired.aligned.bam" \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                        """)
                    else:
                        shell(f"""
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --dta-cufflinks \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                        """)
                else:
                    log_it(logfile, f"Running HISAT2 in Single End mode on RNA data.")
                    shell(f"""
                    hisat2 -p {threads} -x {genome_path} \
                    -U {fastq1} --mm --add-chrname --new-summary --dta-cufflinks \
                    | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                    """)
            else: #if ChIP- or ATAC data
                if config["PAIRED"]:
                    log_it(logfile, f"Running HISAT2 in Paired End mode  on ChIP or ATAC data.")
                    MYNAME = os.path.basename(input.fastq1)
                    logIt(logfile, f"Launching: {MYNAME}...")
                    if keepunpaired == "1":
                        shell(f"""
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --no-spliced-alignment \
                        --un-gz "BAM/{sample}.unpaired.unaligned.bam" --al-gz "BAM/{sample}.unpaired.aligned.bam"\
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                        """)
                    else:
                        shell(f"""
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --no-spliced-alignment \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                        """)
                else:
                    log_it(logfile, f"Running HISAT2 in Single End mode  on ChIP or ATAC data.")
                    shell(f"""
                    hisat2 -p {threads} -x {genome_path} \
                    -U {fastq1} --mm --add-chrname --new-summary --no-spliced-alignment \
                    | samtools view -b - 1> "{outputfolder}"/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                    """)

        # Call the function with parameters
        run_hisat2(params.seq_type, params.paired, threads, params.genome_path, input.fastq1, input.fastq2, params.keepunpaired, params.inputfolder, params.outputfolder, wildcards.sample)