#Rule 3 option HISAT2

## Omnomnomics Snake Rule  ##
import os

# def get_hisat2_input(wildcards):
#     input_folder = master_config['input_folders'][master_config['map_rule_num']-1]
#     sample = wildcards.sample
#     semaphore_path = os.path.join(experiment_dir, "omnomnomics.semaphore")

#     with open(semaphore_path, 'r') as file:
#         lines = file.readlines()
#         if len(lines) < 2 or lines[1].strip() != "1":
#             suffix = ''
            
#         elif config['THEMAPTOOL'] == 'skewer':
#             suffix = '_Skewer'
#         elif config['THEMAPTOOL'] == 'trimmomatic':
#             suffix = '_Trimmomatic'

#     if config['PAIRED']:
#         fastq1 = f"{input_folder}/{sample}_R1{suffix}.trimmed.fastq.gz"
#         fastq2 = f"{input_folder}/{sample}_R2{suffix}.trimmed.fastq.gz"
#         return fastq1, fastq2
#     else:
#         fastq = f"{input_folder}/{sample}{suffix}.trimmed.fastq.gz"
#         return fastq

rule run_hisat2:
    input:
        # get_hisat2_input
        # fastq1=f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R1_Skewer.fastq.gz" if config["PAIRED"] else None,
        # fastq2=f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R2.fastq.gz" if config["PAIRED"] else None,
        # fastq3 = f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample}}_R2.fastq.gz" if not config["PAIRED"] else None
        trimmed_fastq1= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample1}}_R1.trimmed.fastq.gz" if config["PAIRED"] else f"{master_config['output_folders'][master_config['trim_rule_num']-1]}/{{sample}}.trimmed.fastq.gz",
        trimmed_fastq2= f"{master_config['input_folders'][master_config['map_rule_num']-1]}/{{sample1}}_R2.trimmed.fastq.gz" if config['PAIRED'] else None
    output:
        bam=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample1}}.bam",
        stats=f"{master_config['output_folders'][master_config['map_rule_num']-1]}/{{sample1}}.HISAT2_stats.txt"
        #could add the 2 extra  output files if keepunpaired = 1, but not necesarry because will be made then automatically
    params:
        genome_path=os.path.join(f"{config['OMNOM_HOME']}", "genomes", "STAR", f"{config['THEGENOME']}"), ##CHANGE TO HISAT2
        keepunpaired=config.get("KEEPUNPAIRED", "0"),
        seq_type=config["THETYPE"],
        inputfolder = master_config['input_folders'][master_config['map_rule_num']-1],
        outputfolder = master_config['output_folders'][master_config['map_rule_num']-1]
    threads:
        Threads_Per_Rule['3']
    resources:
        mem_mb = Memory_Per_Rule['3']
    run:
        log_it(logfile, "Mapping reads...", f"EXECUTING STEP {master_config['map_rule_num']}")
        log_it(logfile, f"Input folder: {params.inputfolder}")
        log_it(logfile, f"Output folder: {params.outputfolder}")

        hisat2_version = subprocess.check_output(["hisat2", "--version"])
        log_it(logfile, "\n"+hisat2_version.decode("utf-8"), "HISAT2 VERSION")
        print(hisat2_version.decode("utf-8"))
        sanity_check_dir(logfile, params.inputfolder,  master_config['input_file_types'][master_config['map_rule_num']-1])

        def run_hisat2(seq_type, threads, genome_path, fastq1, fastq2, keepunpaired, inputfolder, outputfolder, sample):
            if seq_type == "RNA":
                if config["PAIRED"]:
                    log_it(logfile, f"Running HISAT2 in Paired End mode on RNA data.")
                    MYNAME = os.path.basename(fastq1)
                    log_it(logfile, f"Launching: {MYNAME}...")
                    if keepunpaired == "1":
                        shell(f"""
                        module load samtools && \
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --dta-cufflinks \
                        --un-gz "{outputfolder}/{sample}.unpaired.unaligned.bam" --al-gz "{outputfolder}/{sample}.unpaired.aligned.bam" \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                        """)
                    else:
                        shell(f"""
                        module load samtools && \
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --dta-cufflinks \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                        """)
                else:
                    log_it(logfile, f"Running HISAT2 in Single End mode on RNA data.")
                    shell(f"""
                    module load samtools && \
                    hisat2 -p {threads} -x {genome_path} \
                    -U {fastq1} --mm --add-chrname --new-summary --dta-cufflinks \
                    | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                    """)
            else: #if ChIP- or ATAC data
                if config["PAIRED"]:
                    log_it(logfile, f"Running HISAT2 in Paired End mode  on ChIP or ATAC data.")
                    MYNAME = os.path.basename(fastq1)
                    log_it(logfile, f"Launching: {MYNAME}...")
                    if keepunpaired == "1":
                        shell(f"""
                        module load samtools && \
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --no-spliced-alignment \
                        --un-gz "BAM/{sample}.unpaired.unaligned.bam" --al-gz "BAM/{sample}.unpaired.aligned.bam"\
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                        """)
                    else:
                        shell(f"""
                        module load samtools && \
                        hisat2 -p {threads} -x {genome_path} \
                        -1 {fastq1} -2 {fastq2} --mm --add-chrname --new-summary --no-spliced-alignment \
                        | samtools view -b - 1> "{outputfolder}/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                        """)
                else:
                    log_it(logfile, f"Running HISAT2 in Single End mode  on ChIP or ATAC data.")
                    shell(f"""
                    module load samtools && \
                    hisat2 -p {threads} -x {genome_path} \
                    -U {fastq1} --mm --add-chrname --new-summary --no-spliced-alignment \
                    | samtools view -b - 1> "{outputfolder}"/{sample}.bam" 2> "{outputfolder}/{sample}.HISAT2_stats.txt"
                    """)

        # Call the function with parameters
        run_hisat2(params.seq_type, threads, params.genome_path, input.trimmed_fastq1, input.trimmed_fastq2, params.keepunpaired, params.inputfolder, params.outputfolder, wildcards.sample1)