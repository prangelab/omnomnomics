##############################################
### NGS PIPELINE FOR RNA, ATAC & CHIP DATA ###
##############################################
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
#Main Snakefile
import os
import sys
import time
import datetime
import yaml
import shutil
import glob
import random
import re


# Load the configuration file from command line arguments
configfile: config['config_file']


# Set global variables from the configuration file
workflow_root = config['WORKFLOW_ROOT']
experiment_dir = config["EXPERIMENT_DIR"]
run_date = config["RUNDATE"]
logfile = os.path.join(experiment_dir, "run_logs", f"omnomnomics.run.{run_date}.log")
is_worker_job = "--target-jobs" in sys.argv

# Function to log messages
def log_it(logfile, message, heading=None):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(logfile, 'a') as log:
        if heading: #Check if heading
            log.write("\n{}\n\n".format(heading))
        log.write(f"{timestamp}: {message}\n") #write to log
    print(f"{timestamp}: {message}")

onstart:
    # Upon start, log the start time of the pipeline
    global start_time
    start_time = time.time()
    if not is_worker_job:
        log_it(logfile, "Pipeline started at {}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")), "PIPELINE START TIME")

onerror:
    # Upon error, log the error time of the pipeline and the elapsed time
    end_time = time.time()
    elapsed_time = end_time - start_time
    if not is_worker_job:
        log_it(logfile, "Pipeline failed at {}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")), "PIPELINE RUN TIME")
        log_it(logfile, "Total elapsed time: {:.2f} seconds\n".format(elapsed_time))

# Function for sanity check on directory
def sanity_check_dir(logfile, input_directory, file_ext):
    #check if input directiory exists
    if not os.path.isdir(input_directory):
        log_it(logfile, f"{file_ext} files should be contained in a {input_directory} folder inside your <EXPERIMENT_DIR> ({experiment_dir})! Aborting...", "ERROR")
        log_it(logfile, "Aborting...")
        print(f"{file_ext} files should be contained in a {input_directory} folder inside your <EXPERIMENT_DIR> ({experiment_dir})! Aborting...")
        sys.exit(1)
    #check if input directory contains input files
    if len([f for f in os.listdir(input_directory) if f.endswith(file_ext)]) == 0:
        log_it(logfile, f"The {input_directory} folder inside your <EXPERIMENT_DIR> ({experiment_dir}) does not contain any {file_ext} files! Aborting...", "ERROR")
        log_it(logfile, "Aborting...")
        print(f"The {input_directory} folder inside your <EXPERIMENT_DIR> ({experiment_dir}) does not contain any {file_ext} files! Aborting...")
        sys.exit(1)
    log_it(logfile, "Sanity check completed", "SANITY CHECK")

# Hold our horses for a little while to ensure all files are up to date
time.sleep(0.1)

if not is_worker_job:
    log_it(logfile,"JOB DISPATCHED!", "INITIALIZATION")
##---------------------------------------------------------------------------------------------------------------
## Read defaults from config file or die
##---------------------------------------------------------------------------------------------------------------
# Function to check if the configuration file contains the proper header
def check_config_file_header(logfile, config_file_path, expected_header):
    with open(config_file_path, 'r') as config_file:
        lines = config_file.readlines()
        if len(lines) < 3 or expected_header not in lines[2]:
            log_it(logfile, "Master config file  is malformed! Expected: ## Omnomnomics pipeline config ##. Aborting...", "ERROR")
            print("Master config file  is malformed! Expected: ## Omnomnomics pipeline config ##. Aborting...")
            sys.exit(1)

# Function to load and validate the YAML configuration file
def load_and_validate_yaml(logfile, config_file_path, expected_header):
    # Check for the expected header
    check_config_file_header(logfile, config_file_path, expected_header)

    # Load the YAML content
    with open(config_file_path, 'r') as file:
        try:
            config_content = yaml.safe_load(file)
        except yaml.YAMLError as exc:
            log_it(logfile, f"Error parsing YAML file: {exc}", "ERROR")
            sys.exit(1)
    return config_content

def merge_configs(base_config, override_config):
    merged_config = dict(base_config)
    merged_config.update(override_config)
    return merged_config

workflow_config_file = config['WORKFLOW_CONFIG_FILE']
site_config_file = config['SITE_CONFIG_FILE']

if not os.path.isfile(workflow_config_file):
    log_it(logfile, "Workflow config file does not exist! Aborting...", "ERROR")
    print(f"Workflow config file '{workflow_config_file}' does not exist. Please make sure it exists. Aborting...")
    sys.exit(1)

if not os.path.isfile(site_config_file):
    log_it(logfile, "Site config file does not exist! Aborting...", "ERROR")
    print(f"Site config file '{site_config_file}' does not exist. Please make sure it exists. Aborting...")
    sys.exit(1)

workflow_config = load_and_validate_yaml(logfile, workflow_config_file, "## Omnomnomics pipeline config ##")
site_config = load_and_validate_yaml(logfile, site_config_file, "## Omnomnomics pipeline config ##")
master_config = merge_configs(workflow_config, site_config)

themode = config['THEMODE']

##---------------------------------------------------------------------------------------------------------------
## Final housekeeping
##---------------------------------------------------------------------------------------------------------------
def final_housekeeping(logfile, thecoltable, experiment_dir):
    # Move our camp to the experiment directory
    os.chdir(experiment_dir)
    if not is_worker_job:
        log_it(logfile, f"Changed directory to experiment dir: {experiment_dir}" )
        print(f"Changed directory to experiment dir: {experiment_dir}")
    return

final_housekeeping(logfile, config['THECOLTABLE'], experiment_dir)

##---------------------------------------------------------------------------------------------------------------
## Report some basic stats
##---------------------------------------------------------------------------------------------------------------
if not is_worker_job:
    log_it(logfile, f"Experiment dir: {experiment_dir}")
    log_it(logfile, f"Run date: {run_date}")
    log_it(logfile, f"Log file: {logfile}")
    log_it(logfile, f"Files: {config['NUMFILES']}")
    log_it(logfile, f"Pairs: {config['NUMPAIRS']}")
    log_it(logfile, f"Paired run: {config['PAIRED']}")
    log_it(logfile, f"Type: {config['THETYPE']}")
    log_it(logfile, f"Genome: {config['THEGENOME']}")
    log_it(logfile, f"Mode: {' '.join(map(str, themode))}")
    log_it(logfile, f"Mode min: {min(themode)}")
    log_it(logfile, f"Mode max: {max(themode)}")
    log_it(logfile, f"Trim tool: {config['THETRIMTOOL']}")
    log_it(logfile, f"Map tool: {config['THEMAPTOOL']}")

    log_it(logfile, f"Design formula: {config['MYFORMULA']}", "READ COUNTING SETTINGS")
    log_it(logfile, f"Metadata file: {config['MYMETADATA']}")

    log_it(logfile, f"Input file (MACS3): {config['INPUT']}", "PEAK CALLING SETTINGS")
    log_it(logfile, f"Input file (HOMER): {config['HOMERINPUT']}")
    log_it(logfile, f"Broad peaks: {config['BROAD']}")
    log_it(logfile, f"Peak style: {config['THESTYLE']}")
    log_it(logfile, f"HOMER peaksize: {config['HOMERSIZE']}")
    log_it(logfile, f"HOMER minDist: {config['HOMERMINDIST']}")

    log_it(logfile, f"{config['THEHEAPINIT']} HEAP init", "JAVA MEMORY SETTINGS")
    log_it(logfile, f"{config['THEMEM']} memory per sample")

    log_it(logfile, f"Hub type field: {config['THETYPEFIELD']}", "TRACKHUB SETTINGS")
    log_it(logfile, f"Hub name field: {config['NAMEFIELDS']}")
    log_it(logfile, f"Hub col field: {config['THECOLFIELD']}")
    log_it(logfile, f"Hub separator: {config['THESEPARATOR']}")
    log_it(logfile, f"Hub appendix: {config['THEAPPENDIX']}")
    log_it(logfile, f"Hub overlay: {config['THEOVERLAY']}")
    log_it(logfile, f"Hub coldata folder: {config['THECOLORDATAFOLDER']}")
    log_it(logfile, f"Hub color table: {config['THECOLTABLE']}")
    log_it(logfile, f"Hub mail: {config['THEHUBMAIL']}")

##--------------------------------------------------------------------------------------------------------------
# Obtain Sample Names and Number of Samples
##--------------------------------------------------------------------------------------------------------------
all_outputs = []
routines = []

THEMODERANGEMIN = config['THEMODERANGEMIN']

input_folder = master_config['input_folders'][THEMODERANGEMIN-1]
input_file_type =  master_config['input_file_types'][THEMODERANGEMIN-1]
if THEMODERANGEMIN == 10:
    if config['THETYPE'] == "RNA":
        input_file_type = input_file_type[0]
    else:
        input_file_type = input_file_type[1]
if THEMODERANGEMIN == 11:
    input_file_type = input_file_type[0]
    input_folder = input_folder[0]
if THEMODERANGEMIN == 12:
    input_file_type = input_file_type[0]
    input_folder = input_folder[0]

input_pattern = os.path.join(input_folder, f"*{input_file_type}")
input_files = glob.glob(input_pattern)

FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
FASTQ_READ_SUFFIX_RE = re.compile(r'_(?:R)?[12](?:_[0-9]{3})?$')


def strip_fastq_read_suffix(sample_name):
    return FASTQ_READ_SUFFIX_RE.sub('', sample_name)


def normalize_sample_name(file_path, file_type):
    sample_name = os.path.basename(file_path).replace(file_type, "")
    if "fastq" in file_type or file_type.endswith(".fq.gz") or file_type.endswith(".fq"):
        sample_name = strip_fastq_read_suffix(sample_name)
    sample_name = sample_name.replace(".filtered", "")
    sample_name = sample_name.replace(".sorted.dups_marked", "")
    return sample_name


def fastq_candidate_names(sample, read_label):
    read_number = read_label.replace("R", "")
    suffixes = (
        f"_{read_label}_001",
        f"_{read_label}",
        f"_{read_number}_001",
        f"_{read_number}",
    )
    for suffix in suffixes:
        for extension in FASTQ_EXTENSIONS:
            yield f"{sample}{suffix}{extension}"


def resolve_fastq_input(sample, read_label, input_subdir):
    input_dir = os.path.join(experiment_dir, input_subdir)
    matches = []
    for candidate in fastq_candidate_names(sample, read_label):
        candidate_path = os.path.join(input_dir, candidate)
        if os.path.exists(candidate_path):
            matches.append(candidate_path)
    matches = sorted(set(matches))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise FileNotFoundError(
            f"No FASTQ found for sample '{sample}' read '{read_label}' in {input_dir}. "
            f"Tried common suffixes for Illumina-style FASTQ names."
        )
    raise ValueError(
        f"Multiple FASTQs found for sample '{sample}' read '{read_label}' in {input_dir}: "
        + ", ".join(matches)
    )


# Obtain all the sample names
samples = [normalize_sample_name(f, input_file_type) for f in input_files]
samples2 = [re.sub(r'_L00.', '', string) for string in samples] #From step 4 on, the lane number is not in the sample name anymore

samples = list(set(samples))
samples2 = list(set(samples2))

if config['PAIRED'] == 1 and THEMODERANGEMIN < 4: 
    num_samples = len(samples) / 2
else: 
    num_samples = len(samples)
print(f"NUMBER OF SAMPLES = {num_samples}")

max_nodes = master_config.get('max_nodes', f"{master_config['nodes_in_partition']}") if master_config.get('max_nodes', f"{master_config['nodes_in_partition']}") <= master_config['nodes_in_partition'] else master_config['nodes_in_partition'] 

##--------------------------------------------------------------------------------------------------------------
# Obtain Threads and Memory per rule
##--------------------------------------------------------------------------------------------------------------
def get_rule_core_limits(rule_num):
    if rule_num <= 9:
        min_cores = master_config['mincores_single_sample_step1_9'][rule_num - 1]
        max_cores = master_config['maxcores_single_sample_step1_9'][rule_num - 1]
    else:
        offset = rule_num - 10
        min_cores = master_config['mincores_per_rule_run_step10_13'][offset]
        max_cores = master_config['maxcores_per_rule_run_step10_13'][offset]

    min_cores = (
        min_cores
        if isinstance(min_cores, int) and min_cores > 0
        else master_config['min_slice_cores']
    )
    max_cores = (
        max_cores
        if isinstance(max_cores, int) and max_cores > 0
        else master_config['cores_per_node']
    )

    min_cores = max(min_cores, master_config['min_slice_cores'])
    max_cores = min(max_cores, master_config['cores_per_node'])
    return min_cores, max(max_cores, min_cores)


Threads_Per_Rule = {}
for rule_num in range(1, master_config['max_step'] + 1):
    min_cores, max_cores = get_rule_core_limits(rule_num)
    Threads_Per_Rule[f'{rule_num}'] = min_cores if rule_num <= 9 else max_cores

Memory_Per_Rule = {}
for rule_num in range(1, master_config['max_step'] + 1):
    min_mem_mb = None
    if (
        'min_mem_mb' in master_config
        and isinstance(master_config['min_mem_mb'], list)
        and len(master_config['min_mem_mb']) > (rule_num - 1)
        and master_config['min_mem_mb'][rule_num - 1] is not None
        and isinstance(master_config['min_mem_mb'][rule_num - 1], int)
    ):
        min_mem_mb = master_config['min_mem_mb'][rule_num - 1]

    default_mem_mb = max(
        master_config['min_slice_mem'],
        Threads_Per_Rule[f'{rule_num}'] * master_config['max_mem_per_core_mb'],
        min_mem_mb if min_mem_mb is not None else 0,
    )
    Memory_Per_Rule[f'{rule_num}'] = default_mem_mb

Runtime_Per_Rule = {}
for rule_num in range(1, master_config['max_step'] + 1):
    runtime = master_config.get('default_runtime', 120)
    if (
        'rule_runtime' in master_config
        and isinstance(master_config['rule_runtime'], list)
        and len(master_config['rule_runtime']) > (rule_num - 1)
        and isinstance(master_config['rule_runtime'][rule_num - 1], int)
        and master_config['rule_runtime'][rule_num - 1] > 0
    ):
        runtime = master_config['rule_runtime'][rule_num - 1]
    Runtime_Per_Rule[f'{rule_num}'] = runtime

##--------------------------------------------------------------------------------------------------------------
# Include Snakemake rules for your actual data processing pipeline
##--------------------------------------------------------------------------------------------------------------

def check_and_include_rules(logfile, omnom_home, experiment_dir):
    if not is_worker_job:
        log_it(logfile, "Checking snake rule functions...", "PREFLIGHT")

    rules_dir = os.path.join(omnom_home, "rules")

    # Check if the rules directory exists
    if not os.path.isdir(rules_dir):
        log_it(logfile, f"Snake rules dir {rules_dir} does not exist!", "WORKLFLOW SNAKE RULE CHECK")
        log_it(logfile, "Aborting...")
        print(f"Snake rules dir {rules_dir} does not exist!")
        print("Aborting...")
        sys.exit(1)

    # Check if the directory contains any .smk files
    smk_files = glob.glob(os.path.join(rules_dir, "*.smk"))
    if len(smk_files) == 0:
        log_it(logfile, f"Snake Rules dir {rules_dir} does not contain any .smk files!", "WORKLFLOW SNAKE RULE CHECK")
        log_it(logfile, "Aborting...")
        print(f"Snake Rules dir {rules_dir} does not contain any .smk files!")
        print("Aborting...")
        sys.exit(1)

    # Sanity check if they contain the right header
    valid_smk_files = []
    for smk_file in smk_files:
        with open(smk_file, 'r') as file:
            lines = file.readlines()
            if len(lines) < 3 or lines[2].strip() != "## Omnomnomics Snake Rule ##":
                log_it(logfile, f"Workflow step dir contains a malformed .smk file: {smk_file}", "WORKLFLOW SNAKE RULE CHECK")
                log_it(logfile, f"Aborting...")
                print(f"Workflow step dir contains a malformed .smk file: {smk_file}")
                print(f"Aborting...")
                sys.exit(1)
            valid_smk_files.append(smk_file)


    if not is_worker_job:
        log_it(logfile, "Workflow Snake rules Good!", "WORKLFLOW SNAKE RULES CHECK")
    return valid_smk_files

valid_smk_files = check_and_include_rules(logfile, workflow_root, experiment_dir)
# Include snake rules in the main Snakefile
for smk_file in valid_smk_files:
    include: smk_file

##--------------------------------------------------------------------------------------------------------------
# Execute the desired rules
##--------------------------------------------------------------------------------------------------------------
for rule_num in themode:

    output_folder = master_config['output_folders'][rule_num-1]
    output_file_type =  master_config['output_file_types'][rule_num-1]
    # For all of the specified rules, add its output to input of rule all so that the rule is run
    if rule_num == 1:
        if config['THETRIMTOOL'] == 'skewer':
            if config['PAIRED']: 
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}_R1{output_file_type}", sample = samples)
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}_R2{output_file_type}", sample = samples)
            else: 
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}{output_file_type}", sample = samples)
        else:
            if config['PAIRED']: 
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}_R1{output_file_type}", sample = samples)
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}_R2{output_file_type}", sample = samples)
            else:
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}{output_file_type}", sample = samples)
    if rule_num == 2:
        if config['PAIRED']: 
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}_R1{output_file_type[0]}", sample = samples)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}_R2{output_file_type[0]}", sample = samples)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}_R1{output_file_type[1]}", sample = samples)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}_R2{output_file_type[1]}", sample = samples)
        else:
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}{output_file_type[0]}", sample = samples)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}{output_file_type[1]}", sample = samples)
    if rule_num == 3:
        if config['THEMAPTOOL'] == 'star':
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.bam", sample = samples)
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_3.tmp", sample = samples)
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.STAR_stats.txt", sample = samples)
        if config['THEMAPTOOL'] == 'hisat2':
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.bam", sample = samples)
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_3.tmp", sample = samples)
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.HISAT2_stats.txt", sample = samples)
        if config['THEMAPTOOL'] == 'star_te':
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.bam", sample = samples)
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_3.tmp", sample = samples)
                all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.STAR_TE_stats.txt", sample = samples)
    if rule_num == 4:
        all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.bam", sample = samples2)
        all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_4.tmp",  sample = samples2)
    if rule_num == 5:
        if config['THETYPE'] != "CHIP":
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.sorted.dups_marked.filtered.bam", sample = samples2)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_5.tmp",  sample = samples2)
        else: 
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.filtered.bam", sample = samples2)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_5.tmp",  sample = samples2)
    if rule_num == 6:
        if config['THETYPE'] != "CHIP":
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.sorted.dups_marked.filtered.bam.bai", sample = samples2)
        else: 
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.filtered.bam.bai", sample = samples2)
    if rule_num == 7:
        if config['THETYPE'] != "CHIP":
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.sorted.dups_marked.filtered.bam.stats.txt", sample = samples2)
        else: 
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.filtered.bam.stats.txt", sample = samples2)
    if rule_num == 8:
        if config['THETYPE'] != "CHIP":
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz", sample = samples2)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_8.tmp",  sample = samples2)
        else:
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.filtered.HOMER_tagDir.tar.gz", sample = samples2)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_8.tmp",  sample = samples2)

    if rule_num == 9:
        all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_9.tmp",  sample = samples2)
    if rule_num == 10:
        all_outputs.append( f"{experiment_dir}/{output_folder}/extra_10.tmp")
    if rule_num == 11:
        all_outputs.append( f"{experiment_dir}/{output_folder}/extra_11.tmp")
    if rule_num == 12:
        all_outputs.append( f"{experiment_dir}/{output_folder}/extra_12.tmp")
    if rule_num == 13:
        all_outputs.append( f"{experiment_dir}/{output_folder}/{os.path.basename(config['EXPERIMENT_DIR'])}.results.zip")
    

if not is_worker_job:
    log_it(logfile, '\n'.join(all_outputs), 'ALL OUTPUTS')

# Determine rule priority to resolve any rule ambuigity
if config['THETRIMTOOL'] == "skewer" and config['THEMAPTOOL'] == "star":
    ruleorder: run_skewer > run_trimmomatic > run_star > run_star_te > run_hisat2 > merge_bam 
elif config['THETRIMTOOL'] == "skewer" and config['THEMAPTOOL'] == "star_te":
    ruleorder: run_skewer > run_trimmomatic > run_star_te > run_star > run_hisat2 > merge_bam 
elif config['THETRIMTOOL'] == "skewer" and config['THEMAPTOOL'] == "hisat2":
    ruleorder: run_skewer > run_trimmomatic > run_hisat2 > run_star_te > run_star > merge_bam 
elif config['THETRIMTOOL'] == "trimmomatic" and config['THEMAPTOOL'] == "star":
    ruleorder: run_trimmomatic > run_skewer > run_star > run_star_te > run_hisat2 > merge_bam 
elif config['THETRIMTOOL'] == "trimmomatic" and config['THEMAPTOOL'] == "star_te":
    ruleorder: run_trimmomatic > run_skewer > run_star_te > run_star > run_hisat2 > merge_bam 
elif config['THETRIMTOOL'] == "trimmomatic" and config['THEMAPTOOL'] == "hisat2":
    ruleorder: run_trimmomatic > run_skewer > run_hisat2 > run_star_te > run_star > merge_bam 

#---------------------------------------------------------------------------------------------------------------
# Three line heart of the pipeline to set up the workflow
#---------------------------------------------------------------------------------------------------------------
rule all:
    input:
        all_outputs


onsuccess:
    if not is_worker_job:
        ##Final things:
        #---------------------------------------------------------------------------------------------------------------
        # Run multiqc to gather all stats
        #---------------------------------------------------------------------------------------------------------------
        if not config['NO_MULTIQC']:
            log_it(logfile, "Running multiQC...", "STATS")
            shell(f"""
                eval "$(micromamba shell hook --shell=bash)" && micromamba activate multiqc && \
                multiqc --filename MultiQC/omnomnomics.run.{run_date}.multiqc_report.html --dirs --export .
            """)
        #---------------------------------------------------------------------------------------------------------------
        # Clean up tmp files for workflow
        #---------------------------------------------------------------------------------------------------------------
        if 3 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['merge_rule_num']-1]}/*.extra_3.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 4 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['merge_rule_num']-1]}/*.extra_4.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 5 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['touchup_rule_num']-1]}/*.extra_5.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 8 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['tagdir_rule_num']-1]}/*.extra_8.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 9 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['wig_rule_num']-1]}/*.extra_9.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 10 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['mergewig_rule_num']-1]}/extra_10.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 11 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/extra_11.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 12 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/extra_11.tmp")
            for file in list_of_extra_files:
                os.remove(file)

        #---------------------------------------------------------------------------------------------------------------
        # Remove old BAM files with lane info
        #---------------------------------------------------------------------------------------------------------------
        ## Note that I had to do that here since if I do it before completion rule_all would error that not all of it input files are present
        if 4 in themode:
            for bam_file in glob.glob(f"{experiment_dir}/{master_config['output_folders'][master_config['merge_rule_num']-1]}/*_L00*.bam"):
                os.remove(bam_file)

        #---------------------------------------------------------------------------------------------------------------
        # Log elapsed time and completion
        #---------------------------------------------------------------------------------------------------------------
        end_time = time.time()
        elapsed_time = end_time - start_time
        log_it(logfile, "Pipeline finished at {}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")), "PIPELINE RUN TIME")
        log_it(logfile, "Total elapsed time: {:.2f} seconds\n".format(elapsed_time))

        log_it(logfile, "All done!" ,"FINAL REMARKS")
        log_it(logfile, "Good luck with your downstream analyses!")
