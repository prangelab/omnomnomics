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
import csv
import fcntl
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
log_marker_dir = os.path.join(experiment_dir, "run_logs", f"omnomnomics.run.{run_date}.markers")
step_log_dir = os.path.join(experiment_dir, "run_logs", "steps")
is_worker_job = "--target-jobs" in sys.argv
os.makedirs(log_marker_dir, exist_ok=True)
os.makedirs(step_log_dir, exist_ok=True)

# Function to log messages
def log_it(logfile, message, heading=None):
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(logfile, 'a') as log:
        if heading: #Check if heading
            log.write("\n{}\n\n".format(heading))
        log.write(f"{timestamp}: {message}\n") #write to log
        log.flush()
        os.fsync(log.fileno())
    print(f"{timestamp}: {message}")

def log_once(logfile, marker_name, message, heading=None):
    marker_path = os.path.join(log_marker_dir, marker_name)
    try:
        fd = os.open(marker_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
    except FileExistsError:
        return
    log_it(logfile, message, heading)


def expected_sample_count_for_step(step_num):
    if step_num in {1, 2, 3}:
        return len(lane_samples)
    if step_num in {4, 5, 6, 7, 8, 13}:
        return len(samples2)
    return 1


def step_tracking_paths(step_num):
    step_prefix = f"step{step_num:02d}"
    step_state_dir = os.path.join(step_log_dir, f".{step_prefix}.state")
    os.makedirs(step_state_dir, exist_ok=True)
    return {
        "summary_tsv": os.path.join(step_log_dir, f"{step_prefix}.summary.tsv"),
        "commands_txt": os.path.join(step_log_dir, f"{step_prefix}.commands.txt"),
        "notes_txt": os.path.join(step_log_dir, f"{step_prefix}.notes.txt"),
        "started_dir": os.path.join(step_state_dir, "started"),
        "completed_dir": os.path.join(step_state_dir, "completed"),
        "failed_dir": os.path.join(step_state_dir, "failed"),
        "entered_marker": os.path.join(step_state_dir, "entered.marker"),
        "finished_marker": os.path.join(step_state_dir, "finished.marker"),
    }


def ensure_step_tracking_dirs(step_num):
    paths = step_tracking_paths(step_num)
    for key in ("started_dir", "completed_dir", "failed_dir"):
        os.makedirs(paths[key], exist_ok=True)
    return paths


def append_locked_text(path, text):
    with open(path, "a") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def update_step_summary_row(summary_tsv, row):
    fieldnames = [
        "sample",
        "job_id",
        "status",
        "start_time",
        "end_time",
        "elapsed_seconds",
        "worker_log",
    ]
    rows = []
    if os.path.exists(summary_tsv):
        with open(summary_tsv, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            rows = list(reader)
    updated = False
    for existing in rows:
        if existing["sample"] == row["sample"]:
            existing.update(row)
            updated = True
            break
    if not updated:
        rows.append(row)
    with open(summary_tsv, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def worker_log_path(rule_log_subdir):
    job_id = os.environ.get("SLURM_JOB_ID", "NA")
    if not rule_log_subdir or job_id == "NA":
        return "NA"
    return os.path.join(experiment_dir, "slurm_logs", rule_log_subdir, f"{rule_log_subdir}.{job_id}.out")


def begin_step_sample(step_num, sample, rule_log_subdir):
    paths = ensure_step_tracking_dirs(step_num)
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    started_marker = os.path.join(paths["started_dir"], sample)
    entered_first = False
    try:
        fd = os.open(paths["entered_marker"], os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
        entered_first = True
    except FileExistsError:
        pass
    try:
        fd = os.open(started_marker, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
    except FileExistsError:
        pass
    update_step_summary_row(
        paths["summary_tsv"],
        {
            "sample": sample,
            "job_id": os.environ.get("SLURM_JOB_ID", "NA"),
            "status": "RUNNING",
            "start_time": start_time,
            "end_time": "",
            "elapsed_seconds": "",
            "worker_log": worker_log_path(rule_log_subdir),
        },
    )
    if entered_first:
        log_it(logfile, f"Step {step_num}: first sample entered ({sample})")
    return {"start_time": start_time, "paths": paths}


def record_step_command(step_num, sample, command):
    paths = ensure_step_tracking_dirs(step_num)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    append_locked_text(
        paths["commands_txt"],
        f"[{timestamp}] sample={sample}\n{command}\n\n",
    )


def record_step_note(step_num, sample, message):
    paths = ensure_step_tracking_dirs(step_num)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    append_locked_text(
        paths["notes_txt"],
        f"{timestamp}\t{sample}\t{message}\n",
    )


def finish_step_sample(step_num, sample, rule_log_subdir, start_time, status):
    paths = ensure_step_tracking_dirs(step_num)
    end_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    started_dt = datetime.datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
    ended_dt = datetime.datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")
    elapsed_seconds = f"{(ended_dt - started_dt).total_seconds():.2f}"
    marker_dir = paths["completed_dir"] if status == "OK" else paths["failed_dir"]
    marker_path = os.path.join(marker_dir, sample)
    try:
        fd = os.open(marker_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.close(fd)
    except FileExistsError:
        pass
    update_step_summary_row(
        paths["summary_tsv"],
        {
            "sample": sample,
            "job_id": os.environ.get("SLURM_JOB_ID", "NA"),
            "status": status,
            "start_time": start_time,
            "end_time": end_time,
            "elapsed_seconds": elapsed_seconds,
            "worker_log": worker_log_path(rule_log_subdir),
        },
    )
    expected = expected_sample_count_for_step(step_num)
    completed = len(os.listdir(paths["completed_dir"]))
    failed = len(os.listdir(paths["failed_dir"]))
    if completed + failed >= expected:
        try:
            fd = os.open(paths["finished_marker"], os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.close(fd)
            log_it(logfile, f"Step {step_num}: last sample finished ({sample}); {completed}/{expected} OK, {failed}/{expected} failed")
        except FileExistsError:
            pass

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
def sanity_check_dir(logfile, input_directory, file_ext, marker_name=None):
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
    return

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
        header_lines = [line.strip() for line in config_file if line.strip()]
        if len(header_lines) < 1 or expected_header not in header_lines[0]:
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
    log_it(logfile, f"Duplicate handling: {config['DUPLICATE_HANDLING']}")

    log_it(logfile, f"Design formula: {config['MYFORMULA']}", "READ COUNTING SETTINGS")
    log_it(logfile, f"Metadata file: {config['MYMETADATA']}")

    log_it(logfile, f"Input file (MACS3): {config['INPUT']}", "PEAK CALLING SETTINGS")
    log_it(logfile, f"Broad peaks: {config['BROAD']}")

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
create_homer_tagdirs = config.get('CREATE_HOMER_TAGDIRS', False)

THEMODERANGEMIN = config['THEMODERANGEMIN']

input_folder = master_config['input_folders'][THEMODERANGEMIN-1]
input_file_type =  master_config['input_file_types'][THEMODERANGEMIN-1]
if THEMODERANGEMIN == 11:
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

lane_sample_files = []
lane_sample_sources = [
    (master_config['input_folders'][master_config['trim_rule_num'] - 1], FASTQ_EXTENSIONS),
    (master_config['output_folders'][master_config['trim_rule_num'] - 1], (".trimmed.fastq.gz", ".trimmed.fastq", ".trimmed.fq.gz", ".trimmed.fq")),
    (master_config['output_folders'][master_config['map_rule_num'] - 1], (".bam",)),
]
for folder_name, extensions in lane_sample_sources:
    candidate_files = []
    for extension in extensions:
        candidate_files.extend(glob.glob(os.path.join(experiment_dir, folder_name, f"*{extension}")))
    if folder_name == master_config['output_folders'][master_config['map_rule_num'] - 1]:
        candidate_files = [path for path in candidate_files if re.search(r'_L0\d+\.bam$', os.path.basename(path))]
    if candidate_files:
        lane_sample_files = candidate_files
        break

lane_samples = []
for path in lane_sample_files:
    extension = next((ext for ext in FASTQ_EXTENSIONS if path.endswith(ext)), None)
    if extension is None:
        for ext in (".trimmed.fastq.gz", ".trimmed.fastq", ".trimmed.fq.gz", ".trimmed.fq", ".bam"):
            if path.endswith(ext):
                extension = ext
                break
    if extension is None:
        continue
    lane_samples.append(normalize_sample_name(path, extension))

lane_samples = sorted(set(lane_samples))
if lane_samples:
    lane_sample_wildcard_pattern = "|".join(re.escape(sample_name) for sample_name in lane_samples)
else:
    lane_sample_wildcard_pattern = r"$.^"

merged_sample_wildcard_pattern = "|".join(re.escape(sample_name) for sample_name in sorted(samples2))

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
internal_max_rule_num = max(master_config['max_step'], master_config.get('homer_tagdir_rule_num', master_config['max_step']))
for rule_num in range(1, internal_max_rule_num + 1):
    min_cores, max_cores = get_rule_core_limits(rule_num)
    Threads_Per_Rule[f'{rule_num}'] = min_cores if rule_num <= 9 else max_cores

Memory_Per_Rule = {}
for rule_num in range(1, internal_max_rule_num + 1):
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
for rule_num in range(1, internal_max_rule_num + 1):
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
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.sorted.dups_marked.filtered.bam.qc_summary.pdf", sample = samples2)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.sorted.dups_marked.filtered.bam.qc_summary.svg", sample = samples2)
        else: 
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.filtered.bam.stats.txt", sample = samples2)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.filtered.bam.qc_summary.pdf", sample = samples2)
            all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.filtered.bam.qc_summary.svg", sample = samples2)
    if rule_num == 8:
        all_outputs += expand(f"{experiment_dir}/{output_folder}/{{sample}}.extra_8.tmp",  sample = samples2)
    if rule_num == 9:
        all_outputs.append( f"{experiment_dir}/{output_folder}/extra_9.tmp")
    if rule_num == 10:
        all_outputs.append( f"{experiment_dir}/{output_folder}/extra_10.tmp")
    if rule_num == 11:
        all_outputs.append(f"{experiment_dir}/{output_folder}/{os.path.basename(config['EXPERIMENT_DIR'])}.raw_read_quant.table.txt")
        if config['THETYPE'] == "RNA":
            all_outputs.append(f"{experiment_dir}/{output_folder}/{os.path.basename(config['EXPERIMENT_DIR'])}.featureCounts.summary.txt")
        all_outputs.append(f"{experiment_dir}/{output_folder}/extra_11.tmp")
    if rule_num == 12:
        all_outputs.append( f"{experiment_dir}/{output_folder}/{os.path.basename(config['EXPERIMENT_DIR'])}.results.zip")

if create_homer_tagdirs:
    homer_output_folder = master_config['output_folders'][master_config['homer_tagdir_rule_num'] - 1]
    if config['THETYPE'] != "CHIP":
        all_outputs += expand(f"{experiment_dir}/{homer_output_folder}/{{sample}}.sorted.dups_marked.filtered.HOMER_tagDir.tar.gz", sample=samples2)
        all_outputs += expand(f"{experiment_dir}/{homer_output_folder}/{{sample}}.extra_13.tmp", sample=samples2)
    else:
        all_outputs += expand(f"{experiment_dir}/{homer_output_folder}/{{sample}}.filtered.HOMER_tagDir.tar.gz", sample=samples2)
        all_outputs += expand(f"{experiment_dir}/{homer_output_folder}/{{sample}}.extra_13.tmp", sample=samples2)

if not is_worker_job:
    target_manifest = os.path.join(experiment_dir, "run_logs", f"omnomnomics.run.{run_date}.targets.txt")
    with open(target_manifest, "w") as handle:
        handle.write("\n".join(all_outputs) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    log_it(logfile, f"Planned targets: {len(all_outputs)}. Full target manifest: {target_manifest}", "ALL OUTPUTS")

# Determine rule priority to resolve any rule ambuigity
if config['THETRIMTOOL'] == "fastp" and config['THEMAPTOOL'] == "star":
    ruleorder: run_fastp > run_skewer > run_star > run_star_te > run_hisat2 > merge_bam
elif config['THETRIMTOOL'] == "fastp" and config['THEMAPTOOL'] == "star_te":
    ruleorder: run_fastp > run_skewer > run_star_te > run_star > run_hisat2 > merge_bam
elif config['THETRIMTOOL'] == "fastp" and config['THEMAPTOOL'] == "hisat2":
    ruleorder: run_fastp > run_skewer > run_hisat2 > run_star_te > run_star > merge_bam
elif config['THETRIMTOOL'] == "skewer" and config['THEMAPTOOL'] == "star":
    ruleorder: run_skewer > run_fastp > run_star > run_star_te > run_hisat2 > merge_bam 
elif config['THETRIMTOOL'] == "skewer" and config['THEMAPTOOL'] == "star_te":
    ruleorder: run_skewer > run_fastp > run_star_te > run_star > run_hisat2 > merge_bam 
elif config['THETRIMTOOL'] == "skewer" and config['THEMAPTOOL'] == "hisat2":
    ruleorder: run_skewer > run_fastp > run_hisat2 > run_star_te > run_star > merge_bam 

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
            import csv
            import gzip
            import statistics
            import subprocess
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt

            def normalize_numeric_string(value):
                return value.replace(",", "").replace("%", "").strip()

            def parse_star_mapper_stats(stats_path):
                metrics = {}
                with open(stats_path) as handle:
                    for line in handle:
                        if "|" not in line:
                            continue
                        key, value = [part.strip() for part in line.split("|", 1)]
                        metrics[key] = value
                total_reads = int(normalize_numeric_string(metrics.get("Number of input reads", "0")))
                unique_reads = int(normalize_numeric_string(metrics.get("Uniquely mapped reads number", "0")))
                multi_reads = (
                    int(normalize_numeric_string(metrics.get("Number of reads mapped to multiple loci", "0")))
                    + int(normalize_numeric_string(metrics.get("Number of reads mapped to too many loci", "0")))
                )
                aligned_reads = unique_reads + multi_reads
                aligned_pct = round((aligned_reads / total_reads) * 100, 4) if total_reads else None
                return {
                    "mapper_input_reads": total_reads,
                    "mapper_reported_aligned_reads": aligned_reads,
                    "mapper_uniquely_mapped_reads": unique_reads,
                    "mapper_multimapped_reads": multi_reads,
                    "mapper_reported_alignment_pct": aligned_pct,
                }

            def parse_hisat2_mapper_stats(stats_path):
                total_reads = None
                alignment_pct = None
                with open(stats_path) as handle:
                    for line in handle:
                        total_match = re.match(r"^\s*([\d,]+)\s+reads; of these:", line)
                        if total_match:
                            total_reads = int(normalize_numeric_string(total_match.group(1)))
                        pct_match = re.match(r"^\s*([\d.]+)% overall alignment rate", line)
                        if pct_match:
                            alignment_pct = float(pct_match.group(1))
                aligned_reads = round((total_reads * alignment_pct) / 100) if total_reads is not None and alignment_pct is not None else None
                return {
                    "mapper_input_reads": total_reads,
                    "mapper_reported_aligned_reads": aligned_reads,
                    "mapper_uniquely_mapped_reads": None,
                    "mapper_multimapped_reads": None,
                    "mapper_reported_alignment_pct": alignment_pct,
                }

            def count_fastq_reads(path):
                open_func = gzip.open if path.endswith(".gz") else open
                line_count = 0
                with open_func(path, "rt") as handle:
                    for _ in handle:
                        line_count += 1
                return line_count // 4

            def sample_qc_stats_path(sample_name):
                stats_output_folder = master_config['output_folders'][master_config['stats_rule_num'] - 1]
                if config['THETYPE'] != "CHIP":
                    return f"{experiment_dir}/{stats_output_folder}/{sample_name}.sorted.dups_marked.filtered.bam.stats.txt"
                return f"{experiment_dir}/{stats_output_folder}/{sample_name}.filtered.bam.stats.txt"

            def mapper_stats_path(sample_name):
                map_output_folder = master_config['output_folders'][master_config['map_rule_num'] - 1]
                if config['THEMAPTOOL'] == "hisat2":
                    return f"{experiment_dir}/{map_output_folder}/{sample_name}.HISAT2_stats.txt"
                if config['THEMAPTOOL'] == "star_te":
                    return f"{experiment_dir}/{map_output_folder}/{sample_name}.STAR_TE_stats.txt"
                return f"{experiment_dir}/{map_output_folder}/{sample_name}.STAR_stats.txt"

            sample_to_raw_inputs = {}
            for path in input_files:
                normalized_sample = normalize_sample_name(path, input_file_type)
                sample_to_raw_inputs.setdefault(normalized_sample, []).append(path)

            sample2_to_lane_samples = {}
            for sample_name in samples:
                aggregated_sample = re.sub(r'_L00.', '', sample_name)
                sample2_to_lane_samples.setdefault(aggregated_sample, []).append(sample_name)

            def load_alignment_qc_rows(tsv_paths):
                rows = []
                for tsv_path in tsv_paths:
                    sample_name = None
                    with open(tsv_path, newline="") as handle:
                        for line in handle:
                            if line.startswith("# sample\t"):
                                sample_name = line.strip().split("\t", 1)[1]
                                break
                    with open(tsv_path, newline="") as handle:
                        reader = csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t")
                        for row in reader:
                            value = row["value"]
                            if value == "NA":
                                parsed_value = None
                            else:
                                try:
                                    parsed_value = float(value)
                                except ValueError:
                                    parsed_value = None
                            rows.append({
                                "sample": sample_name or row["sample"],
                                "stage": row["stage"],
                                "metric": row["metric"],
                                "unit": row["unit"],
                                "value": parsed_value,
                            })
                return rows

            def build_flow_qc_rows():
                rows = []
                trim_output_folder = master_config['output_folders'][master_config['trim_rule_num'] - 1]
                mapper_parser = parse_hisat2_mapper_stats if config['THEMAPTOOL'] == "hisat2" else parse_star_mapper_stats

                for sample_name in samples2:
                    lane_samples = sample2_to_lane_samples.get(sample_name, [sample_name])
                    raw_read_total = 0
                    trimmed_read_total = 0

                    for lane_sample in lane_samples:
                        raw_files = sample_to_raw_inputs.get(lane_sample, [])
                        raw_read_total += sum(count_fastq_reads(path) for path in raw_files if os.path.exists(path))

                        if config['PAIRED']:
                            trimmed_paths = [
                                f"{experiment_dir}/{trim_output_folder}/{lane_sample}_R1.trimmed.fastq.gz",
                                f"{experiment_dir}/{trim_output_folder}/{lane_sample}_R2.trimmed.fastq.gz",
                            ]
                        else:
                            trimmed_paths = [f"{experiment_dir}/{trim_output_folder}/{lane_sample}.trimmed.fastq.gz"]
                        trimmed_read_total += sum(count_fastq_reads(path) for path in trimmed_paths if os.path.exists(path))

                    if raw_read_total:
                        rows.append({"sample": sample_name, "stage": "raw_fastq", "metric": "raw_reads", "unit": "reads", "value": raw_read_total})
                    if trimmed_read_total:
                        rows.append({"sample": sample_name, "stage": "trimmed_fastq", "metric": "trimmed_reads", "unit": "reads", "value": trimmed_read_total})
                    if raw_read_total and trimmed_read_total:
                        rows.append({
                            "sample": sample_name,
                            "stage": "trimmed_fastq",
                            "metric": "trimmed_reads_retained_pct",
                            "unit": "percent",
                            "value": round((trimmed_read_total / raw_read_total) * 100, 4),
                        })

                    mapper_totals = {
                        "mapper_input_reads": 0,
                        "mapper_reported_aligned_reads": 0,
                        "mapper_uniquely_mapped_reads": 0,
                        "mapper_multimapped_reads": 0,
                    }
                    mapper_fields_present = {key: False for key in mapper_totals}

                    for lane_sample in lane_samples:
                        stats_path = mapper_stats_path(lane_sample)
                        if not os.path.exists(stats_path):
                            continue
                        lane_metrics = mapper_parser(stats_path)
                        for key in mapper_totals:
                            value = lane_metrics.get(key)
                            if value is not None:
                                mapper_totals[key] += value
                                mapper_fields_present[key] = True

                    if mapper_fields_present["mapper_input_reads"]:
                        rows.append({
                            "sample": sample_name,
                            "stage": "mapper_report",
                            "metric": "mapper_input_reads",
                            "unit": "reads",
                            "value": mapper_totals["mapper_input_reads"],
                        })
                    if mapper_fields_present["mapper_reported_aligned_reads"]:
                        rows.append({
                            "sample": sample_name,
                            "stage": "mapper_report",
                            "metric": "mapper_reported_aligned_reads",
                            "unit": "reads",
                            "value": mapper_totals["mapper_reported_aligned_reads"],
                        })
                    if mapper_fields_present["mapper_uniquely_mapped_reads"]:
                        rows.append({
                            "sample": sample_name,
                            "stage": "mapper_report",
                            "metric": "mapper_uniquely_mapped_reads",
                            "unit": "reads",
                            "value": mapper_totals["mapper_uniquely_mapped_reads"],
                        })
                    if mapper_fields_present["mapper_multimapped_reads"]:
                        rows.append({
                            "sample": sample_name,
                            "stage": "mapper_report",
                            "metric": "mapper_multimapped_reads",
                            "unit": "reads",
                            "value": mapper_totals["mapper_multimapped_reads"],
                        })
                    if mapper_fields_present["mapper_input_reads"] and mapper_fields_present["mapper_reported_aligned_reads"] and mapper_totals["mapper_input_reads"]:
                        rows.append({
                            "sample": sample_name,
                            "stage": "mapper_report",
                            "metric": "mapper_reported_alignment_pct",
                            "unit": "percent",
                            "value": round((mapper_totals["mapper_reported_aligned_reads"] / mapper_totals["mapper_input_reads"]) * 100, 4),
                        })

                return rows

            def metric_values(rows, stage, metric):
                values = []
                for row in rows:
                    if row["stage"] == stage and row["metric"] == metric and row["value"] is not None:
                        values.append(row["value"])
                return values

            def draw_alignment_qc_panel(ax, rows, metric_specs, title, ylabel, log_scale=False):
                point_color = "#355070"
                summary_color = "#e56b6f"
                box_face_color = "#dbe7ee"
                box_edge_color = "#8aa1b1"

                ax.set_title(title)
                if log_scale:
                    ax.set_yscale("symlog", linthresh=1)
                ax.set_ylabel(ylabel)
                ax.grid(axis="y", color="#e6e6e6", linewidth=0.8)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)

                for idx, (label, stage, metric) in enumerate(metric_specs):
                    values = metric_values(rows, stage, metric)
                    if not values:
                        ax.text(idx, 0.5, "NA", ha="center", va="bottom", fontsize=10, color="#777777")
                        continue
                    if len(values) > 1:
                        ax.boxplot(
                            values,
                            positions=[idx],
                            widths=0.5,
                            patch_artist=True,
                            showfliers=False,
                            boxprops={"facecolor": box_face_color, "edgecolor": box_edge_color, "linewidth": 1.2},
                            whiskerprops={"color": box_edge_color, "linewidth": 1.2},
                            capprops={"color": box_edge_color, "linewidth": 1.2},
                            medianprops={"color": box_edge_color, "linewidth": 1.4},
                        )
                    offsets = [idx + ((pos - (len(values) - 1) / 2) * 0.04) for pos in range(len(values))]
                    ax.scatter(offsets, values, color=point_color, alpha=0.85, s=30, zorder=3)
                    mean_value = statistics.mean(values)
                    sd_value = statistics.stdev(values) if len(values) > 1 else 0.0
                    ax.errorbar(
                        idx,
                        mean_value,
                        yerr=sd_value,
                        fmt="o",
                        color=summary_color,
                        ecolor=summary_color,
                        elinewidth=2,
                        capsize=4,
                        markersize=7,
                        zorder=4,
                    )
                    ax.text(idx, mean_value, f"{mean_value:.1f}", ha="center", va="bottom", fontsize=9, color=summary_color)

                ax.set_xticks(range(len(metric_specs)))
                ax.set_xticklabels([label for label, _, _ in metric_specs], rotation=15, ha="right")

            def write_alignment_qc_experiment_summary(tsv_paths):
                os.makedirs(f"{experiment_dir}/MultiQC", exist_ok=True)
                experiment_name = os.path.basename(config["EXPERIMENT_DIR"])
                aggregate_tsv = f"{experiment_dir}/MultiQC/{experiment_name}.alignment_qc_experiment_summary.tsv"
                aggregate_pdf = f"{experiment_dir}/MultiQC/{experiment_name}.alignment_qc_experiment_summary.pdf"
                aggregate_svg = f"{experiment_dir}/MultiQC/{experiment_name}.alignment_qc_experiment_summary.svg"
                rows = load_alignment_qc_rows(tsv_paths) + build_flow_qc_rows()

                with open(aggregate_tsv, "w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["sample", "stage", "metric", "unit", "value"])
                    for row in rows:
                        writer.writerow([
                            row["sample"],
                            row["stage"],
                            row["metric"],
                            row["unit"],
                            "NA" if row["value"] is None else row["value"],
                        ])

                fig, axes = plt.subplots(2, 2, figsize=(14, 10))
                fig.patch.set_facecolor("white")
                fig.suptitle(f"{experiment_name} alignment QC experiment summary", fontsize=16, fontweight="bold")

                draw_alignment_qc_panel(
                    axes[0, 0],
                    rows,
                    [
                        ("Raw reads", "raw_fastq", "raw_reads"),
                        ("Trimmed reads", "trimmed_fastq", "trimmed_reads"),
                        ("Aligned reads", "mapper_report", "mapper_reported_aligned_reads"),
                        ("Filtered mapped reads", "post_filter", "mapped_primary_reads"),
                    ],
                    "Read counts across the pipeline",
                    "Reads",
                    log_scale=True,
                )

                draw_alignment_qc_panel(
                    axes[0, 1],
                    rows,
                    [
                        ("Trim retained %", "trimmed_fastq", "trimmed_reads_retained_pct"),
                        ("Mapper aligned %", "mapper_report", "mapper_reported_alignment_pct"),
                        ("Mapped retained %", "derived", "mapped_primary_reads_retained_pct"),
                        ("Proper pairs retained %", "derived", "properly_paired_templates_retained_pct"),
                    ] if config["PAIRED"] else [
                        ("Trim retained %", "trimmed_fastq", "trimmed_reads_retained_pct"),
                        ("Mapper aligned %", "mapper_report", "mapper_reported_alignment_pct"),
                        ("Mapped retained %", "derived", "mapped_primary_reads_retained_pct"),
                    ],
                    "Pipeline retention metrics",
                    "Percent",
                )

                draw_alignment_qc_panel(
                    axes[1, 0],
                    rows,
                    [
                        ("Mapped reads removed", "derived", "mapped_primary_reads_removed"),
                        ("Duplicate reads removed", "derived", "duplicate_primary_reads_removed"),
                        ("Discordant pairs removed", "derived", "discordant_templates_removed"),
                    ] if config["PAIRED"] else [
                        ("Mapped reads removed", "derived", "mapped_primary_reads_removed"),
                        ("Duplicate reads removed", "derived", "duplicate_primary_reads_removed"),
                    ],
                    "Reads or pairs removed",
                    "Count",
                    log_scale=True,
                )

                draw_alignment_qc_panel(
                    axes[1, 1],
                    rows,
                    [
                        ("Unique mapped reads", "mapper_report", "mapper_uniquely_mapped_reads"),
                        ("Multi mapped reads", "mapper_report", "mapper_multimapped_reads"),
                        ("Pre duplicate reads", "pre_filter", "duplicate_primary_reads"),
                    ] if config["PAIRED"] else [
                        ("Unique mapped reads", "mapper_report", "mapper_uniquely_mapped_reads"),
                        ("Multi mapped reads", "mapper_report", "mapper_multimapped_reads"),
                        ("Pre duplicate reads", "pre_filter", "duplicate_primary_reads"),
                    ],
                    "Mapper and duplication profile",
                    "Reads",
                    log_scale=True,
                )

                fig.tight_layout(rect=[0, 0.02, 1, 0.95])
                fig.savefig(aggregate_pdf)
                fig.savefig(aggregate_svg)
                plt.close(fig)
                return aggregate_tsv, aggregate_pdf, aggregate_svg, rows

            def rows_to_sample_metric_map(rows):
                sample_map = {}
                for row in rows:
                    if row["value"] is None:
                        continue
                    sample_map.setdefault(row["sample"], {})
                    sample_map[row["sample"]][f"{row['stage']}__{row['metric']}"] = row["value"]
                return sample_map

            def write_multiqc_custom_content(rows):
                multiqc_dir = f"{experiment_dir}/MultiQC"
                os.makedirs(multiqc_dir, exist_ok=True)

                generalstats_path = os.path.join(multiqc_dir, "omnomnomics_alignment_flow_generalstats_mqc.yaml")
                table_path = os.path.join(multiqc_dir, "omnomnomics_alignment_flow_table_mqc.yaml")
                sample_metric_map = rows_to_sample_metric_map(rows)

                generalstats_yaml = {
                    "id": "omnomnomics_alignment_flow_generalstats",
                    "plot_type": "generalstats",
                    "headers": {
                        "raw_fastq__raw_reads": {
                            "title": "Raw Reads",
                            "description": "Raw FASTQ reads across all files for this sample.",
                            "scale": "Blues",
                            "format": "{:,.0f}",
                        },
                        "trimmed_fastq__trimmed_reads_retained_pct": {
                            "title": "Trim %",
                            "description": "Trimmed reads retained relative to raw FASTQ reads.",
                            "min": 0,
                            "max": 100,
                            "suffix": "%",
                            "scale": "RdYlGn",
                            "format": "{:,.1f}",
                        },
                        "mapper_report__mapper_reported_alignment_pct": {
                            "title": "Align %",
                            "description": "Mapper-reported aligned reads as a percentage of mapper input reads.",
                            "min": 0,
                            "max": 100,
                            "suffix": "%",
                            "scale": "RdYlGn",
                            "format": "{:,.1f}",
                        },
                        "derived__mapped_primary_reads_retained_pct": {
                            "title": "Filt Retained %",
                            "description": "Post-filter mapped primary reads retained relative to the pre-filter BAM.",
                            "min": 0,
                            "max": 100,
                            "suffix": "%",
                            "scale": "RdYlGn",
                            "format": "{:,.1f}",
                        },
                        "derived__duplicate_primary_reads_removed": {
                            "title": "Dup Removed",
                            "description": "Duplicate primary reads removed between the pre-filter and post-filter BAMs.",
                            "scale": "OrRd",
                            "format": "{:,.0f}",
                        },
                    },
                    "data": sample_metric_map,
                }

                table_headers = {
                    "raw_fastq__raw_reads": {"title": "Raw Reads", "format": "{:,.0f}"},
                    "trimmed_fastq__trimmed_reads": {"title": "Trimmed Reads", "format": "{:,.0f}"},
                    "trimmed_fastq__trimmed_reads_retained_pct": {"title": "Trimmed %", "suffix": "%", "format": "{:,.1f}"},
                    "mapper_report__mapper_input_reads": {"title": "Mapper Input Reads", "format": "{:,.0f}"},
                    "mapper_report__mapper_reported_aligned_reads": {"title": "Mapper Aligned Reads", "format": "{:,.0f}"},
                    "mapper_report__mapper_reported_alignment_pct": {"title": "Mapper Align %", "suffix": "%", "format": "{:,.1f}"},
                    "mapper_report__mapper_uniquely_mapped_reads": {"title": "Unique Mapped Reads", "format": "{:,.0f}"},
                    "mapper_report__mapper_multimapped_reads": {"title": "Multi-mapped Reads", "format": "{:,.0f}"},
                    "pre_filter__mapped_primary_reads": {"title": "Pre-filter Mapped Reads", "format": "{:,.0f}"},
                    "post_filter__mapped_primary_reads": {"title": "Post-filter Mapped Reads", "format": "{:,.0f}"},
                    "derived__mapped_primary_reads_retained_pct": {"title": "Filtered Retained %", "suffix": "%", "format": "{:,.1f}"},
                    "derived__mapped_primary_reads_removed": {"title": "Mapped Reads Removed", "format": "{:,.0f}"},
                    "pre_filter__duplicate_primary_reads": {"title": "Pre-filter Duplicate Reads", "format": "{:,.0f}"},
                    "post_filter__duplicate_primary_reads": {"title": "Post-filter Duplicate Reads", "format": "{:,.0f}"},
                    "derived__duplicate_primary_reads_removed": {"title": "Duplicate Reads Removed", "format": "{:,.0f}"},
                }
                if config["PAIRED"]:
                    table_headers.update({
                        "pre_filter__properly_paired_templates": {"title": "Pre Proper Pairs", "format": "{:,.0f}"},
                        "post_filter__properly_paired_templates": {"title": "Post Proper Pairs", "format": "{:,.0f}"},
                        "derived__properly_paired_templates_retained_pct": {"title": "Proper Pair Retained %", "suffix": "%", "format": "{:,.1f}"},
                        "pre_filter__discordant_templates": {"title": "Pre Discordant Pairs", "format": "{:,.0f}"},
                        "derived__discordant_templates_removed": {"title": "Discordant Pairs Removed", "format": "{:,.0f}"},
                    })

                table_yaml = {
                    "id": "omnomnomics_alignment_flow_table",
                    "section_name": "Omnomnomics Alignment Flow QC",
                    "description": "Aggregated raw, trimmed, mapper, and filtered alignment QC metrics collected across the pipeline.",
                    "plot_type": "table",
                    "pconfig": {
                        "id": "omnomnomics_alignment_flow_table_plot",
                        "title": "Omnomnomics Alignment Flow QC",
                    },
                    "headers": table_headers,
                    "data": sample_metric_map,
                }

                with open(generalstats_path, "w") as handle:
                    yaml.safe_dump(generalstats_yaml, handle, sort_keys=False)
                with open(table_path, "w") as handle:
                    yaml.safe_dump(table_yaml, handle, sort_keys=False)

                return generalstats_path, table_path

            stats_tsv_paths = [sample_qc_stats_path(sample_name) for sample_name in samples2]
            existing_stats_tsvs = [path for path in stats_tsv_paths if os.path.exists(path)]
            if existing_stats_tsvs:
                log_it(logfile, "Generating experiment-level alignment QC summary...", "ALIGNMENT QC SUMMARY")
                aggregate_tsv, aggregate_pdf, aggregate_svg, aggregate_rows = write_alignment_qc_experiment_summary(existing_stats_tsvs)
                log_it(logfile, f"Aggregate alignment QC table: {aggregate_tsv}")
                log_it(logfile, f"Aggregate alignment QC PDF: {aggregate_pdf}")
                log_it(logfile, f"Aggregate alignment QC SVG: {aggregate_svg}")
                multiqc_generalstats, multiqc_table = write_multiqc_custom_content(aggregate_rows)
                log_it(logfile, f"MultiQC custom general stats file: {multiqc_generalstats}")
                log_it(logfile, f"MultiQC custom table file: {multiqc_table}")
            else:
                log_it(logfile, "No sample-level alignment QC tables found. Skipping experiment-level alignment QC summary.", "ALIGNMENT QC SUMMARY")

            log_it(logfile, "Running multiQC...", "STATS")
            multiqc_version = subprocess.check_output(["multiqc", "--version"])
            log_it(logfile, "\n" + multiqc_version.decode("utf-8"), "MULTIQC VERSION")
            multiqc_command = (
                f"multiqc --filename MultiQC/omnomnomics.run.{run_date}.multiqc_report.html --dirs --export ."
            )
            log_it(logfile, multiqc_command, "MULTIQC COMMAND")
            shell(multiqc_command)
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
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['wig_rule_num']-1]}/*.extra_8.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 9 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['mergewig_rule_num']-1]}/extra_9.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 10 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['callpeaks_rule_num']-1]}/extra_10.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if 11 in themode:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['countreads_rule_num']-1]}/extra_11.tmp")
            for file in list_of_extra_files:
                os.remove(file)
        if create_homer_tagdirs:
            list_of_extra_files = glob.glob(f"{master_config['output_folders'][master_config['homer_tagdir_rule_num']-1]}/*.extra_13.tmp")
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
