#Main Snakefile
import os
import sys
import time
import yaml
import shutil
import glob
import random
from datetime import datetime


# Load the configuration file from command line arguments
configfile: config['config_file']


# Set global variables from the configuration file
OMNOM_HOME = config['OMNOM_HOME']
experiment_dir = config["EXPERIMENT_DIR"]
run_date = config["RUNDATE"]
logfile = os.path.join(experiment_dir, f"omnomnomics.run.{run_date}.log")


start_time = time.time()
previous_step_time = 0




# Function to log messages
def log_it(logfile, message, heading=None):
   timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
   with open(logfile, 'a') as log:
       if heading: #Check if heading
           log.write("\n{}\n\n".format(heading))
       log.write(f"{timestamp}: {message}\n") #write to log
   print(f"{timestamp}: {message}")


# Function to log elapsed time
### possibly use benchmarking to log elapsed time.
def log_elapsed(logfile, wildcards = None):
   global previous_step_time
   elapsed = (time.time() - start_time) / 60
   interval = elapsed - previous_step_time
   previous_step_time = elapsed
   log_it(logfile, f"Last step took: {interval:.2f} minutes...", "RUN DURATION")
   log_it(logfile, f"Time elapsed so far: {elapsed:.2f} minutes...")
   # if wildcards:
   #     log_it(f"Completed processing for sample: {wildcards.sample}")
   #return previous_step_time


##Commented out everything regarding omnomnomics.run_in_progress, most likely not going to use it anymore.
# Function to check for already running instance
#def check_instance(logfile):
   #Add a flag to show that we are running snakefile (the job)
   #with open(os.path.join(experiment_dir, f"omnomnomics.run_in_progress.{random.randint(1,1000)}"), 'w') as f:
       #pass
   #Check if semaphore file exists (if we are requed)
   #if not os.path.exists(os.path.join(experiment_dir, "omnomnomics.semaphore")):
       #If not requed, then check if snakefile is already running.
      # if os.path.exists(os.path.join(experiment_dir, "omnomnomics.run_in_progress")):
           #print("Omnomnomics pipeline instance already going to town on this directory! Aborting...")
          # log_it(logfile, "Omnomnomics pipeline instance already going to town on this directory! Aborting...")
           #################################################test if such a slurm job id is accessible here, i dont think so
           #with open(os.path.join(experiment_dir, f"omnomnomics.run_in_progress.job.{os.environ['SLURM_JOB_ID']}.aborted"), 'w') as f:
          # with open(os.path.join(experiment_dir, f"omnomnomics.run_in_progress.job.aborted"), 'w') as f:
          #     pass
          # sys.exit(1)
      # else:
           #Add a flag to show that we are running snakefile (the job)
      #     with open(os.path.join(experiment_dir, "omnomnomics.run_in_progress"), 'w') as f:
      #     pass


##Removed queued flag for now. Most likely not needed.
# Function to remove 'queued' flag
#def remove_queued_flag():
#    queued_flag = os.path.join(experiment_dir, "omnomnomics.run_queued")
#    if os.path.exists(queued_flag):
#        os.remove(queued_flag)


# Function for sanity check on directory
def sanity_check_dir(logfile, input_directory, file_ext):
    #check if input directiory exists
    if not os.path.isdir(input_directory):
        log_it(logfile, f"{file_ext} files should be contained in a {input_directory} folder inside your <EXPERIMENT_DIR> ({experiment_dir})!", "ERROR")
        log_it(logfile, "Aborting...")
        print(f"{file_ext} files should be contained in a {input_directory} folder inside your <EXPERIMENT_DIR> ({experiment_dir})! Aborting...")
    #        os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))
        sys.exit(1)
    #check if input directory contains input files
    if len([f for f in os.listdir(input_directory) if f.endswith(file_ext)]) == 0:
        log_it(logfile, f"The {input_directory} folder inside your <EXPERIMENT_DIR> ({experiment_dir}) does not contain any {file_ext} files!", "ERROR")
        log_it(logfile, "Aborting...")
        print(f"The {input_directory} folder inside your <EXPERIMENT_DIR> ({experiment_dir}) does not contain any {file_ext} files! Aborting...")
    #        os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))
        sys.exit(1)
    log_it(logfile, "Sanity check completed", "SANITY CHECK")


# Initialize micromamba shell and logging micromamba executable
def initialize_micromamba(logfile):
    log_it(logfile, "MAMBA_EXE: {}".format(os.environ.get('MAMBA_EXE', 'Not set')), "MICROMAMBA INIT")
    log_it(logfile, "Initializing shell...")
    os.system('eval "$(micromamba shell hook --shell=bash)"')


# Hold our horses for a little while to let the dispatch script initialise the log file
#############################################################################################necesarry? intialize it above?
time.sleep(0.1)


log_it(logfile,"JOB DISPATCHED!", "INITIALIZATION")
#check_instance(logfile)


##Removed queued flag for now. Most likely not needed.
#remove_queued_flag()




initialize_micromamba(logfile)
sanity_check_dir(logfile, os.path.join(experiment_dir, config["INPUT_FOLDER"]), config["INPUT_FILE_TYPE"])




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
#            os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))
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


CONFIG_FILE = os.path.join(f"{OMNOM_HOME}", "config.yaml")
# Check if CONFIG_FILE exists
if not os.path.isfile(CONFIG_FILE):
   log_it(logfile, "Master config file  does not exist! Aborting...", "ERROR")
   print(f"Master config file '{CONFIG_FILE}' does not exist. Please make sure it exists. Aborting...")
#    os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))
   sys.exit(1)


master_config = load_and_validate_yaml(logfile, CONFIG_FILE, "## Omnomnomics pipeline config ##")


##--------------------------------------------------------------------------------------------------------------
## Check if we are requeued or a fresh run
##---------------------------------------------------------------------------------------------------------------
def check_if_requeued(logfile, experiment_dir, max_step, themode):
   log_it(logfile, "Checking if we are requeued...", "Run START")


   semaphore_path = os.path.join(experiment_dir, "omnomnomics.semaphore")
   if os.path.isfile(semaphore_path):
       log_it(logfile, "Requeued run found!")


       ## Sanity check the semaphore file
       # Check if number of lines is within range of job types
       with open(semaphore_path, 'r') as sem_file:
           sem_lines = sem_file.readlines()


       semlen = len(sem_lines)
       log_it(logfile, f"Testing number of steps in semaphore file: {semlen}")


       if 1 <= semlen <= (max_step + 1):
           log_it(logfile, "Number of lines is within range of job types")
           print("Number of lines is within range of job types")
           log_it(logfile, "Testing job ranges in semaphore file...")
           print("Testing job ranges in semaphore file...")
           malsem = False
           for line in sem_lines:
               # Check if this line is within range of job types
               theline = int(line.strip())
               log_it(logfile, f"Testing job range for line: {theline}")
               print(f"Testing job range for line: {theline}")
               if not (1 <= theline <= max_step):
                   log_it(logfile, "Malformed line found in aborted run log!")
                   print("Malformed line found in aborted run log!")
                   log_it(logfile, semaphore_path)
                   malsem = True
                   break
               log_it(logfile, "Line is within range of job types")
               print("Line is within range of job types")
       else:
           log_it(logfile, "Malformed semaphore file! (too many lines)")
           print("Malformed semaphore file! (too many lines)")
           log_it(logfile, semaphore_path)
           malsem = True


       # Check if final step doesn't overflow the workflow
       if int(sem_lines[-1].strip()) >= max(themode) or int(sem_lines[-1].strip()) >= int(sem_lines[0].strip()):
           log_it(logfile, "Final recorded step in requeued run log implies a finished run, so we cannot resume this run!")
           malsem = True


       # Check the malformed semaphore file flag and act accordingly
       if malsem:
           log_it(logfile, "Cancelling run... Clear the experiment dir of cruft (e.g., the semaphore, unfinished steps, old logs, etc.) and resubmit!")
           print("Cancelling run... Clear the experiment dir of cruft (e.g., the semaphore, unfinished steps, old logs, etc.) and resubmit!")
#           os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))
           sys.exit(1)


       else:
           # Semaphore file is ok: continue using it
           #check if the aborted run ever completed single step
           if semlen == 1:
               log_it(logfile, "The requeued run log shows this run of the job has never actually completed a single step, therefore we will continue using the job mode specified in the config file!")
               log_it(logfile, f"Job mode set to: {themode}")
               log_it(logfile, "Starting a fresh run...")
               log_it(logfile, "Getting ready...")
               return themode
           else:
               # Aborted run completed steps so check where we left off
               finalstep = int(sem_lines[0].strip())
               laststep = int(sem_lines[-1].strip())


               # Resume from the next step
               nextstep = laststep + 1
               log_it(logfile, f"Resuming requeued run from step: {nextstep}")
               log_it(logfile, f"Previously set end point: {finalstep}")
               print(logfile, f"Resuming requeued run from step: {nextstep}")
               print(logfile, f"Previously set end point: {finalstep}")


               # Set the correct mode (need to overwrite previous config)
               themode = list(range(nextstep, finalstep + 1))
               themode_range_min = nextstep
               themode_range_max = finalstep
               log_it(logfile, f"Job mode set to: {themode}")
               log_it(logfile, "Starting a fresh run...")
               log_it(logfile, "Getting ready...")
               return themode


   else:
       # We are not requeued! Let's start fresh!
       log_it(logfile, "Starting a fresh run...")
       #log_it(f"omnomnomics run started on {os.environ.get('SLURMD_NODENAME', 'unknown node')}!")######################### not access this here?
       log_it(logfile, "Getting ready...")
       print("Starting a fresh run...")
       #print(f"omnomnomics run started on {os.environ.get('SLURMD_NODENAME', 'unknown node')}!")##########################can't access this here?
       print("Getting ready...")


       # Initialise semaphore file with the desired end step to keep track of progress
       with open(semaphore_path, 'w') as sem_file:
           sem_file.write(f"{max(themode)}\n")
       return themode




#themode = check_if_requeued(logfile, experiment_dir, master_config['max_step'], config['THEMODE'])
themode = config['THEMODE']


##########nned to return the mode here in any case for if its updated


##---------------------------------------------------------------------------------------------------------------
## Final housekeeping
##---------------------------------------------------------------------------------------------------------------
def final_housekeeping(logfile, thecoltable, experiment_dir):
   # Move color table to experiment dir if it exists
   if os.path.isfile(thecoltable):
       shutil.copy(thecoltable, experiment_dir)
       thecoltable = os.path.basename(thecoltable)
       log_it(logfile, f"Moved color table to experiment dir: {thecoltable}")


   # Move our camp to the experiment directory
   os.chdir(experiment_dir)
   logfile = os.path.basename(logfile)
   log_it(logfile, f"Changed directory to experiment dir: {experiment_dir}")
   print(f"Changed directory to experiment dir: {experiment_dir}")
   print(f"Updated logfile path: {logfile}")


   return thecoltable, logfile


col_table, logfile = final_housekeeping(logfile, config['THECOLTABLE'], experiment_dir)


##---------------------------------------------------------------------------------------------------------------
## Report some basic stats
##---------------------------------------------------------------------------------------------------------------
log_it(logfile, f"Job Start time: {start_time}", "RUN INFO")
log_it(logfile, f"Experiment dir: {experiment_dir}")
log_it(logfile, f"Run date: {run_date}")
log_it(logfile, f"Log file: {logfile}")
#log_it(logfile, f"Cores per sample: {config['MYCORES']}") #needed? kane ruit
#log_it(logfile, f"Available cores: {config['MAXCORES']}") #needed?
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
print(input_folder)
input_file_type =  master_config['input_file_types'][THEMODERANGEMIN-1]
print(input_file_type)


input_pattern = os.path.join(input_folder, f"*{input_file_type}")
print(input_pattern)
input_files = glob.glob(input_pattern)
print(input_files)
samples = [os.path.basename(f).replace(input_file_type, "") for f in input_files]
samples = [f.replace("_R1", "") for f in samples]
samples = [f.replace("_R2", "") for f in samples]
samples = [f.replace("_Skewer", "") for f in samples]
samples = [f.replace("_Trimmomatic", "") for f in samples]
print(samples)
samples = list(set(samples))
print(samples)

if config['PAIRED'] == 1 and THEMODERANGEMIN <= 4: 
    num_samples = len(samples) / 2
else: 
    num_samples = len(samples)

##--------------------------------------------------------------------------------------------------------------
# Obtain Threads and Memory per rule
##--------------------------------------------------------------------------------------------------------------

Threads_Per_Rule = {}
for i in range(1,master_config['max_step']+1):
    rule_num = i
    Threads_Per_Rule[f'{rule_num}'] = (max((master_config['mincores_single_sample_step1_9'][i - 1]
                                    if 'mincores_single_sample_step1_9' in master_config
                                    and isinstance(master_config['mincores_single_sample_step1_9'], list)
                                    and len(master_config['mincores_single_sample_step1_9']) > (i - 1)
                                    and master_config['mincores_single_sample_step1_9'][i - 1] is not None
                                    and isinstance(master_config['mincores_single_sample_step1_9'][i - 1], int)
                                    and master_config['mincores_single_sample_step1_9'][i - 1] > master_config['min_slice_cores']
                                    else master_config['min_slice_cores']),
                                    min((master_config['maxcores_single_sample_step1_9'][i - 1]
                                    if 'maxcores_single_sample_step1_9' in master_config
                                    and isinstance(master_config['maxcores_single_sample_step1_9'], list)
                                    and len(master_config['maxcores_single_sample_step1_9']) > (i - 1)
                                    and master_config['maxcores_single_sample_step1_9'][i - 1] is not None
                                    and isinstance(master_config['maxcores_single_sample_step1_9'][i - 1], int)
                                    and master_config['maxcores_single_sample_step1_9'][i - 1] < master_config['cores_per_node']
                                    else master_config['cores_per_node']), ((master_config['nodes_in_partition']*master_config['cores_per_node'])/num_samples)) ) )
print(Threads_Per_Rule)
Memory_Per_Rule = {}
for i in range(1,master_config['max_step']+1):
    rule_num = i
    Memory_Per_Rule[f'{rule_num}'] =  (master_config['min_mem_mb'][i - 1]
                                    if 'min_mem_mb' in master_config
                                    and isinstance(master_config['min_mem_mb'], list)
                                    and len(master_config['min_mem_mb']) > (i - 1)
                                    and master_config['min_mem_mb'][i - 1] is not None
                                    and isinstance(master_config['min_mem_mb'][i - 1], int)
                                    and master_config['min_mem_mb'][i - 1] > master_config['min_slice_mem']
                                    else (master_config['min_slice_mem'] 
                                        if('min_mem_mb' in master_config
                                            and isinstance(master_config['min_mem_mb'], list)
                                            and len(master_config['min_mem_mb']) > (i - 1)
                                            and master_config['min_mem_mb'][i - 1] is not None
                                            and isinstance(master_config['min_mem_mb'][i - 1], int)
                                            and master_config['min_mem_mb'][i - 1] <= master_config['min_slice_mem'])
                                        else Threads_Per_Rule[f'{rule_num}'] * master_config['max_mem_per_core_mb']))
                                        # else (master_config['mincores_single_sample_step1_9'][i - 1]
                                        # if 'mincores_single_sample_step1_9' in master_config
                                        # and isinstance(master_config['mincores_single_sample_step1_9'], list)
                                        # and len(master_config['mincores_single_sample_step1_9']) > (i - 1)
                                        # and master_config['mincores_single_sample_step1_9'][i - 1] is not None
                                        # and isinstance(master_config['mincores_single_sample_step1_9'][i - 1], int)
                                        # and master_config['mincores_single_sample_step1_9'][i - 1] > master_config['min_slice_cores']
                                        # else master_config['min_slice_cores'])* master_config['max_mem_per_core_mb']))
##--------------------------------------------------------------------------------------------------------------
# Include Snakemake rules for your actual data processing pipeline
##--------------------------------------------------------------------------------------------------------------

def check_and_include_rules(logfile, omnom_home, experiment_dir):
   log_it(logfile, "Checking snake rule functions...", "PREFLIGHT")

   rules_dir = os.path.join(omnom_home, "rules")

   # Check if the rules directory exists
   if not os.path.isdir(rules_dir):
       log_it(logfile, f"Snake rules dir {rules_dir} does not exist!", "WORKLFLOW SNAKE RULE CHECK")
       log_it(logfile, "Aborting...")
       print(f"Snake rules dir {rules_dir} does not exist!")
       print("Aborting...")
#        os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))
       sys.exit(1)

   # Check if the directory contains any .smk files
   smk_files = glob.glob(os.path.join("rules", "*.smk"))
   if len(smk_files) == 0:
       log_it(logfile, f"Snake Rules dir {rules_dir} does not contain any .smk files!", "WORKLFLOW SNAKE RULE CHECK")
       log_it(logfile, "Aborting...")
       print(f"Snake Rules dir {rules_dir} does not contain any .smk files!")
       print("Aborting...")
#       os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))
       sys.exit(1)

   # Sanity check if they contain the right header
   valid_smk_files = []
   for smk_file in smk_files:
       with open(smk_file, 'r') as file:
           lines = file.readlines()
           if len(lines) < 3 or lines[2].strip() != "## Omnomnomics Snake Rule  ##":
               log_it(logfile, f"Workflow step dir contains a malformed .smk file: {smk_file}", "WORKLFLOW SNAKE RULE CHECK")
               log_it(logfile, f"Aborting...")
               print(f"Workflow step dir contains a malformed .smk file: {smk_file}")
               print(f"Aborting...")
#                os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))
               sys.exit(1)
           valid_smk_files.append(smk_file)


   log_it(logfile, "Workflow Snake rules Good!", "WORKLFLOW SNAKE RULES CHECK")
   return valid_smk_files


valid_smk_files = check_and_include_rules(logfile, OMNOM_HOME, experiment_dir)
# Include snake rules in the main Snakefile
#for smk_file in valid_smk_files:
   #include: smk_file
include: "rules/trim_skewer.smk"

include: "rules/trim_trimmomatic.smk"

include: "rules/fastqc.smk"

include: "rules/align_reads_hisat2.smk"

include: "rules/align_reads_STAR.smk"

include: "rules/align_reads_STAR_TE.smk"

include: "rules/merge_lanes_and_clean_names.smk"

include: "rules/touchup_bam.smk"

include: "rules/index_bam.smk"

include: "rules/bam_stats.smk"

##--------------------------------------------------------------------------------------------------------------
# Execute the desired rules
##--------------------------------------------------------------------------------------------------------------

print(themode)
for rule_num in themode:
    # rule_ne = master_config['routines'][rule_num - 1]
    # routine = master_config[rule_name][config['SELECTED_ROUTINE_'+ rule_name.upper()]]
    # routines.append(routine)


    # input_folder = master_config['input_folders'][rule_num-1]
    # print(input_folder)
    output_folder = master_config['output_folders'][rule_num-1]
    print(output_folder)
    # input_file_type =  master_config['input_file_types'][rule_num-1]
    # print(input_file_type)
    output_file_type =  master_config['output_file_types'][rule_num-1]
    print(output_file_type)


    # input_pattern = os.path.join(input_folder, f"*{input_file_type}")
    # print(input_pattern)
    # input_files = glob.glob(input_pattern)
    # print(input_files)
    # samples = [os.path.basename(f).replace(input_file_type, "") for f in input_files]
    # print(samples)
    #all_outputs += expand(f"{output_folder}/rule{rule_num}_output.{{sample}}{output_file_type}", sample = samples)
    if rule_num == 1:
        if config['THETRIMTOOL'] == 'skewer':
            if config['PAIRED']: 
                all_outputs += expand(f"{output_folder}/{{sample}}_R1_Skewer.trimmed{output_file_type}", sample = samples)
                all_outputs += expand(f"{output_folder}/{{sample}}_R2_Skewer.trimmed{output_file_type}", sample = samples)
        else:
            if config['PAIRED']: 
                all_outputs += expand(f"{output_folder}/{{sample}}_R1_Trimmomatic.trimmed{output_file_type}", sample = samples)
                all_outputs += expand(f"{output_folder}/{{sample}}_R2_Trimmomatic.trimmed{output_file_type}", sample = samples)
    if rule_num == 2:
        if config['PAIRED']: 
            all_outputs += expand(f"{output_folder}/{{sample}}_R1{output_file_type[0]}", sample = samples)
            all_outputs += expand(f"{output_folder}/{{sample}}_R2{output_file_type[0]}", sample = samples)
            all_outputs += expand(f"{output_folder}/{{sample}}_R1{output_file_type[1]}", sample = samples)
            all_outputs += expand(f"{output_folder}/{{sample}}_R2{output_file_type[1]}", sample = samples)
        else:
            all_outputs += expand(f"{output_folder}/{{sample}}{output_file_type[0]}", sample = samples)
            all_outputs += expand(f"{output_folder}/{{sample}}{output_file_type[1]}", sample = samples)
    if rule_num == 3:
        if config['THEMAPTOOL'] == 'star':
            if config['PAIRED']: 
                all_outputs += expand(f"{output_folder}/{{sample}}_STAR.bam", sample = samples)
                all_outputs += expand(f"{output_folder}/{{sample}}.STAR_stats.txt", sample = samples)
        if config['THEMAPTOOL'] == 'hisat2':
            if config['PAIRED']: 
                all_outputs += expand(f"{output_folder}/{{sample}}_HISAT2.bam", sample = samples)
                all_outputs += expand(f"{output_folder}/{{sample}}.HISAT2_stats.txt", sample = samples)
        if config['THEMAPTOOL'] == 'star_te':
            if config['PAIRED']: 
                all_outputs += expand(f"{output_folder}/{{sample}}_STAR_TE.bam")
                all_outputs += expand(f"{output_folder}/{{sample}}.STAR_TE_stats.txt", sample = samples)
    #all_output.append("pipeline_completed.txt")
    print(all_outputs)
    #add output of final rule




#add that a log is made and semaphore file is updated after each run.
rule all:
   input:
       all_outputs
##################################################### all code after this, it processes before the snakemake rule all actually runs
##################################################### maybe once i use slurm it doesnt do that




# #Final things:
# #---------------------------------------------------------------------------------------------------------------
# # Run multiqc to gather all stats
# #---------------------------------------------------------------------------------------------------------------
# if config['NO_MULTIQC'] == "0":
#     log_it(logfile, "Running multiQC...", "STATS")
#     subprocess.run(
#         f"micromamba activate multiqc && "
#         f"multiqc --filename omnomnomics.run.{RUN_DATE}.multiqc_report.html --dirs --export . && "
#         f"micromamba deactivate",
#         shell=True,
#         check=True
#     )


# #---------------------------------------------------------------------------------------------------------------
# # Clean up: Compress and package the tag directories if needed
# #---------------------------------------------------------------------------------------------------------------
# if os.path.isdir("HOMER_tagDirs"):
#     log_it(logfile, "Compressing tagDirs...", "CLEANUP")
#     tag_dirs = glob.glob("HOMER_tagDirs/*tagDir")
  
#     for tag_dir in tag_dirs:
#         log_it(logfile, f"Making tar ball of {tag_dir}...")
#         tar_output = f"HOMER_tagDirs/{os.path.basename(tag_dir)}.tar.gz"
#         subprocess.run(["tar", "czf", tar_output, tag_dir], check=True)
  
#     log_it(logfile, "Waiting for tar to complete...")
  
#     for tag_dir in tag_dirs:
#         log_it(logfile, f"Deleting uncompressed tagDir {tag_dir}...")
#         os.rmdir(tag_dir)


#     log_it(logfile, "Cleanup complete")


# #================================================================================================================
# # Log completion
# log_it(logfile, "All done!" "FINAL REMARKS")
# #ELAPSED=(time.time() - start_time)/60s
# #log_it(logfile, f"Final run time: {ELAPSED:.2f} minutes.")
# log_it(logfile, "Good luck with your downstream analyses!")
#os.remove(os.path.join(experiment_dir, "omnomnomics.run_in_progress"))

