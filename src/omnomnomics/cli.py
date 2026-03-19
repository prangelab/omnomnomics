#!/usr/bin/env python3

##############################################
### NGS PIPELINE FOR RNA, ATAC & CHIP DATA ###
##############################################
#=============================================
# Author: Kieran Carroll
# Affiliation: Prangelab AMC / Amsterdam UMC's Core Facility Genomics
# Copyright PrangeLab 2024 ##
#=============================================
import os
import sys
import subprocess
import yaml
import argparse
import re
import glob
import random
import math
import shutil
from pathlib import Path
from datetime import date, datetime

PACKAGE_ROOT = Path(__file__).resolve().parent
WORKFLOW_ROOT = PACKAGE_ROOT / "workflow"
DEFAULT_WORKFLOW_CONFIG = WORKFLOW_ROOT / "config" / "workflow.yaml"
DEFAULT_SITE_CONFIG = WORKFLOW_ROOT / "config" / "site.yaml"

def parse_arguments():
   #Parse command-line arguments
   parser = argparse.ArgumentParser(
       description="Modular HPC pipeline for RNA-seq, ATAC-seq, and ChIP-seq processing.",
       allow_abbrev=False,
   )

   # Define command-line options
   parser.add_argument('-i', '--experiment-dir', help='Path to the experiment directory')
   parser.add_argument('-t', '--type', help='Type of experiment: RNA, ChIP, ATAC')
   parser.add_argument('-g', '--genome', help='Genome version. Avaliable versions: \n \t UCSC/RefSeq/NCBI: mm10, mm39, hg38. \n" "\t ENSMBL/GenBank: GRCh38.p14, GRCm39')
   parser.add_argument('-j', '--mode', help='Job mode. Can be auto, all or a range of jobs. See readme for some examples. \n \t Default: auto')
   parser.add_argument('-T', '--trim-tool', help='Trimming tool choice. can be Skewer or Trimmomatic \n \t Default: Skewer')
   parser.add_argument('-M', '--map-tool', help='Mapping tool choice. Can be HISAT2, STAR, or STAR_TE. STAR(_TE) can only be used for RNA-seq data. \n \t Default: HISAT2')
   parser.add_argument('-f', '--formula', help='RNA: Experimental Design for DE calling. \n \t Default: 1 (just an intercept)')
   parser.add_argument('-I', '--input', help='Input file used for ChIP peak calling. Has to be a .bam file or HOMER tag directory. \n \t Default: do not use input')
   parser.add_argument('-m', '--metadata', help='.txt file with columns of metadata for RNA-seq experiments. \n \t Default: DESeq2 style metadata table describing all samples. Rownames should be samplenames')
   parser.add_argument('-b', '--broad', action='store_true', help='ChIP: Peak calling style. If set, use MACS3 in --broad mode, and use HOMER findPeaks with -size, -minDist and -region settings. Works best with STYLE histone.')
   parser.add_argument('-S', '--style', help='ChIP: Peak calling style for HOMER peak calling. Can be factor or histone. \n \t Default: factor')
   parser.add_argument('-C', '--col-table', help='File specifying which colors to use for the tracks. \n \t Default: gray.tint.color.table. Can be a *txt list file with one color table per line. Different color tables will be used per hub as split by -e. Can be a full (relative) path to a file or a file basename only in conjuction with -P. Use createTrackColorTable.sh to roll your own. Use displayTrackColorTable.sh to visualize existing color tables.')
   parser.add_argument('-P', '--color-data-folder', help='Path to a folder with color tables. \n \t Default: dewintherlab/bin/color_data_for_hubs') ##################change this default??
   parser.add_argument('-o', '--overlay', help='Overlay type (transparentOverlay|stacked|solidOverlay|none) \n \t Default: transparentOverlay')
   parser.add_argument('-L', '--hub-mail', help='Email to use in trackhub \n \t Default: m.dewinther@amsterdamumc.nl')
   parser.add_argument('-X', '--no-multiqc', action='store_true', help='Exclude multiQC stats aggregator. Set if you don not wish to run multiQC.')
   parser.add_argument('-n', '--name-fields', help='Field(s) in filename to use as track name, peak file name, and column header in the count table \n \t Default: 1-3')
   parser.add_argument('-e', '--type-field', help='Field(s) in filename to use as merged hub or peak calling group identifier \n \t Default: 1 (Creates separate merged hubs for each unique entry)')
   parser.add_argument('-c', '--col-field', help='Field(s) in filename to use as color type \n \t Default: 2')
   parser.add_argument('-s', '--separator', help='Separator used in file names. \n \t Default: _')
   parser.add_argument('-a', '--appendix', help='Appendix to add to track name \n \t Default: hub')
   parser.add_argument('-d', '--homer-mindist', help='Minimum distance bewteen peaks. Used for merging regions in HOMER broad peak calling mode. \n \t Default: 2000')
   parser.add_argument('-z', '--homer-size', help='Minimum peak size. Used for defining peaks in HOMER broad peak calling mode. \n \t Default: 500')
   parser.add_argument('-k', '--keepunpaired', action='store_true', help='Keep unpaired or not in HISAT2')
   parser.add_argument('--site-config', help='Optional path to a site-specific config YAML. Default: packaged site config')

   args, unknown = parser.parse_known_args()
   # Check if any unknown arguments are provided
   if unknown:
       parser.error(f"Unrecognized arguments: {' '.join(unknown)}")

   # Check if at least one argument is provided besides the script name
   if len(sys.argv) < 2: 
       print("Not enough input arguments...")
       parser.print_help()
       sys.exit(1)

   return args

# Function to check if the configuration file contains the proper header
def check_config_file_header(config_file_path, expected_header):
   with open(config_file_path, 'r') as config_file:
      
       lines = config_file.readlines()
       if len(lines) < 3 or expected_header not in lines[2]:
           print("Config file does not contain right header. Aborting...")
           sys.exit(1)

# Function to load and validate the YAML configuration file
def load_and_validate_yaml(config_file_path, expected_header):
  
   # Check for the expected header
   check_config_file_header(config_file_path, expected_header)
   # Load the YAML content
   with open(config_file_path, 'r') as file:
       try:
           config_content = yaml.safe_load(file)
       except yaml.YAMLError as exc:
           print(f"Error parsing YAML file: {exc}")
           sys.exit(1)
  
   return config_content

def merge_configs(base_config, override_config):
   merged_config = dict(base_config)
   merged_config.update(override_config)
   return merged_config

##---------------------------------------------------------------------------------------------------------------
## Set required vars or die
##---------------------------------------------------------------------------------------------------------------
def check_required_vars(the_type, experiment_dir, genome, config):
   #Check if required variables are set and valid
   print("Checking experiment dir...")
   if not experiment_dir:
       print("No EXPERIMENT DIR (-i) given! Aborting...", file=sys.stderr)
       sys.exit(1)

   if not os.path.isdir(experiment_dir):
       print("EXPERIMENT DIR (-i) has to be a directory! Aborting...", file=sys.stderr)
       sys.exit(1)

   print("Checking experiment type...")
   if not the_type:
       print("No EXPERIMENT TYPE (-t) specified! Aborting...", file=sys.stderr)
       sys.exit(1)
   else:
       if the_type not in config['TYPES']:
           print(f"Experiment type (-t) should be one of: {', '.join(config['TYPES'])}. Aborting...", file=sys.stderr)
           sys.exit(1)

   print("Checking genome version...")
   if not genome:
       print("No GENOME VERSION (-g) specified! Aborting...", file=sys.stderr)
       sys.exit(1)
   else:
       if genome not in config['GENOMES']:
           print(f"Genome (-g) should be one of: {', '.join(config['GENOMES'])}. Aborting...", file=sys.stderr)
           sys.exit(1)

##---------------------------------------------------------------------------------------------------------------
## Set user subroutine choices
##---------------------------------------------------------------------------------------------------------------
def set_user_subroutine_choices(trim_tool, map_tool, config):
   #Set user subroutine choices
   print("Applying tool selection...")
   if trim_tool:
       trim_tool = trim_tool.lower()
       if trim_tool in config['TRIMMERS']:
           selected_routine_trim = config['TRIMMERS'].index(trim_tool)
           print(f"Trim tool: {trim_tool}")
       else:
           print(f"Trimming tool -T has to be one of: {', '.join(config['TRIMMERS'])}. Aborting...", file=sys.stderr)
           sys.exit(1)
   else:
       selected_routine_trim = config['selected_routine_trim']

   if map_tool:
       map_tool = map_tool.lower()
       if map_tool in config['MAPPERS']:
           selected_routine_map = config['MAPPERS'].index(map_tool)
           print(f"Map tool: {map_tool}")
       else:
           print(f"Mapping tool -M has to be one of: {', '.join(config['MAPPERS'])}. Aborting...", file=sys.stderr)
           sys.exit(1)
   else:
       selected_routine_map = config['selected_routine_map']

   return selected_routine_trim, selected_routine_map
  
##---------------------------------------------------------------------------------------------------------------
## Validate user defined variables
##---------------------------------------------------------------------------------------------------------------
def validate_user_defined_vars(OMNOM_HOME, metadata, experiment_dir, INPUT, the_style, color_data_folder, col_table, overlay, the_type, map_tool, homer_size, homer_mindist, config):
    #Validate user-defined variables
    print("Validating options...")


    # Check metadata file
    if metadata != "NA":
        if not os.path.isfile(metadata):
            print("Metadata file (-m) does not exist! Aborting...", file=sys.stderr)
            sys.exit(1)
        elif not metadata.endswith(".txt"):
            print("Metadata file (-m) is not a .txt file! Aborting...", file=sys.stderr)
            sys.exit(1)
        elif not metadata.startswith(experiment_dir):
            print("Please place your metadata file in your experiment dir!", file=sys.stderr)
            print(f"e.g., -m {os.path.join(f'{OMNOM_HOME}', 'data', 'me', 'my_fantastic_experiment', 'my_mindboggling_metadata.txt')}. Aborting...", file=sys.stderr)
            sys.exit(1)


    # Check input file
    if INPUT != "NA":
        if not os.path.isfile(INPUT):
            print("ChIP/ATAC Input .bam file (-I) does not exist! Aborting...", file=sys.stderr)
            sys.exit(1)
        elif not INPUT.endswith(".bam"):
            print("ChIP/ATAC Input file (-I) is not a .bam file! Aborting....", file=sys.stderr)
            sys.exit(1)
        elif not os.path.isabs(INPUT):
            print("Please provide the full (absolute) path to your input file!", file=sys.stderr)
            print(f"e.g., -I {os.path.join(f'{OMNOM_HOME}', 'genomes', 'input_ChIP', 'my_awesome_input.bam')}. Aborting...", file=sys.stderr)
            sys.exit(1)

        if not os.path.isdir(os.path.join(os.path.dirname(INPUT), f"{os.path.basename(INPUT).replace('.bam', '')}.HOMER_tagDir")):
            print(f"WARNING: No corresponding input HOMER tagDir found! ({homer_input}) does not exist...", file=sys.stderr)
            print("WARNING: Will run HOMER findPeaks without input...", file=sys.stderr)
            homer_input = "NA"
        else:
            homer_input = os.path.join(os.path.dirname(INPUT), f"{os.path.basename(INPUT).replace('.bam', '')}.HOMER_tagDir")
    else:
        homer_input = config['homer_input']

    # Convert peak style to lowercase if given on command line, else take the default from config file
    the_style = the_style.lower()


    # Check peak style
    if the_style not in ["factor", "histone"]:
        print("ChIP peak calling style has to be factor or histone! Aborting...", file=sys.stderr)
        sys.exit(1)


    # Check color data folder
    if not os.path.isdir(color_data_folder):
        print(f"Color table folder (-P) ({color_data_folder}) does not exist! Aborting...", file=sys.stderr)
        sys.exit(1)


    # Remove trailing slash from color data folder path if necessary
    if color_data_folder.endswith('/'):
        color_data_folder = color_data_folder.rstrip('/')


    # Check if a full path to a color table (list) file was provided. If not, prepend the color folder path.
    if not os.path.isfile(col_table):
        col_table = os.path.join(color_data_folder, col_table)
        if not os.path.isfile(col_table):
            print(f"Color table file: {col_table} not found! Aborting...", file=sys.stderr)
            sys.exit(1)
    else:
        # If it's a list file, parse the contents to see if we need to prepend the color folder path
        if col_table.endswith(".txt"):
            with open(col_table, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    line = line.strip()
                    if not os.path.isfile(line):
                        # Check if we have an absolute path already
                        if os.path.isabs(line):
                            print(f"Malformed entry in color table list file: {line} is not a file!", file=sys.stderr)
                            sys.exit(1)
                        # Check if we have a valid entry if we prepend the color data folder
                        elif not os.path.isfile(os.path.join(color_data_folder, line)):
                            print(f"Malformed entry in color table list file: {line} does not exist here or in {color_data_folder}", file=sys.stderr)
                            sys.exit(1)
                        else:
                        #If yes, prepend the color data folder path
                            with open(col_table, 'w') as file:
                                for line in lines:
                                    stripped_line = line.strip()
                                    if re.match(rf'^{re.escape(stripped_line)}$', stripped_line):
                                        modified_line = re.sub(rf'^{re.escape(stripped_line)}$', f'{color_data_folder}/{stripped_line}', stripped_line)
                                        file.write(modified_line + '\n')
                                    else:
                                        file.write(line)

    # Check overlay type
    if overlay and overlay not in config["overlaytypes"]:
        print(f"Overlay type -o has to be one of: {', '.join(config['overlaytypes'])}. Aborting...", file=sys.stderr)
        sys.exit(1)

    # Check valid mapper for the experiment type
    if the_type and the_type != "RNA" and "star" and map_tool and "star" in map_tool :
        print("STAR read aligner can only be used with RNA-seq data! Aborting...", file=sys.stderr)
        sys.exit(1)

    # Check if HOMER peak size and mindist are integers
    if homer_size and not isinstance(homer_size, int):
        print("HOMER peak size (-z) has to be an integer! Aborting...", file=sys.stderr)
        sys.exit(1)
    if homer_mindist and not isinstance(homer_mindist, int):
        print("HOMER minimal distance (-d) has to be an integer! Aborting...", file=sys.stderr)
        sys.exit(1)
    return homer_input

##---------------------------------------------------------------------------------------------------------------
## Finetune some vars
##---------------------------------------------------------------------------------------------------------------
def setup_variables(experiment_dir,config):
    #Setting up variables
    print("Setup variables...")


    # Remove trailing slash from experiment directory path if necessary
    if experiment_dir.endswith('/'):
        experiment_dir = experiment_dir.rstrip('/')

    # Set the date
    run_date = date.today().isoformat()
    return run_date

##---------------------------------------------------------------------------------------------------------------
## Set job mode
##---------------------------------------------------------------------------------------------------------------
def expand_range(mode):
   range_list = []
   ranges = mode.split(',') #Get all the different ranges
   for part in ranges:
       if '-' in part:
           start, end = map(int, part.split('-')) # Get the start and end of a range
           range_list.extend(range(start, end + 1)) # Add the wanted steps to the list of steps
       else:
           range_list.append(int(part)) # Add the wanted steps to the list of steps
   return range_list


def set_job_mode(args, config, experiment_dir, mode):
    #Set job mode
    print("Setup job mode...")
    max_step = config['max_step']

    next_step = None
    final_step = None

    if mode == "auto":
        mode = "all"

    if mode == "all":
    # If we are running the whole pipeline, set the maximum range
        mode_range_min = 1
        mode_range_max = max_step
        mode = f"1-{max_step}"
        mode_steps = list(range(1, max_step + 1))
    else:
    # Handle different job mode patterns (single step, single range, complex patterns)
        #single step pattern
        if re.match(r'^[0-9]+$', mode):
            mode_range_min = int(mode)
            mode_range_max = int(mode)
            
            # Check if it's a valid range
            if mode_range_min < 1 or mode_range_max > max_step:
                print(f"Job mode range should lie between 1 and {max_step}! See help (-h) for job type usage info. Aborting...")
                sys.exit(1)
            
            mode_steps = [mode_range_min]
        
        #single range pattern
        elif re.match(r'^[0-9]+-[0-9]+$', mode):
            mode_range_min, mode_range_max = map(int, mode.split('-'))
            
            # Check if it's a valid range
            if mode_range_min > mode_range_max:
                print("Job mode range has to increase from start to end! Aborting...")
                sys.exit(1)
            
            if mode_range_min < 1 or mode_range_max > max_step:
                print(f"Job mode range should lie between 1 and {max_step}! See help (-h) for job type usage info. Aborting...")
                sys.exit(1)
            
            mode_steps = list(range(mode_range_min, mode_range_max + 1))
        
        #complex pattern
        else:
            # Expand the range for complex patterns
            mode_steps = expand_range(mode)
            
            # Check if it's valid (i.e., only numbers separated by commas). If not, exit!
            if not all(isinstance(step, int) for step in mode_steps):
                print("Job mode should be 'all', 'auto', a single step number, a numeric range of steps separated by a dash (e.g., 1-12), or a number of ranges and steps separated by commas (e.g., 1-3,6,8,10-12,14).")
                print("Aborting...")
                sys.exit(1)
            
            mode_range_min = mode_steps[0]
            mode_range_max = mode_steps[-1]
            
            # Check if it's a valid range
            if mode_range_min < 1 or mode_range_max > max_step:
                print(f"Job mode range should lie between 1 and {max_step}! See help (-h) for job type usage info. Aborting...")
                sys.exit(1)
            
            # Check if the job mode range increases from start to end
            if any(mode_steps[i] >= mode_steps[i+1] for i in range(len(mode_steps) - 1)):
                print("Job mode range has to increase from start to end! Aborting...")
                sys.exit(1)
    mode_steps = sorted(mode_steps)
    return mode_steps, config

##---------------------------------------------------------------------------------------------------------------
## Set some parameters
##---------------------------------------------------------------------------------------------------------------
def validate_input_files(the_type, config, mode_range_min, experiment_dir):
    print("Validating input files...")

    # Check permissions on experiment dir
    if not (os.access(experiment_dir, os.R_OK) and os.access(experiment_dir, os.W_OK)):
        print(f"Permission error! {experiment_dir} is not readable and/or writable! Aborting...", file=sys.stderr)
        sys.exit(1)

    # Check if we have a paired-end run (if we are dealing with FASTQs) and set the flag accordingly
    paired = False
    fastq_files = glob.glob(f"{experiment_dir}/FASTQ/*_R2*fastq.gz") + glob.glob(f"{experiment_dir}/trimmed_FASTQ/*_R2*fastq.gz")
    if len(fastq_files) > 0:
        paired = True

    input_file_types = config['input_file_types']
    input_folders = config['input_folders']

    # Check if the input dir exists for the first step and count files with the correct extension
    num_files = 0
    input_folder_mod_range_min = input_folders[mode_range_min - 1]
    input_file_type_mod_range_min = input_file_types[mode_range_min - 1]

    #For rules with multiple input filetypes, set it to the right one
    if mode_range_min == 10:
        if the_type == "RNA":
            input_file_type_mod_range_min = input_file_type_mod_range_min[0]
        else:
            input_file_type_mod_range_min = input_file_type_mod_range_min[1]
    if mode_range_min == 11:
        input_file_type_mod_range_min = input_file_type_mod_range_min[0]
        input_folder_mod_range_min = input_folder_mod_range_min[0]
    if mode_range_min == 12:
        input_file_type_mod_range_min = input_file_type_mod_range_min[0]
        input_folder_mod_range_min = input_folder_mod_range_min[0]


    if os.path.isdir(f"{experiment_dir}/{input_folder_mod_range_min}"):
        num_files = len(glob.glob(f"{experiment_dir}/{input_folder_mod_range_min}/*{input_file_type_mod_range_min}"))


    # Sanity check file number
    if num_files == 0:
        if mode_range_min == 13 and the_type == "CHIP":
                pass
        elif mode_range_min == 12 and the_type == "CHIP":
                pass
        else: 
            print("No input files detected! Aborting...", file=sys.stderr)
            sys.exit(1)
   
    if (mode_range_min == 13 or mode_range_min == 12) and the_type == "CHIP":
        pass
    else:
        # Check if input files are readable
        input_files = glob.glob(f"{experiment_dir}/{input_folder_mod_range_min}/*{input_file_type_mod_range_min}")
        if not os.access(input_files[0], os.R_OK):
            print(f"Permission error! {input_file_type_mod_range_min} files in {experiment_dir}/{input_folder_mod_range_min} are not readable! Aborting...", file=sys.stderr)
            sys.exit(1)


    # Set the number of pairs (if dealing with FASTQs, else just keep it equal to file number)
    num_pairs = num_files // 2 if paired and mode_range_min < 4 else num_files

    return num_files, num_pairs, paired, input_folder_mod_range_min, input_file_type_mod_range_min

def parse_name_fields(fields):
    #Parse a string of fields like '1,2,4-6' into a list of integers
    field_list = []
    for part in fields.split(','): # Split all the fields
        if '-' in part:
            start, end = part.split('-') # Split a range of steps
            field_list.extend(range(int(start), int(end) + 1)) # Add all the wanted steps to list
        else:
            field_list.append(int(part)) # Add all the wanted steps to list
    return field_list

def run_cut_command(filename, fields, separator):
    #Select specified fields from a filename using the given separator
    try:
        parts = filename.split(separator) # Get all the parts in the filename
        selected_parts = []
        for field in parse_name_fields(fields): # Loop over all the fields
            field_index = field - 1 # -1 since indexing of field names starts at 1
            if field_index < 0 or field_index >= len(parts):
                return f"Error: Field index {field} out of range"
            selected_parts.append(parts[field_index]) # Add to the name
        return ''.join(selected_parts)
    except ValueError:
        return "Error: Invalid field value"

def check_name_field_settings(experiment_dir, separator, name_fields, type_field, col_field, config, input_folder_mod_range_min, input_file_type_mod_range_min):
    try:
        # Get a sample file
        sample_files = glob.glob(f"{experiment_dir}/{input_folder_mod_range_min}/*{input_file_type_mod_range_min}")
        if not sample_files:
            raise ValueError("No sample files found.")
        
        mock_cut = os.path.basename(sample_files[0])
        print(f"Debug: Mock cut - {mock_cut}")  # Debugging line
        
        # Run the cut command simulations
        cut_test = ""
        cut_test += run_cut_command(mock_cut, name_fields, separator)
        cut_test += run_cut_command(mock_cut, type_field, separator)
        cut_test += run_cut_command(mock_cut, col_field, separator)
        
        print(f"Debug: Cut test - {cut_test}")  # Debugging line
        
    except (IndexError, ValueError) as e:
        print("Oops, something is wrong with your field settings! Check your -n, -c, -e, and -s variables!", file=sys.stderr)
        print(f"Error message: {e}")
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

    # If cut_test contains "Error", there was an error
    if "Error" in cut_test:
        print("Oops, something is wrong with your field settings! Check your -n, -c, -e, and -s variables!", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)


def check_unique_sample_names(experiment_dir, input_folder_mod_range_min, input_file_type_mod_range_min, name_fields, separator):
    def parse_fields(field_string):
        #Parse a string representing fields into a list of integers
        fields = []
        for part in field_string.split(','):
            if '-' in part:
                start, end = map(int, part.split('-'))
                fields.extend(range(start, end + 1))
            else:
                fields.append(int(part))
        return fields


    def extract_fields(filename, fields, separator):
        #Extract specified fields from a filename
        parts = filename.split(separator)
        try:
            return separator.join([parts[i - 1] for i in fields])
        except IndexError:
            raise ValueError(f"Invalid field indices: {fields}")
        
    name_fields = parse_fields(name_fields)

    if input_file_type_mod_range_min != ".fastq.gz" and input_file_type_mod_range_min != ".trimmed.fastq.gz" :
        # Get all files with the specified type
        files = glob.glob(f"{experiment_dir}/{input_folder_mod_range_min}/*{input_file_type_mod_range_min}")
    else:
        # Only consider R1 files for FASTQs
        files = glob.glob(f"{experiment_dir}/{input_folder_mod_range_min}/*_R1*{input_file_type_mod_range_min}")

    # Extract names using the specified fields
    sample_names = [extract_fields(os.path.basename(f), name_fields, separator) for f in files]
    print(sample_names)

    # Check for uniqueness
    if len(sample_names) != len(set(sample_names)):
        print(f"Error: using name fields: {name_fields} does not yield unique sample names!", file=sys.stderr)
        print("Please choose a different range of name fields.", file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)

def setup_runtime_parameters(num_pairs, experiment_dir):
    print("Setup runtime parameters...")

    # Make directories for slurm logs, run logs, run configs and MultiQC in experiment directory
    subprocess.run(f"mkdir -p {experiment_dir}/slurm_logs", shell=True, check=True)
    subprocess.run(f"mkdir -p {experiment_dir}/run_logs", shell=True, check=True)
    subprocess.run(f"mkdir -p {experiment_dir}/run_configs", shell=True, check=True)
    subprocess.run(f"mkdir -p {experiment_dir}/MultiQC", shell=True, check=True)

    # Node has 512GiB memory. Keep some margin, use 500 GiB ~ 500 000 MiB. Set Java heap initial size (Xms) to half of the max heap (THEMEM) we calculate per file
    the_mem = 500000 // num_pairs
    the_heap_init = the_mem // 2
    the_mem = f"{the_mem}M"
    the_heap_init = f"{the_heap_init}M"


    return  the_mem, the_heap_init 


def write_run_config(experiment_dir, run_date, config_data):
    #Write the run configuration to a YAML file
    run_config_filename = os.path.join(experiment_dir, "run_configs", f"omnomnomics.run.{run_date}.config.yaml")

    # If the config file already exists, back it up with a random suffix
    if os.path.isfile(run_config_filename):
        backup_filename = f"{run_config_filename}.{random.randint(0, 9999)}.backup"
        os.rename(run_config_filename, backup_filename)

    # Write the new config data to the YAML file
    with open(run_config_filename, 'w') as yaml_file:
        yaml.dump(config_data, yaml_file, default_flow_style=False)

    print(f"Run config file written to {run_config_filename}")

def start_log(experiment_dir, run_date, config):
    #Initialize the log file
    log_file = os.path.join(experiment_dir, "run_logs", f"omnomnomics.run.{run_date}.log")
    if os.path.isfile(log_file):
        backup_log = f"{log_file}.{os.urandom(8).hex()}.backup" #Create a random backup if multiple runs
        os.rename(log_file, backup_log)
        print(f"Existing log file backed up as: {backup_log}")


    with open(log_file, 'w') as log:
        log.write("#################################\n")
        log.write("## Run log for omnomnomics run ##\n")
        log.write("#################################\n\n")
        log.write(f"Pipeline version: {config['omnomnomics']}\n\n")


    return log_file 


##--------------------------------------------------------------------------------------------------------------
# Remove already present outputs of rules that you want to run
##--------------------------------------------------------------------------------------------------------------
def delete_outputs_to_be_updated(mode_steps, config, experiment_dir):
    print("DELETING TO BE UPDATED OUTPUT FILES")
    for num in mode_steps: # Loop over all the to run steps
        outputfolder = config['output_folders'][num-1]
        output_filetype = config['output_file_types'][num-1]
        if isinstance(output_filetype, list): # If multiple output filetypes
            for filetype in output_filetype:
                files = glob.glob(f"{experiment_dir}/{outputfolder}/*{filetype}")
                for file in files:
                    if os.path.exists(file):
                        if filetype == ".hub":
                            shutil.rmtree(file) # Is actually a hub directory and not a file
                        else:
                            os.remove(file)
        else:
            files = glob.glob(f"{experiment_dir}/{outputfolder}/*{output_filetype}")
            for file in files:
                if os.path.exists(file):
                    if num == 4 and output_filetype == ".bam" and re.search(r'L0\d+', os.path.basename(file)):
                        continue
                    elif num == 10:
                        shutil.rmtree(file) # Is actually a hub directory and not a file
                    else:
                        os.remove(file)
##--------------------------------------------------------------------------------------------------------------
# Main function
##--------------------------------------------------------------------------------------------------------------
def main():
    # Define variables for packaged workflow paths
    OMNOM_HOME = str(WORKFLOW_ROOT)
    workflow_config_file = DEFAULT_WORKFLOW_CONFIG

    #Check if packaged workflow root exists
    if not WORKFLOW_ROOT.is_dir():
        print(f"Packaged workflow directory '{WORKFLOW_ROOT}' does not exist. Aborting...")
        sys.exit(1)

    # Check if packaged workflow config exists
    if not workflow_config_file.is_file():
        print(f"Workflow config file '{workflow_config_file}' does not exist. Please make sure it exists. Aborting...")
        sys.exit(1)

    # Capture command line invocation
    the_command = ' '.join(sys.argv)
        
    # Parse command-line arguments
    args = parse_arguments()

    site_config_file = Path(args.site_config).expanduser().resolve() if args.site_config else DEFAULT_SITE_CONFIG
    if not site_config_file.is_file():
        print(f"Site config file '{site_config_file}' does not exist. Please make sure it exists. Aborting...")
        sys.exit(1)

    workflow_config = load_and_validate_yaml(workflow_config_file, "## Omnomnomics pipeline config ##")
    site_config = load_and_validate_yaml(site_config_file, "## Omnomnomics pipeline config ##")
    config = merge_configs(workflow_config, site_config)

    # Access parsed arguments
    experiment_dir = args.experiment_dir
    the_type = args.type.upper()
    genome = args.genome
    trim_tool = args.trim_tool.lower() if args.trim_tool else config.get("trim_tool","skewer").lower()
    map_tool = args.map_tool.lower() if args.map_tool else config.get('map_tool', "hisat2").lower()
    mode = args.mode.lower() if args.mode else config.get('mode', "auto")
    formula = args.formula if args.formula else config.get('formula', "1")
    broad = args.broad if args.broad else config.get('broad', "NA")
    INPUT = args.input if args.input else config.get('input',"NA")
    metadata = args.metadata if args.metadata else config.get('metadata', "NA")
    style = args.style if args.style else config.get('the_style', "factor")
    col_table = args.col_table.replace("{OMNOM_HOME}", OMNOM_HOME) if args.col_table else config.get('color_table', f"{OMNOM_HOME}/bin/color_data_for_hubs/gray.tint.color.table").replace("{OMNOM_HOME}", OMNOM_HOME)
    color_data_folder = args.color_data_folder.replace("{OMNOM_HOME}", OMNOM_HOME) if args.color_data_folder else config.get('color_data_folder', f"{OMNOM_HOME}/bin/color_data_for_hubs").replace("{OMNOM_HOME}", OMNOM_HOME)
    overlay = args.overlay if args.overlay else config.get('overlay', "transparentOverlay")
    hub_mail = args.hub_mail if args.hub_mail else config.get('hub_mail', "m.dewinther@amsterdamumc.nl")
    no_multiqc = args.no_multiqc if args.no_multiqc else config.get('no_multiqc', 0)
    name_fields = args.name_fields if args.name_fields else config.get('name_fields', "1-3")
    type_field = args.type_field if args.type_field else config.get('type_field', "1")
    col_field = args.col_field if args.col_field else config.get('col_field', "2")
    separator = args.separator if args.separator else config.get('separator', "_")
    appendix = args.appendix if args.appendix else config.get('appendix', "hub")
    homer_mindist = args.homer_mindist if args.homer_mindist else config.get('homer_mindist', 2000)
    homer_size = args.homer_size if args.homer_size else config.get('homer_size', 500)
    keep_unpaired = args.keepunpaired if args.keepunpaired else config.get('keep_unpaired', False)

    # Check required variables
    check_required_vars(the_type, experiment_dir, genome, config)

    # Set user subroutine choices
    selected_routine_trim, selected_routine_map = set_user_subroutine_choices(trim_tool, map_tool, config)
    selected_routines = {}
    selected_routines['selected_routine_trim'] = selected_routine_trim
    selected_routines['selected_routine_map'] = selected_routine_map
    selected_routines['selected_routine_qc'] = config['selected_routine_qc']
    selected_routines['selected_routine_merge'] = config['selected_routine_merge']
    selected_routines['selected_routine_touchup'] = config['selected_routine_touchup']   
    selected_routines['selected_routine_index'] = config['selected_routine_index']
    selected_routines['selected_routine_stats'] = config['selected_routine_stats']
    selected_routines['selected_routine_tagdir'] = config['selected_routine_tagdir']
    selected_routines['selected_routine_wig'] = config['selected_routine_wig']
    selected_routines['selected_routine_mergewig'] = config['selected_routine_mergewig']
    selected_routines['selected_routine_callpeaks'] = config['selected_routine_callpeaks']
    selected_routines['selected_routine_countreads'] = config['selected_routine_countreads']
    selected_routines['selected_routine_de'] = config['selected_routine_de']

    # Validate user-defined variables
    homer_input = validate_user_defined_vars(OMNOM_HOME, metadata, experiment_dir, INPUT, style, color_data_folder, col_table, overlay, the_type, map_tool, homer_size, homer_mindist, config)

    # Setup variables
    run_date = setup_variables(experiment_dir, config)

    #getting job mode (which steps)
    mode_steps, config = set_job_mode(args, config, experiment_dir, mode)
    print(f"MODE STEPS = {mode_steps}")

    # Check if pipeline needed or user first has to run separate scripts
    if min(mode_steps) == 12 and the_type == "CHIP":
        print("For ChIP experiments, first determine optimal peak caller settings, then manually run run_quant_peaks.sh and then continue with the next step!")
        return
    elif min(mode_steps) == 13 and the_type == "CHIP":
        print("To call DE peaks for ChIP  data, please first manually determine the best peak calling settings for your experiment and use run_quant_peaks.sh.")
        print("Then execute 'run_call_DE_peaks.sh' on your optimal peak set.")
        return
    elif min(mode_steps) == 11 and the_type == "RNA":
        print("Not a ChIP- or ATAC-seq experiment, skipping step 11 Call Peaks step...")
        mode_steps.pop(0)
        if not mode_steps:
            print("For the rest no steps to run. Done!")
            return
    print(f"MODE STEPS NEW = {mode_steps}")

    #Delete to be updated outputs
    delete_outputs_to_be_updated(mode_steps, config, experiment_dir)

    #checking input files
    num_files, num_pairs, paired, input_folder_mod_range_min, input_file_type_mod_range_min = validate_input_files(the_type, config, min(mode_steps),experiment_dir)

    #checking field settings
    check_name_field_settings(experiment_dir, separator, name_fields, type_field, col_field, config, input_folder_mod_range_min, input_file_type_mod_range_min)

    #checking sample names
    check_unique_sample_names(experiment_dir, input_folder_mod_range_min, input_file_type_mod_range_min, name_fields, separator)

    #setting up runtime parameters (my_cores and max_time commented out for now.)
    the_mem, the_heap_init = setup_runtime_parameters(num_pairs, experiment_dir)

    ## Run configuration for omnomnomics run
    run_config_data = {
        'OMNOM_HOME': OMNOM_HOME,
        'WORKFLOW_CONFIG_FILE': str(workflow_config_file),
        'SITE_CONFIG_FILE': str(site_config_file),
        'EXPERIMENT_DIR': experiment_dir,
        'RUNDATE': run_date,
        'INPUT_FOLDER': input_folder_mod_range_min,
        'INPUT_FILE_TYPE': input_file_type_mod_range_min,
        'INPUT': INPUT,
        'HOMERINPUT': homer_input,
        'BROAD': broad,
        'THESTYLE': style,
        'THEMODE': mode_steps,
        'THETYPE': the_type,
        'NUMFILES': num_files,
        'NUMPAIRS': num_pairs,
        'KEEPUNPAIRED': keep_unpaired,
        'PAIRED': paired,
        'THEMODERANGEMIN': min(mode_steps),
        'THEMODERANGEMAX': max(mode_steps),
        'THEGENOME': genome,
        'MYFORMULA': formula,
        'MYMETADATA': metadata,
        'THETRIMTOOL': trim_tool,
        'THEMAPTOOL': map_tool,
        'NO_MULTIQC': no_multiqc,
        'THETYPEFIELD': type_field,
        'NAMEFIELDS': name_fields,
        'THECOLFIELD': col_field,
        'THESEPARATOR': separator,
        'THEAPPENDIX': appendix,
        'THEOVERLAY': overlay,
        'THECOLORDATAFOLDER': color_data_folder,
        'THEHUBMAIL': hub_mail,
        'THECOLTABLE': col_table,
        'HOMERSIZE': homer_size,
        'HOMERMINDIST': homer_mindist,
        'THEMEM': the_mem,
        'THEHEAPINIT': the_heap_init  
    }
    # Initialize the run config
    write_run_config(experiment_dir, run_date, run_config_data)

    # Initialize the log file
    log_file = start_log(experiment_dir, run_date, config)


    #Print some status info
    print("")
    print("####################################################################################")
    print("Submitting omnomnomics pipeline job...")
    print(f"\tMODE:\t\t{mode_steps}")
    print(f"\tTYPE:\t\t{the_type}")
    print(f"\tFILES:\t\t{num_files}")
    print(f"\tPAIRS:\t\t{num_pairs}")
    print(f"\tExperiment DIR:\t\t{experiment_dir}")
    print(f"\tRun date:\t{run_date}")
    print("")

    #Dispatch the job using Snakemake
    print("Dispatching job...")

    sub_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    print("####################################################################################")
    print(f"Submitted omnomnomics run at {sub_time}")
    print("####################################################################################")

    with open(log_file, 'a') as log:
        log.write(f"Invocation:\t{the_command}\n")
        log.write(f"Snakemake Job submitted at:\t{sub_time}\n")
        log.write(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
        log.write("WAITING FOR SNAKEMAKE & SLURM TO COMPLETE...\n")
        log.write("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")

    # Call to snakemake!!
    cmd = [ "snakemake",  "--profile", os.path.join(OMNOM_HOME, "slurm_profile"), "--snakefile", f"{OMNOM_HOME}/Snakefile.smk",
        "--config", "config_file="+os.path.join(experiment_dir, "run_configs", f'omnomnomics.run.{run_date}.config.yaml'), "--jobs", "1000", 
        "--cores", "1280", "--rerun-triggers", "mtime", "--keep-going"
    ] 

    # For all the specified steps to run, add them to --forcerun so that they are always specified regardless 
    # from if the output is already present or not
    cmd.append("--forcerun")
    for i in mode_steps:
        routine = config['routines'][i-1]
        cmd.append(config[routine][selected_routines[f'selected_routine_{routine}']])

    # Executre the Snakemake command
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print("All done!")

if __name__ == "__main__":
   main()
