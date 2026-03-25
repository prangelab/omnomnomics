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
import shlex
import re
import glob
import random
import math
import shutil
import time
import select
import termios
import tty
import csv
from pathlib import Path
from datetime import date, datetime

from omnomnomics.genomes import genomes_main
from omnomnomics.helpers import create_track_color_table_main, display_track_color_table_main

PACKAGE_ROOT = Path(__file__).resolve().parent
WORKFLOW_ROOT = PACKAGE_ROOT / "workflow"
DEFAULT_WORKFLOW_CONFIG = WORKFLOW_ROOT / "config" / "workflow.yaml"
DEFAULT_SITE_CONFIG = WORKFLOW_ROOT / "config" / "site.yaml"

ANSI_RESET = "\033[0m"
ANSI_BOLD = "\033[1m"
ANSI_COLORS = [
    "\033[38;5;39m",
    "\033[38;5;208m",
    "\033[38;5;112m",
    "\033[38;5;141m",
    "\033[38;5;214m",
    "\033[38;5;45m",
    "\033[38;5;220m",
    "\033[38;5;177m",
    "\033[38;5;81m",
    "\033[38;5;203m",
    "\033[38;5;149m",
    "\033[38;5;69m",
    "\033[38;5;186m",
]
ANSI_STATUS_FIRST = "\033[38;5;81m"
ANSI_STATUS_OK = "\033[38;5;112m"
ANSI_STATUS_FAIL = "\033[38;5;203m"
ANSI_HEADER = "\033[38;5;250m"
ANSI_WARNING = "\033[38;5;214m"
FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
FASTQ_READ_SUFFIX_RE = re.compile(r'_(?:R)?[12](?:_[0-9]{3})?$')


def parse_monitor_arguments(argv):
   parser = argparse.ArgumentParser(
       description="Monitor the latest omnomnomics run log for an experiment directory.",
       allow_abbrev=False,
   )
   parser.add_argument('-i', '--experiment-dir', default='.', help='Path to the experiment directory. Default: current directory')
   parser.add_argument('--refresh', type=float, default=5.0, help='Polling interval in seconds. Default: 5')
   parser.add_argument('--latest-lines', type=int, default=20, help='Number of existing lines to show on startup. Default: 20')
   args, unknown = parser.parse_known_args(argv)
   if unknown:
       parser.error(f"Unrecognized arguments: {' '.join(unknown)}")
   return args


def default_user_site_config():
   xdg_config_home = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")).expanduser()
   return xdg_config_home / "omnomnomics" / "site.yaml"


def resolve_site_config_path(site_config_arg):
   if site_config_arg:
       return Path(site_config_arg).expanduser().resolve()
   user_site_config = default_user_site_config()
   if user_site_config.is_file():
       return user_site_config.resolve()
   return DEFAULT_SITE_CONFIG


def find_latest_run_log(experiment_dir):
   run_logs_dir = Path(experiment_dir) / "run_logs"
   if not run_logs_dir.is_dir():
       print(f"Run log directory '{run_logs_dir}' does not exist. Aborting...", file=sys.stderr)
       sys.exit(1)
   candidates = [
       path for path in run_logs_dir.glob("omnomnomics.run.*.log")
       if path.is_file() and ".backup" not in path.name and not path.name.endswith(".tools.log")
   ]
   if not candidates:
       print(f"No omnomnomics run logs found in '{run_logs_dir}'. Aborting...", file=sys.stderr)
       sys.exit(1)
   return max(candidates, key=lambda path: path.stat().st_mtime)


def colorize_monitor_line(line):
   stripped = line.rstrip("\n")
   if not stripped:
       return stripped

   if "last sample finished" in stripped:
       failed_match = re.search(r"; \d+/\d+ OK, (\d+)/\d+ failed", stripped)
       if failed_match and int(failed_match.group(1)) > 0:
           return f"{ANSI_STATUS_FAIL}{stripped}{ANSI_RESET}"
       return f"{ANSI_STATUS_OK}{stripped}{ANSI_RESET}"

   if "first sample entered" in stripped:
       return f"{ANSI_STATUS_FIRST}{stripped}{ANSI_RESET}"

   if "ERROR" in stripped or "Aborting" in stripped or "Pipeline failed" in stripped:
       return f"{ANSI_STATUS_FAIL}{stripped}{ANSI_RESET}"

   if "Warning" in stripped or "WARNING" in stripped:
       return f"{ANSI_WARNING}{stripped}{ANSI_RESET}"

   step_match = re.search(r"(\[STEP (\d+)\])", stripped)
   if step_match:
       full_marker = step_match.group(1)
       step_num = int(step_match.group(2))
       color = ANSI_COLORS[(step_num - 1) % len(ANSI_COLORS)]
       return stripped.replace(full_marker, f"{ANSI_BOLD}{color}{full_marker}{ANSI_RESET}", 1)

   if stripped.isupper() or stripped.startswith("##"):
       return f"{ANSI_HEADER}{ANSI_BOLD}{stripped}{ANSI_RESET}"

   return stripped


def load_monitor_tail(log_path, latest_lines):
   with log_path.open('r') as handle:
       lines = handle.readlines()
   return [colorize_monitor_line(line) for line in lines[-latest_lines:]]


def strip_fastq_read_suffix(sample_name):
   return FASTQ_READ_SUFFIX_RE.sub('', sample_name)


def normalize_monitor_sample_name(file_path, file_type):
   sample_name = os.path.basename(file_path).replace(file_type, "")
   if "fastq" in file_type or file_type.endswith(".fq.gz") or file_type.endswith(".fq"):
       sample_name = strip_fastq_read_suffix(sample_name)
   sample_name = sample_name.replace(".filtered", "")
   sample_name = sample_name.replace(".sorted.dups_marked", "")
   return sample_name


def collect_monitor_sample_totals(experiment_dir):
   lane_sample_files = []
   lane_sources = [
       ("FASTQ", FASTQ_EXTENSIONS),
       ("trimmed_FASTQ", (".trimmed.fastq.gz", ".trimmed.fastq", ".trimmed.fq.gz", ".trimmed.fq")),
       ("BAM", (".bam",)),
   ]
   for folder_name, extensions in lane_sources:
       folder_path = Path(experiment_dir) / folder_name
       if not folder_path.is_dir():
           continue
       candidate_files = []
       for extension in extensions:
           candidate_files.extend(folder_path.glob(f"*{extension}"))
       if folder_name == "BAM":
           candidate_files = [
               path for path in candidate_files
               if re.search(r'_L0\d+\.bam$', path.name)
           ]
       if candidate_files:
           lane_sample_files = candidate_files
           break

   lane_samples = []
   for path in lane_sample_files:
       extension = next((ext for ext in FASTQ_EXTENSIONS if str(path).endswith(ext)), None)
       if extension is None:
           for ext in (".trimmed.fastq.gz", ".trimmed.fastq", ".trimmed.fq.gz", ".trimmed.fq", ".bam"):
               if str(path).endswith(ext):
                   extension = ext
                   break
       if extension is None:
           continue
       lane_samples.append(normalize_monitor_sample_name(str(path), extension))
   lane_samples = sorted(set(lane_samples))

   if lane_samples:
       merged_samples = sorted(set(re.sub(r'_L00.', '', sample_name) for sample_name in lane_samples))
       return len(lane_samples), len(merged_samples)

   merged_sample_files = []
   merged_sources = [
       ("filtered_BAM", (".sorted.dups_marked.filtered.bam", ".filtered.bam")),
       ("BAM", (".bam",)),
       ("BigWigs", (".plus.bw", ".minus.bw", ".bw")),
   ]
   for folder_name, extensions in merged_sources:
       folder_path = Path(experiment_dir) / folder_name
       if not folder_path.is_dir():
           continue
       candidate_files = []
       for extension in extensions:
           candidate_files.extend(folder_path.glob(f"*{extension}"))
       if folder_name == "BAM":
           candidate_files = [
               path for path in candidate_files
               if not re.search(r'_L0\d+\.bam$', path.name)
           ]
       if candidate_files:
           merged_sample_files = candidate_files
           break

   merged_samples = set()
   for path in merged_sample_files:
       name = path.name
       for extension in (".sorted.dups_marked.filtered.bam", ".filtered.bam", ".plus.bw", ".minus.bw", ".bw", ".bam"):
           if name.endswith(extension):
               merged_samples.add(normalize_monitor_sample_name(str(path), extension))
               break
   return 0, len(merged_samples)


def load_step_monitor_rows(experiment_dir):
   step_log_dir = Path(experiment_dir) / "run_logs" / "steps"
   if not step_log_dir.is_dir():
       return []
   lane_total, merged_total = collect_monitor_sample_totals(experiment_dir)

   step_nums = set()
   for summary_path in step_log_dir.glob("step*.summary.tsv"):
       match = re.search(r"step(\d+)\.summary\.tsv$", summary_path.name)
       if match:
           step_nums.add(int(match.group(1)))
   for state_dir in step_log_dir.glob(".step*.state"):
       match = re.search(r"\.step(\d+)\.state$", state_dir.name)
       if match:
           step_nums.add(int(match.group(1)))

   rows = []
   for step_num in sorted(step_nums):
       state_dir = step_log_dir / f".step{step_num:02d}.state"
       started_dir = state_dir / "started"
       completed_dir = state_dir / "completed"
       failed_dir = state_dir / "failed"
       finished_marker = state_dir / "finished.marker"
       if state_dir.is_dir():
           started = len(list(started_dir.iterdir())) if started_dir.is_dir() else 0
           completed = len(list(completed_dir.iterdir())) if completed_dir.is_dir() else 0
           failed = len(list(failed_dir.iterdir())) if failed_dir.is_dir() else 0
           running = max(0, started - completed - failed)
           finished = finished_marker.exists()
       else:
           running = 0
           completed = 0
           failed = 0
           finished = False
           summary_path = step_log_dir / f"step{step_num:02d}.summary.tsv"
           if summary_path.exists():
               with summary_path.open(newline="") as handle:
                   reader = csv.DictReader(handle, delimiter="\t")
                   for row in reader:
                       status = row.get("status", "")
                       if status == "RUNNING":
                           running += 1
                       elif status == "OK":
                           completed += 1
                       elif status == "FAIL":
                           failed += 1
       if failed:
           state = "FAILED" if running == 0 else "RUNNING/FAIL"
           state_color = ANSI_STATUS_FAIL
       elif running:
           state = "RUNNING"
           state_color = ANSI_STATUS_FIRST
       elif finished:
           state = "DONE"
           state_color = ANSI_STATUS_OK
       elif completed:
           state = "WAITING"
           state_color = ANSI_HEADER
       else:
           state = "PENDING"
           state_color = ANSI_HEADER
       if step_num in {1, 2, 3}:
           total = lane_total
       elif step_num in {4, 5, 6, 7, 8, 13}:
           total = merged_total
       else:
           total = running + completed + failed
       pending = max(0, total - running - completed - failed)
       rows.append({
           "step_num": step_num,
           "state": state,
           "state_color": state_color,
           "running": running,
           "completed": completed,
           "failed": failed,
           "pending": pending,
           "total": total,
       })
   return rows


def render_monitor_screen(log_path, step_rows, tail_lines):
   separator = "-" * 72
   step_width = 6
   state_width = 14
   count_width = 8
   print("\033[2J\033[H", end="")
   print(f"{ANSI_HEADER}{ANSI_BOLD}omnomnomics monitor{ANSI_RESET}")
   print(f"Log: {log_path}")
   print("Press any key to exit.\n")

   if step_rows:
       print(separator)
       print(f"{ANSI_BOLD}Step summary{ANSI_RESET}")
       print(
           f"{'STEP':<{step_width}} "
           f"{'STATE':<{state_width}} "
           f"{'RUNNING':>{count_width}} "
           f"{'PEND':>{count_width}} "
           f"{'DONE':>{count_width}} "
           f"{'FAIL':>{count_width}}"
       )
       for row in step_rows:
           step_label = f"{row['step_num']}"
           padded_state = f"{row['state']:<{state_width}}"
           colored_state = f"{row['state_color']}{padded_state}{ANSI_RESET}"
           print(
               f"{step_label:<{step_width}} "
               f"{colored_state} "
               f"{row['running']:>{count_width}} "
               f"{row['pending']:>{count_width}} "
               f"{row['completed']:>{count_width}} "
               f"{row['failed']:>{count_width}}"
           )
       print(separator)
       print("")

   print(separator)
   print(f"{ANSI_BOLD}Recent log lines{ANSI_RESET}")
   print(separator)
   for line in tail_lines:
       print(line)


class KeypressListener:
   def __init__(self):
       self.enabled = sys.stdin.isatty()
       self.fd = None
       self.original = None

   def __enter__(self):
       if self.enabled:
           self.fd = sys.stdin.fileno()
           self.original = termios.tcgetattr(self.fd)
           tty.setcbreak(self.fd)
       return self

   def __exit__(self, exc_type, exc, tb):
       if self.enabled and self.fd is not None and self.original is not None:
           termios.tcsetattr(self.fd, termios.TCSADRAIN, self.original)

   def wait_for_keypress(self, timeout_seconds):
       if not self.enabled:
           time.sleep(timeout_seconds)
           return False
       end_time = time.time() + timeout_seconds
       while time.time() < end_time:
           remaining = max(0, min(0.2, end_time - time.time()))
           readable, _, _ = select.select([sys.stdin], [], [], remaining)
           if readable:
               os.read(self.fd, 1)
               return True
       return False


def monitor_main(argv):
   args = parse_monitor_arguments(argv)
   experiment_dir = Path(args.experiment_dir).expanduser().resolve()
   latest_log = find_latest_run_log(experiment_dir)

   with KeypressListener() as listener:
       while True:
           current_log = find_latest_run_log(experiment_dir)
           if current_log != latest_log:
               latest_log = current_log
           step_rows = load_step_monitor_rows(experiment_dir)
           tail_lines = load_monitor_tail(latest_log, args.latest_lines)
           render_monitor_screen(latest_log, step_rows, tail_lines)

           if listener.wait_for_keypress(args.refresh):
               print("Monitor stopped.")
               return

def parse_arguments():
   #Parse command-line arguments
   parser = argparse.ArgumentParser(
       description="Modular HPC pipeline for RNA-seq, ATAC-seq, and ChIP-seq processing.",
       allow_abbrev=False,
   )

   # Define command-line options
   parser.add_argument('-i', '--experiment-dir', help='Path to the experiment directory')
   parser.add_argument('-t', '--type', help='Type of experiment: RNA, ChIP, ATAC')
   parser.add_argument('-g', '--genome', help='Genome assembly name. Must match an available directory under the configured genome assembly root')
   parser.add_argument('-j', '--mode', help='Job mode. Can be auto, all or a range of jobs. See readme for some examples. \n \t Default: auto')
   parser.add_argument('-T', '--trim-tool', help='Trimming tool choice. Can be fastp or skewer. \n \t Default: fastp')
   parser.add_argument('-M', '--map-tool', help='Mapping tool choice. Can be HISAT2, STAR, or STAR_TE. STAR(_TE) can only be used for RNA-seq data. \n \t Default: HISAT2')
   parser.add_argument('-f', '--formula', help='RNA: Experimental Design for DE calling. \n \t Default: 1 (just an intercept)')
   parser.add_argument('-I', '--input', help='Input BAM file used for ChIP peak calling with MACS3. \n \t Default: do not use input')
   parser.add_argument('-m', '--metadata', help='.txt file with columns of metadata for RNA-seq experiments. \n \t Default: DESeq2 style metadata table describing all samples. Rownames should be samplenames')
   parser.add_argument('-b', '--broad', action='store_true', help='ChIP: Call broad histone marks with MACS3 --broad mode. Default is TF / narrow peaks.')
   parser.add_argument('-C', '--col-table', help='File specifying which colors to use for the tracks. \n \t Default: gray.tint.color.table from the packaged palette directory. Can be a *txt list file with one color table per line. Different color tables will be used per hub as split by -e. Can be a full path or a file basename in combination with -P. Use `omnomnomics create-track-color-table` to build custom palettes and `omnomnomics display-track-color-table` to preview them.')
   parser.add_argument('-P', '--color-data-folder', help='Path to a folder with color tables. \n \t Default: packaged color_data_for_hubs directory. Use this to point the workflow at custom palettes.') ##################change this default??
   parser.add_argument('-o', '--overlay', help='Overlay type (transparentOverlay|stacked|solidOverlay|none) \n \t Default: transparentOverlay')
   parser.add_argument('-L', '--hub-mail', help='Email to use in trackhub \n \t Default: your@email.com')
   parser.add_argument('-X', '--no-multiqc', action='store_true', help='Exclude multiQC stats aggregator. Set if you don not wish to run multiQC.')
   parser.add_argument('--create-homer-tagdirs', action='store_true', help='Create optional HOMER tag directory exports in addition to the main pipeline outputs.')
   parser.add_argument('--rerun-selected-steps', action='store_true', help='Force recomputation of the selected workflow steps by deleting their existing outputs before submission. Default behavior is to reuse existing outputs when Snakemake considers them up to date.')
   parser.add_argument('--remove-duplicates', action='store_true', help='Remove duplicate reads in step 5. Default is assay-aware: keep for RNA, remove for ATAC and ChIP.')
   parser.add_argument('--keep-duplicates', action='store_true', help='Keep duplicate reads in step 5. Default is assay-aware: keep for RNA, remove for ATAC and ChIP.')
   parser.add_argument('-n', '--name-fields', help='Field(s) in filename to use as track name, peak file name, and column header in the count table \n \t Default: 1-3')
   parser.add_argument('-e', '--type-field', help='Field(s) in filename to use as merged hub or peak calling group identifier \n \t Default: 1 (Creates separate merged hubs for each unique entry)')
   parser.add_argument('-c', '--col-field', help='Field(s) in filename to use as color type \n \t Default: 2')
   parser.add_argument('-s', '--separator', help='Separator used in file names. \n \t Default: _')
   parser.add_argument('-a', '--appendix', help='Appendix to add to track name \n \t Default: hub')
   parser.add_argument('-k', '--keepunpaired', action='store_true', help='Keep unpaired or not in HISAT2')
   parser.add_argument('--dry-run', action='store_true', help='Validate the workflow and build the Snakemake DAG without executing jobs')
   parser.add_argument('--site-config', help='Optional path to a site-specific config YAML. Default: $XDG_CONFIG_HOME/omnomnomics/site.yaml or ~/.config/omnomnomics/site.yaml, then packaged site config')

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
       header_lines = [line.strip() for line in config_file if line.strip()]
       if len(header_lines) < 1 or expected_header not in header_lines[0]:
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

def resolve_config_path(path_value, workflow_root):
   resolved = str(path_value)
   resolved = resolved.replace("{WORKFLOW_ROOT}", str(workflow_root))
   resolved = resolved.replace("{HOME}", str(Path.home()))
   return str(Path(resolved).expanduser().resolve())

def genome_subdir(assembly_root, genome_name, subdir_name):
   return str(Path(assembly_root) / genome_name / subdir_name)

def write_controller_script(script_path, experiment_dir, snakemake_cmd):
   script_lines = [
       "#!/bin/bash",
       "set -euo pipefail",
       f"cd {shlex.quote(str(experiment_dir))}",
       shlex.join(snakemake_cmd),
   ]
   script_path.write_text("\n".join(script_lines) + "\n")
   script_path.chmod(0o755)

def list_available_genomes(assembly_root):
   assembly_root_path = Path(assembly_root)
   if not assembly_root_path.is_dir():
       print(f"Genome assembly root '{assembly_root_path}' does not exist. Aborting...", file=sys.stderr)
       sys.exit(1)
   return sorted(path.name for path in assembly_root_path.iterdir() if path.is_dir())

##---------------------------------------------------------------------------------------------------------------
## Set required vars or die
##---------------------------------------------------------------------------------------------------------------
def check_required_vars(the_type, experiment_dir, genome, available_genomes, config):
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
       if genome not in available_genomes:
           print(f"Genome (-g) should be one of: {', '.join(available_genomes)}. Aborting...", file=sys.stderr)
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
def validate_user_defined_vars(workflow_root, metadata, experiment_dir, INPUT, color_data_folder, col_table, overlay, the_type, map_tool, config):
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
            print(f"e.g., -m {os.path.join(str(workflow_root), 'data', 'me', 'my_fantastic_experiment', 'my_mindboggling_metadata.txt')}. Aborting...", file=sys.stderr)
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
            print(f"e.g., -I {os.path.join(str(Path.home()), 'genomes', 'input_ChIP', 'my_awesome_input.bam')}. Aborting...", file=sys.stderr)
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

    return

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
    if mode_range_min == 11:
        input_file_type_mod_range_min = input_file_type_mod_range_min[0]
        input_folder_mod_range_min = input_folder_mod_range_min[0]


    if os.path.isdir(f"{experiment_dir}/{input_folder_mod_range_min}"):
        num_files = len(glob.glob(f"{experiment_dir}/{input_folder_mod_range_min}/*{input_file_type_mod_range_min}"))


    # Sanity check file number
    if num_files == 0:
        if mode_range_min == 12 and the_type == "CHIP":
                pass
        else: 
            print("No input files detected! Aborting...", file=sys.stderr)
            sys.exit(1)
   
    if mode_range_min == 12 and the_type == "CHIP":
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

def normalize_field_selection_name(filename, file_type):
    normalized = os.path.basename(filename)
    if normalized.endswith(file_type):
        normalized = normalized[:-len(file_type)]
    if any(normalized.endswith(suffix) for suffix in (".plus", ".minus")):
        normalized = re.sub(r'\.(plus|minus)$', '', normalized)
    if "fastq" in file_type or file_type.endswith(".fq.gz") or file_type.endswith(".fq"):
        normalized = strip_fastq_read_suffix(normalized)
    normalized = normalized.replace(".filtered", "")
    normalized = normalized.replace(".sorted.dups_marked", "")
    return normalized

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
        
        mock_cut = normalize_field_selection_name(sample_files[0], input_file_type_mod_range_min)
        
        # Run the cut command simulations
        cut_test = ""
        cut_test += run_cut_command(mock_cut, name_fields, separator)
        cut_test += run_cut_command(mock_cut, type_field, separator)
        cut_test += run_cut_command(mock_cut, col_field, separator)

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
    sample_names = [extract_fields(normalize_field_selection_name(f, input_file_type_mod_range_min), name_fields, separator) for f in files]

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
    subprocess.run(f"mkdir -p {experiment_dir}/slurm_logs/controller", shell=True, check=True)
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
    tools_log_file = os.path.join(experiment_dir, "run_logs", f"omnomnomics.run.{run_date}.tools.log")
    marker_dir = os.path.join(experiment_dir, "run_logs", f"omnomnomics.run.{run_date}.markers")
    tools_marker_dir = os.path.join(experiment_dir, "run_logs", f"omnomnomics.run.{run_date}.tools.markers")
    if os.path.isfile(log_file):
        backup_log = f"{log_file}.{os.urandom(8).hex()}.backup" #Create a random backup if multiple runs
        os.rename(log_file, backup_log)
        print(f"Existing log file backed up as: {backup_log}")
    if os.path.isfile(tools_log_file):
        backup_tools_log = f"{tools_log_file}.{os.urandom(8).hex()}.backup"
        os.rename(tools_log_file, backup_tools_log)
        print(f"Existing tools log backed up as: {backup_tools_log}")
    if os.path.isdir(marker_dir):
        shutil.rmtree(marker_dir)
    if os.path.isdir(tools_marker_dir):
        shutil.rmtree(tools_marker_dir)


    with open(log_file, 'w') as log:
        log.write("#################################\n")
        log.write("## Run log for omnomnomics run ##\n")
        log.write("#################################\n\n")
        log.write(f"Pipeline version: {config['omnomnomics']}\n\n")

    with open(tools_log_file, 'w') as log:
        log.write("###################################\n")
        log.write("## Tool log for omnomnomics run ##\n")
        log.write("###################################\n\n")
        log.write(f"Pipeline version: {config['omnomnomics']}\n\n")


    return log_file 


##--------------------------------------------------------------------------------------------------------------
# Reset selected step bookkeeping for a new run
##--------------------------------------------------------------------------------------------------------------
def reset_step_tracking(mode_steps, experiment_dir):
    step_log_dir = os.path.join(experiment_dir, "run_logs", "steps")
    for num in mode_steps:
        step_prefix = f"step{num:02d}"
        step_state_dir = os.path.join(step_log_dir, f".{step_prefix}.state")
        if os.path.isdir(step_state_dir):
            shutil.rmtree(step_state_dir)
        for suffix in ("summary.tsv", "commands.txt", "notes.txt"):
            step_log_file = os.path.join(step_log_dir, f"{step_prefix}.{suffix}")
            if os.path.exists(step_log_file):
                os.remove(step_log_file)

##--------------------------------------------------------------------------------------------------------------
# Remove already present outputs of rules that you want to run
##--------------------------------------------------------------------------------------------------------------
def delete_outputs_to_be_updated(mode_steps, config, experiment_dir):
    print("FORCING RECOMPUTATION OF SELECTED STEP OUTPUTS")
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
                    elif num == 9:
                        shutil.rmtree(file) # Is actually a hub directory and not a file
                    else:
                        os.remove(file)
        if num == 1:
            trim_metric_files = glob.glob(f"{experiment_dir}/{outputfolder}/*.trim_metrics.tsv")
            for file in trim_metric_files:
                if os.path.exists(file):
                    os.remove(file)
    if config.get('create_homer_tagdirs', False):
        homer_outputfolder = config['output_folders'][config['homer_tagdir_rule_num'] - 1]
        homer_output_filetype = config['output_file_types'][config['homer_tagdir_rule_num'] - 1]
        files = glob.glob(f"{experiment_dir}/{homer_outputfolder}/*{homer_output_filetype}")
        for file in files:
            if os.path.exists(file):
                os.remove(file)
##--------------------------------------------------------------------------------------------------------------
# Main function
##--------------------------------------------------------------------------------------------------------------
def main():
    if len(sys.argv) > 1 and sys.argv[1] == "genomes":
        genomes_main(sys.argv[2:], WORKFLOW_ROOT, DEFAULT_WORKFLOW_CONFIG, DEFAULT_SITE_CONFIG)
        return
    if len(sys.argv) > 1 and sys.argv[1] == "monitor":
        monitor_main(sys.argv[2:])
        return
    if len(sys.argv) > 1 and sys.argv[1] == "create-track-color-table":
        original_argv = sys.argv[:]
        try:
            sys.argv = [original_argv[0], *original_argv[2:]]
            create_track_color_table_main()
        finally:
            sys.argv = original_argv
        return
    if len(sys.argv) > 1 and sys.argv[1] == "display-track-color-table":
        original_argv = sys.argv[:]
        try:
            sys.argv = [original_argv[0], *original_argv[2:]]
            display_track_color_table_main()
        finally:
            sys.argv = original_argv
        return

    # Define variables for packaged workflow paths
    workflow_root = WORKFLOW_ROOT
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

    site_config_file = resolve_site_config_path(args.site_config)
    if not site_config_file.is_file():
        print(f"Site config file '{site_config_file}' does not exist. Please make sure it exists. Aborting...")
        sys.exit(1)

    workflow_config = load_and_validate_yaml(workflow_config_file, "## Omnomnomics pipeline config ##")
    site_config = load_and_validate_yaml(site_config_file, "## Omnomnomics pipeline config ##")
    config = merge_configs(workflow_config, site_config)
    config['WORKFLOW_ROOT'] = str(workflow_root)
    config['genome_assembly_dir'] = resolve_config_path(config['genome_assembly_dir'], workflow_root)
    config['cellranger_reference_dir'] = resolve_config_path(config['cellranger_reference_dir'], workflow_root)
    available_genomes = list_available_genomes(config['genome_assembly_dir'])

    # Access parsed arguments
    experiment_dir = str(Path(args.experiment_dir).expanduser().resolve()) if args.experiment_dir else args.experiment_dir
    the_type = args.type.upper()
    genome = args.genome
    trim_tool = args.trim_tool.lower() if args.trim_tool else config.get("trim_tool","fastp").lower()
    map_tool = args.map_tool.lower() if args.map_tool else config.get('map_tool', "hisat2").lower()
    mode = args.mode.lower() if args.mode else config.get('mode', "auto")
    formula = args.formula if args.formula else config.get('formula', "1")
    broad = args.broad if args.broad else config.get('broad', "NA")
    INPUT = args.input if args.input else config.get('input',"NA")
    metadata = args.metadata if args.metadata else config.get('metadata', "NA")
    col_table = resolve_config_path(args.col_table, workflow_root) if args.col_table else resolve_config_path(config.get('color_table', f"{workflow_root}/bin/color_data_for_hubs/gray.tint.color.table"), workflow_root)
    color_data_folder = resolve_config_path(args.color_data_folder, workflow_root) if args.color_data_folder else resolve_config_path(config.get('color_data_folder', f"{workflow_root}/bin/color_data_for_hubs"), workflow_root)
    overlay = args.overlay if args.overlay else config.get('overlay', "transparentOverlay")
    hub_mail = args.hub_mail if args.hub_mail else config.get('hub_mail', "your@email.com")
    no_multiqc = args.no_multiqc if args.no_multiqc else config.get('no_multiqc', 0)
    create_homer_tagdirs = args.create_homer_tagdirs or config.get('create_homer_tagdirs', False)
    rerun_selected_steps = args.rerun_selected_steps
    if args.remove_duplicates and args.keep_duplicates:
        print("Please specify only one of --remove-duplicates or --keep-duplicates. Aborting...", file=sys.stderr)
        sys.exit(1)
    duplicate_handling = config.get('duplicate_handling', 'auto')
    if args.remove_duplicates:
        duplicate_handling = "remove"
    elif args.keep_duplicates:
        duplicate_handling = "keep"
    elif duplicate_handling == "auto":
        duplicate_handling = "keep" if the_type == "RNA" else "remove"
    name_fields = args.name_fields if args.name_fields else config.get('name_fields', "1-3")
    type_field = args.type_field if args.type_field else config.get('type_field', "1")
    col_field = args.col_field if args.col_field else config.get('col_field', "2")
    separator = args.separator if args.separator else config.get('separator', "_")
    appendix = args.appendix if args.appendix else config.get('appendix', "hub")
    keep_unpaired = args.keepunpaired if args.keepunpaired else config.get('keep_unpaired', False)
    dry_run = args.dry_run

    # Check required variables
    check_required_vars(the_type, experiment_dir, genome, available_genomes, config)
    config['create_homer_tagdirs'] = create_homer_tagdirs
    config['duplicate_handling'] = duplicate_handling

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
    selected_routines['selected_routine_wig'] = config['selected_routine_wig']
    selected_routines['selected_routine_mergewig'] = config['selected_routine_mergewig']
    selected_routines['selected_routine_callpeaks'] = config['selected_routine_callpeaks']
    selected_routines['selected_routine_countreads'] = config['selected_routine_countreads']
    selected_routines['selected_routine_de'] = config['selected_routine_de']

    # Validate user-defined variables
    validate_user_defined_vars(workflow_root, metadata, experiment_dir, INPUT, color_data_folder, col_table, overlay, the_type, map_tool, config)

    # Setup variables
    run_date = setup_variables(experiment_dir, config)

    #getting job mode (which steps)
    mode_steps, config = set_job_mode(args, config, experiment_dir, mode)
    print(f"MODE STEPS = {mode_steps}")

    # Check if pipeline needed or user first has to run separate scripts
    if min(mode_steps) == 11 and the_type == "CHIP":
        print("For ChIP experiments, first determine optimal peak caller settings and quantify peaks with your chosen downstream workflow before continuing.")
        return
    if 10 in mode_steps and the_type == "RNA":
        print("Not a ChIP- or ATAC-seq experiment, skipping step 10 Call Peaks step...")
        mode_steps = [step for step in mode_steps if step != 10]
        if not mode_steps:
            print("For the rest no steps to run. Done!")
            return
    if 12 in mode_steps and the_type != "RNA":
        print("Step 12 DE calling is RNA-only. Skipping it for this assay type...")
        mode_steps = [step for step in mode_steps if step != 12]
        if not mode_steps:
            print("For the rest no steps to run. Done!")
            return
    if 12 in mode_steps and the_type == "RNA":
        print("Step 12 DE calling is not implemented yet for RNA, skipping it for now...")
        mode_steps = [step for step in mode_steps if step != 12]
        if not mode_steps:
            print("For the rest no steps to run. Done!")
            return
    print(f"MODE STEPS NEW = {mode_steps}")

    # Reset selected step bookkeeping for the new run
    reset_step_tracking(mode_steps, experiment_dir)

    # Optionally force recomputation of selected step outputs
    if rerun_selected_steps:
        delete_outputs_to_be_updated(mode_steps, config, experiment_dir)
    else:
        print("REUSING EXISTING OUTPUTS WHEN POSSIBLE")

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
        'WORKFLOW_CONFIG_FILE': str(workflow_config_file),
        'SITE_CONFIG_FILE': str(site_config_file),
        'WORKFLOW_ROOT': str(workflow_root),
        'SCRIPT_DIR': str(workflow_root / "bin" / "scripts"),
        'COLOR_DATA_FOLDER_DEFAULT': str(workflow_root / "bin" / "color_data_for_hubs"),
        'GENOME_ASSEMBLY_DIR': config['genome_assembly_dir'],
        'CELLRANGER_REFERENCE_DIR': config['cellranger_reference_dir'],
        'HISAT2_GENOME_DIR': genome_subdir(config['genome_assembly_dir'], genome, "hisat2"),
        'STAR_GENOME_DIR': genome_subdir(config['genome_assembly_dir'], genome, "star"),
        'GENOME_AUX_DIR': genome_subdir(config['genome_assembly_dir'], genome, "aux"),
        'EXPERIMENT_DIR': experiment_dir,
        'RUNDATE': run_date,
        'INPUT_FOLDER': input_folder_mod_range_min,
        'INPUT_FILE_TYPE': input_file_type_mod_range_min,
        'INPUT': INPUT,
        'BROAD': broad,
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
        'CREATE_HOMER_TAGDIRS': create_homer_tagdirs,
        'RERUN_SELECTED_STEPS': rerun_selected_steps,
        'DUPLICATE_HANDLING': duplicate_handling,
        'THETYPEFIELD': type_field,
        'NAMEFIELDS': name_fields,
        'THECOLFIELD': col_field,
        'THESEPARATOR': separator,
        'THEAPPENDIX': appendix,
        'THEOVERLAY': overlay,
        'THECOLORDATAFOLDER': color_data_folder,
        'THEHUBMAIL': hub_mail,
        'THECOLTABLE': col_table,
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
    print("Dispatching controller job...")

    sub_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    print("####################################################################################")
    print(f"Submitted omnomnomics run at {sub_time}")
    print("####################################################################################")

    with open(log_file, 'a') as log:
        log.write(f"Invocation:\t{the_command}\n")
        log.write(f"Controller submission at:\t{sub_time}\n")

    # Call to snakemake!!
    max_cores = str(config.get("max_nodes", 1) * config.get("cores_per_node", 1))
    max_jobs = str(config.get("max_jobs", 100))
    max_jobs_per_second = str(config.get("max_jobs_per_second", 2))
    max_status_checks_per_second = str(config.get("max_status_checks_per_second", 2))
    cmd = [ "snakemake",  "--profile", os.path.join(str(workflow_root), "slurm_profile"), "--snakefile", f"{workflow_root}/Snakefile.smk",
        "--config", "config_file="+os.path.join(experiment_dir, "run_configs", f'omnomnomics.run.{run_date}.config.yaml'), "--jobs", max_jobs,
        "--cores", max_cores, "--max-jobs-per-second", max_jobs_per_second, "--max-status-checks-per-second", max_status_checks_per_second,
        "--default-resources", f"partition={config['partition']}", f"runtime={config['default_runtime']}", "--rerun-triggers", "mtime", "--keep-going"
    ]

    if dry_run:
        cmd.append("--dry-run")

    # Only force selected routines when the user explicitly asks to recompute them
    if rerun_selected_steps:
        cmd.append("--forcerun")
        for i in mode_steps:
            routine = config['routines'][i-1]
            cmd.append(config[routine][selected_routines[f'selected_routine_{routine}']])

    controller_script = Path(experiment_dir) / "run_configs" / f"omnomnomics.run.{run_date}.controller.sh"
    write_controller_script(controller_script, experiment_dir, cmd)

    controller_job_name = f"omnomnomics-{Path(experiment_dir).name}"
    controller_log = Path(experiment_dir) / "slurm_logs" / "controller" / "controller.%j.out"
    sbatch_cmd = [
        "sbatch",
        "--parsable",
        "--export=ALL",
        "--chdir",
        str(experiment_dir),
        "--job-name",
        controller_job_name,
        "--partition",
        str(config.get("controller_partition", config["partition"])),
        "--cpus-per-task",
        str(config.get("controller_cores", 1)),
        "--mem",
        str(config.get("controller_mem_mb", 8192)),
        "--time",
        str(config.get("controller_runtime", 1440)),
        "--output",
        str(controller_log),
        "--error",
        str(controller_log),
    ]
    if config.get("controller_constraint"):
        sbatch_cmd.extend(["--constraint", str(config["controller_constraint"])])
    sbatch_cmd.append(str(controller_script))

    try:
        submission = subprocess.run(
            sbatch_cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}", file=sys.stderr)
        if e.stderr:
            print(e.stderr, file=sys.stderr)
        sys.exit(1)

    controller_job_id = submission.stdout.strip().split(";")[0]
    controller_log_path = str(controller_log).replace("%j", controller_job_id)

    with open(log_file, 'a') as log:
        log.write(f"Controller job ID:\t{controller_job_id}\n")
        log.write(f"Controller script:\t{controller_script}\n")
        log.write(f"Controller log:\t{controller_log_path}\n")
        log.write(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
        log.write("CONTROLLER JOB SUBMITTED TO SLURM.\n")
        log.write("CHECK THE CONTROLLER LOG AND RUN LOG FOR PROGRESS.\n")
        log.write("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n")

    print(f"Controller job ID:\t{controller_job_id}")
    print(f"Controller log:\t{controller_log_path}")
    print("CLI finished after successful controller submission.")

if __name__ == "__main__":
   main()
