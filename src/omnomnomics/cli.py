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
import json
from pathlib import Path
from datetime import date, datetime

from omnomnomics.de_config import DEConfigError, resolve_de_config
from omnomnomics.genomes import genomes_main
from omnomnomics.helpers import create_track_color_table_main, display_track_color_table_main
from omnomnomics.metadata import (
    MetadataError,
    derive_metadata_rows,
    normalize_metadata_filename,
    read_metadata_table,
    resolve_de_metadata,
    write_metadata_table,
)

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
FASTQ_READ_SUFFIX_RE = re.compile(r'(?:_R[12]|_[12]_[0-9]{3})$')
MERGED_LANE_SUFFIX_RE = re.compile(r'_L00.')
METADATA_REQUIRED_STEPS = {9, 10, 11, 12, 13, 14, 15, 16}
ASSAY_PUBLIC_INTERNAL_STEP_MAP = {
   "RNA": {i: i for i in range(1, 13)},
   "ATAC": {
      1: 1,
      2: 2,
      3: 3,
      4: 4,
      5: 5,
      6: 6,
      7: 7,
      8: 8,
      9: 9,
      10: 10,
      11: 13,
      12: 14,
      13: 11,
      14: 15,
      15: 16,
   },
   "CHIP": {
      1: 1,
      2: 2,
      3: 3,
      4: 4,
      5: 5,
      6: 6,
      7: 7,
      8: 8,
      9: 9,
      10: 10,
      11: 13,
      12: 14,
      13: 11,
      14: 15,
      15: 16,
   },
}


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


def parse_de_app_arguments(argv):
   parser = argparse.ArgumentParser(
       description="Launch the omnomnomics DE Shiny app.",
       allow_abbrev=False,
   )
   parser.add_argument(
       "-i",
       "--project-dir",
       default=".",
       help="Project directory or DE_calling directory to pre-load in the app. Default: current directory",
   )
   parser.add_argument("--host", default="127.0.0.1", help="Host interface for Shiny. Default: 127.0.0.1")
   parser.add_argument("--port", type=int, default=3838, help="Port for Shiny. Default: 3838")
   parser.add_argument(
       "--no-browser",
       action="store_true",
       help="Do not auto-open a browser window.",
   )
   args, unknown = parser.parse_known_args(argv)
   if unknown:
       parser.error(f"Unrecognized arguments: {' '.join(unknown)}")
   return args


def print_top_level_help(exit_code=0):
   help_text = """
usage: omnomnomics <command> [options]

Assay workflow commands:
  omnomnomics rna   ...   Run RNA-seq workflow
  omnomnomics atac  ...   Run ATAC-seq workflow
  omnomnomics chip  ...   Run ChIP-seq workflow

Utility commands:
  omnomnomics monitor                   Monitor latest run log
  omnomnomics de-app                    Launch DE Shiny app
  omnomnomics genomes ...               Genome helper subcommands
  omnomnomics create-track-color-table  Build custom track color table
  omnomnomics display-track-color-table Preview an existing track color table

Examples:
  omnomnomics rna -i /path/to/experiment -g GRCh38
  omnomnomics monitor -i /path/to/experiment
  omnomnomics de-app --project-dir /path/to/experiment/DE_calling

Use assay-specific help:
  omnomnomics rna --help
  omnomnomics atac --help
  omnomnomics chip --help
"""
   print(help_text.strip())
   sys.exit(exit_code)


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


def de_app_main(argv):
   args = parse_de_app_arguments(argv)
   app_dir = WORKFLOW_ROOT / "R" / "shiny_app"
   app_entry = app_dir / "app.R"
   if not app_entry.is_file():
       print(f"Shiny app entry file '{app_entry}' does not exist. Aborting...", file=sys.stderr)
       sys.exit(1)

   project_dir = str(Path(args.project_dir).expanduser().resolve())
   env = os.environ.copy()
   env["OMNOMNOMICS_DE_APP_PROJECT"] = project_dir

   launch_browser = "FALSE" if args.no_browser else "TRUE"
   app_dir_r = json.dumps(str(app_dir))
   host_r = json.dumps(str(args.host))
   r_expr = (
       f"shiny::runApp({app_dir_r}, "
       f"host={host_r}, "
       f"port={int(args.port)}, "
       f"launch.browser={launch_browser})"
   )
   cmd = ["Rscript", "-e", r_expr]
   try:
       completed = subprocess.run(cmd, env=env)
   except FileNotFoundError:
       print("Rscript was not found in PATH. Activate an environment with R installed and retry.", file=sys.stderr)
       sys.exit(1)
   sys.exit(completed.returncode)

def parse_arguments(argv=None):
   #Parse command-line arguments
   parser = argparse.ArgumentParser(
       description="Modular HPC pipeline for RNA-seq, ATAC-seq, and ChIP-seq processing.",
       allow_abbrev=False,
   )

   # Define command-line options
   parser.add_argument('-i', '--experiment-dir', help='Path to the experiment directory')
   parser.add_argument('-g', '--genome', help='Genome assembly name. Must match an available directory under the configured genome assembly root')
   parser.add_argument('-j', '--mode', help='Job mode. Can be auto, all or a range of jobs. See readme for some examples. \n \t Default: auto')
   parser.add_argument('-T', '--trim-tool', help='Trimming tool choice. Can be fastp or skewer. \n \t Default: fastp')
   parser.add_argument('-M', '--map-tool', help='Mapping tool choice. Can be HISAT2, STAR, or STAR_TE. STAR(_TE) can only be used for RNA-seq data. \n \t Default: HISAT2')
   parser.add_argument('--de-formula', help='Explicit DESeq2 design formula for step 12. If provided, it overrides --de-columns and --de-block.')
   parser.add_argument('--de-config', action='append', help='Optional YAML file with DE analysis settings for step 12. Can be passed multiple times to run multiple DE analyses sequentially.')
   parser.add_argument('--de-out-dir', help='Optional DE output subdirectory inside the DE_calling folder. Overrides de_config io.out_dir.')
   parser.add_argument('--de-enable-custom-modules', action='store_true', help='Enable optional custom module enrichment (phase 3) in step 12. Requires a GMT file from --de-custom-modules-gmt or de_config enrichment.custom_modules.gmt_file.')
   parser.add_argument('--de-custom-modules-gmt', help='Path to a custom GMT file for optional custom module enrichment in step 12.')
   parser.add_argument('-I', '--input', help='Input BAM file used for ChIP peak calling with MACS3. \n \t Default: do not use input')
   parser.add_argument('-m', '--metadata', help='Tabular metadata file. The first column must be named filename. Metadata drives sample naming, peak grouping, trackhub grouping, and DE design.')
   parser.add_argument('--broad-mode', choices=['off', 'domain', 'genebody', 'diffuse'], help='ChIP broad-mark handling mode. off keeps the narrow/TF-like path, domain uses MACS3 broad domains, genebody will use gene-body features, and diffuse will use tiled windows. Default: off')
   parser.add_argument('--chip-broad-qvalue', type=float, help='Relaxed MACS3 q-value used for ChIP broad domain mode pooled/replicate/pseudoreplicate calls. Default: 0.05')
   parser.add_argument('--chip-broad-cutoff', type=float, help='MACS3 --broad-cutoff value for ChIP broad domain mode. Default: 0.1')
   parser.add_argument('--chip-broad-min-length', type=int, help='Optional MACS3 --min-length for ChIP broad domain mode. Default: unset')
   parser.add_argument('--chip-broad-max-gap', type=int, help='Optional MACS3 --max-gap for ChIP broad domain mode. Default: unset')
   parser.add_argument('--chip-broad-replicate-fraction', type=float, help='Minimum fraction of true replicates that must support a pooled broad domain in ChIP domain mode. Value in [0,1]. Default: 1.0')
   parser.add_argument('--chip-broad-overlap-fraction', type=float, help='Minimum reciprocal overlap fraction used when evaluating support for pooled broad domains in ChIP domain mode. Value in (0,1]. Default: 0.5')
   parser.add_argument('--chip-diffuse-bin-size', type=int, help='Fixed genomic bin size in bp for ChIP diffuse mode. Default: 10000')
   parser.add_argument('--chip-diffuse-merge-gap', type=int, help='Maximum gap in bp used to merge adjacent significant diffuse bins into domains. Default: one bin width')
   parser.add_argument('--chip-diffuse-keep-nonstandard-chroms', action='store_true', help='Keep nonstandard chromosomes in ChIP diffuse mode. Default filters to chr1-22,X,Y style chromosomes only.')
   parser.add_argument('--chip-diffuse-keep-chrm', action='store_true', help='Keep chrM/MT bins in ChIP diffuse mode. Default excludes chrM/MT bins.')
   parser.add_argument('-C', '--col-table', help='File specifying which colors to use for the tracks. \n \t Default: gray.tint.color.table from the packaged palette directory. Can be a *txt list file with one color table per line. Different color tables will be used per hub when multiple sample_type groups are present. Can be a full path or a file basename in combination with -P. Use `omnomnomics create-track-color-table` to build custom palettes and `omnomnomics display-track-color-table` to preview them.')
   parser.add_argument('-P', '--color-data-folder', help='Path to a folder with color tables. \n \t Default: packaged color_data_for_hubs directory. Use this to point the workflow at custom palettes.') ##################change this default??
   parser.add_argument('-o', '--overlay', help='Overlay type (transparentOverlay|stacked|solidOverlay|none) \n \t Default: transparentOverlay')
   parser.add_argument('-L', '--hub-mail', help='Email to use in trackhub \n \t Default: your@email.com')
   parser.add_argument('-X', '--no-multiqc', action='store_true', help='Exclude multiQC stats aggregator. Set if you don not wish to run multiQC.')
   parser.add_argument('--create-homer-tagdirs', action='store_true', help='Create optional HOMER tag directory exports in addition to the main pipeline outputs.')
   parser.add_argument('--rerun-selected-steps', action='store_true', help='Force recomputation of the selected workflow steps by deleting their existing outputs before submission. Default behavior is to reuse existing outputs when Snakemake considers them up to date.')
   parser.add_argument('--remove-duplicates', action='store_true', help='Remove duplicate reads in step 5. Default is assay-aware: keep for RNA, remove for ATAC and ChIP.')
   parser.add_argument('--keep-duplicates', action='store_true', help='Keep duplicate reads in step 5. Default is assay-aware: keep for RNA, remove for ATAC and ChIP.')
   parser.add_argument('--sample-name', help='Metadata columns used to derive sample_id. Accepts comma-separated column names or 1-based indices. Required when a metadata-driven step is selected.')
   parser.add_argument('--sample-type', help='Metadata columns used to derive sample_type. Accepts comma-separated column names or 1-based indices. Default: all_samples')
   parser.add_argument('--sample-color', help='Metadata columns used to derive sample_color categories. Accepts comma-separated column names or 1-based indices. Default: sample_type, else all_samples')
   parser.add_argument('--de-columns', help='Metadata columns of interest for auto-building the DESeq2 design when --de-formula is not given.')
   parser.add_argument('--de-block', help='Metadata columns to include as blocking terms when auto-building the DESeq2 design.')
   parser.add_argument('--de-interactions', action='store_true', help='Include interaction terms when auto-building the DESeq2 design from exactly two --de-columns values.')
   parser.add_argument('--atac-peak-opt-mode', choices=['none', 'fast', 'full'], help='ATAC step-10 peak optimization mode. none runs one fixed MACS3 call per group, fast runs a reduced candidate grid, full runs the full candidate grid.')
   parser.add_argument('--chip-peak-opt-mode', choices=['none', 'fast', 'full'], help='ChIP narrow (TF-like) step-10 peak optimization mode. none runs one fixed MACS3 call per group, fast runs a reduced candidate grid, full runs the full candidate grid.')
   parser.add_argument('--narrow-peak-strategy', choices=['idr', 'macs3'], help='Narrow-peak strategy for ATAC and narrow ChIP in step 10. idr uses replicate-consensus IDR (default), macs3 uses MACS3-only optimization.')
   parser.add_argument('--idr-mode', choices=['basic', 'encode'], help='IDR mode when --narrow-peak-strategy idr is active. basic runs true-replicate pairwise IDR consensus; encode additionally runs pseudo-replicate IDR diagnostics.')
   parser.add_argument('--idr-pair-fraction', type=float, help='Minimum fraction of replicate pairs that must support a peak in IDR basic/encode consensus. Value in [0,1]. Default: 0.5')
   parser.add_argument('--idr-pairing-policy', choices=['all_pairs', 'anchor_vs_all'], help='Replicate pairing policy for IDR when groups have >2 replicates. all_pairs uses all pairwise combinations; anchor_vs_all pairs the first replicate with every other replicate.')
   parser.add_argument('--spp-gate', choices=['none', 'warn', 'drop', 'strict'], help='SPP QC gate mode for ATAC/ChIP peak QC. none disables SPP gating, warn reports flags only, drop excludes flagged samples from downstream count/DE, strict aborts if any sample fails thresholds.')
   parser.add_argument('-a', '--appendix', help='Appendix to add to track name \n \t Default: hub')
   parser.add_argument('-k', '--keepunpaired', action='store_true', help='Keep unpaired or not in HISAT2')
   parser.add_argument('--dry-run', action='store_true', help='Validate the workflow and build the Snakemake DAG without executing jobs')
   parser.add_argument('--site-config', help='Optional path to a site-specific config YAML. Default: $XDG_CONFIG_HOME/omnomnomics/site.yaml or ~/.config/omnomnomics/site.yaml, then packaged site config')
   parser.add_argument('--retention-policy', choices=['all', 'pruned', 'minimal'], help='Post-run output retention policy. all keeps everything, pruned keeps FASTQ plus reusable downstream outputs, minimal keeps FASTQ plus only the requested terminal outputs. Default: all')
   parser.add_argument('--max-project-size', help='Soft project-size cap such as 200G or 800GB. Omnomnomics may delete safe intermediates and skip BigWig or trackhub creation when the cap would be exceeded.')

   args, unknown = parser.parse_known_args(argv)
   # Check if any unknown arguments are provided
   if unknown:
       parser.error(f"Unrecognized arguments: {' '.join(unknown)}")

   # Check if at least one argument is provided besides the script name
   argv_to_check = argv if argv is not None else sys.argv[1:]
   if len(argv_to_check) < 1:
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


def positive_int_or_none(value):
   try:
       parsed = int(value)
   except (TypeError, ValueError):
       return None
   return parsed if parsed > 0 else None


def apply_assay_runtime_defaults(config, assay_type):
   if assay_type not in {"ATAC", "CHIP"}:
       return

   chromatin_default_runtime = positive_int_or_none(config.get("chromatin_default_runtime"))
   if chromatin_default_runtime is not None:
       current_default_runtime = positive_int_or_none(config.get("default_runtime")) or 0
       config["default_runtime"] = max(current_default_runtime, chromatin_default_runtime)

   chromatin_controller_runtime = positive_int_or_none(config.get("chromatin_controller_runtime"))
   if chromatin_controller_runtime is not None:
       current_controller_runtime = positive_int_or_none(config.get("controller_runtime")) or 0
       config["controller_runtime"] = max(current_controller_runtime, chromatin_controller_runtime)

   chromatin_rule_runtime = config.get("chromatin_rule_runtime")
   if not isinstance(chromatin_rule_runtime, list):
       return

   rule_runtime = list(config.get("rule_runtime", []))
   for index, chromatin_runtime in enumerate(chromatin_rule_runtime):
       parsed_chromatin_runtime = positive_int_or_none(chromatin_runtime)
       if parsed_chromatin_runtime is None:
           continue
       while len(rule_runtime) <= index:
           rule_runtime.append(None)
       current_rule_runtime = positive_int_or_none(rule_runtime[index]) or 0
       rule_runtime[index] = max(current_rule_runtime, parsed_chromatin_runtime)
   config["rule_runtime"] = rule_runtime


def normalize_de_config_paths(args_de_config, config_de_config):
    if args_de_config:
        raw_items = args_de_config
    else:
        raw_items = config_de_config

    if raw_items in (None, "", "NA"):
        return []
    if isinstance(raw_items, str):
        raw_items = [raw_items]
    if not isinstance(raw_items, list):
        raw_items = [str(raw_items)]

    normalized = []
    for item in raw_items:
        if item in (None, "", "NA"):
            continue
        normalized.append(str(Path(str(item)).expanduser().resolve()))
    return normalized


def de_config_has_explicit_out_dir(config_path):
    try:
        loaded = yaml.safe_load(Path(config_path).read_text()) or {}
    except Exception:
        return False
    if not isinstance(loaded, dict):
        return False
    io_cfg = loaded.get("io")
    if not isinstance(io_cfg, dict):
        return False
    out_dir_value = io_cfg.get("out_dir")
    return out_dir_value is not None and str(out_dir_value).strip() not in ("", "NA")

def resolve_config_path(path_value, workflow_root):
   resolved = str(path_value)
   resolved = resolved.replace("{WORKFLOW_ROOT}", str(workflow_root))
   resolved = resolved.replace("{HOME}", str(Path.home()))
   return str(Path(resolved).expanduser().resolve())

def parse_size_to_bytes(size_value):
   if size_value is None:
       return 0
   if isinstance(size_value, (int, float)):
       return int(size_value)

   size_text = str(size_value).strip().upper()
   if not size_text or size_text == "NA":
       return 0

   match = re.fullmatch(r"(\d+(?:\.\d+)?)\s*([KMGT]?B?)?", size_text)
   if not match:
       raise ValueError(f"Invalid size value: {size_value}")

   number = float(match.group(1))
   suffix = match.group(2) or "B"
   multipliers = {
       "B": 1,
       "K": 1024,
       "KB": 1024,
       "M": 1024 ** 2,
       "MB": 1024 ** 2,
       "G": 1024 ** 3,
       "GB": 1024 ** 3,
       "T": 1024 ** 4,
       "TB": 1024 ** 4,
   }
   if suffix not in multipliers:
       raise ValueError(f"Invalid size suffix in value: {size_value}")
   return int(number * multipliers[suffix])

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
def validate_user_defined_vars(workflow_root, metadata, experiment_dir, INPUT, color_data_folder, col_table, overlay, the_type, map_tool, config, de_config):
    #Validate user-defined variables
    print("Validating options...")


    # Check metadata file
    if metadata != "NA":
        if not os.path.isfile(metadata):
            print("Metadata file (-m) does not exist! Aborting...", file=sys.stderr)
            sys.exit(1)
        elif not metadata.endswith((".txt", ".tsv", ".csv")):
            print("Metadata file (-m) must be a .txt, .tsv, or .csv file. Aborting...", file=sys.stderr)
            sys.exit(1)

    de_config_paths = []
    if isinstance(de_config, list):
        de_config_paths = [item for item in de_config if item and item != "NA"]
    elif de_config and de_config != "NA":
        de_config_paths = [de_config]
    for de_cfg in de_config_paths:
        if not os.path.isfile(de_cfg):
            print(f"DE config file (--de-config) does not exist: {de_cfg}. Aborting...", file=sys.stderr)
            sys.exit(1)
        if not de_cfg.endswith((".yaml", ".yml")):
            print(f"DE config file (--de-config) must be a .yaml or .yml file: {de_cfg}. Aborting...", file=sys.stderr)
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


def resolve_public_mode_steps(mode, public_max_step):
    if mode == "auto":
        mode = "all"
    if mode == "all":
        return list(range(1, public_max_step + 1))

    if re.match(r'^[0-9]+$', mode):
        steps = [int(mode)]
    elif re.match(r'^[0-9]+-[0-9]+$', mode):
        start, end = map(int, mode.split('-'))
        if start > end:
            print("Job mode range has to increase from start to end! Aborting...")
            sys.exit(1)
        steps = list(range(start, end + 1))
    else:
        steps = expand_range(mode)
        if not all(isinstance(step, int) for step in steps):
            print("Job mode should be 'all', 'auto', a single step number, a numeric range of steps separated by a dash (e.g., 1-12), or a number of ranges and steps separated by commas (e.g., 1-3,6,8,10-12).")
            print("Aborting...")
            sys.exit(1)
        if any(steps[i] >= steps[i + 1] for i in range(len(steps) - 1)):
            print("Job mode range has to increase from start to end! Aborting...")
            sys.exit(1)

    if any(step < 1 or step > public_max_step for step in steps):
        print(f"Job mode range should lie between 1 and {public_max_step} for this assay. Aborting...")
        sys.exit(1)
    return sorted(steps)


def map_public_to_internal_steps(assay_type, public_steps):
    mapping = ASSAY_PUBLIC_INTERNAL_STEP_MAP.get(assay_type, {})
    internal = []
    for step in public_steps:
        if step not in mapping:
            print(f"Public step {step} is not defined for assay {assay_type}. Aborting...", file=sys.stderr)
            sys.exit(1)
        mapped = mapping[step]
        if mapped not in internal:
            internal.append(mapped)
    return sorted(internal)

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

    # For rules with multiple input filetypes, use the first entry as primary sanity input
    if isinstance(input_file_type_mod_range_min, list):
        input_file_type_mod_range_min = input_file_type_mod_range_min[0]
    if isinstance(input_folder_mod_range_min, list):
        input_folder_mod_range_min = input_folder_mod_range_min[0]


    if os.path.isdir(f"{experiment_dir}/{input_folder_mod_range_min}"):
        num_files = len(glob.glob(f"{experiment_dir}/{input_folder_mod_range_min}/*{input_file_type_mod_range_min}"))


    # Sanity check file number
    if num_files == 0:
        print("No input files detected! Aborting...", file=sys.stderr)
        sys.exit(1)

    # Check if input files are readable
    input_files = glob.glob(f"{experiment_dir}/{input_folder_mod_range_min}/*{input_file_type_mod_range_min}")
    if not os.access(input_files[0], os.R_OK):
        print(f"Permission error! {input_file_type_mod_range_min} files in {experiment_dir}/{input_folder_mod_range_min} are not readable! Aborting...", file=sys.stderr)
        sys.exit(1)


    # Set the number of pairs (if dealing with FASTQs, else just keep it equal to file number)
    num_pairs = num_files // 2 if paired and mode_range_min < 4 else num_files

    return num_files, num_pairs, paired, input_folder_mod_range_min, input_file_type_mod_range_min

def normalize_field_selection_name(filename, file_type):
    _ = file_type
    return normalize_metadata_filename(filename)

def merged_sample_roots_for_mode(experiment_dir, input_folder, input_file_type):
    if input_file_type not in FASTQ_EXTENSIONS and input_file_type not in (".trimmed.fastq.gz", ".trimmed.fastq", ".trimmed.fq.gz", ".trimmed.fq"):
        files = glob.glob(f"{experiment_dir}/{input_folder}/*{input_file_type}")
    else:
        files = []
        for extension in FASTQ_EXTENSIONS if input_file_type in FASTQ_EXTENSIONS else (input_file_type,):
            files.extend(glob.glob(f"{experiment_dir}/{input_folder}/*{extension}"))

    normalized_files = sorted(
        set(normalize_field_selection_name(file_path, input_file_type) for file_path in files)
    )
    return sorted(set(MERGED_LANE_SUFFIX_RE.sub("", file_name) for file_name in normalized_files))


def metadata_required_for_mode(mode_steps):
    return bool(METADATA_REQUIRED_STEPS.intersection(mode_steps))


def validate_metadata_sample_roots(derived_rows, expected_sample_roots):
    metadata_roots = []
    for row in derived_rows:
        filename_key = row.get("filename_key", "").strip()
        if filename_key:
            metadata_roots.append(filename_key)
            continue
        metadata_roots.append(normalize_metadata_filename(row["filename"]))
    metadata_roots = sorted(metadata_roots)
    expected_roots = sorted(expected_sample_roots)
    missing = sorted(set(expected_roots) - set(metadata_roots))
    extra = sorted(set(metadata_roots) - set(expected_roots))
    if missing or extra:
        problems = []
        if missing:
            problems.append("missing metadata rows for: " + ", ".join(missing))
        if extra:
            problems.append("metadata rows without matching samples: " + ", ".join(extra))
        raise MetadataError("Metadata sample matching failed: " + "; ".join(problems))


def count_table_sample_ids(count_table_path):
    with open(count_table_path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
    if len(header) < 2:
        raise MetadataError(f"Count table '{count_table_path}' is missing sample columns.")
    return header[1:]


def validate_metadata_sample_ids(derived_rows, expected_sample_ids):
    expected_id_set = set(expected_sample_ids)
    metadata_sample_id_set = {row["sample_id"] for row in derived_rows}
    metadata_filename_key_set = {
        (row.get("filename_key", "").strip() or normalize_metadata_filename(row["filename"]))
        for row in derived_rows
    }

    if expected_id_set == metadata_sample_id_set:
        return
    if expected_id_set == metadata_filename_key_set:
        return

    missing = sorted(expected_id_set - metadata_sample_id_set)
    extra = sorted(metadata_sample_id_set - expected_id_set)
    problems = []
    if missing:
        problems.append("missing metadata rows for count-table samples: " + ", ".join(missing))
    if extra:
        problems.append("metadata rows without matching count-table samples: " + ", ".join(extra))
    raise MetadataError(
        "Metadata sample matching failed: "
        + "; ".join(problems)
        + ". Count table headers did not match either metadata sample_id or normalized metadata filename keys."
    )


def validate_deseq_design_full_rank(metadata_path, formula):
    r_expression = (
        "args <- commandArgs(trailingOnly = TRUE); "
        "metadata_path <- args[[1]]; "
        "design_formula <- args[[2]]; "
        "metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE); "
        "metadata[] <- lapply(metadata, function(column) if (is.character(column)) factor(column) else column); "
        "design <- as.formula(design_formula); "
        "model_matrix <- model.matrix(design, data = metadata); "
        "if (nrow(model_matrix) == 0) stop('Design matrix is empty.'); "
        "if (qr(model_matrix)$rank < ncol(model_matrix)) stop('Design matrix is not full rank.'); "
        "cat('Design formula OK\\n')"
    )
    try:
        completed = subprocess.run(
            ["Rscript", "-e", r_expression, str(metadata_path), formula],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise MetadataError("Rscript is required to validate the DESeq2 design but was not found.") from exc
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.strip() or exc.stdout.strip() or "Unknown R validation error."
        raise MetadataError(f"DESeq2 design validation failed: {stderr}") from exc
    return completed.stdout.strip()

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
    if len(sys.argv) == 1:
        print_top_level_help(exit_code=1)
    if len(sys.argv) > 1 and sys.argv[1] in {"-h", "--help"}:
        print_top_level_help(exit_code=0)

    forwarded_argv = None
    assay = None
    assay_verbs = {"rna", "atac", "chip"}
    if len(sys.argv) > 1 and sys.argv[1].lower() in assay_verbs:
        assay = sys.argv[1].lower()
        forwarded_argv = list(sys.argv[2:])

    if len(sys.argv) > 1 and sys.argv[1] == "genomes":
        genomes_main(sys.argv[2:], WORKFLOW_ROOT, DEFAULT_WORKFLOW_CONFIG, DEFAULT_SITE_CONFIG)
        return
    if len(sys.argv) > 1 and sys.argv[1] == "monitor":
        monitor_main(sys.argv[2:])
        return
    if len(sys.argv) > 1 and sys.argv[1] == "de-app":
        de_app_main(sys.argv[2:])
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

    if assay is None:
        print(
            "Please select an assay verb: 'omnomnomics rna ...', 'omnomnomics atac ...', or 'omnomnomics chip ...'. Aborting...",
            file=sys.stderr,
        )
        sys.exit(1)

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
    args = parse_arguments(forwarded_argv)

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
    the_type = assay.upper()
    apply_assay_runtime_defaults(config, the_type)
    genome = args.genome
    trim_tool = args.trim_tool.lower() if args.trim_tool else config.get("trim_tool","fastp").lower()
    map_tool = args.map_tool.lower() if args.map_tool else config.get('map_tool', "hisat2").lower()
    mode = args.mode.lower() if args.mode else config.get('mode', "auto")
    de_formula = args.de_formula if args.de_formula else config.get('de_formula', "NA")
    de_config_list = normalize_de_config_paths(args.de_config, config.get('de_config', "NA"))
    de_config = de_config_list[0] if de_config_list else "NA"
    de_out_dir = args.de_out_dir if args.de_out_dir else config.get('de_out_dir', "")
    de_enable_custom_modules = args.de_enable_custom_modules or config.get('de_enable_custom_modules', False)
    de_custom_modules_gmt = (
        str(Path(args.de_custom_modules_gmt).expanduser().resolve())
        if args.de_custom_modules_gmt
        else str(config.get('de_custom_modules_gmt', "")).strip()
    )
    if de_custom_modules_gmt and de_custom_modules_gmt != "NA":
        de_custom_modules_gmt = str(Path(de_custom_modules_gmt).expanduser().resolve())
    else:
        de_custom_modules_gmt = ""
    broad_mode = (args.broad_mode if args.broad_mode else config.get('broad_mode', 'off')).strip().lower()
    chip_broad_qvalue = args.chip_broad_qvalue if args.chip_broad_qvalue is not None else config.get('chip_broad_qvalue', 0.05)
    chip_broad_cutoff = args.chip_broad_cutoff if args.chip_broad_cutoff is not None else config.get('chip_broad_cutoff', 0.1)
    chip_broad_min_length = args.chip_broad_min_length if args.chip_broad_min_length is not None else config.get('chip_broad_min_length', None)
    chip_broad_max_gap = args.chip_broad_max_gap if args.chip_broad_max_gap is not None else config.get('chip_broad_max_gap', None)
    chip_broad_replicate_fraction = args.chip_broad_replicate_fraction if args.chip_broad_replicate_fraction is not None else config.get('chip_broad_replicate_fraction', 1.0)
    chip_broad_overlap_fraction = args.chip_broad_overlap_fraction if args.chip_broad_overlap_fraction is not None else config.get('chip_broad_overlap_fraction', 0.5)
    chip_diffuse_bin_size = args.chip_diffuse_bin_size if args.chip_diffuse_bin_size is not None else config.get('chip_diffuse_bin_size', 10000)
    chip_diffuse_merge_gap = args.chip_diffuse_merge_gap if args.chip_diffuse_merge_gap is not None else config.get('chip_diffuse_merge_gap', "auto")
    chip_diffuse_standard_chroms_only = not args.chip_diffuse_keep_nonstandard_chroms
    if 'chip_diffuse_standard_chroms_only' in config and not args.chip_diffuse_keep_nonstandard_chroms:
        chip_diffuse_standard_chroms_only = bool(config.get('chip_diffuse_standard_chroms_only', True))
    chip_diffuse_exclude_chrm = not args.chip_diffuse_keep_chrm
    if 'chip_diffuse_exclude_chrm' in config and not args.chip_diffuse_keep_chrm:
        chip_diffuse_exclude_chrm = bool(config.get('chip_diffuse_exclude_chrm', True))
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
    sample_name = args.sample_name if args.sample_name else config.get('sample_name', "")
    sample_type = args.sample_type if args.sample_type else config.get('sample_type', "")
    sample_color = args.sample_color if args.sample_color else config.get('sample_color', "")
    de_columns = args.de_columns if args.de_columns else config.get('de_columns', "")
    de_block = args.de_block if args.de_block else config.get('de_block', "")
    de_interactions = args.de_interactions or config.get('de_interactions', False)
    atac_peak_opt_mode = (args.atac_peak_opt_mode if args.atac_peak_opt_mode else config.get('atac_peak_opt_mode', 'fast')).strip().lower()
    chip_peak_opt_mode = (args.chip_peak_opt_mode if args.chip_peak_opt_mode else config.get('chip_peak_opt_mode', 'fast')).strip().lower()
    narrow_peak_strategy = (args.narrow_peak_strategy if args.narrow_peak_strategy else config.get('narrow_peak_strategy', 'idr')).strip().lower()
    idr_mode = (args.idr_mode if args.idr_mode else config.get('idr_mode', 'encode')).strip().lower()
    idr_pair_fraction = args.idr_pair_fraction if args.idr_pair_fraction is not None else config.get('idr_pair_fraction', 0.5)
    idr_pairing_policy = (args.idr_pairing_policy if args.idr_pairing_policy else config.get('idr_pairing_policy', 'all_pairs')).strip().lower()
    try:
        idr_pair_fraction = float(idr_pair_fraction)
    except (TypeError, ValueError):
        print("--idr-pair-fraction must be a numeric value in [0,1]. Aborting...", file=sys.stderr)
        sys.exit(1)
    if not (0.0 <= idr_pair_fraction <= 1.0):
        print("--idr-pair-fraction must be within [0,1]. Aborting...", file=sys.stderr)
        sys.exit(1)
    try:
        chip_broad_qvalue = float(chip_broad_qvalue)
        chip_broad_cutoff = float(chip_broad_cutoff)
        chip_broad_replicate_fraction = float(chip_broad_replicate_fraction)
        chip_broad_overlap_fraction = float(chip_broad_overlap_fraction)
    except (TypeError, ValueError):
        print("ChIP broad-domain numeric parameters must be valid numbers. Aborting...", file=sys.stderr)
        sys.exit(1)
    if not (0.0 < chip_broad_qvalue <= 1.0):
        print("--chip-broad-qvalue must be within (0,1]. Aborting...", file=sys.stderr)
        sys.exit(1)
    if not (0.0 < chip_broad_cutoff <= 1.0):
        print("--chip-broad-cutoff must be within (0,1]. Aborting...", file=sys.stderr)
        sys.exit(1)
    if not (0.0 <= chip_broad_replicate_fraction <= 1.0):
        print("--chip-broad-replicate-fraction must be within [0,1]. Aborting...", file=sys.stderr)
        sys.exit(1)
    if not (0.0 < chip_broad_overlap_fraction <= 1.0):
        print("--chip-broad-overlap-fraction must be within (0,1]. Aborting...", file=sys.stderr)
        sys.exit(1)
    try:
        chip_diffuse_bin_size = int(chip_diffuse_bin_size)
    except (TypeError, ValueError):
        print("--chip-diffuse-bin-size must be an integer. Aborting...", file=sys.stderr)
        sys.exit(1)
    if chip_diffuse_bin_size <= 0:
        print("--chip-diffuse-bin-size must be > 0. Aborting...", file=sys.stderr)
        sys.exit(1)
    if chip_diffuse_merge_gap in ("", "NA", None, "auto"):
        chip_diffuse_merge_gap = chip_diffuse_bin_size
    try:
        chip_diffuse_merge_gap = int(chip_diffuse_merge_gap)
    except (TypeError, ValueError):
        print("--chip-diffuse-merge-gap must be an integer or omitted. Aborting...", file=sys.stderr)
        sys.exit(1)
    if chip_diffuse_merge_gap < 0:
        print("--chip-diffuse-merge-gap must be >= 0. Aborting...", file=sys.stderr)
        sys.exit(1)
    if chip_broad_min_length in ("", "NA"):
        chip_broad_min_length = None
    if chip_broad_max_gap in ("", "NA"):
        chip_broad_max_gap = None
    if chip_broad_min_length is not None:
        try:
            chip_broad_min_length = int(chip_broad_min_length)
        except (TypeError, ValueError):
            print("--chip-broad-min-length must be an integer. Aborting...", file=sys.stderr)
            sys.exit(1)
        if chip_broad_min_length <= 0:
            print("--chip-broad-min-length must be > 0. Aborting...", file=sys.stderr)
            sys.exit(1)
    if chip_broad_max_gap is not None:
        try:
            chip_broad_max_gap = int(chip_broad_max_gap)
        except (TypeError, ValueError):
            print("--chip-broad-max-gap must be an integer. Aborting...", file=sys.stderr)
            sys.exit(1)
        if chip_broad_max_gap < 0:
            print("--chip-broad-max-gap must be >= 0. Aborting...", file=sys.stderr)
            sys.exit(1)
    spp_gate = (args.spp_gate if args.spp_gate else config.get('spp_gate', 'warn')).strip().lower()
    appendix = args.appendix if args.appendix else config.get('appendix', "hub")
    keep_unpaired = args.keepunpaired if args.keepunpaired else config.get('keep_unpaired', False)
    dry_run = args.dry_run
    retention_policy = args.retention_policy if args.retention_policy else config.get('retention_policy', 'all')
    max_project_size_raw = args.max_project_size if args.max_project_size else config.get('max_project_size', "NA")
    try:
        max_project_size_bytes = parse_size_to_bytes(max_project_size_raw)
    except ValueError as exc:
        print(f"{exc}. Aborting...", file=sys.stderr)
        sys.exit(1)

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
    selected_routines['selected_routine_peakqc'] = config['selected_routine_peakqc']
    selected_routines['selected_routine_analyzepeaks'] = config['selected_routine_analyzepeaks']
    selected_routines['selected_routine_dechrom'] = config['selected_routine_dechrom']
    selected_routines['selected_routine_analyzepeaksde'] = config['selected_routine_analyzepeaksde']

    # Validate user-defined variables
    validate_user_defined_vars(workflow_root, metadata, experiment_dir, INPUT, color_data_folder, col_table, overlay, the_type, map_tool, config, de_config_list)

    # Setup variables
    run_date = setup_variables(experiment_dir, config)

    # Resolve assay-specific public job mode and map it to internal workflow steps
    public_step_map = ASSAY_PUBLIC_INTERNAL_STEP_MAP.get(the_type)
    if not public_step_map:
        print(f"No public step map is defined for assay type {the_type}. Aborting...", file=sys.stderr)
        sys.exit(1)
    public_max_step = max(public_step_map)
    public_mode_steps = resolve_public_mode_steps(mode, public_max_step)
    mode_steps = map_public_to_internal_steps(the_type, public_mode_steps)
    print(f"PUBLIC MODE STEPS = {public_mode_steps}")
    print(f"INTERNAL MODE STEPS = {mode_steps}")

    # Check if pipeline needed or user first has to run separate scripts
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
    print(f"INTERNAL MODE STEPS NEW = {mode_steps}")

    # Reset selected step bookkeeping for the new run
    reset_step_tracking(mode_steps, experiment_dir)

    # Optionally force recomputation of selected step outputs
    if rerun_selected_steps:
        delete_outputs_to_be_updated(mode_steps, config, experiment_dir)
    else:
        print("REUSING EXISTING OUTPUTS WHEN POSSIBLE")

    #checking input files
    num_files, num_pairs, paired, input_folder_mod_range_min, input_file_type_mod_range_min = validate_input_files(the_type, config, min(mode_steps),experiment_dir)

    metadata_required = metadata_required_for_mode(mode_steps)
    derived_metadata_path = "NA"
    resolved_de_formula = "NA"
    de_design_mode = "NA"
    sample_name_columns = []
    sample_type_columns = []
    sample_color_columns = []
    de_columns_resolved = []
    de_block_resolved = []
    resolved_de_config_path = "NA"
    resolved_de_config = {}
    de_config_file_path = "NA"
    resolved_de_config_paths = []
    resolved_de_configs = []
    de_config_file_paths = []

    if metadata_required:
        if metadata == "NA":
            print(
                "A metadata file is required for the selected workflow steps. Provide -m/--metadata. Aborting...",
                file=sys.stderr,
            )
            sys.exit(1)
        if not sample_name:
            print(
                "--sample-name is required when running metadata-driven workflow steps. Aborting...",
                file=sys.stderr,
            )
            sys.exit(1)

        try:
            run_configs_dir = os.path.join(experiment_dir, "run_configs")
            os.makedirs(run_configs_dir, exist_ok=True)
            metadata_fieldnames, metadata_rows = read_metadata_table(metadata)
            derived_fieldnames, derived_rows, selector_map = derive_metadata_rows(
                metadata_fieldnames,
                metadata_rows,
                sample_name_selector=sample_name,
                sample_type_selector=sample_type,
                sample_color_selector=sample_color,
            )
            de_steps = {
                config.get('de_rule_num'),
                config.get('dechrom_rule_num'),
            }
            de_steps = {step for step in de_steps if isinstance(step, int)}
            if de_steps and min(mode_steps) in de_steps:
                count_table_path = os.path.join(
                    experiment_dir,
                    config['input_folders'][config['de_rule_num'] - 1],
                    f"{os.path.basename(experiment_dir)}.raw_read_quant.table.txt",
                )
                validate_metadata_sample_ids(
                    derived_rows,
                    count_table_sample_ids(count_table_path),
                )
            else:
                expected_sample_roots = merged_sample_roots_for_mode(
                    experiment_dir,
                    input_folder_mod_range_min,
                    input_file_type_mod_range_min,
                )
                validate_metadata_sample_roots(derived_rows, expected_sample_roots)
            sample_name_columns = selector_map["sample_name_columns"]
            sample_type_columns = selector_map["sample_type_columns"]
            sample_color_columns = selector_map["sample_color_columns"]

            if (config.get('de_rule_num') in mode_steps) or (config.get('dechrom_rule_num') in mode_steps):
                derived_fieldnames, derived_rows, resolved_de_formula, de_context = resolve_de_metadata(
                    metadata_fieldnames,
                    derived_fieldnames,
                    derived_rows,
                    None if de_formula == "NA" else de_formula,
                    de_columns,
                    de_block,
                    bool(de_interactions),
                )
                de_design_mode = str(de_context["mode"])
                de_columns_resolved = list(de_context["de_columns"])
                de_block_resolved = list(de_context["de_block"])
                if len(de_config_list) > 1 and de_out_dir:
                    raise MetadataError(
                        "Global --de-out-dir cannot be combined with multiple --de-config files. "
                        "Set io.out_dir inside each DE config YAML."
                    )

                de_configs_to_resolve = de_config_list if de_config_list else ["NA"]
                resolved_out_dirs = []
                for idx, one_de_config in enumerate(de_configs_to_resolve, start=1):
                    if one_de_config != "NA" and len(de_configs_to_resolve) > 1:
                        if not de_config_has_explicit_out_dir(one_de_config):
                            raise MetadataError(
                                "When multiple --de-config files are provided, each config must explicitly set io.out_dir. "
                                f"Missing io.out_dir in: {one_de_config}"
                            )
                    try:
                        one_file_path, one_resolved = resolve_de_config(
                            None if one_de_config == "NA" else one_de_config,
                            resolved_de_formula,
                            de_out_dir_override=(de_out_dir or None) if len(de_configs_to_resolve) == 1 else None,
                            custom_modules_gmt=de_custom_modules_gmt or None,
                            enable_custom_modules=bool(de_enable_custom_modules),
                        )
                    except DEConfigError as exc:
                        raise MetadataError(str(exc)) from exc

                    one_out_dir = str(one_resolved.get("io", {}).get("out_dir", "")).strip()
                    if one_out_dir in resolved_out_dirs:
                        raise MetadataError(
                            "Multiple DE configs resolved to the same io.out_dir. "
                            f"Use unique out_dir values per config. Duplicate: {one_out_dir}"
                        )
                    resolved_out_dirs.append(one_out_dir)
                    de_config_file_paths.append(one_file_path)
                    resolved_de_configs.append(one_resolved)

                    one_resolved_path = os.path.join(
                        run_configs_dir,
                        f"omnomnomics.run.{run_date}.de_config.resolved.{idx:02d}.yaml",
                    )
                    with open(one_resolved_path, "w") as handle:
                        yaml.safe_dump(one_resolved, handle, sort_keys=False)
                    resolved_de_config_paths.append(one_resolved_path)

                resolved_de_config = resolved_de_configs[0]
                de_config_file_path = de_config_file_paths[0]
                resolved_de_config_path = resolved_de_config_paths[0]
                resolved_de_formula = str(resolved_de_config["design"]["formula"])
                if resolved_de_config["design"].get("formula") and de_config_file_path != "NA":
                    de_design_mode = "config_formula"

            derived_metadata_path = os.path.join(
                run_configs_dir,
                f"omnomnomics.run.{run_date}.metadata_derived.tsv",
            )
            write_metadata_table(derived_metadata_path, derived_fieldnames, derived_rows)

            if 12 in mode_steps:
                validation_message = validate_deseq_design_full_rank(
                    derived_metadata_path,
                    resolved_de_formula,
                )
                print(validation_message)
        except MetadataError as exc:
            print(f"{exc} Aborting...", file=sys.stderr)
            sys.exit(1)

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
        'default_runtime': config.get('default_runtime'),
        'controller_runtime': config.get('controller_runtime'),
        'rule_runtime': config.get('rule_runtime'),
        'HISAT2_GENOME_DIR': genome_subdir(config['genome_assembly_dir'], genome, "hisat2"),
        'STAR_GENOME_DIR': genome_subdir(config['genome_assembly_dir'], genome, "star"),
        'GENOME_AUX_DIR': genome_subdir(config['genome_assembly_dir'], genome, "aux"),
        'EXPERIMENT_DIR': experiment_dir,
        'RUNDATE': run_date,
        'INPUT_FOLDER': input_folder_mod_range_min,
        'INPUT_FILE_TYPE': input_file_type_mod_range_min,
        'INPUT': INPUT,
        'BROAD_MODE': broad_mode,
        'CHIP_BROAD_QVALUE': chip_broad_qvalue,
        'CHIP_BROAD_CUTOFF': chip_broad_cutoff,
        'CHIP_BROAD_MIN_LENGTH': chip_broad_min_length if chip_broad_min_length is not None else "NA",
        'CHIP_BROAD_MAX_GAP': chip_broad_max_gap if chip_broad_max_gap is not None else "NA",
        'CHIP_BROAD_REPLICATE_FRACTION': chip_broad_replicate_fraction,
        'CHIP_BROAD_OVERLAP_FRACTION': chip_broad_overlap_fraction,
        'CHIP_DIFFUSE_BIN_SIZE': chip_diffuse_bin_size,
        'CHIP_DIFFUSE_MERGE_GAP': chip_diffuse_merge_gap,
        'CHIP_DIFFUSE_STANDARD_CHROMS_ONLY': bool(chip_diffuse_standard_chroms_only),
        'CHIP_DIFFUSE_EXCLUDE_CHRM': bool(chip_diffuse_exclude_chrm),
        'THEMODE': mode_steps,
        'PUBLIC_MODE_STEPS': public_mode_steps,
        'THETYPE': the_type,
        'NUMFILES': num_files,
        'NUMPAIRS': num_pairs,
        'KEEPUNPAIRED': keep_unpaired,
        'PAIRED': paired,
        'THEMODERANGEMIN': min(mode_steps),
        'THEMODERANGEMAX': max(mode_steps),
        'THEGENOME': genome,
        'DE_FORMULA': de_formula,
        'RESOLVED_DE_FORMULA': resolved_de_formula,
        'DE_COLUMNS': de_columns,
        'DE_BLOCK': de_block,
        'DE_INTERACTIONS': bool(de_interactions),
        'ATAC_PEAK_OPT_MODE': atac_peak_opt_mode,
        'CHIP_PEAK_OPT_MODE': chip_peak_opt_mode,
        'NARROW_PEAK_STRATEGY': narrow_peak_strategy,
        'IDR_MODE': idr_mode,
        'IDR_PAIR_FRACTION': idr_pair_fraction,
        'IDR_PAIRING_POLICY': idr_pairing_policy,
        'SPP_GATE': spp_gate,
        'DE_DESIGN_MODE': de_design_mode,
        'DE_CONFIG_FILE': de_config_file_path,
        'DE_CONFIG_RESOLVED_FILE': resolved_de_config_path,
        'DE_CONFIG_RESOLVED_JSON': json.dumps(resolved_de_config, sort_keys=True) if resolved_de_config else "{}",
        'DE_CONFIG_FILES': de_config_file_paths,
        'DE_CONFIG_RESOLVED_FILES': resolved_de_config_paths,
        'DE_CONFIG_RESOLVED_LIST_JSON': json.dumps(resolved_de_configs, sort_keys=True) if resolved_de_configs else "[]",
        'MYMETADATA': metadata,
        'DERIVED_METADATA_FILE': derived_metadata_path,
        'METADATA_REQUIRED': metadata_required,
        'SAMPLE_NAME_SELECTOR': sample_name,
        'SAMPLE_TYPE_SELECTOR': sample_type,
        'SAMPLE_COLOR_SELECTOR': sample_color,
        'SAMPLE_NAME_COLUMNS': sample_name_columns,
        'SAMPLE_TYPE_COLUMNS': sample_type_columns,
        'SAMPLE_COLOR_COLUMNS': sample_color_columns,
        'DE_COLUMNS_RESOLVED': de_columns_resolved,
        'DE_BLOCK_RESOLVED': de_block_resolved,
        'THETRIMTOOL': trim_tool,
        'THEMAPTOOL': map_tool,
        'NO_MULTIQC': no_multiqc,
        'CREATE_HOMER_TAGDIRS': create_homer_tagdirs,
        'RERUN_SELECTED_STEPS': rerun_selected_steps,
        'DUPLICATE_HANDLING': duplicate_handling,
        'THEAPPENDIX': appendix,
        'THEOVERLAY': overlay,
        'THECOLORDATAFOLDER': color_data_folder,
        'THEHUBMAIL': hub_mail,
        'THECOLTABLE': col_table,
        'THEMEM': the_mem,
        'THEHEAPINIT': the_heap_init,
        'RETENTION_POLICY': retention_policy,
        'MAX_PROJECT_SIZE': str(max_project_size_raw),
        'MAX_PROJECT_SIZE_BYTES': int(max_project_size_bytes)
    }
    # Initialize the run config
    write_run_config(experiment_dir, run_date, run_config_data)

    # Initialize the log file
    log_file = start_log(experiment_dir, run_date, config)


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

    with open(log_file, 'a') as log:
        log.write(f"Invocation:\t{the_command}\n")

    if dry_run:
        print("")
        print("####################################################################################")
        print("Running omnomnomics dry-run...")
        print(f"\tMODE:\t\t{mode_steps}")
        print(f"\tTYPE:\t\t{the_type}")
        print(f"\tFILES:\t\t{num_files}")
        print(f"\tPAIRS:\t\t{num_pairs}")
        print(f"\tExperiment DIR:\t\t{experiment_dir}")
        print(f"\tRun date:\t{run_date}")
        print("####################################################################################")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Dry-run failed with exit code {e.returncode}.", file=sys.stderr)
            sys.exit(e.returncode)
        with open(log_file, 'a') as log:
            log.write("Dry-run completed successfully without controller submission.\n")
        print("Dry-run completed. No controller job was submitted.")
        return

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
        log.write(f"Controller submission at:\t{sub_time}\n")

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
