#!/bin/bash
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --job-name=chip-input-bench

set -euo pipefail
: "${BENCHMARK_SOURCE:?}"
: "${BENCHMARK_ROOT:?}"
eval "$(micromamba shell hook --shell bash)"
micromamba activate omnomnomics
cd "$BENCHMARK_SOURCE"
export PYTHONPATH="$BENCHMARK_SOURCE/src"
export NUMBA_CACHE_DIR="$TMPDIR/chip_benchmark_numba_${SLURM_JOB_ID}"
export MPLCONFIGDIR="$TMPDIR/chip_benchmark_mpl_${SLURM_JOB_ID}"
export XDG_CACHE_HOME="$TMPDIR/chip_benchmark_xdg_${SLURM_JOB_ID}"
export BENCHMARK_PARALLEL_IDR=1
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
git rev-parse HEAD
python --version
macs3 --version
samtools --version
featureCounts -v
Rscript --version
python scripts/validation/chip_input_benchmark.py --root "$BENCHMARK_ROOT" --index "${SLURM_ARRAY_TASK_ID:-0}" --threads "$SLURM_CPUS_PER_TASK"
