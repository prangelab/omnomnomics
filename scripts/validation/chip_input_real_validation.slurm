#!/bin/bash
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=01:30:00
#SBATCH --job-name=chip-input-real

set -euo pipefail
: "${BENCHMARK_SOURCE:?}"
: "${BENCHMARK_ROOT:?}"
eval "$(micromamba shell hook --shell bash)"
micromamba activate omnomnomics
cd "$BENCHMARK_SOURCE"
export PYTHONPATH="$BENCHMARK_SOURCE/src"
export NUMBA_CACHE_DIR="$TMPDIR/chip_real_numba_${SLURM_JOB_ID}"
export MPLCONFIGDIR="$TMPDIR/chip_real_mpl_${SLURM_JOB_ID}"
export XDG_CACHE_HOME="$TMPDIR/chip_real_xdg_${SLURM_JOB_ID}"
export BENCHMARK_PARALLEL_IDR=1
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
git rev-parse HEAD
python scripts/validation/chip_input_real_validation.py --root "$BENCHMARK_ROOT" --threads "$SLURM_CPUS_PER_TASK"
