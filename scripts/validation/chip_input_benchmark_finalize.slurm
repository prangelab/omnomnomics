#!/bin/bash
#SBATCH --partition=rome
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=00:30:00
#SBATCH --job-name=chip-input-report

set -euo pipefail
: "${BENCHMARK_SOURCE:?}"
: "${BENCHMARK_ROOT:?}"
eval "$(micromamba shell hook --shell bash)"
micromamba activate omnomnomics
cd "$BENCHMARK_SOURCE"
export PYTHONPATH="$BENCHMARK_SOURCE/src:$BENCHMARK_SOURCE/scripts/validation"
export NUMBA_CACHE_DIR="$TMPDIR/chip_report_numba_${SLURM_JOB_ID}"
export MPLCONFIGDIR="$TMPDIR/chip_report_mpl_${SLURM_JOB_ID}"
export XDG_CACHE_HOME="$TMPDIR/chip_report_xdg_${SLURM_JOB_ID}"
mkdir -p "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR" "$XDG_CACHE_HOME"
git rev-parse HEAD
python - <<'PY'
import os
from pathlib import Path
import shutil
import subprocess
from chip_input_benchmark import LENGTH, SCENARIOS
from omnomnomics.chip_analysis import read_regions, support_annotations, write_tsv

root = Path(os.environ["BENCHMARK_ROOT"])
for scenario in SCENARIOS:
    for seed in range(4):
        if not (root / f"{scenario}_seed{seed}" / "complete.json").exists():
            raise RuntimeError(f"Incomplete realization: {scenario} seed {seed}")
smoke = root / "B1_seed0"
catalogues = {condition: smoke / "J/peak_calling" / f"{condition}.MACS3.optimized.bed" for condition in ("WT", "KO")}
for arm in ("J", "U", "A"):
    project, export = smoke / arm, smoke / "metrics" / arm
    cache, = (project / "DE_calling").rglob("chip_model_cache.rds")
    subprocess.run(["Rscript", "scripts/validation/chip_input_benchmark_export.R", str(cache), str(export), "export"], check=True)
    regions = read_regions(project / "DE_calling/chip_testing_regions.bed", {"chr1": LENGTH})
    header, rows = support_annotations(regions, catalogues)
    write_tsv(export / "discovery_support.tsv", header, rows)
    shutil.copy2(smoke / "J/peak_calling/chip_idr_peak_calling/idr_selected_peaks.tsv", export / "group_idr_selected_peaks.tsv")
PY
python scripts/validation/chip_input_benchmark_profiles.py --root "$BENCHMARK_ROOT" --output "$BENCHMARK_ROOT/../profiles"
python scripts/validation/chip_input_benchmark_report.py --root "$BENCHMARK_ROOT" --output "$BENCHMARK_ROOT/../report"
