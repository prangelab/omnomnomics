#!/usr/bin/env bash
#SBATCH --job-name=omnom_chr_smoke
#SBATCH --partition=rome
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=16
#SBATCH --output=/scratch-shared/kprange/omnomtest/chip_input_validation_orchestration/chromosome_smoke.%j.out
#SBATCH --error=/scratch-shared/kprange/omnomtest/chip_input_validation_orchestration/chromosome_smoke.%j.out

set -euo pipefail

: "${PINNED_COMMIT:?PINNED_COMMIT is required}"

cd "${HOME}/omnomnomics"
test "$(git rev-parse HEAD)" = "${PINNED_COMMIT}"
eval "$(micromamba shell hook --shell bash)"
micromamba activate omnomnomics
export PYTHONPATH="${HOME}/omnomnomics/src"

python - <<'PY'
import json
from pathlib import Path

from omnomnomics.chromosomes import read_reference_fai
from omnomnomics.genomes import resolve_blacklist_bed

assembly_root = Path('/scratch-shared/kprange/omnomtest/genomes/assemblies')
name = 'GRCh38.p14'
fai = assembly_root / name / 'fasta/genome.fa.fai'
reference = read_reference_fai(fai)
blacklist = resolve_blacklist_bed(name, assembly_root, download_if_missing=False)
provenance = json.loads(Path(f'{blacklist}.source.json').read_text())
with blacklist.open() as handle:
    lines = [line.rstrip('\n') for line in handle if line.strip()]
assert lines
assert all(line.split('\t', 1)[0] in reference for line in lines)
assert provenance['source_records'] == provenance['written_records'] + provenance['skipped_records']
assert provenance['chromosome_mapping'].get('chr1') == '1'
print(f"blacklist={blacklist}")
print(f"source_records={provenance['source_records']}")
print(f"written_records={provenance['written_records']}")
print(f"skipped_records={provenance['skipped_records']}")
print(f"chr1_to_1={provenance['chromosome_mapping'].get('chr1')}")
PY

run_root=/scratch-shared/kprange/omnomtest/chip_input_validation_d6ca74c90779
peaks="${run_root}/projects/chip_narrow_human_gm12878_vs_k562_ctcf/peak_calling/GM12878.MACS3.optimized.bed"
blacklist=/scratch-shared/kprange/omnomtest/genomes/assemblies/GRCh38.p14/aux/GRCh38.p14.blacklist.reference.bed
test -s "${peaks}"
test -s "${blacklist}"
printf 'source_peaks=%s\n' "$(wc -l < "${peaks}")"
printf 'blacklist_overlap_peaks=%s\n' "$(bedtools intersect -u -a "${peaks}" -b "${blacklist}" | wc -l)"
