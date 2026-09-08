# Troubleshooting

## Stale Snakemake lock

Confirm that no controller or worker process still targets the experiment. Then
unlock using the exact Snakefile, profile, and run configuration recorded for
the run:

```bash
snakemake --unlock \
  --snakefile /path/to/Snakefile.smk \
  --directory EXPERIMENT \
  --profile /path/to/slurm_profile \
  --configfile EXPERIMENT/run_configs/RUN.config.yaml
```

Do not remove `.snakemake/` while jobs are active.

## Quota or filesystem failures

`Disk quota exceeded`, `No space left on device`, and write failures from tools
such as `samtools sort` can all indicate exhausted project or scratch space.
Inspect the filesystem quota and largest directories, cancel remaining jobs, and
free space before resuming.

```bash
du -sh EXPERIMENT/* | sort -h
```

Use `--retention-policy pruned|minimal` and `--max-project-size` on subsequent
runs. The size guard cannot recover a filesystem that is already full enough to
prevent its own bookkeeping writes.

## Scheduler submission failures

Transient `sbatch` errors may be scheduler-side rather than pipeline failures.
Check the project compute budget and site submission pacing. On Snellius:

```bash
budget-overview -p rome
```

## Missing HOMER genome

Optional HOMER export requires separately installed HOMER genome data:

```bash
configureHomer.pl -install hg38
```

Common aliases are `GRCh38`/`GRCh38.p14` to `hg38` and `GRCm39` to `mm39`.

## Finding the actual error

The final controller lines often only report that a worker failed. Find the rule
and external job ID in `slurm_logs/controller/`, then inspect
`slurm_logs/<rule>/<rule>.<jobid>.out`.
