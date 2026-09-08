# Quickstart

## Prepare the project

Place one assay in an experiment directory. Raw reads belong in `FASTQ/`:

```text
experiment/
└── FASTQ/
    ├── sample_A_L001_R1_001.fastq.gz
    ├── sample_A_L001_R2_001.fastq.gz
    ├── sample_B_L001_R1_001.fastq.gz
    └── sample_B_L001_R2_001.fastq.gz
```

Use one consistent field separator in filenames. Metadata is recommended for
any grouping or differential analysis.

## Validate

```bash
micromamba activate omnomnomics
omnomnomics rna -i /path/to/experiment -g GRCh38 --dry-run
```

Replace `rna` with `atac` or `chip` as appropriate. A dry run validates inputs,
resolves configuration, and builds the Snakemake DAG without submitting work.

## Run

```bash
omnomnomics rna \
  -i /path/to/experiment \
  -g GRCh38 \
  -j all \
  --retention-policy pruned \
  --max-project-size 300G
```

The command submits a controller job with `sbatch`. The controller then submits
the rule jobs, so do not submit the command through a second `sbatch` wrapper.

## Monitor

```bash
omnomnomics monitor -i /path/to/experiment
```

The monitor reports workflow step states and recent pipeline log lines. Rule-level
errors are written under `slurm_logs/<rule>/`.

## Resume or rerun

By default, Snakemake reuses valid outputs. Select a later stage range to resume:

```bash
omnomnomics rna -i /path/to/experiment -g GRCh38 -j 5-12
```

Use `--rerun-selected-steps` only when those selected outputs must be rebuilt.
It is not required after an ordinary interrupted run.
