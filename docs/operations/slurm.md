# Slurm execution

The CLI is a submitter. A normal workflow command performs preflight validation,
writes a resolved run configuration and controller script, and submits that
controller with `sbatch`. The controller runs Snakemake and dispatches worker
jobs through the packaged Slurm profile.

```text
login shell
└── omnomnomics ...
    └── controller Slurm job
        └── Snakemake worker jobs
```

Do not wrap the CLI invocation in another `sbatch` command.

Site configuration defines controller resources, worker defaults, partition,
maximum concurrent jobs, and scheduler pacing. Conservative defaults for
`max_jobs_per_second` and `max_status_checks_per_second` reduce transient Slurm
submission failures on busy clusters.

## Inspect jobs

```bash
squeue -u "$USER"
scontrol show job JOB_ID
```

## Cancel a workflow

Cancel the controller and its worker jobs before modifying or removing active
outputs. Job names and IDs are visible in `squeue`; rule logs are preserved in
`slurm_logs/`.

After an interrupted run, make sure no jobs remain before clearing a stale
Snakemake lock. See [Troubleshooting](troubleshooting.md).
