# Monitoring and logs

Start the terminal monitor with:

```bash
omnomnomics monitor -i EXPERIMENT
```

The step table uses the assay step numbers accepted by `-j` and distinguishes running,
pending/waiting, completed, and failed work. A later step can run for one sample
while an earlier step is still completing other samples because Snakemake
schedules dependencies per output, not as global barriers.

## Provenance locations

| Location | Content |
| --- | --- |
| `run_logs/omnomnomics.run.*.log` | run settings, stage events, cleanup decisions, and final status |
| `run_logs/steps/` | step summaries, commands, and notes |
| `run_logs/*.tools.log` | tool versions and tool-provenance records |
| `run_configs/` | resolved run configuration and controller script |
| `slurm_logs/controller/` | Snakemake controller output |
| `slurm_logs/<rule>/` | rule-specific worker output and exceptions |
| `.snakemake/log/` | native Snakemake execution log |

When the monitor shows a failed step, inspect the corresponding worker log
first. Controller output identifies the failed rule and external Slurm job ID.
