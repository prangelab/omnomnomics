# Storage and retention

Sequencing depth can make raw and intermediate data substantially larger than a
sample count alone suggests. Omnomnomics provides explicit retention policies
and a soft project-size guard.

## Retention policies

| Policy | Behavior |
| --- | --- |
| `all` | keep all outputs and intermediates |
| `pruned` | keep FASTQ and reusable downstream products; remove bulky upstream intermediates after successful use |
| `minimal` | keep FASTQ and only terminal outputs required by the requested branches |

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  --retention-policy pruned
```

Cleanup is stage-aware. `trimmed_FASTQ/` is not eligible until lane merging has
completed, and `BAM/` is not eligible until filtered BAM creation has completed.
Lightweight trim and mapper metrics are copied to `run_logs/flow_qc_cache/`
before their source directories are removed, preserving later reporting.

## Soft size guard

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  --retention-policy minimal \
  --max-project-size 300G
```

The guard measures project usage at storage checkpoints and attempts only safe
cleanup. If essential FASTQs and filtered BAMs already exceed the cap, the run
records a warning. Space-heavy BigWig/track-hub branches may be skipped while
requested count-table work continues.

This is a soft workflow policy, not a reservation. It cannot prevent failures
caused by another project consuming the same filesystem quota or by temporary
job-local scratch exhaustion.
