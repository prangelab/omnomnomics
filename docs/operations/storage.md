# Storage and retention

Sequencing depth can make raw and intermediate data substantially larger than a
sample count alone suggests. Omnomnomics provides explicit retention policies
and a soft project-size guard.

## Retention policies

| Policy | Behavior |
| --- | --- |
| `all` | keep all outputs and intermediates |
| `pruned` | keep source inputs and reusable downstream products; remove eligible bulky intermediates after successful use |
| `minimal` | keep source inputs and terminal outputs required by the requested branches |

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  --retention-policy pruned
```

Cleanup is stage-aware. `trimmed_FASTQ/` is not eligible until lane merging has
completed, and `BAM/` is not eligible until filtered BAM creation has completed.
These completion boundaries do not make source inputs eligible for deletion.
Lightweight trim and mapper metrics are copied to `run_logs/flow_qc_cache/`
before their source directories are removed, preserving later reporting.

## Source input protection

Before cleanup, omnomnomics records the earliest available sequence file type
for each filename identity: FASTQ, trimmed FASTQ, BAM, filtered BAM, or BigWig.
For example, when a sequencing facility supplies BAMs and no earlier reads are
available for those samples, their BAMs are source inputs and are retained.
FASTQs for another sample do not remove that protection. Table/feature-only
projects retain their earliest available input files instead.

Protection is stored in `run_configs/source_protection.json` and included in each
run configuration. Existing protection survives later runs, even if additional
upstream files are added. BAM indexes, explicitly supplied ChIP control BAMs and
the metadata file are protected too. Renamed sample identities are treated
conservatively; the pipeline does not assume differently named files are
regenerable from one another.

Retention cleanup, size-pressure cleanup, lane-BAM cleanup and forced-output
cleanup all respect this protection. A directory containing protected sources is
retained as a whole, so it can also retain generated files stored alongside the
sources. The run log records skipped deletions. Original sources take precedence
over the soft size cap. Keep the protection manifest when moving or archiving a
project; old run configurations without a source snapshot skip these cleanup
operations rather than infer source ownership after processing has started.

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

For metadata-declared ChIP preparation, the library manifest also records known
generated paths across selective restarts. This prevents a renamed intermediate
from being mistaken for a new supplied source. Existing protections are never
removed. Preflight rejects a planned output that would overwrite a protected
source. An existing index beside a source filtered BAM is reused as a protected
input; missing indexes can be generated.

Custom ChIP testing BEDs are protected too. ChIP counting archives canonical
testing coordinates and source provenance next to the matrices; these remain
terminal products for DE-only runs. This permits later DE entry without
generated BAMs or the peak tree. Existing group-support annotations retain
their recorded provenance; fresh annotation marks missing catalogues unassessed.
Generated MACS fragment files currently follow retention of their containing
peak tree, which is conservative and can retain substantial intermediates.

On this route, cleanup waits for all selected consumers of the folder, including
alignment QC of pre-filter BAMs. Failed consumers do not qualify as completion.
Completion markers and generated lane BAMs remain reusable under `all`; eligible
intermediates can still be removed under `pruned`, `minimal` or size pressure.
