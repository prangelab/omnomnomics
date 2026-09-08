# ATAC-seq

ATAC runs remove duplicates by default, filter at MAPQ 30, report mitochondrial
read statistics, and use replicate-aware IDR consensus for narrow peaks.

```bash
omnomnomics atac \
  -i EXPERIMENT \
  -g GRCh38 \
  -m metadata.tsv \
  --sample-name condition,donor,replicate \
  --sample-type condition
```

The workflow continues after common step 9 with peak calling, peak QC, pre-DE
interpretation, counting, differential chromatin analysis, and post-DE
interpretation (steps 10-15).

## Peak construction

`--narrow-peak-strategy idr` is the default. MACS3 calls each biological
replicate separately, and IDR derives reproducible group-level peaks from the
replicate calls. `--narrow-peak-strategy macs3` instead pools each group and can
choose MACS3 parameters using the optimizer below.

| Control | Behavior |
| --- | --- |
| `--idr-mode basic` | run pairwise IDR on the true biological replicates |
| `--idr-mode encode` | add pooled and per-replicate pseudoreplicate diagnostics to the true-replicate analysis |
| `--idr-pairing-policy all_pairs` | compare every replicate pair; most complete for groups with more than two replicates |
| `--idr-pairing-policy anchor_vs_all` | compare the first replicate with each remaining replicate; fewer jobs for large groups |

With the MACS3 strategy, `--atac-peak-opt-mode none` makes one fixed call at
q=0.01, shift=-100, and extension size=200. `fast` compares that call with a
shift=-75/extension=150 alternative. `full` evaluates six combinations spanning
q=0.01 or 0.05 and shift/extension pairs of -100/200, -75/150, and -50/100.

For each group, the optimizer favors high FRiP and agreement between replicate
peak sets, penalizes blacklist overlap, and weakly favors a median peak width
near 250 bp. It writes all candidate scores and the selected parameters under
`atac_peak_call_optimization/`.

## SPP quality gate

`--spp-gate none` disables decisions based on NSC/RSC. `warn` records threshold
failures without changing downstream samples. `drop` excludes flagged samples
from counting and differential analysis. `strict` stops the workflow if any
sample fails the configured thresholds.

Install the IDR helper for the default strategy and the SPP helper when NSC/RSC
metrics or SPP gate decisions are required.
