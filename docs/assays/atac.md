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

The public workflow continues after common step 9 with peak calling, peak QC,
pre-DE interpretation, counting, differential chromatin analysis, and post-DE
interpretation (steps 10-15).

Important controls include:

- `--narrow-peak-strategy idr|macs3`
- `--idr-mode basic|encode`
- `--idr-pairing-policy all_pairs|anchor_vs_all`
- `--atac-peak-opt-mode none|fast|full`
- `--spp-gate none|warn|drop|strict`

Install the IDR helper for the default strategy and the SPP helper when NSC/RSC
metrics or SPP gate decisions are required.
