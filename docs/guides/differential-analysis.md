# Differential analysis

RNA differential expression is public step 12. ATAC/ChIP differential chromatin
analysis is public step 14. Both use metadata-driven DESeq2 designs.

## Automatic design

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  -j 12 \
  -m metadata.tsv \
  --de-columns genotype,stim \
  --de-block donor \
  --de-interactions
```

Use `--de-formula` when the model must be specified directly. It overrides
`--de-columns` and `--de-block`.

## YAML configuration

`src/omnomnomics/workflow/config/de_config.example.yaml` documents filtering,
contrasts, shrinkage, QC, plotting, enrichment, and runtime controls.

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  -j 12 \
  -m metadata.tsv \
  --de-config de_grouped.yaml
```

Repeat `--de-config` to run several analyses sequentially. Each file must set a
unique `io.out_dir`; do not combine repeated config files with global
`--de-out-dir`.

The workflow writes resolved metadata and rendered R scripts so the statistical
analysis can be audited and customized later.

## Chromatin interpretation

ATAC/ChIP post-DE analysis writes set manifests, signal-run statuses, and motif-
run statuses. Signal plotting can be controlled with
`--post-de-signal-policy auto|require|skip`. Motif `TIMEOUT` and `FAIL` statuses
are report-level outcomes rather than silent omissions.
