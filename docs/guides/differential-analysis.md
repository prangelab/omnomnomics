# Differential analysis

RNA differential expression runs at step 12. ATAC/ChIP differential chromatin
analysis runs at step 14. Both use metadata-driven DESeq2 designs.

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

Use a YAML file when an analysis needs more control than the design options on
the command line. A compact genotype-by-stimulation configuration can look like
this:

```yaml
version: 1
io:
  out_dir: "genotype_by_stim"

design:
  formula: "~ donor + genotype * stim"
  reference_levels:
    genotype: "NT"
    stim: "C"

contrasts:
  mode: "explicit"
  explicit:
    items:
      - [stim, L, C]

filtering:
  enabled: true
  method: "min_count_samples"
  min_count: 10
  min_samples: 2

deseq2:
  lfc_shrink:
    enabled: true
    type: "apeglm"
    use_for_tables: true

thresholds:
  alpha: 0.05
  lfc_for_sig: 1.0

enrichment:
  enabled: true
```

`design.formula` defines the DESeq2 model, while `reference_levels` fixes the
baseline used to interpret coefficients. Explicit contrasts name the comparison
to report. In this interaction model, `[stim, L, C]` is the stimulation effect
at the reference genotype (`NT`); genotype-specific response differences require
an additional interaction-coefficient contrast. Filtering removes weakly
detected features before fitting;
`thresholds` controls significance classification and plotting rather than the
model formula. Each analysis writes to its own `io.out_dir`.

The
[complete example configuration](https://github.com/prangelab/omnomnomics/blob/main/src/omnomnomics/workflow/config/de_config.example.yaml)
also documents automatic contrasts, QC and plotting controls, latent factors,
clusterProfiler and decoupler enrichment, custom gene sets, and runtime policy.
Values omitted from a user configuration retain the packaged defaults.

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

After differential chromatin testing, omnomnomics creates BED sets for
significant, increased, and decreased regions. It can then add signal profiles,
heatmaps, and motif reports for those sets.

| `--post-de-signal-policy` | Behavior |
| --- | --- |
| `auto` | create signal plots when the required BigWigs are available |
| `require` | fail if the requested signal plots cannot be produced |
| `skip` | omit signal plotting while retaining the differential results |

Every signal and motif attempt is recorded in status tables. A motif status of
`TIMEOUT` or `FAIL` therefore identifies an optional report that did not finish;
the differential result tables and region sets remain available.
