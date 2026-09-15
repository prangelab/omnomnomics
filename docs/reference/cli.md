# CLI reference

The installed command is the authoritative reference:

```bash
omnomnomics --help
omnomnomics rna --help
omnomnomics atac --help
omnomnomics chip --help
omnomnomics monitor --help
omnomnomics genomes --help
```

## Common workflow options

| Option | Purpose |
| --- | --- |
| `-i`, `--experiment-dir` | experiment directory |
| `-g`, `--genome` | normalized assembly name |
| `-j`, `--mode` | `auto`, `all`, a range, or a comma-separated selection |
| `-m`, `--metadata` | tabular metadata whose first column is `filename` |
| `-T`, `--trim-tool` | `skewer` or `fastp`; default is Skewer |
| `-M`, `--map-tool` | RNA mapper: HISAT2, STAR, or the lab-specific, multimapping-permissive STAR-TE preset |
| `--dry-run` | validate and construct the DAG without submission |
| `--rerun-selected-steps` | delete and recompute selected stage outputs |
| `--site-config` | override the discovered user site configuration |
| `--retention-policy` | `all`, `pruned`, or `minimal` |
| `--max-project-size` | soft cap such as `300G` or `800GB` |
| `--keep-duplicates` | override assay default and retain duplicates |
| `--remove-duplicates` | override assay default and remove duplicates |
| `-X`, `--no-multiqc` | disable final MultiQC aggregation |

## Metadata and DE options

| Option | Purpose |
| --- | --- |
| `--sample-name` | columns used to derive unique sample IDs; required whenever metadata is supplied |
| `--input-match` | ChIP input matching columns; requires role metadata, used in preparation and analytical control resolution |
| `-I`, `--input` | legacy single shared ChIP control BAM; cannot be combined with metadata roles |
| `--chip-regions-bed` | protected custom non-overlapping ChIP testing BED; bypasses joint discovery |
| `--chip-joint-peak-q` | permissive joint testing-region MACS q threshold; default 0.1 |
| `--chip-se-fragment-length` | inferred SE fragment extension for ChIP MACS calls; default 200 bp |
| `--chip-effective-genome-size` | effective size used consistently by ChIP MACS calls; default approximate total reference length |
| `--chip-input-tracks` | `fold_enrichment` (default), `log2_ratio`, `qpois`, or `none`; visualization only |
| `--sample-type` | columns used to group peaks and hubs |
| `--sample-color` | columns used to derive palette categories |
| `--de-columns` | biological variables in an automatic DE design |
| `--de-block` | blocking variables in an automatic DE design |
| `--de-interactions` | include the interaction of two DE columns |
| `--de-formula` | explicit formula overriding automatic design options |
| `--de-config` | repeatable DE analysis YAML |
| `--de-out-dir` | output subtree for one DE analysis |

See the assay pages for ATAC/ChIP peak and broad-mode options.
