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
| `-M`, `--map-tool` | HISAT2, STAR, or STAR-TE where supported |
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
| `--sample-name` | columns used to derive unique sample IDs |
| `--sample-type` | columns used to group peaks and hubs |
| `--sample-color` | columns used to derive palette categories |
| `--de-columns` | biological variables in an automatic DE design |
| `--de-block` | blocking variables in an automatic DE design |
| `--de-interactions` | include the interaction of two DE columns |
| `--de-formula` | explicit formula overriding automatic design options |
| `--de-config` | repeatable DE analysis YAML |
| `--de-out-dir` | output subtree for one DE analysis |

See the assay pages for ATAC/ChIP peak and broad-mode options.
