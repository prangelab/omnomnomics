# ChIP-seq

ChIP runs remove duplicates by default and supports feature definitions matched
to distinct signal geometries.

| `--broad-mode` | Intended signal | Feature strategy |
| --- | --- | --- |
| `off` | transcription factors and narrow marks | narrow MACS3 plus IDR by default |
| `domain` | broad domains | MACS3 broad peaks with replicate support |
| `genebody` | gene-body-associated marks | annotation-derived gene bodies |
| `diffuse` | very diffuse chromatin signal | fixed genomic bins |

```bash
# Narrow / TF-like
omnomnomics chip -i EXPERIMENT -g GRCh38 -m metadata.tsv --broad-mode off

# Broad domains
omnomnomics chip -i EXPERIMENT -g GRCh38 -m metadata.tsv --broad-mode domain

# Gene-body signal
omnomnomics chip -i EXPERIMENT -g GRCh38 -m metadata.tsv --broad-mode genebody

# Diffuse signal
omnomnomics chip -i EXPERIMENT -g GRCh38 -m metadata.tsv --broad-mode diffuse
```

Use `-I` to provide a matching input/control BAM where appropriate. Domain mode
exposes MACS3 broad cutoff and replicate-overlap controls. Diffuse mode exposes
bin-size and merge-gap controls.

The public ChIP stage numbering matches ATAC: peak/feature construction at step
10 through post-DE interpretation at step 15.
