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

Use `-I` to provide the matching input/control BAM where appropriate. ChIP
feature construction runs at step 10 and post-DE interpretation at step 15.

## Narrow and transcription-factor ChIP

`--narrow-peak-strategy idr` is the default and derives reproducible peaks from
separate biological-replicate calls. `--idr-mode encode` adds pseudoreplicate
diagnostics; `basic` uses true replicates only. For groups with more than two
replicates, `--idr-pairing-policy all_pairs` evaluates every pair, while
`anchor_vs_all` compares the first replicate with each remaining replicate.

The alternative `--narrow-peak-strategy macs3` enables
`--chip-peak-opt-mode`:

| Mode | Candidate calls per group |
| --- | --- |
| `none` | one model-based MACS3 call at q=0.01 |
| `fast` | model-based calls at q=0.01 and q=0.001 |
| `full` | the two model-based calls plus fixed shift/extension calls at -75/150 and -100/200 for both q-values |

For `fast` and `full`, omnomnomics calls each candidate on both pooled and
replicate BAMs. It selects the candidate using a documented composite score:
45% replicate peak-set agreement, 35% FRiP, 15% low blacklist overlap, and 5%
preference for a median peak width near 250 bp. Full candidate tables, selected
parameters, and diagnostic plots are retained under
`chip_narrow_peak_call_optimization/`. This optimization is useful when a
MACS3-only peak set is required; IDR remains the preferred default for
replicated TF ChIP-seq.

## Broad-domain ChIP

`--broad-mode domain` calls broad MACS3 domains and filters them by replicate
support. The principal controls are:

| Control | Purpose |
| --- | --- |
| `--chip-broad-qvalue` | relaxed q-value used for pooled, replicate, and pooled-pseudoreplicate MACS3 calls |
| `--chip-broad-cutoff` | MACS3 significance cutoff used to link nearby enriched regions into broad domains |
| `--chip-broad-min-length` | optional minimum domain length passed to MACS3 |
| `--chip-broad-max-gap` | optional maximum gap that MACS3 may bridge within a domain |
| `--chip-broad-replicate-fraction` | minimum fraction of biological replicates that must support a pooled domain |
| `--chip-broad-overlap-fraction` | minimum reciprocal overlap used to count replicate support |

The defaults require support from every replicate and at least 50% reciprocal
overlap. Relaxing these values increases sensitivity but also admits less
reproducible domains.

## Gene-body and diffuse ChIP

`--broad-mode genebody` quantifies annotation-derived gene bodies instead of
calling peaks. `--broad-mode diffuse` partitions the genome into fixed bins,
with bin-size and merge-gap controls for signals that do not form discrete
domains. These modes change the feature definition used for counting and
differential analysis; they are not merely plotting choices.

The SPP gate has the same behavior as for ATAC: `none` disables gating, `warn`
reports flags, `drop` removes flagged samples from downstream analysis, and
`strict` aborts on any threshold failure.
