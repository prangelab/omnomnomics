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

Use metadata roles for shared or separately matched inputs, or retain `-I` for
one shared control BAM. ChIP feature construction runs at public step 10,
counting at 13, DE at 14 and post-DE interpretation at 15.

## Metadata-declared inputs

On `input_dev`, declare `role=chip` and `role=input` libraries together in the
normal metadata file. `--input-match` selects the fields used to associate them.
FASTQ and supplied BAM libraries can be mixed, including SE/PE layouts.
Matched controls now feed peak calling, joint discovery, regional enrichment
and input-relative tracks. See the [metadata guide](../guides/metadata.md#chip-input-libraries)
for shared-input, replicate-matched and BAM-entry examples.

Every biological-replicate call uses its matched controls. A pooled ChIP group
uses the union of distinct control libraries; shared controls contribute once.
Multiple inputs contribute all reads without depth equalization. Pseudoreplicates
use their parent's full controls. MACS receives three-column BEDPE fragments:
observed proper-pair inserts for PE and inferred extensions for SE, controlled
by `--chip-se-fragment-length` (default 200 bp). This changes the former SE
model/shift candidate behavior. The inferred length needs to suit the library.
Legacy `-I` supplies a prepared control BAM; it does not preprocess that control.
Use compatible filtering and duplicate handling for the control and ChIP libraries.
Optional HOMER tag-directory mode is not enabled for role-based inputs.
All ChIP MACS calls use the same declared genome-size convention.
`--chip-effective-genome-size` overrides the default total reference length,
which is an approximation and includes bases that may not be uniquely mappable.
Use an appropriate effective size when available. This also removes the former
implicit human-genome default from group calls on other organisms.

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
| `none` | one fragment-based MACS3 call at q=0.01 |
| `fast` | fragment-based calls at q=0.01 and q=0.001 |
| `full` | the same two calls; redundant shift variants are omitted for fixed fragment coordinates |

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

## Differential testing regions

Peak modes use a separate permissive joint MACS call over all ChIP libraries
in the run, with the unique union of their inputs. `--chip-joint-peak-q` defaults
to 0.1; in domain mode the joint broad cutoff also uses this value. Run one
compatible mark/antibody and genome per discovery cohort. All reads contribute
to discovery, so deeper libraries and conditions with more replicates contribute
more. Discovery pooling defines coordinates; biological replicates are counted
separately for inference.

Existing condition/group peak catalogues remain available. Their overlap is
annotation only: any positive overlap is reported, with overlap bases/fraction,
support categories and available IDR selection/fallback status. Missing catalogues
are unassessed. Overlap does not admit or remove testing regions.

Use `--chip-regions-bed atlas.bed` to bypass joint discovery. BEDs must be
tab-separated, nonempty, zero-based half-open, non-overlapping and free of
duplicate coordinates, with valid reference contigs/bounds. Original names and
checksums are archived; the supplied file is protected. This option conflicts
with gene-body/bin selectors, which already define coordinates. A supplied BED
is independent only if its source is independent of the measurements being
tested. BED-only counting and DE do not require peak calling or IDR.

## Counts and enrichment

The raw DE matrix contains ChIP libraries only. Separate input counts, mapped
fragment depths and regional enrichment share the frozen coordinates. Counts
use unique feature assignment: SE aligned reads and PE proper same-contig
fragments. Reads ambiguous between overlapping annotation-derived gene bodies
are excluded. Regional enrichment sums input counts/depths over distinct matched
libraries, then divides ChIP CPM by pooled input CPM. No pseudocount is added;
input counts below the site-configured `chip_enrichment_min_input` (default 5)
are labelled `low_input` with undefined enrichment. Zero ChIP counts have fold
enrichment zero and undefined log2 enrichment. No-input results are labelled
`no_input`. These ratios are not ChIP-qPCR percent input.

DESeq2 uses raw ChIP counts without input subtraction/division and requires
residual biological replication. The initial filter uses total counts across
ChIP samples, rather than requiring a peak/count in both conditions. Reports
record the filter and normalization. Results describe relative abundance under
that normalization; input-aware discovery does not remove background from raw
counts or recover absolute/global ChIP changes. Support subsets retain the
adjusted p-values of the full tested family.

Counts archive coordinates and provenance, enabling DE-only restarts after
generated BAM/peak cleanup or moving the project. Project-local source paths
are compared relative to the project root; external source paths retain their
identities. Changed region/control, regional-enrichment or QC sample-selection
settings require public step 13 to refresh counts/summaries. SPP exclusions
apply only with `--spp-gate drop`; other modes ignore an existing drop list.
The automatic DE design is selected from ChIP metadata alone, excluding input
libraries. Per-library count caches and DE model caches reuse unchanged work;
`--rerun-selected-steps` bypasses those caches.

## Input-relative tracks

Ordinary ChIP CPM coverage remains available. `--chip-input-tracks` selects:

| Mode | Output scale |
| --- | --- |
| `fold_enrichment` (default with inputs) | MACS fragment pileup divided by depth-scaled modelled local background |
| `log2_ratio` | direct 10-bp CPM read-coverage log2 ChIP/input ratio, pseudocount 1 CPM |
| `qpois` | MACS local-background significance, −log10(q), using unscaled fragment pileups |
| `none` | no additional input-relative tracks |

Tracks and JSON command/scale provenance are stored in `BigWigs/input_relative`.
Track hubs place them in separately labelled containers. Read-coverage ratios,
modelled-background FE and significance are different quantities. Average
profiles remain descriptive; enrichment tracks are not substituted into DE
counts, and loci/bins are not biological replicates for treatment significance.
With a project-size cap, these optional additional tracks are omitted; ordinary
track size guarding remains in force.
