# Metadata and sample identity

Metadata provides stable sample identity independently of FASTQ suffixes and
drives grouping, color categories, and differential designs.

The first column must be named `filename` and map to input filenames after
normalization. Additional columns describe biological and technical structure.

```text
filename	genotype	stim	donor	replicate
NT_C_D1_A_S50_L003	NT	C	D1	A
NT_L_D1_A_S58_L003	NT	L	D1	A
KD24_C_D1_A_S54_L003	KD24	C	D1	A
KD24_L_D1_A_S62_L003	KD24	L	D1	A
```

Selectors accept comma-separated column names or one-based indices:

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  -m metadata.tsv \
  --sample-name genotype,stim,donor,replicate \
  --sample-type genotype,stim \
  --sample-color genotype
```

- `sample_id` must uniquely identify each biological sample.
- `sample_type` defines analysis and track-hub groups.
- `sample_color` defines palette categories for non-stranded tracks.

Technical lanes are processed independently and merged at BAM level. Do not
remove lane identity from `filename`; use `sample_id` derivation to identify the
biological sample to which lanes belong.

For stranded RNA hubs, plus tracks are red and minus tracks are blue regardless
of `sample_color`.

## ChIP input libraries

On `input_dev`, ChIP and input libraries can share one metadata file. Add a
`role` column containing `chip` or `input`, and choose unique library identities
with `--sample-name`. The same selectors are used during preprocessing.

```text
filename	role	library	cell_type	condition	replicate	read_layout
chip_1	chip	chip_1	macrophage	untreated	1	PE
chip_2	chip	chip_2	macrophage	untreated	2	PE
input_1	input	input_1	macrophage	untreated	1	PE
facility_input.bam	input	input_2	macrophage	untreated	2	SE
```

```bash
omnomnomics chip -i EXPERIMENT -g GRCh38 -j 1-7 \
  -m metadata.tsv --sample-name library \
  --input-match cell_type,condition,replicate
```

`--input-match` accepts column names or one-based indices. All selected values
must match exactly (apart from surrounding whitespace):

| Selector | Association |
| --- | --- |
| `cell_type` | share inputs within a cell type |
| `cell_type,condition` | match each condition within a cell type |
| `cell_type,condition,replicate` | match individual replicates |
| `input_group` | match explicit user-defined sharing groups |

A constant `input_group` value provides a single shared-input group. Matching
several inputs records a set of distinct library IDs. MACS pools their fragment
files internally; shared controls are included once per ChIP pool. Technical
units sharing a library ID must follow the existing
`technical_replicate` convention and agree on role, matching fields and layout.
A missing match is an error. Unassigned inputs are listed and prepared separately.
Input-only preparation is allowed without `--input-match`. Metadata roles cannot
be combined with legacy `-I`.

Put raw reads in `FASTQ/` and supplied aligned BAMs in `BAM/`. Either role can
start from either source. Use an explicit `.bam` filename to select a supplied
BAM; use a FASTQ filename when unknown BAM and FASTQ sources compete for a bare
filename root. Existing pipeline outputs are tracked in the manifest. Filenames
may use `.fastq`, `.fq`, compressed variants, lane suffixes and `_R1_001` /
`_R2_001`. Paired reads must be complete. An isolated R1 requires an explicit
`read_layout=SE`; otherwise it could be a missing mate. Unsuffixed reads are
single-end. Optional `read_layout` values `SE`/`PE` are checked against discovered
FASTQs or BAM flags. Mixed layouts are supported across libraries, not within a
library's technical units.

BAMs are checked for readability, basic integrity, contig names/lengths and a
consistent layout among the first 10,000 alignments inspected. This bounded
check cannot establish the BAM's full processing history or prove its assembly
provenance. The reference must have `fasta/genome.fa.fai` or the normalized
chromosome-size file. Supplied BAMs enter before normal ChIP filtering; record
any filtering already performed by the sequencing facility.

For BAM-only entry, use `-j 4-7` (including technical merging), or `-j 5-7` for
existing per-library BAMs. Existing filtered BAMs in `filtered_BAM/` support
`-j 6-7`; noncanonical names are copied to the derived library name. A missing
intermediate requires explicitly selecting its producing step.

`run_configs/library_manifest.json` records sources, layout, roles, work and
outputs. `input_associations.tsv` records the resolved matches. Hash-named
snapshots remain available for the run configuration; a separate
`metadata_experiments.tsv` run artifact contains only ChIP rows. Original sources
remain protected across restarts and retention policies.

Analytical stages use these associations for every peak call, joint discovery,
input-relative tracks and regional enrichment on peaks, bins or gene bodies.
Input libraries remain outside the raw ChIP DE matrix. After preprocessing,
select public `-j 10,13-14` for group peaks, joint testing regions, counting and
DE, or supply `--chip-regions-bed atlas.bed` with `-j 13-14` for predefined
coordinates. DE-only `-j 14` can reuse archived counts without BAMs; changed
region/control settings require recounting. See the [ChIP guide](../assays/chip.md)
for scaling, coordinate rules and normalization limitations. Input-only runs
support preparation only; optional HOMER export is not enabled for role metadata.
Metadata-free processing and legacy `-I` remain available.
