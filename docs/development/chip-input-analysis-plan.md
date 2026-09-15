# ChIP input use and differential-analysis implementation plan

Date: 2026-09-11. Branch: `input_dev`.
Status: local analytical implementation completed on 2026-09-14; broader
statistical and visual validation remains an agreed review checkpoint.
The detailed contract below records the agreed design. See the
[implementation and validation record](chip-input-analysis-implementation.md)
for implemented choices, local checks and remaining statistical review.

## Starting point and scope

Input preparation through public steps 1–7 is implemented locally, including
metadata matching, FASTQ/BAM entry, mixed library layouts and source protection.
See [the preparation record](chip-input-implementation.md) and
[the preparation plan](chip-input-plan.md). Its reported 65 unittest tests and
45 standalone checks validate preparation, not the analytical changes below.
Existing branch changes remain uncommitted at the time of this record.

The preparation milestone rejected downstream analytical requests. This phase
enables the consumers with explicit dependencies and control use. Preserve
legacy `-I` as a single shared control and the existing
prohibition on combining it with metadata-declared inputs.

The manuscript ledger remains paused. No manuscript edits or validation reruns
are part of recording this plan. Later manuscript ChIP validations must use
inputs, after analytical support has been built and validated.

## Agreed analytical contract

### Control resolution and peak calling

- A single matched input is used directly. Multiple matched input libraries
  contribute all their reads to a pooled control; do not equalize their depths
  or downsample them. Preserve the individual libraries and their identities.
- A pooled ChIP group uses the union of distinct input IDs associated with its
  members. Count a shared input once, regardless of how many ChIP samples use it.
- Individual ChIP replicates use their own matched control sets. Pseudoreplicate
  calls use the full control set of their parent ChIP library or pooled group.
- Apply controls consistently to replicate, pooled, pseudoreplicate, narrow,
  broad and optimization candidate calls, covering split and monolithic paths.
  IDR operates on the resulting peaks; it does not itself consume input BAMs.
- Missing or incompatible declared inputs are errors. Never silently fall back
  to another control or to an uncontrolled call. Explicit no-input ChIP remains
  a supported, clearly recorded choice outside controlled validation runs.
- Handle caller SE/PE representations explicitly. Mixed-layout preparation
  support is not proof that a particular MACS invocation supports the same mix.

### Differential abundance and regional enrichment

- Count ChIP samples and each distinct input library on exactly the same region
  coordinates, using consistent fragment-counting conventions. Keep separate
  raw ChIP and input matrices plus the association table.
- Routine differential analysis uses raw ChIP fragment counts, biological
  replicates and the relevant experimental design. Inputs are not extra
  treatment replicates or automatic DESeq2 subtraction/division factors.
- Interpret results as relative ChIP abundance under the chosen between-ChIP
  normalization assumptions. Input-aware peak calls do not remove background
  from counts, and input correction does not solve global scaling by itself.
- Retain input counts for enrichment summaries and diagnostics, including
  background and possible copy-number differences. An input-relative
  differential model is future optional work, not a release prerequisite.
- Peak, gene-body and bin summaries report raw counts, depth-adjusted ChIP and
  input abundances, and input-relative enrichment. Aggregate counts across each
  region before forming enrichment; do not average unstable tiny-bin ratios.
- Define and record scaling, denominator handling and reliability flags. Use
  undefined/low-information indicators where background information is
  inadequate. Do not label sequencing ratios as ChIP-qPCR percent input.
- Keep regional enrichment separate from the raw-count DE input. Final numeric
  thresholds and normalization reference choices require implementation review
  and validation; they are not established by this planning record.

### Tracks and profiles

- Retain ordinary ChIP coverage/CPM tracks. Add modelled MACS-background fold
  enrichment as the principal enrichment track, with optional significance
  tracks and optional direct ChIP/input log-ratio output.
- Make input-relative visualization configurable and record its scale. Do not
  average or compare independently scaled tracks without a justified convention.
- Average profiles are descriptive by default. Do not add treatment-significance
  stars by treating peaks or adjacent profile bins as biological replicates.
  Any later inferential profile mode needs an explicit question and a model
  respecting biological replication. Locus-resampling bands must be labelled
  as such and not presented as uncertainty across biological experiments.

## Region selection for differential testing

### Default: permissive joint discovery

Keep condition-specific peak catalogues and IDR/reproducibility outputs. Add a
separate candidate-region catalogue for differential testing, discovered jointly
across the relevant ChIP samples without selecting on the treatment contrast.
Pool only compatible measurements of the same target/mark and genome within the
intended analysis scope; do not combine unrelated antibodies or assays.

Start with a configurable permissive joint-discovery threshold of approximately
MACS q <= 0.1. This is a candidate discovery setting, not the differential FDR
threshold or a guarantee of exhaustive sensitivity. Verify actual MACS narrow
and broad parameter semantics before exposing settings; do not conflate MACS
q-values, p-values and IDR thresholds.

Use the union of distinct matched controls for the discovery pool. Pooling is
for coordinate discovery only: count each biological replicate separately and
retain replicate-level inference. Do not require a peak in both conditions.

Apply any additional background/low-information filtering using suitable
aggregate abundance or enrichment across samples, without using the tested
condition contrast. Document the filter and assess its null behavior. Merely
omitting condition labels does not prove independence from the test statistic.
Freeze the final testing coordinates and assign stable region IDs before DE.

### Existing peaks provide annotation, not eligibility

Do not intersect joint candidates with condition-specific peak sets as a
mandatory admission rule. Such a gate would restore condition-dependent
selection even if joint coordinates were retained.

For each candidate, report overlap with each existing group peak catalogue,
the number/list of supported groups, and whether it has joint-discovery-only
support. For two-group contrasts, provide readable categories: group A only,
group B only, both, or joint-only. Preserve group-specific reproducibility
status, including pooled fallbacks, instead of calling every overlap IDR support.

Keep candidate coordinates unchanged by annotation. Record overlap bases and
fractions so users can distinguish marginal from substantial overlap. The exact
display classification threshold remains to be specified and logged; it must
not silently become a DE filter. Missing catalogues mean support not assessed,
not no overlap.

Calculate differential adjusted p-values over the full testing family retained
after the documented filters, separately as appropriate for each contrast.
Support-based table/plot subsets retain those values. Do not recompute FDR on a
post-selected support subset or imply that its subset FDR is newly guaranteed.

### Alternative: user-supplied BED

Add a custom region BED option that bypasses joint discovery for testing.
Proposed CLI names, subject to existing interface conventions:
`--chip-regions-bed` and `--chip-joint-peak-q` (initial default 0.1).
Both options are now implemented under these names.

- Validate readability, nonempty BED records, zero-based half-open integer
  coordinates, positive widths, known contigs and coordinate bounds against
  the selected reference. Header compatibility cannot establish assembly
  provenance; record the declared genome and supplied path.
- Preserve the original BED byte-for-byte as a protected user source. Record
  its checksum and archive a lightweight canonical region manifest with stable
  IDs and original names where present.
- Explicitly resolve duplicate and overlapping intervals and their fragment
  assignment behavior before implementation. Do not silently merge, clip or
  rename biological features without an auditable mapping. Prefer validation
  errors for unsupported ambiguity over arbitrary changes.
- A supplied BED may come from an independent atlas or annotation, but file
  import does not establish statistical independence. Document that a BED
  derived from these same condition-specific ChIP measurements retains its
  selection history.
- Permit BED-only differential entry without requiring fresh peak calling or
  IDR outputs. Annotate existing peak support when available; otherwise mark
  it unassessed. No caller dependencies merely to populate optional annotations.
- Make BED, gene-body and fixed-bin region modes explicit and unambiguous.
  Reject conflicting selectors. Gene-body and fixed-bin modes retain their
  predefined coordinates and do not require the joint discovery call.

## Implementation map and order

Internal rule filenames below do not redefine public step numbering.

| Unit | Main locations | Deliverable |
| --- | --- | --- |
| 1. Review existing state and contracts | `cli.py`, `chip_inputs.py`, `library_manifest.py`, `workflow/Snakefile.smk`, `workflow/config/workflow.yaml` | Reconcile the existing preparation manifest, downstream guards, sample groups, public steps and normalization code with this plan. Define analysis scope, region modes and CLI validation. |
| 2. Control-set products | Shared helpers and dedicated workflow rules | Deterministic unique-ID control sets, reusable pooled products, explicit BAM/index/association dependencies, complete provenance and correct layout handling. |
| 3. Existing peak consumers | `workflow/rules/10.call_peaks.smk`, split IDR jobs and their companion rules | Correct controls on every applicable call while preserving group peak outputs and reproducibility behavior. |
| 4. Testing-region construction | New focused region helper/rules; rule 10; CLI/config | Permissive joint discovery and custom BED routes, stable coordinate manifest, explicit discovery settings and input identities. Keep the group-union catalogue distinct from the DE catalogue. |
| 5. Support and filtering | Region helper, `peak_annotation.py`, rules `13.peak_qc.smk` and `14.analyze_peaks.smk` as applicable | Aggregate eligibility filters with recorded reasons; per-group overlap annotations that cannot accidentally select the testing set. |
| 6. Counts and enrichment | `workflow/rules/11.count_reads.smk`, focused summary helpers | Separate raw ChIP/input counts aligned by stable region ID, associations, regional enrichment and reliability indicators for peaks, gene bodies and bins. |
| 7. Differential results and explorer | `15.call_DE_chrom.smk`, `templates/de_core_chrom.R.tmpl`, shared R helpers and Shiny modules | Raw ChIP DE only, audited normalization/design, correct testing-family adjustment, retained support annotations and descriptive profile behavior. Audit `16.analyze_peaks_de.smk` for downstream assumptions. |
| 8. Tracks and reporting | `8.create_wiggles.smk`, `9.merge_wiggles.smk`, profile/QC consumers | Clearly labelled coverage and input-relative tracks with role-safe grouping and documented scaling. |
| 9. Workflow lifecycle | CLI, Snakefile targets/monitoring, `storage.py`, source manifests | Selective entry, unchanged-run reuse, changed-region/control invalidation, cleanup after successful consumption and source protection. Enable new-route downstream stages only when their full requested path is supported. |
| 10. Documentation and acceptance | Assay/CLI/metadata/storage guides and tests | User examples, normalization limitations, source provenance, meaningful functional/regression checks and a reviewable validation report. |

Current source anchors: `10.call_peaks.smk:2312` merges final split-IDR group
BEDs, and `11.count_reads.smk:368` selects the shared ChIP counting BED, with a
QC-filtered fallback. Audit the monolithic merge paths too. Replacing only one
filename would leave inconsistent consumers. Trace annotations, feature IDs,
counts, profiles and DE together before changing interfaces.

Keep lightweight semantic dependency manifests stable across unchanged runs.
The existing mtime-based rerun policy must detect changes in BED content,
discovery threshold, sample membership and control mapping through actual
dependencies; a path present only in parameters is insufficient. Unchanged
support annotations should not rerun peak calling or statistical fitting.

Retain original FASTQs, imported BAMs/indexes and custom BEDs. Generated pools
are eligible for cleanup only after all successful consumers, including track
and discovery jobs, no longer need them. Keep checksums, commands, versions,
association tables, region manifests, filter summaries and normalization records
for long-term traceability. Do not add dependencies or environment changes
incidentally; use the established managed environment and review any actual
dependency need under project policy.

## Validation and acceptance

### Functional and workflow checks

- Single, shared, per-condition, per-replicate and multiple inputs; duplicate
  shared controls included exactly once in every relevant pool. Legacy `-I`
  and deliberate no-input paths retain their explicit semantics.
- Real small MACS calls with known controls: replicate/group/pseudoreplicate
  paths, optimization candidates, narrow and broad behavior, and each supported
  library layout. Missing/incompatible controls fail before substantive work.
- Joint discovery includes a strong condition-specific region even where the
  opposite condition has no called peak. Joint-only regions survive unless
  removed by a documented aggregate filter, regardless of support annotations.
- Changing group labels alone within a fixed discovery cohort must not alter
  its candidate coordinates or aggregate eligibility. Null-calibration tests
  are still required; this invariance check is necessary but insufficient.
- Custom BED validation covers malformed/out-of-bounds/unknown-contig/empty
  records, duplicates, overlaps, names and reference mismatch. Verify exact
  expected counts, stable IDs, original-file checksums and no discovery jobs
  during BED-only entry.
- ChIP/input matrices share coordinates; inputs never become treatment
  replicates. Raw counts remain unchanged by enrichment settings. Gene-body
  and fixed-bin modes bypass joint discovery.
- Support annotations handle one/both/no groups, multi-group cohorts, partial
  overlap, missing catalogues and IDR fallbacks. Annotation/view filtering must
  not change fitted DE results or adjusted p-values.
- Real local DAG checks cover partial runs, downstream-only restarts, reuse,
  changed BED/controls/q settings, failed jobs, all retention policies and size
  pressure. Protected supplied sources remain byte-identical.
- Run appropriate existing unittest/integration and standalone regressions,
  including RNA/ATAC compatibility; the previous preparation pass is not a
  substitute for acceptance of these changes.

### Statistical and visual review after implementation

Use controlled simulations with known truth and small suitable biological
examples. Simulations should reproduce read placement, discovery and selection,
not merely feed ideal count matrices directly into DESeq2.

Evaluate no-change data, sparse changes, condition-specific sites, broad/local
changes, unequal sequencing depths, unequal replicate numbers, variable
background and shared versus distinct controls. Include global-shift examples
to expose normalization limitations rather than claiming universal recovery.

Compare permissive joint candidates, the current group-union route and an
independent region set where available. Report candidate counts, tested counts,
support categories, widths, count/background distributions, effect sizes,
empirical false-positive/FDR behavior and sensitivity. Inspect representative
browser regions and replicate-aware summaries, especially joint-only hits.
Evaluate normalization-reference choices independently of discovery threshold.

The initial joint q threshold is 0.1. Explore a small prespecified sensitivity
comparison with tighter settings after the baseline is built. Do not choose a
threshold because it produces preferred biological findings or more stars.
Record all evaluated settings and their consequences.

After building and testing, present the resulting tables, plots and diagnostics
for discussion. Then decide whether discovery should be tightened, aggregate
filtering adjusted, or overlap annotations used to simplify presentation.
Do not silently convert support annotations into an inferential admission gate.
Any proposed gate requires explicit reconsideration of selection and error
control. This post-validation review is an agreed checkpoint, not work to do now.

Substantial computation must use the project's Slurm/environment procedures.
Manuscript validations are a later milestone; archive lightweight figures,
reports and provenance only under the prescribed artifact rules.

## Implementation details still to resolve

These details do not reopen the agreed analytical direction:

- Whether a compatible multi-condition experiment uses one discovery cohort
  or explicit per-analysis cohorts, and how that is represented in metadata.
- How unequal library depths/replicate counts affect discovery contributions.
  The agreed all-reads pooling rule for matched inputs remains unchanged;
  any different ChIP discovery weighting requires a documented decision.
- Exact aggregate filter, reliability thresholds and between-ChIP normalization
  reference; audit current behavior and validate assumptions before fixing them.
- BED overlap/duplicate semantics, default support display threshold, and
  final CLI names consistent with the existing configuration system.
- Caller layout representation and narrow/broad settings; physical pool files
  versus supported multiple-BAM arguments, with equivalent unique membership.

## Methodological references

- [Lun and Smyth, 2014: region selection and error control](https://academic.oup.com/nar/article/42/11/e95/1442937).
- [csaw book: filtering and aggregate enrichment](https://bioconductor.org/books/3.17/csawBook/chap-filter.html).
- [Bonhoure et al., 2014: experimental reference normalization](https://pmc.ncbi.nlm.nih.gov/articles/PMC4079971/).
- [Orlando et al., 2014: ChIP-Rx and global changes](https://pubmed.ncbi.nlm.nih.gov/25437568/).

The literature motivates the design; it does not establish calibration of the
exact Omnomnomics implementation before the planned validation is performed.
