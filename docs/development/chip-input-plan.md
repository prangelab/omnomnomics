# ChIP input declarations and preparation: source map and implementation plan

Date: 2026-09-11. Status: preparation milestone implemented on `input_dev`; validation details in [the implementation record](chip-input-implementation.md).

Development branch: `input_dev`, created from the existing `codex/DE_analyses`
checkout at `579730ee9b0f92c163e1f1bda4189add9f5152c5`.
The tracked working tree was clean before branching. Manuscript material and
the paused manuscript review ledger are outside this change.

Immediate prerequisite: the existing cleanup paths could delete supplied BAMs.
`storage.py` now records each filename identity's earliest available source stage
before cleanup, persists that protection, and guards retention, size-pressure,
lane-BAM and forced-output cleanup. Mixed source/output directories are retained
conservatively. The library identity/metadata manifest now refines this protection with durable
source and generated-output provenance, without discarding existing protections.

## Scope and completion boundary

Support ChIP and input libraries in the same normal metadata file, originating
from FASTQ or supplied BAM files. Resolve their relationships using selected
metadata fields. Prepare controls through the shared preprocessing stages and
record their identities, sources, layout, outputs and associations.

This phase does not choose control pooling, change MACS3 control selection,
create input-relative BigWigs, normalize gene-body/bin counts, change statistical
models, or rerun manuscript validations. These decisions follow separately.
Technical lane/replicate merging within a library remains distinct from pooling
different input libraries for a ChIP analysis.

The first acceptance milestone is preparation through public steps 1–7 plus a
resolved control association manifest. A downstream rollout must not imply that
newly declared controls are already being used analytically. Implemented interim
behavior: reject downstream requests using the new control declarations until
their consumers are enabled in the next phase, with a preparation-only example
in the message. Preserve existing metadata-free preprocessing and legacy `-I`
analysis routes. This temporary development boundary is not the final product UI.

## How the current application works

Paths and line numbers below refer to the source revision recorded above.

| Layer | Current behavior and implementation anchor | Consequence for inputs |
| --- | --- | --- |
| CLI and configuration | `src/omnomnomics/cli.py:1531` dispatches assay verbs, merges packaged/site defaults, resolves resources and genome paths. `workflow/config/workflow.yaml:5` defines internal rule numbers, input/output folders and suffixes. | Extend the existing interface/configuration rather than adding an independent control pipeline. |
| Step selection | `cli.py:1107` resolves public steps; `cli.py:1137` maps ChIP/ATAC public steps 11–15 to internal 13,14,11,15,16. | Keep public numbering stable; refer to internal numbers only in implementation. |
| File discovery | `cli.py:1161` selects one folder/suffix from the earliest internal step. It infers one global paired flag from R2 FASTQs in the project. | Mixed FASTQ/BAM entry and library-specific layouts cannot be represented reliably. BAM-only entry does not determine pairing from BAM contents. |
| Metadata | `metadata.py:34` reads TSV/CSV with `filename` first. `:106` normalizes names; `:132` derives sample IDs, types and colors, permitting repeated IDs only for declared technical replicates. | Add role/matching without conflating control identity with biological grouping or technical merging. |
| Metadata activation | `cli.py:80` and `:1900` restrict normal metadata derivation to selected later stages; `:1953` builds DE design from the derived rows. | Supplied metadata must be interpreted in preprocessing-only runs too. Controls must not affect DE formula inference or replicate structure. |
| Runtime handoff | `cli.py:1354` writes run YAML; `:2017` writes derived metadata. `:2149` constructs Snakemake invocation, using mtime rerun triggers and a Slurm profile, then dry-runs or submits a controller. | Freeze one resolved library/association description for controller and workers. New mapping changes need explicit dependency tracking, not merely config parameters. |
| Workflow identity | `workflow/Snakefile.smk:435` loads metadata; `:497` handles merged names; `:950` independently discovers files; `:1027` creates `samples` and `samples2`. Later entry can replace them with all metadata sample IDs at `:1043`. | CLI and workflow currently duplicate discovery. Input rows would enter ordinary analysis sample lists unless separated explicitly. |
| Lanes and technical merges | `Snakefile.smk:1062` searches FASTQ, trimmed FASTQ and BAM folders, stopping at the first populated source. `:1114` resolves merge units. `rules/4.merge_lanes_and_clean_names.smk:15` merges, renames or marks a passthrough BAM. | Discovery must operate per library; BAM-supplied controls must not inherit FASTQ alignment-marker dependencies. |
| Preprocessing | Rules 1–3 trim, FastQC and align. Rule 4 merges; rule 5 sorts, marks/removes duplicates and filters; rules 6–7 index and report alignment QC. ChIP filtered output is `<sample>.filtered.bam`. | Reuse these stages for input libraries, preserving source-specific entry and layout. |
| Target construction | `Snakefile.smk:1298` expands requested outputs over common sample lists. All rules are included at `:1290`; dependencies may bring in upstream work. | A BAM-supplied library must not receive trim/alignment targets. Selective entry needs explicit dependency rules. |
| Existing control use | `cli.py:564` accepts one absolute BAM via `-I`; `rules/10.call_peaks.smk:303`, `:329`, `:404`, `:475`, `:638` pass the global `INPUT` through parameters. | This is an external path, not an input sample or a complete DAG dependency. Later control integration must cover split IDR and monolithic peak-calling routes. |
| Downstream analyses | Peak grouping, peak QC, feature counting and profile generation iterate `samples2`. Gene-body/diffuse feature construction bypasses control-aware peak calling. | Keep controls out of treatment groups, IDR replicate lists, ordinary count columns and treatment profiles. Analytical use remains separate work. |
| Tracks | Rule 8 generates coverage; `rules/9.merge_wiggles.smk:151` independently enumerates every `.bw` in the folder. | Merely filtering `samples2` is insufficient: directory scans can reintroduce controls or stale tracks. |
| DE and explorer | Rules 12/15 render R templates and copy derived metadata into analysis outputs; `de_app.py` serves the resulting explorer. | Maintain experiment-only DE metadata and unchanged result interfaces in this phase. |
| Monitoring and reporting | `Snakefile.smk:149` predicts step totals from shared lists; `:291` marks completion and triggers size cleanup. `:1614` aggregates flow QC; `cli.py:297` separately rediscovers monitor totals. | Totals must reflect actual planned work per stage and include input QC with role labels. |
| Storage and recomputation | `Snakefile.smk:576`, `:615`, `:842` can delete whole intermediate folders. `:2118` deletes lane-named BAMs after merging. `cli.py:1436` deletes selected output suffixes before metadata processing. | Supplied BAMs need persistent source ownership and protection before any deletion can occur. |

## Proposed metadata contract

Keep `filename` as the first column. Add an explicit `role` column for new ChIP
control declarations, with values `chip` and `input`. For backward compatibility,
absence of this column means ordinary experimental samples. Do not reserve or
reinterpret unrelated RNA/ATAC metadata columns globally. In a ChIP sheet that
declares roles, reject empty/unknown roles rather than guessing.

Reuse existing column-selector syntax for a proposed `--input-match` option;
names or one-based indices are accepted. Match exact tuples of trimmed original
metadata values, before identifier sanitization. Match only against `input`
rows. Empty matching fields are invalid, not wildcards; matching is independent
of `--sample-type`, track colors and DE columns.

Example metadata (one row per existing metadata input unit, not one row per read):

```text
filename        role    cell_type    condition    replicate    library
chip_A_1        chip    macrophage   untreated    1            chip_A_1
chip_A_2        chip    macrophage   untreated    2            chip_A_2
input_A_1       input   macrophage   untreated    1            input_A_1
input_A_2.bam   input   macrophage   untreated    2            input_A_2
```

Here `--sample-name library` keeps identities unique. Matching on `cell_type`
shares eligible inputs by cell type; adding `condition` narrows the association;
adding `replicate` gives replicate matching. An `input_group` column can express
irregular sharing using the same matcher (`--input-match input_group`), without a
separate mapping file or a second matching algorithm. A common label provides
the single shared-control case.

Multiple matches produce a deterministic set of distinct input library IDs.
They do not merge files or create extra biological replicates. Resolve technical
units to their library identity before deduplicating associations. Technical
units that merge must agree on role and matching metadata. Reject chip/input
sample-ID collisions before existing technical-replicate logic can merge them.

For declared matching, report zero-match ChIP samples as errors, list unmatched
input libraries separately, and never fall back to a different group. Extra
unassigned input libraries may still be prepared and retained with that status.
Allow input-only preprocessing without requiring a treatment association.
Reject new metadata-based control declarations together with legacy `-I`;
continue supporting `-I` unchanged on its existing route.

## Source resolution and runtime records

The implementation uses small testable Python modules, including
`src/omnomnomics/library_manifest.py` and `src/omnomnomics/chip_inputs.py`.
Reuse `metadata.py` selectors and normalization; reconcile the existing CLI and
Snakefile suffix handling against compatibility tests rather than creating a
third filename parser.

Each resolved library record should include original metadata/file key, runtime
sample ID, role, technical units, original source paths/type, SE/PE layout,
selected starting stage, expected BAM/index paths and source ownership. Keep
original source and current reusable intermediate separate. Matching records
contain ChIP ID, input ID, selector names and the values that established the
match. Save a readable association TSV alongside the runtime manifest and archive
both with the run configuration. Preserve source ownership across later runs.

FASTQs remain in `FASTQ/`; externally supplied aligned BAMs remain in `BAM/`.
The BAM contract is an aligned library entering before ordinary ChIP filtering;
it is not an arbitrary file assumed already filtered. Existing filtered-BAM
partial entry remains supported. Already removed reads cannot be recovered or
preprocessing inferred reliably from a filename; document supplied BAM history.

Explicit `.bam` metadata filenames can identify supplied BAM sources. Bare roots
remain supported where the source is unambiguous. When FASTQ and BAM coexist,
reuse known pipeline intermediates from recorded provenance; unknown competing
sources must not be chosen merely by timestamp or folder search order. Protect
supplied paths, indexes and resolved symlink targets from deletion/overwrite.

Validate readability, complete read pairs, identity collisions and BAM integrity.
Check BAM reference names/lengths against the selected genome and record the
declared assembly; compatible headers alone cannot prove assembly provenance.
Infer BAM layout from alignment flags with a bounded inspection and reject
unresolved/inconsistent layouts rather than defaulting silently to single-end.

Represent layout per library from the start. Support preparation of SE inputs
alongside PE ChIP libraries, and vice versa. Do not broaden downstream mixed-layout
ChIP analysis support incidentally; retain and validate its existing constraints.
Because some Snakemake output paths are selected at parse time, use constrained
SE/PE rule variants where necessary, not a global flag changed between jobs.

## Ordered implementation units

| Order | What changes | Where | Reviewable result |
| --- | --- | --- | --- |
| 1 | Define roles, exact matching, source records, technical-unit consistency and backward compatibility as pure functions. | `metadata.py`; new `chip_inputs.py`, `library_manifest.py`; focused tests. | Metadata resolves deterministically into libraries and associations without running tools or writing project data. |
| 2 | Parse supplied metadata before discovery/recomputation; resolve source ownership; protect original data during every cleanup route. | `cli.py:1161`, `:1436`, `:1870`, `:1900`; `Snakefile.smk:576`, `:615`, `:842`, `:2118`; shared cleanup helper if needed. | Invalid metadata cannot delete files; imported BAMs survive retention, size pressure and forced reruns. Protect sources before enabling mixed-source jobs. |
| 3 | Replace duplicate discovery for the new route with a serialized manifest; resolve stage-specific work and preserve ordinary/legacy runs. | `cli.py` arguments/config writing/controller invocation; `workflow/config/workflow.yaml`; `Snakefile.smk:435`, `:950`, `:1027`, `:1062`, `:1298`. | Mixed FASTQ/BAM dry-run has exactly the intended jobs and source paths. |
| 4 | Adapt preprocessing dependencies, layout selection and merge passthrough; keep filtering defaults. | `rules/1.trim_fastp.smk`, `1.trim_skewer.smk`, `2.fastqc.smk`, `3.align_reads_hisat2.smk`, `3.align_reads_STAR.smk`, `3.align_reads_STAR_TE.smk` where shared wiring applies, `4.merge_lanes_and_clean_names.smk`, `5.touchup_bam.smk`, `6.index_bam.smk`, `7.bam_stats.smk`. | Both roles reach filtered/indexed BAMs; supplied BAMs never request absent FASTQs, mapper stats or alignment markers. No mapper algorithm changes. |
| 5 | Make preprocessing, experimental and input collections explicit; separate all-library metadata from experiment-only DE metadata. | `Snakefile.smk` collections; `cli.py` DE derivation; rules 8–16 and optional HOMER consumers; track directory scans. | Input libraries cannot enter treatment groups, treatment count columns or profiles through either lists or folder scans. Later analytical consumers remain disabled for the new route pending their design. |
| 6 | Record matching/source provenance and stage-aware progress; make reruns reproducible. | `cli.py:297`, run configuration/log writing and allowed-rule lists; `Snakefile.smk:149`, `:291`, target manifest, flow/alignment QC aggregation. | Logs show prepared inputs and associations; imported BAMs show unavailable upstream metrics rather than fabricated zero counts. Reused jobs do not prevent completion or trigger cleanup prematurely. |
| 7 | Document the interface and validate the preparation milestone. | `docs/guides/metadata.md`, `docs/getting-started/project-layout.md`, `docs/assays/chip.md`, `docs/reference/cli.md`, `docs/operations/storage.md`, `README.md`, tests. | Examples cover shared, condition-matched, replicate-matched and mixed-source preparation; existing RNA/ATAC behavior is preserved. |

Avoid replacing every occurrence of `samples2` mechanically. Rules 1–7 need
applicable libraries of both roles; experimental consumers need ChIP samples;
QC aggregation needs all prepared libraries; technical merging needs units
within one library. Preserve filename-key versus derived sample-ID behavior on
each existing entry route.

Limit scheduling changes to manifest-resolved inputs and targets. An omitted
upstream stage may use an existing valid product under the existing partial-entry
contract; it must not unexpectedly schedule an entire upstream analysis. Keep
original public step numbering and update allowed-rule companion lists when
new SE/PE variants are introduced. A downstream-only metadata/table restart
must not suddenly require retained raw FASTQs or imported BAMs it does not use.

Cleanup must use known ownership and actual successful consumers, including
reused products, rather than deleting a whole folder because a counter finished.
Failure must never qualify as successful consumption. Preserve lightweight flow
QC caches before deleting eligible generated files. If selective deletion cannot
prove safety, retain the folder and log the reason. Protect imported lane-named
BAMs under the existing end-of-run lane cleanup too.

Keep semantic manifests stable across unchanged reruns, with immutable archived
copies per run. Current `--rerun-triggers mtime` does not notice a control path
changed only in `params`; future consumers must depend on their actual BAMs and
an association artifact that changes when relevant matching changes. Do not
invalidate every output merely because a new timestamped run YAML exists.

## Validation before calling preparation support complete

Use existing managed environments; no new dependency is required by this plan.
The source review itself required no test execution. The subsequent source
protection fix has dedicated behavioral tests in `tests/test_source_protection.py`.

| Test level | Required cases |
| --- | --- |
| Metadata behavior | Shared, condition, replicate and explicit-group matching; multiple inputs retained as a set; empty/unknown roles and matching values; zero matches; unused inputs; input-only preparation; role/ID collisions; conflicting technical units; values that sanitize to the same identifier; legacy metadata and `-I` conflict. |
| Source behavior | All FASTQ; all supplied BAM; ChIP FASTQ plus input BAM; ChIP BAM plus input FASTQ; FASTQ suffix/lane conventions; technical replicates; known intermediates coexisting with sources; ambiguous sources; missing pairs; BAM layout/reference mismatch; mixed SE/PE preparation. |
| Real DAG behavior | Use tiny local fixtures and an existing Snakemake environment with a local executor. Inspect scheduled jobs and dependencies for public `1–3`, `1–7`, `4–7`, `5–7`, index/QC-only restarts and existing later-entry routes. Imported BAMs must have no alignment producer, no alignment sentinel dependency and no source/output collision. |
| Ownership and failure | Temporary-project tests for every deletion route, including force-recompute and lane cleanup; imported sources/indexes remain byte-identical; failed/partial stages retain required inputs; changing metadata cannot erase source ownership. |
| Reuse and reporting | An unchanged successful rerun schedules no substantive work; mixed reused/new jobs report correct totals; source provenance survives cleanup; role-specific QC remains visible; no-mapper BAM entry does not invent alignment/trim statistics. |
| Consumer isolation | A fixture with input BAMs and stale input BigWigs cannot add them to treatment groups, IDR pairs, ordinary count columns, DE metadata or track groups. New-route downstream requests cannot silently run without applying the declared controls. |
| Existing regressions | Run the established standalone regression runner and unittest discovery: the standalone runner skips unittest classes and functions requiring fixtures. Include metadata partial entry, BAM layout, step mapping, config layering, count staging, peak layout and chromatin dependency suites. Add behavioral/DAG tests instead of relying solely on source-string assertions. |

Actual tiny BAM preprocessing checks are appropriate after implementation.
Substantial validation datasets and Snellius jobs are not part of this phase.
Any later cluster test must use the documented environment and Slurm workflow.

## Analytical follow-up after this milestone

The subsequent design discussion is recorded in
[the ChIP input use and differential-analysis plan](chip-input-analysis-plan.md).
It covers matched controls throughout peak calling, permissive joint discovery,
custom BED regions, overlap annotations, separate input enrichment and raw-ChIP
differential analysis, plus tracks, validation and a post-validation review.
That phase is planned, not implemented. Work is paused after recording the plan.
Only after analytical implementation and validation should manuscript ChIP
validation runs be replaced.
