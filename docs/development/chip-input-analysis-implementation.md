# ChIP analytical input implementation

Branch: `input_dev`. Date: 2026-09-14. Changes are local and uncommitted.

## Implemented behavior

Metadata-declared controls now feed analytical ChIP stages while preserving
legacy `-I` and deliberate no-input runs. This extends the separately tested
[preparation milestone](chip-input-implementation.md). The manuscript review
ledger remains paused; no manuscript files or manuscript validation runs were
changed, and no Snellius jobs were submitted.

- Single controls are used directly; MACS pools multiple fragment files with all
  reads. Group controls are a unique union of member associations. Replicate
  calls use their own inputs and pseudoreplicates inherit full parent controls.
  These rules cover split IDR, MACS-only optimization and broad-domain calls.
- ChIP MACS uses BEDPE fragments: PE proper-pair inserts and inferred SE
  extensions (default 200 bp). Mixed layouts across libraries are explicit.
  Redundant optimizer shift variants are omitted; `fast` and `full` evaluate
  q=0.01 and q=0.001. All ChIP calls use the same genome-size convention, with
  `--chip-effective-genome-size` available to override approximate total
  reference length. This removes implicit human-size group calls on other
  organisms. Actual or inferred fragment coordinates replace former SE model
  fitting, so historical peak sets need not be identical.
- Peak-backed DE uses a separate joint q=0.1 catalogue (joint broad cutoff 0.1
  for domain mode), or `--chip-regions-bed`. All compatible ChIP libraries in
  one run contribute all reads to discovery, with unique matched controls.
  One mark/antibody and genome per run is the supported cohort convention;
  unequal depths/replicate numbers therefore weight discovery unequally.
- Custom BEDs require non-overlapping, unique, reference-bounded coordinates.
  Original names and checksums are retained; the user file stays protected.
  Gene-body and bin modes retain their predefined coordinates and bypass
  joint MACS. Overlapping annotation-derived gene bodies are permitted;
  ambiguous feature assignments are excluded. Identical gene coordinates
  already coalesce labels in the annotation builder.
- Raw ChIP and input count matrices share coordinates but remain separate.
  Counting excludes unmapped, secondary, supplementary and QC-failed reads;
  PE additionally requires proper same-contig pairs. Counts use unique feature
  assignment based on aligned-read/mate overlap, rather than every position
  inside an inferred insert. Library mapped depths include eligible fragments
  outside regions; mixed-layout enrichment needs this geometry caveat.
- Regional enrichment is ChIP CPM / summed matched-input CPM after summing
  input counts and mapped depths. Input contributions are depth weighted.
  No pseudocount; input counts below 5 are `low_input` with undefined ratios
  (site configuration can change that threshold). Zero ChIP has fold zero and
  undefined log2 ratio. No controls produces `no_input`. This is not percent input.
- Routine DE uses raw ChIP counts with biological replication and its chosen
  experimental design. No input subtraction or division. The ChIP prefilter
  uses aggregate raw counts, including the configured minimum total count;
  it does not require support in both conditions. Normalization/design/filter
  choices remain in existing DE configuration and reports. No normalization
  default was changed to match the pilot below.
- Group peak overlap is annotation only. Tables retain union overlap bases and
  fractions, supported/unassessed groups, original names and available IDR
  consensus/fallback status. Display support requires any positive overlap.
  Factor contrasts add numerator-only/denominator-only/both/joint-only/unassessed
  categories according to the run's actual peak grouping. Full and significant
  result tables retain support metadata and the full-family adjusted p-values.
- Ordinary coverage stays available. Additional tracks select MACS modelled
  background fold enrichment, direct CPM read-coverage log2 ratio (1 CPM
  pseudocount, 10 bp bins), MACS −log10(q), or none. JSON records commands,
  identities and scale. Track hubs use separate labelled containers.
  Optional extra track generation is omitted and logged when a project-size
  cap is set.
- Post-DE profiles remain descriptive means across regions with one line per
  library. ChIP labels/provenance state this; no treatment test uses loci/bins
  as biological replicates. These profiles use ordinary coverage tracks.
- Semantic dependency manifests separate control, group-caller and joint-region
  settings. BED/controls/q changes invalidate their consumers. Group relabelling
  does not change joint coordinates. Pooled/pseudo fragment conversion checks
  source BAM size/mtime and layout/extension before reuse.
- Counts archive coordinates, names and source/settings provenance, allowing
  DE-only restarts without generated BAMs or the peak tree. Changed region or
  control/enrichment settings require public step 13 to refresh counts/summaries. Per-library count caches
  and ChIP DE model caches reuse unchanged work; explicit recomputation bypasses
  them. Support changes do not refit the model or alter its p-values.
- Retention/size/force source protection remains active for original reads,
  imported BAMs/indexes, metadata and custom BEDs. Generated MACS fragment files
  currently follow retention of the containing peak tree; this is conservative
  and can retain bulky intermediates. The terminal-directory mapping now uses
  internal step IDs correctly for count/DE/HOMER products.

Optional HOMER tag-directory mode remains explicitly unavailable for role-based
inputs. It is not silently routed to an uncontrolled analytical path.

## Local functional validation

Acceptance checks: **78 tests passed in the complete unittest suite**. Three
additional configuration/storage tests and focused reruns cover the final
threshold forwarding/invalidation and capped-track logging changes (81 distinct
unittest tests in total). **45 standalone regression checks passed** and
`git diff --check` passed. Both workflow includes are present in package data. Archived outputs:
[full suite](chip-input-analysis-validation/unittest.txt),
[focused counts/configuration](chip-input-analysis-validation/focused_counts.txt),
[storage cap](chip-input-analysis-validation/storage_cap.txt) and
[standalone regressions](chip-input-analysis-validation/standalone.txt).

Commands:

```text
python -m unittest discover -s tests
PYTHONPATH=src:tests python -m unittest test_chip_analysis_config
PYTHONPATH=src:tests python -m unittest test_chip_analysis.AnalysisWorkflowTests.test_custom_bed_counting_mixed_layout_and_reuse test_chip_analysis.AnalysisWorkflowTests.test_gene_body_and_bins_have_input_enrichment_without_joint_macs test_chip_analysis.AnalysisWorkflowTests.test_legacy_shared_input_and_explicit_no_input_counting
python tests/run_regression_tests.py
```

The final full suite took 364 seconds; the storage-cap DAG check passed in 2 seconds; the focused count/configuration rerun
passed 5 tests in 79 seconds. An earlier run overlapped source edits and had a
controller/worker settings-manifest mismatch; the subsequent clean full run
passed. These results are local software checks, not a large biological or
statistical calibration.

Tests use temporary datasets in the existing `omnomnomics-macos-test`
environment. No dependencies were installed or updated. Workflow jobs run
outside the sandbox only to access Snakemake's existing macOS cache. Other
checks use writable task-specific Python cache directories.

Real workflow checks cover mixed SE/PE exact counts, unique PE fragments,
custom BED protection/reuse/invalidation, controlled joint/narrow/broad and
split-IDR pseudo calls, predefined gene bodies/bins, all three enrichment track
modes and separate track-hub containers. Full small BED DE reports run without
peak QC, retain annotations, reuse statistical fitting after support changes,
and restart after generated BAM/peak removal. Legacy shared control/no-input
counting and existing preparation/source-protection/RNA/ATAC regressions are
included. Joint catalogue tests check a strong condition-specific signal,
coordinate invariance under grouping changes and q-setting invalidation.

The Shiny table reader was inspected: it retains arbitrary result columns.
An interactive Shiny browser session and full post-DE motif/signal jobs were
not executed. Cluster scheduling, production-scale resource use and every
external mapper variant remain beyond these local checks.

Existing tools include MACS3 3.0.4, samtools 1.23.1, featureCounts 2.1.1 and
R 4.4.1. MACS bedGraphs are converted with the existing deepTools-associated
pyBigWig; no ad hoc installation was used for the locally absent UCSC converter.
The new workflow include is declared in package data.

## Read-level scientific pilot

The reproducible [pilot script](../../scripts/validation/chip_input_analysis_pilot.py)
uses the implementation's BEDPE conversion, joint q=0.1 discovery and raw
featureCounts counting. Each case has 80 seeded sites, 3 WT + 3 KO libraries
and uniform shared input. Region abundance is rounded gamma noise with simple
read placement, not a fully specified negative-binomial/biological model.
DESeq2 uses the pilot's explicit poscounts reference and local fit. See the
[results TSV](chip-input-analysis-validation/pilot.tsv) and
[provenance JSON](chip-input-analysis-validation/pilot.json).

| Scenario | Findings at adjusted p < 0.05 |
| --- | --- |
| Five null seeds | seeds 0 and 1 each have one false positive; seeds 2–4 have none |
| Sparse loss: 12 of 80 sites | 12 true positives and 2 false positives; realized FDP 2/14 ≈ 0.143 |
| Global loss: 80 of 80 sites | no significant sites; relative normalization absorbs the global shift |

These results do **not** establish FDR calibration. In two of five null runs
there was at least one false discovery; the sparse example also exceeds 5%
realized FDP. Five simplified null seeds cannot attribute this to discovery,
model misspecification or sampling variability, and a single run's realized
FDP is not the method's expected FDR. The global case demonstrates reference
dependence, not failure to recover an absolute quantity the method measures.
Do not turn these pilot results into a validated accuracy or calibration claim.

## Agreed review checkpoint

The software baseline is built. Keep q=0.1 and overlap annotation as agreed;
no support admission gate or outcome-selected tightening was added. Before
manuscript validation/release claims, perform the broader statistical/visual
review in the [agreed plan](chip-input-analysis-plan.md): more realistic nulls,
condition-specific and broad changes, unequal depths/replication/background,
shared/distinct inputs, joint-vs-group-union-vs-independent regions and a
prespecified threshold sensitivity comparison. Biological examples/browser
inspection and normalization-reference evaluation are still required. This
record distinguishes passing software checks from that remaining review.

## Changed-file review map

The map includes the existing preparation/source-protection changes in this
working branch. Line anchors identify the first changed line of tracked files,
or line 1 for newly added files; links are absolute workspace paths.

| File | Change starts |
| --- | --- |
| [README.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/README.md:386>) | 386 |
| [docs/assays/chip.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/assays/chip.md:27>) | 27 |
| [docs/development/chip-input-analysis-implementation.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-implementation.md:1>) | 1 |
| [docs/development/chip-input-analysis-plan.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-plan.md:1>) | 1 |
| [docs/development/chip-input-analysis-validation/focused_counts.txt](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/focused_counts.txt:1>) | 1 |
| [docs/development/chip-input-analysis-validation/pilot.json](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/pilot.json:1>) | 1 |
| [docs/development/chip-input-analysis-validation/pilot.tsv](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/pilot.tsv:1>) | 1 |
| [docs/development/chip-input-analysis-validation/standalone.txt](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/standalone.txt:1>) | 1 |
| [docs/development/chip-input-analysis-validation/storage_cap.txt](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/storage_cap.txt:1>) | 1 |
| [docs/development/chip-input-analysis-validation/unittest.txt](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/unittest.txt:1>) | 1 |
| [docs/development/chip-input-implementation.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-implementation.md:1>) | 1 |
| [docs/development/chip-input-plan.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-plan.md:1>) | 1 |
| [docs/getting-started/project-layout.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/getting-started/project-layout.md:29>) | 29 |
| [docs/guides/metadata.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/guides/metadata.md:39>) | 39 |
| [docs/operations/storage.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/operations/storage.md:12>) | 12 |
| [docs/reference/cli.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/reference/cli.md:37>) | 37 |
| [pyproject.toml](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/pyproject.toml:41>) | 41 |
| [scripts/validation/chip_input_analysis_pilot.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/scripts/validation/chip_input_analysis_pilot.py:1>) | 1 |
| [src/omnomnomics/chip_analysis.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/chip_analysis.py:1>) | 1 |
| [src/omnomnomics/chip_inputs.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/chip_inputs.py:1>) | 1 |
| [src/omnomnomics/cli.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/cli.py:36>) | 36 |
| [src/omnomnomics/library_manifest.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/library_manifest.py:1>) | 1 |
| [src/omnomnomics/metadata.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/metadata.py:26>) | 26 |
| [src/omnomnomics/preprocessing.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/preprocessing.py:1>) | 1 |
| [src/omnomnomics/storage.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/storage.py:1>) | 1 |
| [src/omnomnomics/workflow/Snakefile.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/Snakefile.smk:23>) | 23 |
| [src/omnomnomics/workflow/chip_analysis.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/chip_analysis.smk:1>) | 1 |
| [src/omnomnomics/workflow/chip_preparation.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/chip_preparation.smk:1>) | 1 |
| [src/omnomnomics/workflow/config/workflow.yaml](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/config/workflow.yaml:26>) | 26 |
| [src/omnomnomics/workflow/rules/1.trim_fastp.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/1.trim_fastp.smk:20>) | 20 |
| [src/omnomnomics/workflow/rules/1.trim_skewer.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/1.trim_skewer.smk:18>) | 18 |
| [src/omnomnomics/workflow/rules/10.call_peaks.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/10.call_peaks.smk:118>) | 118 |
| [src/omnomnomics/workflow/rules/11.count_reads.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/11.count_reads.smk:19>) | 19 |
| [src/omnomnomics/workflow/rules/15.call_DE_chrom.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/15.call_DE_chrom.smk:217>) | 217 |
| [src/omnomnomics/workflow/rules/16.analyze_peaks_de.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/16.analyze_peaks_de.smk:885>) | 885 |
| [src/omnomnomics/workflow/rules/2.fastqc.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/2.fastqc.smk:14>) | 14 |
| [src/omnomnomics/workflow/rules/3.align_reads_STAR.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/3.align_reads_STAR.smk:17>) | 17 |
| [src/omnomnomics/workflow/rules/3.align_reads_STAR_TE.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/3.align_reads_STAR_TE.smk:17>) | 17 |
| [src/omnomnomics/workflow/rules/3.align_reads_hisat2.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/3.align_reads_hisat2.smk:17>) | 17 |
| [src/omnomnomics/workflow/rules/4.merge_lanes_and_clean_names.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/4.merge_lanes_and_clean_names.smk:25>) | 25 |
| [src/omnomnomics/workflow/rules/5.touchup_bam.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/5.touchup_bam.smk:27>) | 27 |
| [src/omnomnomics/workflow/rules/6.index_bam.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/6.index_bam.smk:13>) | 13 |
| [src/omnomnomics/workflow/rules/7.bam_stats.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/7.bam_stats.smk:33>) | 33 |
| [src/omnomnomics/workflow/rules/9.merge_wiggles.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/9.merge_wiggles.smk:17>) | 17 |
| [src/omnomnomics/workflow/templates/de_core_chrom.R.tmpl](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/templates/de_core_chrom.R.tmpl:23>) | 23 |
| [tests/test_bam_touchup_layout.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_bam_touchup_layout.py:11>) | 11 |
| [tests/test_chip_analysis.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_chip_analysis.py:1>) | 1 |
| [tests/test_chip_analysis_config.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_chip_analysis_config.py:1>) | 1 |
| [tests/test_chip_inputs.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_chip_inputs.py:1>) | 1 |
| [tests/test_chromatin_partial_dependencies.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_chromatin_partial_dependencies.py:23>) | 23 |
| [tests/test_source_protection.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_source_protection.py:1>) | 1 |

## Review repairs (2026-09-14)

- Legacy `-I` controls receive an internal ID that cannot collide with any
  analytical or declared library ID. Existing `legacy_input` output labels are
  retained when that name is available. A ChIP library named `legacy_input`
  keeps its own BAM and counts; its control gets a distinct ID.
- SPP drop lists affect counting only in `drop` mode. The count settings track
  the gate, and a controller-written sample-selection dependency invalidates
  counts when a drop list is removed or changes the selected libraries, even
  with mtime-only reuse. Workers do not rewrite this dependency. Count
  provenance records the gate and applied exclusions. DE-only entry checks
  archived count columns against the current selection and requires public
  step 13 when they differ.
- Automatic DE-column selection considers ChIP rows alone. Full library
  metadata remains archived, and explicit or assisted designs retain their
  existing behavior.
- Count provenance records the project root. Region/control comparison uses
  relative project-local BAM/BED paths while retaining external path
  identities and all other settings, including the BED checksum. Older
  archives infer their root from canonical ChIP BAM paths. Relocating a project
  therefore permits DE-only entry without BAMs; changed controls, layouts,
  coordinates or discovery settings still require refreshed counts.

Validation: six new regression tests passed, including real SE/PE counting,
gate transitions, drop-list removal, old/new archive relocation without BAMs,
changed BED rejection and automatic design selection. All 45 standalone
regression checks passed. The full unittest suite passed: **87 tests in
448.001 seconds**. `git diff --check` passed.

Archived test logs:

- [Focused repairs](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/review_fixes.txt:1>)
- [Full unittest suite](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/review_full.txt:1>)
- [Standalone regressions](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-analysis-validation/review_standalone.txt:1>)

Changed source and documentation anchors for these repairs:

| File | Lines | Change |
| --- | --- | --- |
| [chip_analysis.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/chip_analysis.py:28>) | 28, 60 | Collision-free legacy identity; portable region-spec comparison |
| [chip_analysis.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/chip_analysis.smk:23>) | 4, 23, 37, 118, 185, 220 | Control identity, gate-aware selection and provenance |
| [11.count_reads.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/11.count_reads.smk:18>) | 18 | Archive validation and sample-selection dependency |
| [metadata.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/metadata.py:257>) | 257, 323 | Separate rows for automatic design selection |
| [cli.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/cli.py:2076>) | 2076 | Select automatic design using ChIP rows |
| [test_chip_review_regressions.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_chip_review_regressions.py:15>) | 15, 37, 128 | Six regression tests |
| [chip.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/assays/chip.md:156>) | 156 | Restart, QC and automatic-design behavior |

These repairs remain local on `input_dev`. No statistical benchmark, Snellius
validation, manuscript edits or manuscript-ledger work was performed.
