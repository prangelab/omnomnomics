# ChIP input preparation implementation

Branch: `input_dev`. Date: 2026-09-11.

The preparation milestone implements metadata-declared ChIP/input libraries,
exact metadata matching, mixed FASTQ/BAM entry and library-specific SE/PE
processing through public steps 1–7. It preserves source files and records
associations without making analytical choices about how controls are used.

## Behavior

- Declare libraries in one metadata file with `role=chip` or `role=input`.
  Select matching fields with `--input-match`; all selected fields must agree.
  Multiple matches remain a deterministic set of input library IDs.
- FASTQs enter from `FASTQ/`, supplied aligned BAMs from `BAM/`, and existing
  filtered BAMs from `filtered_BAM/`. Explicit filenames resolve source
  ambiguity. Technical merging remains within a library.
- The CLI validates metadata, layout, reference compatibility, source collisions
  and selective-stage prerequisites before output deletion or job submission.
  Role declarations cannot be combined with legacy `-I`.
- The shared preprocessing rules prepare both roles. Stage totals count only
  applicable units/libraries, distinguish reused jobs, and retain successful
  completion products for unchanged reruns.
- Source protection covers retention, forced deletion, size-pressure cleanup and
  lane-BAM removal. Recorded ownership survives later stage selections; failed
  or unfinished consumers cannot authorize cleanup. Protected source indexes
  are reused. Mixed source/output directories are retained conservatively.
- Run configurations reference immutable hash-named library and association
  snapshots. Readable current manifests, role-labelled step/QC records and
  experiment-only metadata are also retained.

## Boundary at the preparation milestone

At this milestone, role-based runs explicitly rejected downstream steps beyond 7 and optional
HOMER output. MACS3/IDR control integration, pooling, input-relative BigWigs,
bin/gene-body normalization and manuscript validation reruns remain for the
next phase. Existing metadata-free preprocessing and legacy `-I` analysis
remain on their original routes. The manuscript and review ledger were not edited.

The subsequent [analytical implementation](chip-input-analysis-implementation.md)
supersedes the steps 1–7 boundary. The acceptance figures below describe the
preparation milestone only.

BAM validation uses quickcheck, reference header comparison and a bounded
inspection of 10,000 alignments; it cannot prove processing history or assembly
provenance. SE reads named only R1 need `read_layout=SE` to distinguish them from
a missing mate. Runtime manifests contain absolute paths; a relocated project
needs a fresh CLI run to resolve paths at its new location.

## Validation

Acceptance results: **65 unittest tests passed** (including local Snakemake
integration jobs), **45 standalone regression checks passed**, and
`git diff --check` passed.
Commands: `python -m unittest discover -s tests` and
`python tests/run_regression_tests.py`.

The existing `omnomnomics-macos-test` environment was used, with writable Python
cache directories and its existing Java runtime. No dependencies were installed
or updated. Tests use temporary local fixtures; no Snellius jobs were submitted.

Real job checks exercise mixed-layout trimming/FastQC, noncanonical supplied
trimmed filenames, supplied BAM filtering/indexing/QC, filtered-BAM import,
unchanged reruns and mixed reused/new work. Original files are compared byte for
byte. DAG checks cover mixed source entry, forced SE/PE rule selection and
legacy RNA/ATAC BAM entry. Pure tests cover matching, collision/error cases,
partial restarts, ownership and cleanup after failures.

Skewer and STAR variants are parsed with the workflow; their external programs
were not executed in the local acceptance tests. Large datasets and cluster
scheduling remain untested in this phase.

## Changed-file review map

Links point to the first changed line in each tracked file, or line 1 of a new
file. Source-protection work from the prerequisite is included in this branch's
working diff.

| File | Change starts |
| --- | --- |
| [README.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/README.md:386>) | 386 |
| [docs/assays/chip.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/assays/chip.md:30>) | 30 |
| [docs/development/chip-input-plan.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/development/chip-input-plan.md:1>) | 1 |
| [docs/getting-started/project-layout.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/getting-started/project-layout.md:29>) | 29 |
| [docs/guides/metadata.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/guides/metadata.md:39>) | 39 |
| [docs/operations/storage.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/operations/storage.md:12>) | 12 |
| [docs/reference/cli.md](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/docs/reference/cli.md:37>) | 37 |
| [pyproject.toml](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/pyproject.toml:41>) | 41 |
| [src/omnomnomics/chip_inputs.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/chip_inputs.py:1>) | 1 |
| [src/omnomnomics/cli.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/cli.py:36>) | 36 |
| [src/omnomnomics/library_manifest.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/library_manifest.py:1>) | 1 |
| [src/omnomnomics/metadata.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/metadata.py:26>) | 26 |
| [src/omnomnomics/preprocessing.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/preprocessing.py:1>) | 1 |
| [src/omnomnomics/storage.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/storage.py:1>) | 1 |
| [src/omnomnomics/workflow/Snakefile.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/Snakefile.smk:23>) | 23 |
| [src/omnomnomics/workflow/chip_preparation.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/chip_preparation.smk:1>) | 1 |
| [src/omnomnomics/workflow/config/workflow.yaml](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/config/workflow.yaml:26>) | 26 |
| [src/omnomnomics/workflow/rules/1.trim_fastp.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/1.trim_fastp.smk:20>) | 20 |
| [src/omnomnomics/workflow/rules/1.trim_skewer.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/1.trim_skewer.smk:18>) | 18 |
| [src/omnomnomics/workflow/rules/2.fastqc.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/2.fastqc.smk:14>) | 14 |
| [src/omnomnomics/workflow/rules/3.align_reads_STAR.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/3.align_reads_STAR.smk:17>) | 17 |
| [src/omnomnomics/workflow/rules/3.align_reads_STAR_TE.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/3.align_reads_STAR_TE.smk:17>) | 17 |
| [src/omnomnomics/workflow/rules/3.align_reads_hisat2.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/3.align_reads_hisat2.smk:17>) | 17 |
| [src/omnomnomics/workflow/rules/4.merge_lanes_and_clean_names.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/4.merge_lanes_and_clean_names.smk:25>) | 25 |
| [src/omnomnomics/workflow/rules/5.touchup_bam.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/5.touchup_bam.smk:27>) | 27 |
| [src/omnomnomics/workflow/rules/6.index_bam.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/6.index_bam.smk:13>) | 13 |
| [src/omnomnomics/workflow/rules/7.bam_stats.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/7.bam_stats.smk:33>) | 33 |
| [src/omnomnomics/workflow/rules/9.merge_wiggles.smk](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/src/omnomnomics/workflow/rules/9.merge_wiggles.smk:153>) | 153 |
| [tests/test_bam_touchup_layout.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_bam_touchup_layout.py:11>) | 11 |
| [tests/test_chip_inputs.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_chip_inputs.py:1>) | 1 |
| [tests/test_source_protection.py](</Users/k.h.prange/Library/CloudStorage/OneDrive-AmsterdamUMC/Documenten/Tech/omnomnomics/tests/test_source_protection.py:1>) | 1 |
