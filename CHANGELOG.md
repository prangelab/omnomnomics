# Changelog

All notable changes to this project are documented in this file.

## Unreleased

### Added
- ATAC/ChIP post-DE interpretation reports now write explicit `signal_runs.tsv` and `motif_runs.tsv` status tables covering computed, reused, skipped, no-motif, timeout, and failure outcomes.
- ChIP `--broad-mode genebody` and `--broad-mode diffuse` documentation for feature definitions, DE behavior, and post-DE interpretation outputs.
- Peak annotations and differential-chromatin result tables now include explicit `nearest_promoter_gene` and `distance_to_nearest_promoter_bp` fields.

### Changed
- Post-DE profile documentation now reflects the custom matplotlib renderer used to keep legends outside the signal axes.
- Motif analysis documentation now describes the permanent MEME-format motif database cache under the configured genome assembly root.
- Peak annotation derives gene spans and strand-aware promoters from transcript or exon records when a GTF has no `gene` features, while retaining alternative transcript start sites.
- DE summary plots suppress redundant metadata annotations, enrichment discovery includes nested clusterProfiler result directories, and enrichment plot margins accommodate long titles.
- QC distance plots use sample identifiers, equivalent metadata partitions collapse to one annotation, and combined-summary pathway and motif labels are compacted without altering result tables.
- Post-DE signal profiles use fixed figure margins so labels and external legends remain inside the PDF canvas.
- Chromatin partial reruns enter peak QC and pre-DE analysis through durable peak and annotation BED files instead of disposable completion markers.

### Notes
- Pre-DE deepTools profile plots may still place legends inside the plot panel; post-DE profile plots use the newer renderer.

## 0.5.0 - 2026-04-14

### Added
- Installable Python package layout under `src/` with CLI entry point and packaged workflow assets.
- Environment specifications for reproducible installs, including dedicated macOS test environment support.
- Site-config resolution with user-level config discovery and packaged fallback.
- Genome helper commands (`omnomnomics genomes ...`) for assembly discovery and installation workflows.
- Track color helper subcommands (`create-track-color-table`, `display-track-color-table`).
- Controller-job submission model for HPC orchestration, with controller script generation and controller logging.
- Monitor mode with step-state summary and recent log tail display.
- Retention policy support (`all`, `pruned`, `minimal`) with post-run cleanup behavior.
- Soft max project size guard with stage-aware cleanup and optional skipping of space-heavy branches.
- Flow-QC metric caching to preserve downstream reporting when cleanup removes intermediate files.
- Step-level tracking artifacts (started/completed/failed markers, per-step summaries, step command/notes logs).
- Aggregate alignment QC outputs and MultiQC custom content generation.
- Metadata-driven sample derivation for naming, grouping, and coloring:
  - `sample_id`
  - `sample_type`
  - `sample_color`
- Metadata-derived table output per run (`metadata_derived.tsv`) in run config output.
- Differential-expression design preflight for RNA step 12:
  - auto-build mode from selected metadata columns
  - explicit formula mode
  - optional interaction-term toggle for auto-built designs
  - full-rank design validation via `Rscript`/DESeq2-compatible model matrix checks

### Changed
- Pipeline internals moved from legacy repository-root script assumptions to packaged path resolution.
- Configuration responsibilities clarified across workflow defaults, site config, and per-run resolved config.
- BigWig generation modernized to `bamCoverage`-based outputs, including stranded RNA plus/minus tracks.
- Trackhub construction decoupled from HOMER defaults and driven directly from BigWig outputs.
- Count-table generation moved to modern assay-specific paths (`featureCounts` for RNA; assay-aware count logic for downstream modes).
- Step selection and resume behavior improved with clearer mode filtering and selected-step bookkeeping.
- Logging improved for command provenance, tool version capture, and storage-guard decisions.
- Metadata requirements standardized:
  - first metadata column is `filename`
  - reserved derived columns are blocked in user metadata
  - selectors accept named columns or 1-based indices
- Filename normalization standardized across CLI preflight and workflow runtime so metadata `filename` values can map from FASTQ/BAM/BigWig-style names to internal sample keys.

### Fixed
- Multiple same-day run logging and backup handling issues for run and tools logs.
- Snakemake submission pacing and status-refresh behavior under HPC scheduler load.
- Merge and parsing bugs in sample-name normalization and grouping separators.
- Step 10/11 and related DAG edge cases that could block valid execution paths.
- RNA-specific track coloring and grouped hub generation edge cases.
- Legacy sample-matching mismatch in step 12 by supporting both:
  - derived `sample_id` headers
  - normalized filename-key style headers from existing count tables
- `--dry-run` behavior to prevent unintended controller `sbatch` submission and run Snakemake dry-run directly.
- Step 9 failure mode when no `.bw` inputs are present by skipping trackhub creation gracefully instead of aborting the whole run.

### Notes
- Optional HOMER tag directory export remains available but is no longer required on the default path.
- Planned next hardening area: resumable recovery after partial cleanup/failure when later-stage artifacts exist but earlier sentinels/intermediates are missing.
