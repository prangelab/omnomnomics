# analyze_peaks_de (post-DE) specification

## Scope
- Assay: ATAC and ChIP chromatin DE modes.
- Position in flow: after DE/count-model outputs are available.
- Purpose: run peak-set interpretation on DE-significant subsets, top-ranked DE fallback subsets, and optional pre-DE unique peaks.

## Inputs
- Peak-level quantification and DE results tables.
- Pre-DE set outputs from `analyze_peaks` (union/shared/unique beds).
- Peak annotation beds from peak QC (`*.annotated.bed`) for promoter/distal split.
- Existing per-sample BigWig signal tracks from the BigWig creation stage.

## Outputs
- `peak_calling/analyze_peaks_de/sets/`
  - `de_significant/*.bed`
  - `de_up/*.bed`
  - `de_down/*.bed`
  - `de_significant_promoter/*.bed`
  - `de_significant_distal/*.bed`
  - `top_de_ranked*.bed`
  - `unique_from_prede/*.bed` (optional extra layer from step 14)
- `peak_calling/analyze_peaks_de/signal/`
  - deepTools matrices and heatmaps for DE subsets.
  - matplotlib-rendered profile PDFs built from the deepTools matrices.
  - `signal_runs.tsv` with per-set `OK`, `REUSE`, or `SKIP` status.
- `peak_calling/analyze_peaks_de/motifs/`
  - promoter, distal, and top-ranked motif runs per DE subset.
  - optional motif runs for pre-DE unique sets.
  - motif summary table.
- `peak_calling/analyze_peaks_de/summary/`
  - `set_manifest.tsv`
  - `analyze_peaks_de_report.tsv`

## Processing blocks
1. DE subset creation:
   - read contrast-level DE peak tables.
   - generate combined/up/down DE sets.
   - generate top-ranked DE fallback sets from the full DE table where available.
2. Context split:
   - split DE sets into promoter and distal subsets using peak metadata.
3. Optional unique layer:
   - import pre-DE unique peak BED files from `analyze_peaks/intersections/unique`.
4. Signal summaries:
   - run deepTools on sufficiently sized sets using existing BigWigs.
   - with `--post-de-signal-policy auto`, schedule missing BigWigs through step 8 before plotting.
   - with `--post-de-signal-policy require`, fail clearly if required BigWigs are absent.
   - with `--post-de-signal-policy skip`, disable post-DE signal plotting and record explicit summary-table skips.
   - render heatmaps with deepTools `plotHeatmap`.
   - render profile PDFs from the deepTools matrix with matplotlib, with the sample legend placed outside the signal axes.
   - reuse equivalent signal outputs when another set has the same capped BED content, BigWig list, sample labels, region label, and matrix settings.
   - mark copied signal outputs as `REUSE` in `signal_runs.tsv`, including the source set name and any region-cap reason.
5. Motif analysis:
   - run MEME SEA on sufficiently sized promoter/distal DE sets and top-ranked fallback sets.
   - run MEME AME on top-ranked DE sets using rank scores derived from BED order.
   - write MEME text-mode TSV outputs to avoid HTML-template failures on minimal HPC MEME installations.
   - optionally run on pre-DE unique sets.
   - prioritize promoter/promoter-region sets before distal, unique, and all-peak fallback sets.
   - cap motif set count, peaks per set, centered FASTA window size, per-set runtime, and motif threads via workflow defaults so motif failures remain report-level outcomes.

## Config notes
- Current implementation runs with practical defaults and tool-availability checks.
- Motif defaults are `post_de_motif_max_sets: 6`, `post_de_motif_max_peaks: 100`, `post_de_motif_window_bp: 200`, `post_de_motif_timeout_seconds: 1200`, `post_de_motif_threads: 1`, and `post_de_motif_database: auto`.
- `post_de_motif_database: auto` first uses the permanent MEME-format motif database cache under the configured genome assembly root, then falls back to the active environment or `OMNOMNOMICS_MEME_MOTIF_DATABASE`, `MEME_MOTIF_DATABASE`, or `JASPAR_MOTIF_DATABASE`. If no database is available, omnomnomics downloads the default JASPAR vertebrate database into `motif_databases/` for reuse. Run `omnomnomics genomes motifs` to prefetch or refresh this cache, or set `post_de_motif_database` to an explicit `.meme` motif file for reproducible database pinning.
- Fine-grained motif background matching controls are reserved for a follow-up refinement.

## Status
- Implemented as workflow step 16 (`analyze_peaks_de`) for ATAC and ChIP.
- Fail-soft behavior:
  - DE set generation and manifests are always attempted.
  - deepTools and MEME motif blocks are skipped with clear logs if executables are unavailable.
  - Post-DE signal plotting does not regenerate BigWigs inside the post-DE worker; missing BigWigs are either scheduled through step 8 or required according to `--post-de-signal-policy`, while `skip` disables signal plotting regardless of BigWig availability.
  - Post-DE signal plotting records per-set outcomes in `signal_runs.tsv`: `OK` for newly computed signal artifacts, `REUSE` for copied equivalent artifacts, and `SKIP` for ineligible or unavailable signal sets.
  - Motif runs write progress incrementally and record per-set `TIMEOUT`, `FAIL`, `NO_MOTIFS`, `SKIP`, or `OK` statuses.
