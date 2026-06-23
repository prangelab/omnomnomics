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
  - deepTools matrices, heatmaps, and profiles for DE subsets.
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
   - with `--post-de-signal-policy skip`, skip signal plots with a summary-table reason if required BigWigs are absent.
5. Motif analysis:
   - run GimmeMotifs on sufficiently sized promoter/distal DE sets and top-ranked fallback sets.
   - optionally run on pre-DE unique sets.
   - prioritize promoter/promoter-region sets before distal, unique, and all-peak fallback sets.
   - cap motif set count, peaks per set, per-set runtime, and motif threads via workflow defaults so motif failures remain report-level outcomes.

## Config notes
- Current implementation runs with practical defaults and tool-availability checks.
- Motif defaults are `post_de_motif_max_sets: 6`, `post_de_motif_max_peaks: 100`, `post_de_motif_timeout_seconds: 1200`, and `post_de_motif_threads: 1`.
- Fine-grained motif background matching controls are reserved for a follow-up refinement.

## Status
- Implemented as workflow step 16 (`analyze_peaks_de`) for ATAC and ChIP.
- Fail-soft behavior:
  - DE set generation and manifests are always attempted.
  - deepTools and GimmeMotifs blocks are skipped with clear logs if executables are unavailable.
  - Post-DE signal plotting does not regenerate BigWigs inside the post-DE worker; missing BigWigs are either scheduled through step 8, required, or skipped according to `--post-de-signal-policy`.
  - Motif runs write progress incrementally and record per-set `TIMEOUT`, `FAIL`, `NO_MOTIFS`, `SKIP`, or `OK` statuses.
