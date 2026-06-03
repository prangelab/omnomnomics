# analyze_peaks_de (post-DE) specification

## Scope
- Assay: ATAC (and future ChIP if DE-like quantification is enabled).
- Position in flow: after DE/count-model outputs are available.
- Purpose: run peak-set interpretation specifically on DE-significant subsets, with an optional layer for pre-DE unique peaks.

## Inputs
- Peak-level quantification and DE results tables.
- Pre-DE set outputs from `analyze_peaks` (union/shared/unique beds).
- Peak annotation beds from peak QC (`*.annotated.bed`) for promoter/distal split.
- BAM or BigWig signal tracks.

## Outputs
- `peak_calling/analyze_peaks_de/sets/`
  - `de_significant/*.bed`
  - `de_up/*.bed`
  - `de_down/*.bed`
  - `de_significant_promoter/*.bed`
  - `de_significant_distal/*.bed`
  - `unique_from_prede/*.bed` (optional extra layer from step 14)
- `peak_calling/analyze_peaks_de/signal/`
  - deepTools matrices, heatmaps, and profiles for DE subsets.
- `peak_calling/analyze_peaks_de/motifs/`
  - promoter and distal motif runs per DE subset.
  - optional motif runs for pre-DE unique sets.
  - motif summary table.
- `peak_calling/analyze_peaks_de/summary/`
  - `set_manifest.tsv`
  - `analyze_peaks_de_report.tsv`

## Processing blocks
1. DE subset creation:
   - read contrast-level DE peak tables.
   - generate combined/up/down DE sets.
2. Context split:
   - split DE sets into promoter and distal subsets using peak metadata.
3. Optional unique layer:
   - import pre-DE unique peak BED files from `analyze_peaks/intersections/unique`.
4. Signal summaries:
   - run deepTools on sufficiently sized sets.
5. Motif analysis:
   - run GimmeMotifs on sufficiently sized promoter/distal DE sets.
   - optionally run on pre-DE unique sets.

## Config notes
- Current implementation runs with practical defaults and tool-availability checks.
- Fine-grained motif background matching controls are reserved for a follow-up refinement.

## Status
- Implemented as workflow step 16 (`analyze_peaks_de`) for ATAC.
- Fail-soft behavior:
  - DE set generation and manifests are always attempted.
  - deepTools and GimmeMotifs blocks are skipped with clear logs if executables are unavailable.
