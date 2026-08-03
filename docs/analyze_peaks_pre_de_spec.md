# analyze_peaks (pre-DE) specification

## Scope
- Assay: ATAC and ChIP chromatin modes.
- Position in flow: after peak QC and before count table generation in the ATAC/ChIP public step flow.
- Purpose: create reusable peak, region, or bin set structure plus signal-level summaries independent of DE significance.

## Inputs
- `peak_calling/*.bed` from peak calling or feature generation.
- `filtered_BAM/*.sorted.dups_marked.filtered.bam` and matching `.bai` for ATAC.
- `filtered_BAM/*.filtered.bam` and matching `.bai` for ChIP.
- Derived metadata (for group labels via `sample_type` / `sample_id`).

## Outputs
- `peak_calling/analyze_peaks/sets/group_unions/*.union.bed` for peak-like ATAC and narrow/domain ChIP modes.
- `peak_calling/analyze_peaks/intersections/unique/*.bed` for peak-like modes with multiple groups.
- `peak_calling/analyze_peaks/intersections/shared/*.bed` for peak-like modes with multiple groups.
- gene-body or diffuse-bin metadata and summaries for ChIP `--broad-mode genebody` and `--broad-mode diffuse`.
- `peak_calling/analyze_peaks/summary/set_manifest.tsv`
- `peak_calling/analyze_peaks/summary/overlap_matrix.tsv`
- `peak_calling/analyze_peaks/signal/matrices/*.matrix.gz`
- `peak_calling/analyze_peaks/signal/heatmaps/*.pdf`
- `peak_calling/analyze_peaks/signal/profiles/*.pdf`

## Processing blocks
1. Build group unions:
   - merge per-group peak beds into one sorted merged union bed per biological group.
2. Build pairwise intersections:
   - shared peaks (`intersect -u`) and unique peaks (`intersect -v`) per group pair.
   - overlap table with shared/unique counts and Jaccard index.
3. ChIP feature-mode summaries:
   - for `--broad-mode genebody`, summarize annotation-derived gene-body features.
   - for `--broad-mode diffuse`, summarize filtered fixed genomic bins.
4. Signal visualization with deepTools:
   - create temporary sample bigWigs from BAM (`bamCoverage`, CPM normalization).
   - `computeMatrix` on union sets (reference-point mode, center ±3 kb by default).
   - `plotHeatmap` and `plotProfile` per group set.

## Guardrails
- If assay/mode has no meaningful pre-DE peak, region, or bin summaries, the step writes a sentinel and logs the skip reason.
- If fewer than two peak groups are detected for intersection-style summaries, intersection outputs are skipped with an explicit reason.
- Required executables: `bedtools`, `bamCoverage`, `computeMatrix`, `plotHeatmap`, `plotProfile`.

## Defaults
- `bamCoverage`: `--binSize 25 --normalizeUsing CPM`
- `computeMatrix`: `reference-point`, `center`, `-b 3000`, `-a 3000`, `--binSize 50`, `--skipZeros`
- Heatmap sorting: descending mean.

## Notes
- This step is intentionally pre-DE and does not filter by differential significance.
- DE-aware peak subset analysis is specified separately in `analyze_peaks_de_spec.md`.
- Pre-DE profile plots are still generated through the deepTools `plotProfile` path; the custom matplotlib profile renderer is used for post-DE signal profiles.
