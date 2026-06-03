# analyze_peaks (pre-DE) specification

## Scope
- Assay: ATAC only.
- Position in flow: after peak QC and before count table generation in ATAC public step flow.
- Purpose: create reusable peak set structure and signal-level summaries independent of DE significance.

## Inputs
- `peak_calling/*.bed` from peak calling.
- `filtered_BAM/*.sorted.dups_marked.filtered.bam` and matching `.bai`.
- Derived metadata (for group labels via `sample_type` / `sample_id`).

## Outputs
- `peak_calling/analyze_peaks/sets/group_unions/*.union.bed`
- `peak_calling/analyze_peaks/intersections/unique/*.bed`
- `peak_calling/analyze_peaks/intersections/shared/*.bed`
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
3. Signal visualization with deepTools:
   - create temporary sample bigWigs from BAM (`bamCoverage`, CPM normalization).
   - `computeMatrix` on union sets (reference-point mode, center ±3 kb by default).
   - `plotHeatmap` and `plotProfile` per group set.

## Guardrails
- If assay is not ATAC, step is skipped and sentinel is written.
- If fewer than two peak groups are detected, step writes sentinel and logs skip reason.
- Required executables: `bedtools`, `bamCoverage`, `computeMatrix`, `plotHeatmap`, `plotProfile`.

## Defaults
- `bamCoverage`: `--binSize 25 --normalizeUsing CPM`
- `computeMatrix`: `reference-point`, `center`, `-b 3000`, `-a 3000`, `--binSize 50`, `--skipZeros`
- Heatmap sorting: descending mean.

## Notes
- This step is intentionally pre-DE and does not filter by differential significance.
- DE-aware peak subset analysis is specified separately in `analyze_peaks_de_spec.md`.
