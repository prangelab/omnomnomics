# Outputs

## Quality control

Omnomnomics writes per-sample BAM statistics and PDF/SVG summaries, an
experiment-level QC table and figures, and a MultiQC report unless disabled.
Flow summaries combine raw/trimmed read counts, mapper alignment statistics, and
post-filter BAM metrics. Cached lightweight metrics allow reporting after safe
intermediate cleanup.

## RNA

- `DE_calling/*.raw_read_quant.table.txt`: featureCounts matrix
- `DE_calling/<analysis>/`: differential tables, plots, enrichment, and scripts
- `BigWigs/*.plus.bw` and `*.minus.bw`: stranded signal
- `merged_hubs/*.hub/`: grouped UCSC track hubs

## ATAC and ChIP

- `peak_calling/`: peaks or assay-specific feature definitions
- `peak_calling/analyze_peaks/`: unions, intersections, annotations, and signal summaries
- `DE_calling/`: count matrix and differential chromatin results
- `DE_calling/analyze_peaks_de/`: significant sets, signal plots, motifs, and status tables

Important post-DE audit tables include:

- `summary/set_manifest.tsv`
- `summary/analyze_peaks_de_report.tsv`
- `signal/signal_runs.tsv`
- `motifs/motif_runs.tsv`

Statuses explicitly distinguish new results, reused results, intentional skips,
timeouts, failures, and valid no-enrichment outcomes.
