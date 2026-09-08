# Workflow stages

Steps 1 through 9 are shared preprocessing and signal-track stages.

| Step | Stage | Result |
| --- | --- | --- |
| 1 | Trim reads | `trimmed_FASTQ/` |
| 2 | FastQC | `fastqc_reports/` |
| 3 | Align reads | lane-level files in `BAM/` |
| 4 | Merge technical lanes | sample-level BAMs in `BAM/` |
| 5 | Touch up BAMs | filtered BAMs in `filtered_BAM/` |
| 6 | Index BAMs | `.bai` files |
| 7 | Alignment QC | tabular and PDF/SVG summaries |
| 8 | Create BigWigs | `BigWigs/` |
| 9 | Build grouped track hubs | `merged_hubs/` |

Step 5 runs `samtools collate`, `fixmate`, `sort`, `markdup`, and `view`.
Duplicates are kept by default for RNA and removed by default for ATAC and ChIP.
Use `--keep-duplicates` or `--remove-duplicates` to override the assay default.

## RNA public stages

| Step | Stage |
| --- | --- |
| 10 | Track-hub completion boundary |
| 11 | Create count table |
| 12 | Differential expression |

## ATAC and ChIP public stages

| Step | Stage |
| --- | --- |
| 10 | Call peaks or construct assay features |
| 11 | Peak QC |
| 12 | Pre-DE peak analysis |
| 13 | Create count table |
| 14 | Differential chromatin analysis |
| 15 | Post-DE peak interpretation |

Public stage numbers are stable user interfaces. They are translated to the
appropriate internal Snakemake rules for each assay.

## Selecting stages

`-j auto` and `-j all` select the full public assay workflow. Ranges and comma-
separated selections are accepted:

```bash
# Map only
omnomnomics rna -i EXPERIMENT -g GRCh38 -j 1-3

# Resume from lane merging through RNA counting
omnomnomics rna -i EXPERIMENT -g GRCh38 -j 4-11

# Skip track-hub generation
omnomnomics rna -i EXPERIMENT -g GRCh38 -j 1-8,11-12
```

Every selected step must either receive an existing valid input or have that
input produced by another selected/dependency step. Prefer contiguous ranges
unless you understand the relevant file contracts.
