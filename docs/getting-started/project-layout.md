# Project layout

Omnomnomics keeps inputs, intermediate files, results, and provenance in one
experiment directory.

```text
EXPERIMENT_DIR/
├── FASTQ/
├── trimmed_FASTQ/
├── fastqc_reports/
├── BAM/
├── filtered_BAM/
├── BigWigs/
├── merged_hubs/
├── peak_calling/
├── DE_calling/
├── MultiQC/
├── run_configs/
├── run_logs/
└── slurm_logs/
```

Not every assay or selected stage range creates every directory.

| Directory | Purpose |
| --- | --- |
| `FASTQ` | original compressed reads; retained by all retention policies |
| `trimmed_FASTQ` | adapter-trimmed reads and trim metrics |
| `BAM` | supplied aligned BAMs, aligner output, lane-merged BAMs, and mapper statistics |
| `filtered_BAM` | sorted, filtered BAMs, indexes, and BAM QC |
| `BigWigs` | per-sample browser signal tracks |
| `merged_hubs` | grouped UCSC track hubs |
| `peak_calling` | ATAC/ChIP feature definitions and pre-DE reports |
| `DE_calling` | count matrices and differential-analysis results |
| `run_configs` | fully resolved per-run configuration |
| `run_logs` | run provenance, stage summaries, and cached lightweight metrics |
| `slurm_logs` | controller and worker standard output grouped by rule |

The reference root is separate from experiment data. A normalized assembly has
`fasta/genome.fa`, `annotation/genes.gtf`, aligner index directories, and `aux/`.

For metadata-declared ChIP preparation, both experiment and input FASTQs belong
in `FASTQ/`; externally supplied BAMs belong in `BAM/`. The earliest available
source for each library is retained under every retention policy. Resolved
library and input-association manifests are written to `run_configs/`. See the
[metadata guide](../guides/metadata.md#chip-input-libraries).
