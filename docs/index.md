# omnomnomics

**omnomnomics** is an installable Snakemake application for modular RNA-seq,
ATAC-seq, and ChIP-seq processing on Slurm-based HPC systems. It supports
complete FASTQ-to-analysis workflows and selective stage ranges for resuming or
reusing existing data.

The pipeline provides:

- FASTQ trimming, FastQC, alignment, lane merging, and assay-aware BAM filtering
- alignment and experiment-level quality-control reports
- stranded RNA BigWigs and grouped UCSC track hubs
- RNA count tables and DESeq2 differential-expression analysis
- replicate-aware ATAC/ChIP peak calling, peak QC, differential chromatin
  analysis, motif analysis, and signal summaries
- auditable run configuration, tool-version, command, scheduler, and stage logs
- retention policies and a soft project-size guard for quota-constrained HPC work

## Start here

1. [Install omnomnomics](getting-started/installation.md) and configure the HPC site.
2. Prepare an experiment directory with a `FASTQ/` folder.
3. Read the [quickstart](getting-started/quickstart.md) and the page for your assay.
4. Run a dry run before dispatching a full analysis.

```bash
omnomnomics rna \
  -i /path/to/experiment \
  -g GRCh38 \
  --dry-run
```

Remove `--dry-run` to submit the controller job. The controller runs Snakemake
on one Slurm allocation and dispatches worker jobs for individual rules. You do
not need to wrap the `omnomnomics` command in `sbatch`.

## Choose an assay

| Assay | Command | Principal terminal outputs |
| --- | --- | --- |
| RNA-seq | `omnomnomics rna` | count tables, DE results, stranded BigWigs, track hubs |
| ATAC-seq | `omnomnomics atac` | consensus peaks, peak QC, counts, differential regions |
| ChIP-seq | `omnomnomics chip` | narrow or broad features, QC, counts, differential regions |

Use `omnomnomics --version` to report the installed release and
`omnomnomics <assay> --help` for the authoritative options of that release.
