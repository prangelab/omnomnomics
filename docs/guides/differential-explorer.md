# Differential Explorer

The Differential Explorer is a packaged Shiny application for browsing completed
RNA, peak, gene-body, and genomic-bin differential result trees.

From the main environment:

```bash
omnomnomics de-app --project-dir /path/to/project_or_DE_calling
```

From the lightweight workstation environment:

```bash
micromamba activate omnomnomics-explorer
omnomnomics-de-app --project-dir /path/to/project_or_DE_calling
```

The launcher imports only Python standard-library modules and starts the packaged
Shiny application. The Explorer environment intentionally excludes Snakemake,
Slurm, aligners, and peak callers.

Keep the Explorer environment and application code on the same omnomnomics
release as the result project to avoid schema mismatches.
