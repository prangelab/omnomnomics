Omnomnomics 0.6.1 repairs forced chromatin DE/post-DE reruns and improves the
Differential Explorer. The Linux environment locks are unchanged from 0.6.0.

## Fixes

- Selecting and forcing chromatin DE and post-DE now executes both stages
  when previous outputs already exist. Normal selective reuse is preserved.
- PCA grouping supports metadata columns containing spaces and single groups.

## Differential Explorer

- Responsive layout with project/contrast sidebar, clearer controls and live
  result counters.
- Larger round PCA points with adjustable size.
- Settings for sample grouping, palette presets, custom gradients and group
  colours shared by PCA, boxplots, heatmap annotations and their PDF exports.
- Separate ORA up/down dotplots and GSEA signed-NES plots per collection,
  with readable pathway labels and every panel included in PDF export.
- Optional use of an already installed coloRspacer package for gradients;
  base R provides the fallback without adding a mandatory dependency.

## Validation and installation

See the attached validation JSON for exact-candidate artifact installation,
regression checks and hashes. The J14 repair was also exercised on Snellius,
and Explorer controls and exports were tested in a local browser.

The complete clean protocol revalidation against this published release is
planned separately. Do not interpret repaired or selective-run durations as
full-pipeline performance.

Review the bundled lock provenance before creating environments. Install the
wheel inside the matching environment with:

```bash
python -m pip install --no-deps --no-index /path/to/omnomnomics-0.6.1-py3-none-any.whl
```

Reference genomes, experimental data, Slurm configuration and external tools
remain separate from the Python package. The Linux lock bundle contains the
matching creation, verification and helper-linking scripts.
