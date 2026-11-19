# 0.6.1 release preparation

This patch release contains the forced chromatin DE/post-DE dependency repair
and the Differential Explorer interface, palette and enrichment-plot updates.
Publication is manual. Do not move the v0.6.0 tag or replace its assets.

## Validation scope

- J14 passed on repair commit `2bf1c0f8ced97201f3fc22e94a26935111b24089`
  on Snellius: actual DE and post-DE workers ran, the log2-fold-change cutoff
  changed from 1 to 2, and 5,570 significant regions satisfied the new filters.
  Raw counts, all eight ChIP tracks and nine source FASTQs were unchanged.
- Explorer layout, filtering, palette controls, larger PCA points, separate
  ORA/GSEA panels and six-page enrichment PDF export were tested locally.
  The optional coloRspacer 0.1.0 integration and base R fallback both passed.
- Final versioned wheel/sdist installation and full regression acceptance
  are recorded in the accompanying `omnomnomics-0.6.1-validation.json` asset.
  Read that record for the exact tested commit, checks and artifact hashes.
- The full clean protocol run against the published 0.6.1 assets, including
  Data S1, fresh environments/references, all assay routes, recovery and J14,
  is a separate planned check after publication. Earlier repaired or v0.6.0
  runs do not establish that acceptance.

## Dependencies

The four Linux x86-64 locks and readable environment specifications are
unchanged from v0.6.0 (729 main, 140 IDR, 154 SPP and 146 Explorer packages).
No new mandatory dependency is introduced. coloRspacer is optional; custom
gradients also work with base R. Existing environments are preserved during
artifact testing by installing the project into separate validation venvs.

## Manual publication

1. Publish from the prepared `main` after its checks pass. Do not tag an older
   main commit. Release preparation is pushed directly to `main`; no pull
   request is required. The validation asset identifies the tested code and
   any later documentation-only commit without changing artifact identity.
2. Create a GitHub release named **Omnomnomics 0.6.1**, with new tag `v0.6.1`
   targeting the merged `main`. Do not publish it as a prerelease.
3. Use `RELEASE_NOTES.md` from this directory as the release description.
4. Attach the prepared wheel, source distribution, Linux lock bundle,
   validation JSON and `SHA256SUMS`. GitHub supplies its own source archives;
   these do not replace the tested wheel and source distribution.
5. Publish manually. Verify the resulting tag and downloaded asset hashes
   before updating stable-installation instructions or the protocol.

The citation metadata carries the prepared version. Add the actual publication
date after release; no publication date is asserted during preparation.
