# 0.6.1 release validation

This patch release contains the forced chromatin DE/post-DE dependency repair
and the Differential Explorer interface, palette and enrichment-plot updates.
Published manually on 19 November 2026 as
[0.6.1](https://github.com/prangelab/omnomnomics/releases/tag/0.6.1), pointing
to `61d47d66b828837747c7e445563db8a54da4b5f9`. The published tag has no `v`
prefix. Its software and version metadata match tested commit
`b7f783b3fe0c41bdeb2cc4f2baef5765a79d061c`; the later commit changes only
release instructions. Published tags and assets remain unchanged.

## Validation scope

- J14 passed on repair commit `2bf1c0f8ced97201f3fc22e94a26935111b24089`
  on Snellius: actual DE and post-DE workers ran, the log2-fold-change cutoff
  changed from 1 to 2, and 5,570 significant regions satisfied the new filters.
  Raw counts, all eight ChIP tracks and nine source FASTQs were unchanged.
- Explorer layout, filtering, palette controls, larger PCA points, separate
  ORA/GSEA panels and six-page enrichment PDF export were tested locally.
  The optional coloRspacer 0.1.0 integration and base R fallback both passed.
- Final wheel/sdist and standalone Explorer installations passed, with 4,002
  byte-verified workflow assets per installation. All 122 unittests passed
  without skips and all 49 standalone regressions passed. The accompanying
  `omnomnomics-0.6.1-validation.json` records exact identities and scope.
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

## Publication record

The author published the release and edited its public notes to retain only
user-facing changes. `RELEASE_NOTES.md` mirrors that published description.
Development and validation details remain in this record.

Five uploaded assets contain the wheel, source distribution, Linux lock bundle,
validation JSON and `SHA256SUMS`. The tag and asset identities are checked
separately from the planned complete clean protocol run.
