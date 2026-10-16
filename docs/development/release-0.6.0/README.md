# 0.6.0 release validation

Started 16 October 2026 on `input_dev`. This is release-candidate validation,
not a declaration that 0.6.0 has passed every acceptance check.

## Candidate and environment provenance

- Established analysis source: `5bffbd83d5c3e979d7ec69edc8637663223c8cd1`.
- Initial exact-lock commit: `3b84a7b3acccae9ed186c0aa587ded0566b9a722`.
- Artifact-test source: `04f4c53a75856f430c053b7f1b2ff690165942bf`.
- Runtime work root: `/scratch-shared/kprange/omnomtest/release_060_20261016`.
- Inventory job: `28136339`; constrained dependency-solve job: `28136688`.
  Both completed. Micromamba version: 2.3.3.

The exact Linux x86-64 locks include 695 main, 140 IDR, 154 SPP and 143
Explorer Conda packages. Every artifact comes from the public Conda-forge or
Bioconda publishing community. Exact versions, builds, URLs and MD5/SHA-256
hashes are in `locks/linux-64/package_provenance.tsv`. No lock-generation
package was installed: project scripts serialize the checked Micromamba plans.

The initial dependency closure preserves the established versions/builds.
The solve also includes the declared but previously absent `bioconductor-sva`
3.58.0 (`r45h01b2380_0`, Bioconda) and its dependency closure. That adds
`bioconductor-edger` 4.8.2, `bioconductor-genefilter` 1.92.0,
`bioconductor-limma` 3.66.0 and `r-statmod` 1.5.2. It selects
`r-survival` 3.8_12 where this dependency was outside the original selected
closure. These changes must pass the new-environment tests; they do not prove
numerical equivalence with every previous analysis.

SVA is the long-established Bioconductor package maintained by Jeffrey Leek,
John Storey and Evan Johnson, distributed here through Bioconda. It is already
declared by `environment.yml`; no new analysis feature is being introduced.
See the [upstream package record](https://bioconductor.org/packages/sva/).

## Acceptance checks

| Check | Evidence/status |
|---|---|
| Exact dependency solves | All four succeeded, job 28136688 |
| New prefixes and fresh package caches | Array job 28137199; main running at last check; IDR, SPP and Explorer completed and passed exact inventory verification |
| Existing analysis environments retained | Creation uses separate scratch prefixes and cache roots; no updates to the established environments |
| Installed wheel assets and CLI | Job 28137628 scheduled after successful completion of the creation array |
| Source-distribution installation | Same job; separate Python venv inherits only the newly reconstructed main dependency environment |
| R namespaces, helper wrappers and Explorer loading | Same job; pending |
| Complete unittest and standalone regression runners | Same job; pending; temporary-path standalone tests are now executed rather than silently skipped |
| Small real-tool workflow checks | Included in the bioinformatics unittest tier; pending |
| Three full written tutorial routes | Not submitted; staged source FASTQs/reference files require restoration after scratch purging |
| Final release/tag/archive | Not created; depends on successful validation of the final source and locks |

The installed-package smoke check runs outside the checkout, verifies that
Python imports resolve into the new prefix, compares every bundled workflow
asset byte-for-byte against the wheel, checks version consistency, loads the
required R namespaces and Explorer code, and runs `pip check`. The full suites
include real HISAT2, MACS3, input-relative track and mixed-entry/retention tests.
Static source assertions and installed-code execution retain their distinct
roles; passing one is not a substitute for the other.

## Installation behavior changes

- `scripts/create_locked_environment.sh` refuses existing prefixes and verifies
  lock checksums plus the entire installed Conda package set.
- `scripts/link_locked_helpers.sh` binds separately locked IDR/SPP prefixes
  after checking their inventories.
- Named IDR/SPP installers consume the release locks and reject a mismatched
  existing helper environment instead of silently solving an update.
- The source distribution includes the environment specifications, locks,
  installation scripts and test sources through `MANIFEST.in`.

Large environments, package caches and tutorial analysis data remain on
scratch. Archive compact logs, checksums, manifests and reports locally before
the 14-day scratch lifetime expires. Selective test duration is not a
full-pipeline performance measurement.
