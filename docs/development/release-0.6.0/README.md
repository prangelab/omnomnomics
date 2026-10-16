# 0.6.0 release validation

Started 16 October 2026 on `input_dev`. This is release-candidate validation,
not a declaration that 0.6.0 has passed every acceptance check.

## Candidate and environment provenance

- Established analysis source: `5bffbd83d5c3e979d7ec69edc8637663223c8cd1`.
- Initial exact-lock commit: `3b84a7b3acccae9ed186c0aa587ded0566b9a722`.
- Initial artifact-test source: `04f4c53a75856f430c053b7f1b2ff690165942bf`.
- Runtime work root: `/scratch-shared/kprange/omnomtest/release_060_20261016`.
- Inventory job: `28136339`; constrained dependency-solve job: `28136688`.
  Both completed. Micromamba version: 2.3.3.

The corrected Linux x86-64 locks contain 729 main, 140 IDR, 154 SPP and 146
Explorer Conda packages (initially 695/140/154/143). IDR/SPP package identities
remain unchanged. In main, MultiQC is the only existing package whose version
or build changes; PyArrow and its previously absent closure are added, and the
deprecated Polars LTS transition package is removed. Every artifact comes from the public Conda-forge or
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
| New prefixes and fresh package caches | Array 28137199 completed all four exact Conda inventories, but main failed runtime acceptance; corrected clean main/Explorer reconstruction pending |
| Existing analysis environments retained | Creation uses separate scratch prefixes and cache roots; no updates to the established environments |
| Installed wheel assets and CLI | Job 28137628 built wheel/sdist and passed installed asset/CLI checks before failing the R runtime check; corrected artifact attempt pending |
| Source-distribution installation | Not reached in failed job 28137628; separate Python venv will inherit only the reconstructed main dependency environment |
| R namespaces, helper wrappers and Explorer loading | Helpers passed in 28137628; main R failed because annotation payload hooks ran before their dependencies; full repeat required |
| Complete unittest and standalone regression runners | Pending corrected artifact attempt; temporary-path standalone tests are executed rather than silently skipped |
| Small real-tool workflow checks | Included in the bioinformatics unittest tier; pending |
| Three full written tutorial routes | Not submitted; checksummed FASTQ restoration job 28137955 is running; reference restoration remains pending |
| Final release/tag/archive | Not created; depends on successful validation of the final source and locks |

The installed-package smoke check runs outside the checkout, verifies that
Python imports resolve into the new prefix, compares every bundled workflow
asset byte-for-byte against the wheel, checks version consistency, loads the
required R namespaces and Explorer code, and runs `pip check`. The full suites
include real HISAT2, MACS3, input-relative track and mixed-entry/retention tests.
Static source assertions and installed-code execution retain their distinct
roles; passing one is not a substitute for the other.

## Failures found by clean reconstruction

The first creation array returned success and matched all Conda records, but
artifact job 28137628 failed its R runtime check. Diagnosis job 28141650 traced
missing GO.db and organism annotation payloads to alphabetic serialization of
the lock: Bioconda post-link hooks ran before `yq` and R were linked. Locks now
preserve the solver's dependency order. Creation additionally loads the required
R namespaces and runs `pip check`, so a recorded wrapper package cannot hide a
missing annotation payload. The separate annotation sources/checksums are
recorded in `locks/bioconductor-data-payloads.json`.

The R smoke check also used an incorrect namespace spelling (`decoupler`);
the actual namespace is `decoupleR`. This check is corrected without changing
the installed package. Diagnosis identified unsatisfied upstream MultiQC 1.33
Python requirements (`pyarrow` and `polars-lts-cpu`); Conda's deprecated LTS
transition package does not supply the latter Python distribution.

The reviewed replacement is Bioconda `multiqc` 1.35, build `pyhdfd78af_2`,
maintained upstream by Phil Ewels and Vlad Savelyev. Its
[versioned dependency specification](https://github.com/MultiQC/MultiQC/blob/v1.35/pyproject.toml)
requires PyArrow and Polars with its compatibility runtime. The manifest now
declares PyArrow and retains Polars/runtime-compat 1.39.3. Changed package
archives and transitive dependencies are recorded in the provenance TSV.
Exploratory solve 28142306 selected Polars 2.0.0 and was not accepted. Solve
28142798 retained Polars 1.39.3 but allowed unnecessary transitive updates;
constrained solve 28143746 retains the existing dependency pins. Its first
submission 28143572 failed on an orchestration inventory path, before solving;
the corrected attempt uses the preserved version-3 inventory.

Explorer's lightweight environment is intentional. Package metadata now keeps
PyYAML in the core and moves genomepy, pysam and deepTools into the `pipeline`
extra. Explorer's manifest declares PyYAML explicitly. Fresh Explorer creation,
wheel installation, CLI, app loading and `pip check` must pass before acceptance.
Earlier reports and failed install logs are preserved; no failed attempt is
reported as a passed release gate.

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
