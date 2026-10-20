# 0.6.0 release validation

Started 16 October and completed 19 October 2026 on `input_dev`. All planned
software, reference, full tutorial and actual recovery checks passed on exact
`54eeba92fa4b2417b451307783d4492fc36f4c01`. The author published
[v0.6.0](https://github.com/prangelab/omnomnomics/releases/tag/v0.6.0) on
19 October. Tag identity and all five downloaded release assets were verified
against the accepted source and prepared checksums. Earlier attempts below
retain their dated context.

## Candidate and environment provenance

- Established analysis source: `5bffbd83d5c3e979d7ec69edc8637663223c8cd1`.
- Initial exact-lock commit: `3b84a7b3acccae9ed186c0aa587ded0566b9a722`.
- Initial artifact-test source: `04f4c53a75856f430c053b7f1b2ff690165942bf`.
- Current accepted installation/artifact candidate: `54eeba92fa4b2417b451307783d4492fc36f4c01`.
- Previous accepted candidate: `96abebe9ff4bf978bbcf9ffb305fe37f4a340daa`.
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
| Exact dependency solves | Initial job 28136688 and corrected constrained job 28143746 succeeded |
| New prefixes and fresh package caches | Corrected main and Explorer passed exact inventories, R namespace loads and pip checks in array 28144176; unchanged IDR/SPP verified and executed |
| Existing analysis environments retained | Creation uses separate scratch prefixes and cache roots; no updates to the established environments |
| Installed wheel assets and CLI | Passed job 28280899, including byte-for-byte verification of 3,998 bundled workflow assets |
| Source-distribution installation | Passed job 28280899 in a separate Python venv inheriting only the fresh main dependency environment |
| R namespaces, helper wrappers and Explorer loading | Passed job 28280899; IDR 2.0.4.2, SPP 1.16.0, main R, Shiny discovery and lightweight Explorer verified |
| Complete unittest and standalone regression runners | Passed job 28280899: 116 unittest cases, none skipped, and 49 standalone functions with zero failures |
| Small real-tool workflow checks | Bioinformatics unittest tier passed within job 28280899 |
| Three full written tutorial routes | Passed all three complete runs and output audit 28282108; actual selective recovery and unchanged count/track hashes passed28288688 |
| Public release/tag/assets | v0.6.0 published 19 October; exact tested tag commit and all five asset hashes verified; separate Zenodo archive/DOI remains open |

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

## Accepted artifacts and next gate

Corrected artifact/regression job 28153260 completed successfully on 96abebe. The
unchanged verified dependency locks were reused; both installed artifacts were
rebuilt and the complete suites rerun after the GTF compatibility fix. Wheel SHA-256:
`a12e66540dbcae383d0099a77971b5c12a7ce78ebfe1823414d60c59636f566c`. Source-distribution SHA-256:
`ca870f6ba3034aecd311864ffd3b5efea8a44c22d2ef450b4a43d1f29f33bd70`. Both artifacts and compact reports
are archived locally under `MANUSCRIPT/release_validation/20261016/attempt_96abebe/`.

Reference preparation first stopped before any reference download because the
preparatory tutorial audit treated a normalized contrast mapping as a list.
The audit now verifies its factor, numerator and denominator explicitly. This
was a tutorial-audit script correction, not a change to the accepted package.
Reference/command-audit acceptance and the three full tutorial routes remain
required before the final release.

## UCSC annotation compatibility gate

Real reference preparation 28152153 found a source defect after the successful
software gates above: genomepy's UCSC mm39 annotation contains transcript/exon
records with `gene_id`, but no explicit `gene` records. Its provider README
identifies the actual annotation table as `ncbiRefSeq`, despite the requested
`refGene` style. The installer previously rejected this valid GTF and would
produce empty gene BEDs from it. The three dependent tutorial launchers
28152268/28152269/28152271 were canceled without starting controllers.

The correction accepts gene-associated transcript/exon GTFs and derives gene
BED bounds from their gene_id spans, separately for each contig and strand.
Existing explicit gene features take precedence. The original GTF, identifiers
and chromosome names are retained. Reference-size, primary-contig, coordinate
and missing-gene-identifier checks remain enforced. This source change passed
new wheel/sdist installation checks, 112 unittests and 49 standalone regressions
in job 28153260. Real reference retry 28153327 then successfully produced
nonempty gene BEDs and all HISAT2 parts; its final acceptance check stopped
because genomepy has no native mm39 blacklist entry.

## Native mm39 exclusion reference

Genomepy 0.16.3's blacklist plugin covers mm10 but has no mm39 entry; repeating
its download cannot supply an mm39 blacklist. The release tutorial therefore
uses the native mm39 exclusion set from
[Ogata et al., excluderanges](https://doi.org/10.1093/bioinformatics/btad198),
distributed by Mikhail Dozmorov through immutable
[Zenodo record 21480958](https://doi.org/10.5281/zenodo.21480958).
This is a separate reference dataset, not an ENCODE mm10 list or a software
dependency. No dependency installation or package-code change is required.

The 21,010-byte GRanges RDS has upstream MD5
`31fb7836575fd41f4f5e213f2fbebf1e` and independently recorded SHA-256
`ba7c625d422e4930b589cbff7b6a09fd6cd48058e9eea71b310c6031c23eff11`.
The tutorial helper verifies both before reading, asserts the mm39 assembly and
3,147 regions, and exports one-based inclusive GRanges as zero-based half-open
BED3 using start minus one and unchanged end. Existing chromosome harmonization
then checks the converted coordinates against this FASTA. Full source metadata
and the conversion rule are frozen in the tutorial reference manifest.

Reference completion reuses the successfully installed genome, annotation and
index. Earlier failed reference checks and canceled launchers remain recorded;
none count as completed tutorial runs. Historical SoftwareX blacklist provenance
is distinct from this newly pinned tutorial resource.


Reference completion 28157855 passed: 61 reference contigs, complete inspected
HISAT2 index, nonempty annotation BEDs, 3,147 normalized exclusion intervals,
JASPAR2024 motifs and 20 SHA-256 receipts. Lightweight reference reports and
source-preparation code/checksums are archived locally. No package dependencies
or pipeline code changed during this reference-data completion.

| Tutorial | Launcher | Controller | Collector | Confirmed public steps |
|---|---|---|---|---|
| RNA | 28157969 | 28157972 | 28157973 | 1–12 |
| Matched-input ChIP | 28157970 | 28158003 | 28158018 | 1–15; 8 libraries and 8 associations |
| ATAC | 28157971 | 28158005 | 28158017 | 1–15 |

These fresh full tutorial runs and all collectors completed successfully.
Independent output audit 28178821 accepted all three output sets: 24 source
FASTQ SHA-256 values remained unchanged, BAM headers/indexes and BigWigs were
readable against mm39, source protection and pruned retention were correct,
integer count tables and DE/post-DE outputs were present, and all worker
accounting records reported successful completion. Filtered peak sets had zero
remaining overlaps with the pinned exclusion reference. ChIP had eight matched
associations, four separately documented input-relative tracks and raw IP counts
for DE without input subtraction or division.

| Route | Full controller duration | Count-table features | Features with adjusted p value | Significant features |
|---|---|---|---|---|
| RNA | 00:41:20 | 39,998 | 19,235 | 3,013 |
| ChIP | 01:28:13 | 104,201 | 54,845 | 18,053 |
| ATAC | 02:19:52 | 29,701 | 29,701 | 2,847 |

Significance uses adjusted p < 0.05 and absolute log2 fold change >= 1.
For ChIP, 81,799 features had a p value; independent filtering left adjusted
p values for 54,845. These categories must not be conflated. Full controller
durations include scheduler waiting during execution. They are distinct from
the log's internal elapsed time, CPU time and billed core hours. Actual Slurm
allocation is retained in compact resource reports, including allocations above
worker thread counts.

The same audit job then stopped at the RNA 5–12 recovery example: step 5 needs
unfiltered BAM inputs already removed by successful `pruned` retention. The
output audit itself passed. A separate dry-run job 28275680 checks RNA 6–12
from retained filtered BAMs and chromatin 11–15. Actual recovery execution is
still pending; no selective duration is accepted as full-run performance.
Compact full-run evidence is archived in
`MANUSCRIPT/release_validation/20261016/tutorial_acceptance/`.

## Recovery defect found after full output acceptance

Retained-input dry-run job 28275680 completed for all three routes. ChIP/ATAC
11–15 reused all outputs. RNA 6–12 scheduled four coverage jobs solely because
step-8 completion markers had been deleted at the end of the full run, despite
both strand tracks being present. This is a real reuse defect and blocks final
release acceptance; the earlier full-run output audit remains valid for 96abebe.

The correction makes coverage targets and track-hub dependencies use the real
BigWig files when no size cap is configured, declares both RNA strand outputs,
and retains marker-based completion for optional size-capped branches. Recovery
must reuse intact tracks after marker cleanup and rebuild a missing strand.
Two real Snakemake regression cases cover RNA and ATAC, including an actual
small RNA track rebuild. A new exact artifact candidate must pass the complete
software gates and installed recovery checks before release. Dependencies and
analysis algorithms are unchanged.


## Earlier e473f1c software acceptance, 19 October

Exact candidate `e473f1cacf9ade8db74734ea6ba7ca5ed9d750e1` passed artifact job
28275721: installed wheel and sdist outside the checkout, byte verification of
3,998 workflow assets per installation, CLI/R/helper/Explorer/pip checks,
114 unittests without skips, 49 standalone regressions without failures and
Shiny discovery. Both real coverage-recovery regressions also passed separately
in job 28275719. No dependencies changed.

Wheel SHA-256:
`ba69be5b0374efa36be513fc9dd078d0cefcab360d6e9ad6522d4137b84279c6`.
Source-distribution SHA-256:
`1c3c28aa43f6b1304821b711e809033f788768ad6efcead9a1029d743517f0f7`.
Both artifacts and compact reports are archived locally at
`MANUSCRIPT/release_validation/20261016/attempt_e473f1c/`; copied artifacts were
independently checked against the SHA-256 receipts. Earlier accepted artifacts
and failures remain separate.

| Fresh tutorial | Launcher | Controller | Collector | Confirmed public steps |
|---|---|---|---|---|
| RNA | 28275724 | 28276659 | 28276660 | 1–12 |
| Matched-input ChIP | 28275725 | 28276662 | 28276663 | 1–15; 8 libraries, 8 associations |
| ATAC | 28275726 | 28276661 | 28276664 | 1–15 |

All three controllers are running; launch checks and DAGs were inspected.
Output audit 28277242 waits for all three resource collectors. Actual selective
recovery is scheduled after audit acceptance, with explicit entry ranges,
unchanged count/track SHA-256 checks and zero new analysis workers required.
Neither these pending runs nor recovery duration is accepted as final full-run
performance. Original 96abebe measurements above remain identified as historical
candidate evidence. Final release/tag/publication remains open.


## Selective retention defect found by actual recovery

All three e473f1c full tutorials and independent output audit 28277242 passed.
Actual recovery controllers 28280729/28280736/28280741 completed with no analysis
workers and unchanged raw count tables. Verification 28280742 nevertheless
failed: pruned retention deleted existing ChIP/ATAC BigWigs and track hubs when
steps 8–9 were absent from the selective range. Source FASTQs remain intact.
The original accepted full-run evidence is archived separately before recovery.
This is a real retention defect; final acceptance remains open.

The correction preserves existing reusable analysis branches under `pruned`,
independent of the current selection: filtered BAMs, BigWigs, track hubs, peak
and DE outputs, and HOMER tag directories. Heavy early intermediates remain
eligible for cleanup. `minimal` retains its explicit requested-terminal-output
behavior and source protection. Two regressions verify preserved branches with
pruned intermediate cleanup and unchanged minimal behavior. Exact-candidate
software and full tutorial/recovery gates will be rerun without dependency
changes. The failed recovery evidence is retained and no check is suppressed.


## Retention-fix software acceptance, 19 October

Exact candidate `54eeba92fa4b2417b451307783d4492fc36f4c01` passed artifact job
28280899: 116 unittests without skips, 49 standalone regressions without
failures, wheel/sdist installations outside the checkout, byte verification of
3,998 workflow assets per installation, CLI/R/helpers/Explorer/pip checks and
Shiny discovery. All 19 targeted source-protection/retention cases passed in
28280839. Dependencies and their verified prefixes are unchanged.

Wheel SHA-256:
`ba5db4c8bb811b2fd100a3da93805dcdd91867d584d9cfd3e6435018f252d837`.
Source-distribution SHA-256:
`c1e1862e8d0484783474910f8cf1c77b7e0e5f89e160cd296c363c45f636e280`.
Copied artifacts and reports under
`MANUSCRIPT/release_validation/20261016/attempt_54eeba9/` were independently
verified against their SHA-256 receipts.

| Fresh tutorial | Launcher | Controller | Collector | Confirmed public steps |
|---|---|---|---|---|
| RNA | 28280910 | 28281924 | 28281925 | 1–12 |
| Matched-input ChIP | 28280913 | 28281927 | 28281928 | 1–15; 8 libraries, 8 associations |
| ATAC | 28280914 | 28281926 | 28281930 | 1–15 |

All three full controllers and collectors completed. Output audit 28282108
passed every route. Actual recovery controllers 28288656/28288658/28288687 and
verification 28288688 passed: zero analysis workers, three unchanged raw count
tables, all 20 unchanged ordinary/input-relative tracks, all 24 source FASTQs and
existing track hubs retained. Earlier failed evidence and checks remain intact.

| Accepted full tutorial | Controller duration | Significant features |
|---|---|---|
| RNA | 00:42:35 | 3,013 |
| Matched-input ChIP | 01:27:25 | 18,053 |
| ATAC | 02:14:43 | 2,847 |

These are complete-controller tutorial measurements, distinct from recovery
and historical SoftwareX performance. Resource reports preserve allocation,
billing, CPU and memory separately. Compact final evidence is archived locally
in tutorial_acceptance_54eeba9 and tutorial_recovery_54eeba9, with the complete
acceptance record at MANUSCRIPT/release_validation/20261016/FINAL_ACCEPTANCE.md.
Validation completed before publication. The author subsequently published
v0.6.0; the identity verification is recorded below.

## Published release verification, 19 October

The public [v0.6.0 release](https://github.com/prangelab/omnomnomics/releases/tag/v0.6.0)
was published at `2026-10-19T15:32:59Z`. Its tag points directly to
`54eeba92fa4b2417b451307783d4492fc36f4c01`, the exact source accepted above.
The tag is lightweight and has no cryptographic tag signature; this check
establishes source and artifact identity, not signed publisher attestation.

All five uploaded assets were downloaded independently. Their sizes and
SHA-256 values match the prepared upload manifest and GitHub's asset digests:
the tested wheel and source distribution, exact Linux lock bundle, compact
validation JSON and SHA256SUMS. The source and artifacts were not rebuilt or
changed after validation, so no additional analysis run was needed.

The GitHub release is public, with `draft=false` and `prerelease=false`.
GitHub reports `immutable=false`; the pinned source commit and archived asset
hashes provide the identity record. A separate Zenodo version DOI is not yet
recorded. Local publication receipts are under
`MANUSCRIPT/release_validation/20261016/published_v0.6.0/`.
