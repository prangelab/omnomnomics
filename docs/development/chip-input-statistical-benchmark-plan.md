# ChIP input benchmark plan

Date: 2026-09-14. Branch: `input_dev`. Status: all 24 simulated datasets,
72 analyses and the one practical validation complete; practical assertions
passed after bounded input quality selection. Revised scope: **24 simulated datasets**, followed by
one small real-data Snellius validation. This replaces the larger calibration
proposal. Pipeline defaults and the manuscript ledger remain unchanged.

## Purpose and fixed scope

Check the implemented ChIP analysis against known simulated changes, identify
failure modes, and demonstrate its normalization limits. This is focused
software validation, not a formal FDR-calibration study or a general comparison
of differential-binding methods.

Each dataset contains three WT and three KO ChIP libraries plus its controls.
Six scenarios use four independent seeds each: **24 datasets in total**. Three
comparison arms reuse those reads, giving 72 primary analyses. There is no
separate confirmation round or additional scenario grid. Expand only a scenario
that exposes a concrete unresolved problem, documenting the reason and proposed
extent before doing so.

## Six scenarios

| ID | Scenario | Perturbation and question | Seeds |
| --- | --- | --- | --- |
| B1 | No biological difference | Unchanged specific signal/background; inspect false calls | 0–3 |
| B2 | No difference, unequal depth | Same biology, KO exposure 4x WT; check depth handling | 0–3 |
| B3 | Local signal changes | 20% of sites change, half fourfold up and half fourfold down; assess discovery and correct-direction recovery | 0–3 |
| B4 | Signal in one condition only | Specific component lost in KO at 20% of sites; background remains; check condition-specific candidates | 0–3 |
| B5 | Changing input background | Specific binding unchanged; background doubles in KO in 20% of blocks; use condition-matched inputs; distinguish binding from total ChIP changes | 0–3 |
| B6 | Global signal loss | All specific components halve, background stays fixed, sequencing totals are fixed; demonstrate normalization-reference dependence | 0–3 |

Use stable scenario/library-specific random streams derived from a recorded
master seed. A repeated seed number in different scenarios does not share a
random stream. Freeze effects and coordinates before sampling reads.

## Simulation and analysis configuration

- Synthetic mapped BAMs: one 30 Mb contig, 1,000 predetermined signal sites
  with weak/medium/strong classes and 1,000 independent decoy intervals.
- Paired-end: 50 bp mates, proper pairs, insert lengths centred on 200 bp.
  Existing functional tests cover SE/mixed-layout correctness; these do not
  add another simulation grid.
- Generate fragment numbers using Gamma biological variation (mean one,
  dispersion 0.05) followed by Poisson sequencing. Signal/background mixtures
  need not share one NB dispersion. Use nonuniform spatial background and
  roughly 20/100/400 reference-depth counts across strength classes.
- Baseline expected depths: 400,000 fragments per ChIP and 600,000 per input.
  Record actual depths. Use one shared cell-type input except B5, which has
  condition-matched inputs. B6 uses exact sequencing budgets with multinomial
  allocation from perturbed molecular weights; label its model separately.
- Primary settings: explicit `~ condition`, KO versus WT, joint MACS q=0.1,
  raw ChIP counts, DESeq2 `ratio/parametric`, resolved aggregate filtering and
  configured shrinkage. Record effective fit and fallbacks.
- Regional enrichment: unique controls pooled by depth, minimum observed
  input count 5, no pseudocount. Inputs never become DE treatment replicates.

The pilot used `poscounts/local`; it does not validate these defaults. There
are no normalization or q sweeps. A concrete problem with primary settings
can motivate a specific follow-up comparison.

## Three paired comparison arms

1. **Joint:** implemented permissive joint q=0.1 discovery; support is annotation.
2. **Group-union:** union of actual current group peak catalogues with recorded
   default IDR/optimization settings. Preserve pooled-fallback provenance.
3. **Independent BED:** fixed, disjoint, reference-bounded 1 kb atlas intervals
   from latent sites/decoys, defined before condition/read generation. This
   establishes simulation independence, not the provenance of real user BEDs.

All arms use the same reads, raw-count DE method and documented filters. Report
widths and testing-family sizes; different coordinates affect sensitivity and
effect dilution. Support subsets retain their full-family adjusted p-values.

## Truth, metrics and interpretation

Calculate expected total-ChIP/input assignment for each interval using the
implementation's aligned-mate overlap and unique-assignment rules. Keep
specific-signal truth separate from total-ChIP truth. Merged intervals may
contain changed/unchanged sites; classify their expected net effect, not any
overlap with a changed site. Join results by stable region ID, never by row
position after filtering or sorting.

Use known generator exposure and predetermined invariant sites as truth
references, not fitted size factors or significant hits. In B5 a total-ChIP
change caused by background can be real while specific binding is unchanged.
In B6 report normalization-reference dependence and background dilution;
ordinary ChIP results do not establish an absolute-occupancy scale.

For each dataset/arm report:

- Candidate/tested counts, widths, support categories and missed planted sites.
- True/false calls and correct-direction site recovery, before and after
  candidate discovery/aggregate filtering.
- Effects versus expected total-ChIP effects, size factors, input reliability,
  warnings and failed/unavailable results.
- Statistical calls at `padj < 0.05` separately from displayed hits with the
  configured fold-change cutoff and shrinkage.
- Runtime/storage and representative browser/profile plots.

Report dataset `FDP = false_calls / max(all_calls, 1)` and the number of null
datasets with any rejection. Four seeds per scenario give descriptive evidence;
do not claim nominal 5% FDR is calibrated. Under a complete null any rejection
gives FDP one; zero observed errors would not establish safety. Unavailable
p/padj are non-rejections and stay visible. Retain failed/empty datasets in
the report instead of discarding them.

Review concrete behavior: correct associations, raw-count invariance under
support/enrichment-only changes, retention of condition-only candidates,
reasonable correct-direction recovery in B3/B4, low-input flags and honest
interpretation of B5/B6. Repeated false calls, wrong directions or failures
require inspection of the affected seeds. No formal confidence-bound pass/fail
gate or automatic change to q, normalization or support filtering applies.

## Execution and reporting

1. Build a small harness under `scripts/validation/` using existing fragment,
   MACS/control and count helpers. Match resolved pipeline filtering, fitting
   and results. Verify parity with a full workflow on the same realization:
   exact counts; unshrunk LFC/p/padj within `1e-10` numerical tolerance.
2. Use B1 seed 0 as the first smoke run, counted within the 24. Check generator
   truth, parity and cost, then run the remaining 23.
3. Produce one compact report and inspect predetermined windows from B3–B6.
   Profiles remain descriptive; regions are not biological replicates. Use
   existing functional tests for the full track-mode matrix.
4. Discuss concrete issues and follow-ups. Then validate one small real
   ChIP/input project on Snellius: end-to-end outputs, selective restarts,
   pooling, tracks and retention. Manuscript reruns follow later.

Use existing managed environments and declared dependencies. Record commands,
versions, configs, seeds, truth, checksums and source revision/snapshot; the
current uncommitted branch name alone is insufficient provenance.

Submit substantial compute with `sbatch` on Snellius, using
`/scratch-shared/kprange/omnomnomics/` for shared data and `$TMPDIR` for temporary
files. Synchronize the remote branch/revision and inspect its documented
environment first. Profile allocation/concurrency; respect rome's 16-core
minimum and allocation-coupled memory. IDR/MACS subjobs can reserve two threads
within each 16-core allocation; statistical settings remain unchanged. No installs in jobs. The benchmark
harness is separate from the CLI's own controller; do not wrap ordinary CLI
submission in another `sbatch` job.

Archive lightweight metrics, figures, manifests and logs under
`docs/development/chip-input-statistical-benchmark-results/`. Keep large BAMs
and count products out of the repository/manuscript figures. Archive essentials
before shared scratch's 14-day expiry. See the execution record for jobs and
source snapshots.

## Methodological references

[DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html):
raw counts, normalization references, independent filtering and unavailable
adjusted p-values.

[Lun and Smyth's read-to-region workflow](https://pmc.ncbi.nlm.nih.gov/articles/PMC4706055/):
region construction/filtering must be considered with differential testing.

[Bonhoure et al., 2014](https://genome.cshlp.org/content/24/7/1157):
experimental spike-in references support sample-to-sample ChIP normalization;
this simulation does not add such a reference to ordinary ChIP-seq.
