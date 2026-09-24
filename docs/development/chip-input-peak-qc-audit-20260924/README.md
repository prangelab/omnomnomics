# GRCh38 ChIP peak-QC contig audit

The human narrow CTCF and domain H3K27me3 validation runs used the Ensembl-style GRCh38.p14 reference. The source peak BEDs use chromosome names such as `1`. Step 13 currently retains only names matching `^chr([0-9]+|X|Y)$` (`src/omnomnomics/workflow/rules/13.peak_qc.smk`, line 93 at commit `2c53a766310612862ea0e6258edf77d4618cd56c`). This removes standard chromosomes in the run's reference naming scheme.

| Project and group | Optimized peaks before peak-QC | After peak-QC |
|---|---:|---:|
| CTCF GM12878 | 34,002 | 0 |
| CTCF K562 | 34,926 | 0 |
| H3K27me3 GM12878 | 4,014 | 0 |
| H3K27me3 HepG2 | 28,929 | 0 |

The CTCF IDR summary reports 34,002 and 34,926 consensus peaks with status `ok`. The domain concordance summary reports 4,014 and 28,929 final union peaks. The corresponding source BEDs were inspected on Snellius and begin with chromosome `1`; the downstream peak-QC summaries report zero retained peaks. The attached four TSVs were copied from `/scratch-shared/kprange/omnomtest/chip_input_validation_d6ca74c90779/projects/` on 2026-09-24. No full peak BEDs or raw data were copied.

This establishes a peak-QC contig-name bug. It does not assess whether the upstream peaks are biologically optimal or whether the separate joint-discovery DE calls are well calibrated. The `MANUSCRIPT/validation_evidence/chip_input_validation_20260922/README.md` description of the empty sets as "legitimately empty" should be corrected after repair and revalidation.

Next, make the standard-chromosome filter accept both `1` and `chr1` naming schemes, verify blacklist coordinate compatibility, add a regression test, and selectively rerun the affected peak-QC and downstream outputs. Preserve the upstream peak calls for comparison.
