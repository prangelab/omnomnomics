# GRCh38 ChIP peak-QC contig audit

The human CTCF and H3K27me3 validation projects used an Ensembl-style GRCh38.p14 reference, whose primary contigs are named `1`, `2`, and so on. Peak-QC at commit `2c53a766310612862ea0e6258edf77d4618cd56c` accepted only UCSC-style names such as `chr1`, so it discarded every reference-native primary peak. The original four compact TSVs in this directory document the failure; they have been retained unchanged.

| Project and group | Source optimized peaks | Before repair | After repair | Blacklist overlaps removed |
|---|---:|---:|---:|---:|
| CTCF GM12878 | 34,002 | 0 | 33,799 | 187 |
| CTCF K562 | 34,926 | 0 | 34,674 | 226 |
| H3K27me3 GM12878 | 4,014 | 0 | 3,993 | 17 |
| H3K27me3 HepG2 | 28,929 | 0 | 28,755 | 148 |

The difference between source and retained counts also includes nonprimary contigs. Reference-native primary filtering and blacklist removal are separately counted in the revalidated `peak_filtering_summary.tsv` files. Merged sets retained 44,005 narrow peaks and 31,228 domains. Their annotations contain the same number of rows plus a header, and the pre-DE profiles and heatmaps are populated.

The repair on `input_dev` at `d5179a16fe6491053ce04770f27043973c435f6b` makes primary-contig selection reference-aware and maps external BED contig aliases to the installed reference, preserving coordinates and recording conversion provenance. It does not indiscriminately remove `chr`: UCSC-style references retain their native names. Genome installation also checks that primary GTF contigs match the reference. The full source implementation passed 104/104 regression tests in Snellius job `27110546`; the final commit's eight chromosome unit tests and end-to-end custom-BED alias test passed in job `27111143`. A real GRCh38.p14 smoke test (`27110615`) converted all 636 blacklist intervals to reference-native names and confirmed that 187 of 34,002 GM12878 CTCF peaks overlapped the corrected blacklist.

Selective Snellius reruns completed on 2026-09-24/25 for narrow CTCF (controller `27153865`), domain H3K27me3 (`27155208`), and diffuse H3K9me3 (`27155210`). Each controller and collector completed successfully. The diffuse run retained 285,511 of 310,096 initial 10-kb bins; its DE analysis tested 278,817 features and reported 71,079 significant features, compared with 282,847 and 72,757 before harmonization. Narrow and domain joint-discovery DE feature counts remained 67,161 and 171,257 because those regions are defined independently of the optimized peak-QC sets. These result changes should be interpreted as a corrected feature definition, not as evidence that a particular peak caller or biological comparison is optimal.

The reruns used a shared-node site override for count workers after the `scratch-node` queue estimate exceeded the controller time limit. The override changed scheduling only; count provenance still records raw ChIP counts without input subtraction or division. All three selective resource summaries and compact result evidence are archived under `MANUSCRIPT/validation_evidence/chip_chromosome_validation_20260925/`, with PDF figures under `MANUSCRIPT/validation_figures/chip_chromosome_validation_20260925/`. Their wall times must not be reported as full-pipeline performance.
