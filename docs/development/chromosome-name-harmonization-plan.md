# Reference-native chromosome names for BED inputs

## Problem and decision

`omnomnomics genomes install` currently normalizes the installation layout and builds or copies indexes. It copies FASTA, FASTA index, GTF, chromosome sizes, and blacklist without harmonizing contig names (`src/omnomnomics/genomes.py`). The GRCh38.p14 validation reference on Snellius uses `1`, `10`, etc.; its cached UCSC blacklist uses `chr10`, etc. In step 13, the standard-chromosome filter accepts only `chr`-prefixed names, removing every optimized peak in the human CTCF and H3K27me3 validations. The measured losses are preserved in `docs/development/chip-input-peak-qc-audit-20260924/`.

The installed FASTA index (`fasta/genome.fa.fai`) defines the names used by alignments and internal genomic outputs. Never rename the FASTA, BAM headers, aligner indexes, or internally generated peak BEDs to impose UCSC or Ensembl naming. Convert only external interval resources to reference-native names at a controlled boundary, and use reference-aware filtering for internal BEDs. A `chr` prefix change is not a genome-build conversion or liftover.

## Mapping contract

1. Read the exact reference contig names and lengths from `genome.fa.fai`. Prefer exact name matches.
2. Permit only explicit, unambiguous aliases for primary chromosomes (`1`/`chr1`, `X`/`chrX`, `Y`/`chrY`) and recognized mitochondrial names (`MT`/`M`/`chrM`) when the target reference has exactly one matching contig. Do not guess names for patches, decoys, or alternative haplotypes; accept them only by exact match.
3. Reject mappings that collapse two distinct input contigs onto one reference contig, or where more than one reference target is possible. Validate integer, zero-based half-open BED coordinates and target-contig bounds after mapping.
4. Preserve each external source file byte-for-byte. Record its checksum, selected assembly, reference-index checksum, mapping table, transformed and skipped record counts, and helper version with any generated reference-native BED. Log every name change. Keep BED coordinates and other fields unchanged.
5. A name alias cannot establish assembly identity. Require the selected assembly as the user's declared target and reject out-of-bounds records. Do not infer that an arbitrary `chr1` BED is GRCh38 merely because its names can be mapped.

## Implementation order

### 1. Shared reference-name helper

Add a small module under `src/omnomnomics/` for loading the `.fai`, resolving the limited alias set, identifying primary chromosomes in the installed reference, validating coordinates, and writing provenance. Use one implementation from genome installation, CLI validation, and workflow rules. Keep it independent of ChIP-specific code.

### 2. Installation and blacklist boundary

In `src/omnomnomics/genomes.py`, validate that the installed FASTA index and chromosome-sizes file agree and that annotation contigs required for analysis exist in the reference. Preserve the provider's original blacklist. Create a separately named reference-native blacklist with checksum-based provenance after either `genomes install` or `genomes blacklist`; also support existing installations when a blacklist is first resolved. Do not overwrite the provider file.

For a blacklist, retain exact-match contigs and safe primary aliases. Report records on unsupported contigs that are absent from the reference; fail if no usable intervals remain or if mapping is ambiguous. Make the workflow consume the reference-native derivative, including for an Ensembl reference with a UCSC blacklist. Avoid downloading a blacklist unexpectedly during a running analysis job.

### 3. Peak-QC and other interval operations

Replace the `^chr([0-9]+|X|Y)$` filter in `src/omnomnomics/workflow/rules/13.peak_qc.smk` with primary-contig membership based on the selected reference. Apply the same policy to `src/omnomnomics/workflow/rules/10.call_peaks.smk`, which already accepts both prefix styles for diffuse bins but uses blacklists for diffuse filtering and peak-optimization scoring. Use the derived blacklist in every `bedtools intersect` call. Report counts before filtering, after primary-contig selection, and after blacklist exclusion separately. Fail when a source containing reference-primary peaks loses all of them at primary-contig selection; report a truly empty source or an all-blacklisted result as distinct cases.

Include the derived blacklist and mapping/reference specification as explicit rule inputs or stable settings fingerprints wherever they affect results. A new blacklist or changed alias policy must invalidate cached peak-QC, bins, and optimization outputs in full and selective runs.

### 4. User-supplied ChIP region BED

Update `--chip-regions-bed` validation in `src/omnomnomics/cli.py` and region construction in `src/omnomnomics/workflow/chip_analysis.smk` and `src/omnomnomics/chip_analysis.py` to accept safe aliases and emit reference-native testing regions. Continue to protect the original BED as an immutable user source. Keep its original feature names in the region manifest and record the source checksum and all contig mappings. Reject unknown or ambiguous contigs, duplicates created by mapping, invalid bounds, and unsupported overlap patterns before counting. BED-only selective entry must work without requiring peak calling.

### 5. Verification and targeted reanalysis

- Unit tests: Ensembl-style and UCSC-style references; exact-match precedence; both alias directions; mitochondria; alternative contigs; mixed/ambiguous names; out-of-bounds coordinates; source preservation; deterministic derivative/provenance; blacklist with zero usable intervals.
- Workflow tests: nonempty peak sets survive step 13 under both reference styles; blacklist intervals actually intersect reference-native peaks or bins; changed blacklist or helper version triggers selective recomputation; custom BED entry produces stable reference-native region IDs without modifying the source.
- Run the regression suite locally, then synchronize the fixed `input_dev` revision to Snellius. Submit non-trivial validation through `sbatch`. Rerun the affected human CTCF and H3K27me3 peak-QC and dependent outputs while retaining their upstream peak calls. Recheck their peak counts, annotations, plots, provenance, and downstream peak-dependent results. Reassess diffuse-bin outputs because correcting blacklist names can change bin eligibility. Confirm UCSC-style and mouse validations still behave as before.
- Only after the reruns, correct the `MANUSCRIPT/validation_evidence/chip_input_validation_20260922/README.md` claim that the empty peak sets were legitimate, refresh lightweight evidence and figure copies under the repository artifact rules, and update the manuscript review ledger. Do not treat selective rerun times as full pipeline performance measurements.

No parameter thresholds, peak-calling methods, or DE normalization are changed by this naming fix.
