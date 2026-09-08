# Metadata and sample identity

Metadata provides stable sample identity independently of FASTQ suffixes and
drives grouping, color categories, and differential designs.

The first column must be named `filename` and map to input filenames after
normalization. Additional columns describe biological and technical structure.

```text
filename	genotype	stim	donor	replicate
NT_C_D1_A_S50_L003	NT	C	D1	A
NT_L_D1_A_S58_L003	NT	L	D1	A
KD24_C_D1_A_S54_L003	KD24	C	D1	A
KD24_L_D1_A_S62_L003	KD24	L	D1	A
```

Selectors accept comma-separated column names or one-based indices:

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  -m metadata.tsv \
  --sample-name genotype,stim,donor,replicate \
  --sample-type genotype,stim \
  --sample-color genotype
```

- `sample_id` must uniquely identify each biological sample.
- `sample_type` defines analysis and track-hub groups.
- `sample_color` defines palette categories for non-stranded tracks.

Technical lanes are processed independently and merged at BAM level. Do not
remove lane identity from `filename`; use `sample_id` derivation to identify the
biological sample to which lanes belong.

For stranded RNA hubs, plus tracks are red and minus tracks are blue regardless
of `sample_color`.
