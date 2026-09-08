# Track hubs

Step 8 creates per-sample BigWigs and step 9 groups them into UCSC track hubs.
Use metadata-driven `sample_type` to determine the hub grouping.

For samples such as `NT_C_D1_A`, `NT_L_D1_A`, `KD24_C_D1_A`, and
`KD24_L_D1_A`, grouping by `genotype,stim` creates four hubs: `NT_C`, `NT_L`,
`KD24_C`, and `KD24_L`.

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  -j 8-9 \
  -m metadata.tsv \
  --sample-name genotype,stim,donor,replicate \
  --sample-type genotype,stim
```

RNA BigWigs are stranded. Plus-strand tracks use `255,0,0` and minus-strand
tracks use `0,0,255`; palette selection does not override these strand colors.

For ATAC and ChIP, `--sample-color` categories are mapped through the selected
color table. Use the helper commands to create or inspect palettes:

```bash
omnomnomics create-track-color-table --help
omnomnomics display-track-color-table --help
```

The packaged palette directory is read-only. Supply a writable custom palette
directory with `-P` when creating tables.
