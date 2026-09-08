# RNA-seq

RNA runs default to HISAT2 mapping, keep duplicate reads, use MAPQ 15 filtering,
and create separate plus- and minus-strand BigWigs.

```bash
omnomnomics rna \
  -i EXPERIMENT \
  -g GRCh38 \
  -m metadata.tsv \
  --sample-name genotype,stim,donor,replicate \
  --sample-type genotype,stim \
  --de-columns genotype,stim \
  --de-block donor
```

The default full public flow runs through step 12:

1. preprocessing and alignment (steps 1-4)
2. filtered BAMs and alignment QC (steps 5-7)
3. stranded BigWigs and hubs (steps 8-10)
4. featureCounts table (step 11)
5. DESeq2 analysis (step 12)

STAR and STAR-TE are available with `-M`. Use `--remove-duplicates` only when the
experimental design specifically requires duplicate removal; retaining RNA
duplicates is the assay-aware default.

For storage-constrained projects, `--retention-policy pruned` preserves reusable
filtered BAMs while removing bulky upstream intermediates after they are safe to
delete.
