# Genomes and annotations

Reference assemblies are managed outside the code checkout. Omnomnomics uses a
normalized assembly layout:

```text
REFERENCE_ROOT/
├── assemblies/
│   └── GRCh38/
│       ├── fasta/genome.fa
│       ├── annotation/genes.gtf
│       ├── hisat2/
│       ├── star/
│       └── aux/
└── motif_databases/
```

Common operations:

```bash
omnomnomics genomes list --species human
omnomnomics genomes installed
omnomnomics genomes install --species mouse
omnomnomics genomes blacklist --assembly GRCh38
omnomnomics genomes motifs
```

`genomes install` normalizes FASTA and GTF locations and builds requested
HISAT2/STAR indexes directly. STAR indexing can require substantially more
memory than HISAT2 indexing.

When available, ENCODE blacklist files are cached in the assembly `aux/`
directory. The default JASPAR MEME database is cached once under
`motif_databases/`; use `genomes motifs --force` to refresh it. For strict
reproducibility, set an explicit MEME-format motif database in the workflow
configuration.
