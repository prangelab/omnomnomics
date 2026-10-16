# Release environment locks

The `linux-64/` files are **release candidates under validation**, not an
approved 0.6.0 release. They target Linux x86-64. They do not establish support
for macOS, Windows, or other processor architectures.

Four environments are kept separate: `omnomnomics` (main pipeline),
`omnomnomics-idr`, `omnomnomics-spp`, and `omnomnomics-explorer`. The ordinary
`environment*.yml` files remain readable dependency specifications. The
`*-lock.yml` files include exact direct and transitive package versions,
archive URLs, MD5 and SHA-256 checksums in Mamba's unified lock format. The
`.explicit.txt` files provide the equivalent platform-specific explicit
specifications; `.packages.json` files support post-install verification.

The candidate locks were solved against the project manifests while pinning
the dependency closure from the established Snellius environments wherever
available. Unrelated installed packages were excluded. The main environment
also includes the build-system requirements already declared in
`pyproject.toml`. SVA was declared in `environment.yml` but absent from the
established environment; the candidate solve adds it and its dependencies.

Every archive's exact name, version, build, publishing community, source URL
and checksums appears in `linux-64/package_provenance.tsv`. Only the public
Conda-forge and Bioconda channels occur in these candidate locks. No package
from an unrecorded pip install is copied into the locks. Omnomnomics itself
must be built and installed separately from the pinned source artifact.

To create one new environment after reviewing its provenance:

```bash
bash scripts/create_locked_environment.sh main /absolute/path/to/new/main
```

Replace `main` with `idr`, `spp`, or `explorer` for another component. The
script rejects existing prefixes, checks the lock-file checksums, creates the
locked environment without site configuration overrides, and verifies the
installed Conda records. It does not update an existing environment, install
Omnomnomics, or connect helper wrappers automatically.

Regeneration uses `scripts/release/prepare_release_specs.py` for constrained
specifications, Micromamba dry-run plans for dependency solving, and
`scripts/release/build_environment_locks.py` for reviewed lock serialization.
The generator rejects failed solves, unapproved source URLs, missing archive
checksums, and incomplete dependency sets. Before accepting regenerated locks,
review dependency changes and repeat clean creation, artifact installation,
regressions, and relevant workflow validation.

Reference: [Mamba specification-file documentation](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html#specification-files).
