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

The initial clean-install attempt exposed an ordering defect in lock
serialization: sorting packages alphabetically ran Bioconductor data hooks
before their dependencies were linked. Corrected locks preserve the solver's
installation order. An exact Conda inventory alone is insufficient evidence of
a functional environment; creation also checks R namespaces, helper execution
and Python dependency consistency. Failed initial reports remain preserved. Corrected clean creation and
installed-artifact acceptance are recorded separately in
`docs/development/release-0.6.0/README.md`.

Bioconda's annotation packages install separate, versioned Bioconductor data
payloads through a package hook. `bioconductor-data-payloads.json` records their
canonical sources, exact versions and MD5 values from the SHA-256-pinned
Bioconda metadata package. The upstream hook verifies the payload checksum
before invoking R installation. Runtime namespace checks catch a missing
payload even when Conda has already recorded its wrapper package as installed.

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
installed Conda records, required R namespaces and Python dependencies. It does not update an existing environment, install
Omnomnomics, or connect helper wrappers automatically.

Regeneration uses `scripts/release/prepare_release_specs.py` for constrained
specifications, Micromamba dry-run plans for dependency solving, and
`scripts/release/build_environment_locks.py` for reviewed lock serialization.
The generator rejects failed solves, unapproved source URLs, missing archive
checksums, and incomplete dependency sets. Before accepting regenerated locks,
review dependency changes and repeat clean creation, artifact installation,
regressions, and relevant workflow validation.

Reference: [Mamba specification-file documentation](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html#specification-files).
