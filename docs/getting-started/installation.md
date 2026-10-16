# Installation

## Locked Linux environments

The 0.6.0 candidate includes complete Linux x86-64 environment locks. These
locks are undergoing clean-install and release validation; see
[lock status and provenance](https://github.com/prangelab/omnomnomics/blob/input_dev/locks/README.md).
Use a reviewed release checkout or immutable commit for installation.

After cloning the repository, choose a new environment root on storage
appropriate for software installations, then create the isolated components:

```bash
OMNOM_ENV_ROOT=/absolute/path/to/omnomnomics-environments
bash scripts/create_locked_environment.sh main "$OMNOM_ENV_ROOT/main"
bash scripts/create_locked_environment.sh idr "$OMNOM_ENV_ROOT/idr"
bash scripts/create_locked_environment.sh spp "$OMNOM_ENV_ROOT/spp"
micromamba activate "$OMNOM_ENV_ROOT/main"
python -m pip install --no-deps --no-build-isolation --no-index .
bash scripts/link_locked_helpers.sh "$OMNOM_ENV_ROOT/main" \
  "$OMNOM_ENV_ROOT/idr" "$OMNOM_ENV_ROOT/spp"
```

The creation script verifies the lock checksums, refuses to modify an existing
prefix, and checks the complete installed Conda inventory. Helper linking also
verifies all three inventories. Normal pipeline runs activate only the main
environment. The full pipeline requires Slurm; the environment locks alone do
not configure the cluster or install reference data.

The locks are platform-specific and do not establish macOS or Windows support.
The ordinary YAML specifications below are useful for development and other
platform investigations, but re-solving them does not reproduce the locked
release environment.

## Main environment

Clone the repository, create the managed environment, and install the package:

```bash
git clone https://github.com/prangelab/omnomnomics.git
cd omnomnomics
micromamba env create -f environment.yml
micromamba activate omnomnomics
python -m pip install -e .
```

The supplied environment is the supported source of Snakemake and workflow tool
versions. Omnomnomics does not require a `conda` executable for rule execution;
the environment is created and activated with Micromamba before the controller
is submitted.

## Optional companion environments

Default narrow ATAC/ChIP analysis uses IDR. Install its isolated helper
environment once:

```bash
bash scripts/install_idr_helper.sh
```

SPP cross-correlation QC is optional and uses a second helper environment:

```bash
bash scripts/install_spp_helper.sh
```

Both installers use the Linux release locks and expose wrappers to the main
environment. An existing helper environment is reused only when its complete
inventory matches the lock; a mismatch stops installation without updating it.
Normal runs still only require `micromamba activate omnomnomics`. For custom
environment prefixes, use the explicit locked-creation and linking commands
above instead of the named-environment helpers.

## Differential Explorer only

For local inspection of completed differential-analysis outputs without the HPC
toolchain:

```bash
micromamba env create -f environment.explorer.yml
micromamba activate omnomnomics-explorer
python -m pip install --no-deps --no-build-isolation .
omnomnomics-de-app --project-dir /path/to/project_or_DE_calling
```

Use `environment.explorer.yml` from the same release that produced the results.
The candidate Linux Explorer lock can instead be created with
`bash scripts/create_locked_environment.sh explorer /absolute/path/to/new/explorer`.
Activate that prefix and install the application code as shown above. Explorer
installation and application startup are checked separately from the pipeline.

## Site configuration

Copy the packaged template to the user configuration directory and edit the
cluster-specific values:

```bash
mkdir -p "${XDG_CONFIG_HOME:-$HOME/.config}/omnomnomics"
cp src/omnomnomics/workflow/config/site.yaml \
  "${XDG_CONFIG_HOME:-$HOME/.config}/omnomnomics/site.yaml"
```

The search order is:

1. `--site-config PATH`
2. `$XDG_CONFIG_HOME/omnomnomics/site.yaml`
3. `~/.config/omnomnomics/site.yaml`
4. the packaged site defaults

The packaged cluster defaults target the
[SURF Snellius national supercomputer](https://servicedesk.surf.nl/wiki/spaces/WIKI/pages/30660184/Snellius).
Before running on another HPC system, review
[Configuration](../operations/configuration.md) and provide a site configuration
for that environment.
