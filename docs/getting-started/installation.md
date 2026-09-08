# Installation

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

Both installers expose wrappers to the main environment. Normal runs still only
require `micromamba activate omnomnomics`.

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
