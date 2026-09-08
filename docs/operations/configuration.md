# Configuration

Three configuration layers are resolved for each run:

1. packaged workflow defaults define assay behavior
2. site configuration defines cluster resources and scheduler behavior
3. the generated run configuration records resolved CLI and runtime settings

Run configurations are written to `EXPERIMENT_DIR/run_configs/` and should be
retained with the analysis provenance.

## User site configuration

The packaged template is
`src/omnomnomics/workflow/config/site.yaml`. Copy it to:

```text
$XDG_CONFIG_HOME/omnomnomics/site.yaml
```

or, when `XDG_CONFIG_HOME` is unset:

```text
~/.config/omnomnomics/site.yaml
```

Key settings include:

```yaml
partition: "rome"
max_jobs: 100
max_jobs_per_second: 2
max_status_checks_per_second: 2

controller_partition: "rome"
controller_cores: 1
controller_mem_mb: 8192
controller_runtime: 1440

worker_constraint: null
heavy_io_constraint: null
genome_assembly_dir: "{HOME}/genomes/assemblies"
```

Cluster feature names are site-specific. Set `worker_constraint` only when all
worker jobs require that feature; otherwise leave it empty. Use `--site-config`
for a one-run override without modifying the persistent user file.
