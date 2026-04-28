required_packages <- c("shiny")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages for DE explorer app: ", paste(missing_packages, collapse = ", "))
}

app_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
shared_dir <- file.path(app_dir, "..", "de_shared")
modules_dir <- file.path(app_dir, "modules")

source(file.path(shared_dir, "io.R"), local = FALSE)
source(file.path(shared_dir, "de_tables.R"), local = FALSE)
source(file.path(shared_dir, "plots_volcano.R"), local = FALSE)
source(file.path(shared_dir, "plots_heatmap.R"), local = FALSE)
source(file.path(shared_dir, "plots_boxplot.R"), local = FALSE)
source(file.path(shared_dir, "qc_plots.R"), local = FALSE)
source(file.path(shared_dir, "enrichment_plots.R"), local = FALSE)
source(file.path(modules_dir, "mod_project_loader.R"), local = FALSE)
source(file.path(modules_dir, "mod_analysis_selector.R"), local = FALSE)
source(file.path(modules_dir, "mod_de_table.R"), local = FALSE)
source(file.path(modules_dir, "mod_volcano.R"), local = FALSE)
source(file.path(modules_dir, "mod_heatmap.R"), local = FALSE)
source(file.path(modules_dir, "mod_boxplots.R"), local = FALSE)
source(file.path(modules_dir, "mod_qc.R"), local = FALSE)
source(file.path(modules_dir, "mod_enrichment.R"), local = FALSE)
