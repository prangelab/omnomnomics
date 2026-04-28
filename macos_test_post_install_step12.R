#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install(
  c(
    "sva",
    "decoupleR",
    "OmnipathR"
  ),
  ask = FALSE,
  update = FALSE
)

message("Post-install complete for macOS test DE packages.")
