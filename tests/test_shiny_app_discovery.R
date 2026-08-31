script_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
repo_root <- dirname(dirname(script_path))
source(file.path(repo_root, "src", "omnomnomics", "workflow", "R", "de_shared", "io.R"))

suffixes <- c("genes", "peaks", "gene_bodies", "bins")

for (suffix in suffixes) {
  project <- tempfile(paste0("shiny-discovery-", suffix, "-"))
  de_root <- file.path(project, "DE_calling")
  analysis_dir <- file.path(de_root, "results")
  contrast_label <- "condition_A_vs_condition_B"
  contrast_dir <- file.path(analysis_dir, contrast_label)
  dir.create(file.path(de_root, "qc"), recursive = TRUE)
  dir.create(contrast_dir, recursive = TRUE)

  utils::write.table(
    data.frame(sample_id = c("A_rep1", "B_rep1"), condition = c("A", "B")),
    file.path(de_root, "metadata_derived.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  utils::write.table(
    data.frame(
      contrast_label = contrast_label,
      contrast_type = "factor",
      numerator = "condition_A",
      denominator = "condition_B"
    ),
    file.path(analysis_dir, "contrast_plan.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  table_path <- file.path(
    contrast_dir,
    paste0(contrast_label, ".diff_", suffix, ".DESeq2.txt")
  )
  utils::write.table(
    data.frame(gene_id = "feature_1", log2FoldChange = 2, padj = 0.01),
    table_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  discovered <- discover_de_calling(project)
  stopifnot(length(discovered$analyses) == 1)
  contrast <- discovered$analyses[[1]]$contrast_index
  stopifnot(nrow(contrast) == 1)
  stopifnot(identical(contrast$de_table[[1]], normalizePath(table_path, mustWork = TRUE)))
  stopifnot(identical(contrast$contrast_type[[1]], "factor"))
  stopifnot(identical(contrast$numerator[[1]], "condition_A"))
  stopifnot(identical(contrast$denominator[[1]], "condition_B"))

  unlink(project, recursive = TRUE)
}

cat("Shiny project discovery tests passed for:", paste(suffixes, collapse = ", "), "\n")
