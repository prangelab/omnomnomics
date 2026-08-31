normalize_existing_dir <- function(path_value, label = "directory") {
  if (is.null(path_value) || !nzchar(trimws(path_value))) {
    stop(label, " is empty.")
  }
  resolved <- normalizePath(path.expand(path_value), winslash = "/", mustWork = FALSE)
  if (!dir.exists(resolved)) {
    stop(label, " does not exist: ", resolved)
  }
  resolved
}

de_table_suffixes <- c("genes", "peaks", "gene_bodies", "bins")

de_table_pattern <- function() {
  paste0("\\.diff_(", paste(de_table_suffixes, collapse = "|"), ")\\.DESeq2\\.txt$")
}

find_contrast_de_table <- function(contrast_dir, contrast_label) {
  candidates <- file.path(
    contrast_dir,
    paste0(contrast_label, ".diff_", de_table_suffixes, ".DESeq2.txt")
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[[1]]
}

resolve_de_root <- function(path_value) {
  root <- normalize_existing_dir(path_value, "Input directory")
  if (basename(root) == "DE_calling") {
    return(root)
  }
  nested <- file.path(root, "DE_calling")
  if (dir.exists(nested)) {
    return(normalizePath(nested, winslash = "/", mustWork = FALSE))
  }
  stop(
    "Could not resolve DE_calling directory from input: ", root,
    ". Provide either the project directory containing DE_calling or the DE_calling directory itself."
  )
}

list_analysis_dirs <- function(de_root) {
  de_root <- resolve_de_root(de_root)
  subdirs <- list.dirs(de_root, recursive = FALSE, full.names = TRUE)
  if (length(subdirs) == 0) {
    return(character(0))
  }

  keep <- vapply(subdirs, function(one_dir) {
    base_name <- basename(one_dir)
    if (startsWith(base_name, ".")) {
      return(FALSE)
    }
    if (identical(base_name, "qc")) {
      return(FALSE)
    }
    has_manifest <- file.exists(file.path(one_dir, "contrast_plan.tsv")) ||
      file.exists(file.path(one_dir, "contrast_summary.tsv"))
    if (has_manifest) {
      return(TRUE)
    }
    nested <- list.dirs(one_dir, recursive = FALSE, full.names = TRUE)
    any(vapply(nested, function(candidate) {
      length(list.files(candidate, pattern = de_table_pattern(), full.names = TRUE)) > 0
    }, logical(1)))
  }, logical(1))

  sort(subdirs[keep])
}

index_contrasts <- function(analysis_dir) {
  analysis_dir <- normalize_existing_dir(analysis_dir, "analysis directory")
  manifest_path <- file.path(analysis_dir, "contrast_plan.tsv")
  contrast_labels <- character(0)
  manifest_tbl <- NULL

  if (file.exists(manifest_path)) {
    manifest_tbl <- utils::read.table(
      manifest_path,
      sep = "\t",
      header = TRUE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if ("contrast_label" %in% colnames(manifest_tbl)) {
      contrast_labels <- unique(as.character(manifest_tbl$contrast_label))
      contrast_labels <- contrast_labels[!is.na(contrast_labels) & nzchar(contrast_labels)]
    }
  }

  if (length(contrast_labels) == 0) {
    candidate_dirs <- list.dirs(analysis_dir, recursive = FALSE, full.names = TRUE)
    contrast_labels <- basename(candidate_dirs[vapply(candidate_dirs, function(one_dir) {
      length(list.files(one_dir, pattern = de_table_pattern(), full.names = TRUE)) > 0
    }, logical(1))])
  }

  if (length(contrast_labels) == 0) {
    return(data.frame(
      contrast_id = character(0),
      contrast_dir = character(0),
      de_table = character(0),
      sig_table = character(0),
      clusterprofiler_dir = character(0),
      decoupler_dir = character(0),
      custom_modules_dir = character(0),
      combined_summary_plot = character(0),
      stringsAsFactors = FALSE
    ))
  }

  contrast_labels <- sort(unique(contrast_labels))
  rows <- lapply(contrast_labels, function(label) {
    contrast_dir <- file.path(analysis_dir, label)
    de_table <- find_contrast_de_table(contrast_dir, label)
    sig_table <- file.path(contrast_dir, paste0(label, ".sig_only.tsv"))
    combined_plot <- file.path(contrast_dir, paste0(label, ".combined_summary.pdf"))
    row <- data.frame(
      contrast_id = label,
      contrast_dir = contrast_dir,
      de_table = de_table,
      sig_table = if (file.exists(sig_table)) sig_table else NA_character_,
      clusterprofiler_dir = if (dir.exists(file.path(contrast_dir, "clusterProfiler"))) file.path(contrast_dir, "clusterProfiler") else NA_character_,
      decoupler_dir = if (dir.exists(file.path(contrast_dir, "decoupler"))) file.path(contrast_dir, "decoupler") else NA_character_,
      custom_modules_dir = if (dir.exists(file.path(contrast_dir, "custom_modules"))) file.path(contrast_dir, "custom_modules") else NA_character_,
      combined_summary_plot = if (file.exists(combined_plot)) combined_plot else NA_character_,
      stringsAsFactors = FALSE
    )
    if (!is.null(manifest_tbl) && "contrast_label" %in% colnames(manifest_tbl)) {
      manifest_idx <- match(label, as.character(manifest_tbl$contrast_label))
      if (!is.na(manifest_idx)) {
        for (column in setdiff(colnames(manifest_tbl), "contrast_label")) {
          row[[column]] <- manifest_tbl[[column]][[manifest_idx]]
        }
      }
    }
    row
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

read_metadata_derived <- function(de_root) {
  de_root <- resolve_de_root(de_root)
  metadata_path <- file.path(de_root, "metadata_derived.tsv")
  if (!file.exists(metadata_path)) {
    stop("Missing metadata_derived.tsv in ", de_root)
  }
  utils::read.table(
    metadata_path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

read_qc_assets <- function(de_root) {
  de_root <- resolve_de_root(de_root)
  qc_dir <- file.path(de_root, "qc")
  if (!dir.exists(qc_dir)) {
    stop("Missing shared qc directory: ", qc_dir)
  }
  all_files <- list.files(qc_dir, recursive = TRUE, full.names = TRUE)
  qc_prefix <- paste0(normalizePath(qc_dir, winslash = "/", mustWork = FALSE), "/")
  rel_paths <- vapply(all_files, function(one_file) {
    file_norm <- normalizePath(one_file, winslash = "/", mustWork = FALSE)
    if (startsWith(file_norm, qc_prefix)) {
      return(substr(file_norm, nchar(qc_prefix) + 1, nchar(file_norm)))
    }
    basename(one_file)
  }, character(1))
  data.frame(
    file = all_files,
    rel_path = rel_paths,
    stringsAsFactors = FALSE
  )
}

discover_de_calling <- function(de_root) {
  de_root <- resolve_de_root(de_root)
  metadata_path <- file.path(de_root, "metadata_derived.tsv")
  qc_dir <- file.path(de_root, "qc")
  if (!file.exists(metadata_path)) {
    stop("Required file is missing: ", metadata_path)
  }
  if (!dir.exists(qc_dir)) {
    stop("Required directory is missing: ", qc_dir)
  }

  analysis_dirs <- list_analysis_dirs(de_root)
  analyses <- lapply(analysis_dirs, function(one_dir) {
    list(
      analysis_id = basename(one_dir),
      analysis_dir = one_dir,
      contrast_index = index_contrasts(one_dir)
    )
  })
  names(analyses) <- vapply(analyses, function(x) x$analysis_id, character(1))

  list(
    de_root = de_root,
    metadata_path = metadata_path,
    metadata = read_metadata_derived(de_root),
    qc_dir = qc_dir,
    qc_index = read_qc_assets(de_root),
    analyses = analyses
  )
}
