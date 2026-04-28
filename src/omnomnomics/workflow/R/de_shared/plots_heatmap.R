load_vst_matrix <- function(project_index_obj) {
  if (is.null(project_index_obj)) {
    return(NULL)
  }
  vst_path <- file.path(project_index_obj$qc_dir, "vst_matrix.tsv")
  if (!file.exists(vst_path)) {
    return(NULL)
  }
  tbl <- utils::read.table(
    vst_path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (ncol(tbl) < 2) {
    return(NULL)
  }
  gene_col <- colnames(tbl)[[1]]
  mat <- as.matrix(tbl[, -1, drop = FALSE])
  rownames(mat) <- as.character(tbl[[gene_col]])
  mode(mat) <- "numeric"
  mat
}

strip_id_version <- function(x) {
  sub("\\..*$", "", as.character(x))
}

resolve_heatmap_genes <- function(tbl, vst_mat, top_n = 50) {
  if (is.null(tbl) || !nrow(tbl) || is.null(vst_mat) || !nrow(vst_mat)) {
    return(character(0))
  }

  candidate_cols <- c("gene_id", "Geneid", "Geneid_noVersion", "symbol", "gene_symbol", "SYMBOL", "external_gene_name")
  candidate_cols <- candidate_cols[candidate_cols %in% colnames(tbl)]
  if (!length(candidate_cols)) {
    candidate_cols <- colnames(tbl)[[1]]
  }

  vst_ids <- rownames(vst_mat)
  vst_ids_nover <- strip_id_version(vst_ids)

  for (one_col in candidate_cols) {
    vals <- unique(as.character(tbl[[one_col]]))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    vals <- utils::head(vals, top_n * 4)
    if (!length(vals)) {
      next
    }

    direct <- vals[vals %in% vst_ids]
    if (length(direct) >= 2) {
      return(utils::head(unique(direct), top_n))
    }

    vals_nover <- strip_id_version(vals)
    idx <- match(vals_nover, vst_ids_nover)
    mapped <- vst_ids[stats::na.omit(idx)]
    mapped <- unique(mapped)
    if (length(mapped) >= 2) {
      return(utils::head(mapped, top_n))
    }
  }

  character(0)
}

resolve_gene_display_labels <- function(tbl, vst_ids) {
  if (is.null(tbl) || !nrow(tbl) || !length(vst_ids)) {
    return(stats::setNames(as.character(vst_ids), as.character(vst_ids)))
  }

  id_candidates <- c("gene_id", "Geneid", "Geneid_noVersion")
  id_candidates <- id_candidates[id_candidates %in% colnames(tbl)]
  symbol_candidates <- c("gene_symbol", "symbol", "SYMBOL", "external_gene_name")
  symbol_candidates <- symbol_candidates[symbol_candidates %in% colnames(tbl)]

  labels <- stats::setNames(as.character(vst_ids), as.character(vst_ids))
  if (!length(symbol_candidates)) {
    return(labels)
  }

  symbol_col <- symbol_candidates[[1]]
  work_tbl <- tbl
  work_tbl$.symbol <- as.character(work_tbl[[symbol_col]])

  if (length(id_candidates)) {
    for (id_col in id_candidates) {
      work_tbl$.id <- as.character(work_tbl[[id_col]])
      work_tbl <- work_tbl[!is.na(work_tbl$.id) & nzchar(work_tbl$.id), , drop = FALSE]
      if (!nrow(work_tbl)) {
        next
      }

      id_map <- stats::setNames(work_tbl$.symbol, work_tbl$.id)
      nover_map <- stats::setNames(work_tbl$.symbol, strip_id_version(work_tbl$.id))
      for (gid in vst_ids) {
        sym <- id_map[[gid]]
        if (is.null(sym) || is.na(sym) || !nzchar(sym)) {
          sym <- nover_map[[strip_id_version(gid)]]
        }
        if (!is.null(sym) && !is.na(sym) && nzchar(sym)) {
          labels[[gid]] <- sym
        }
      }
      return(labels)
    }
  }

  fallback_symbols <- unique(as.character(work_tbl$.symbol))
  fallback_symbols <- fallback_symbols[!is.na(fallback_symbols) & nzchar(fallback_symbols)]
  if (length(fallback_symbols)) {
    for (i in seq_along(vst_ids)) {
      if (i <= length(fallback_symbols)) {
        labels[[vst_ids[[i]]]] <- fallback_symbols[[i]]
      }
    }
  }
  labels
}

row_zscore <- function(mat) {
  scaled <- t(scale(t(mat)))
  scaled[is.na(scaled)] <- 0
  scaled
}
