build_pca_scores <- function(vst_mat, metadata_df, top_n = 1000) {
  if (is.null(vst_mat) || nrow(vst_mat) < 2 || ncol(vst_mat) < 2 || is.null(metadata_df) || !nrow(metadata_df)) {
    return(NULL)
  }
  sample_col <- if ("sample_id" %in% colnames(metadata_df)) "sample_id" else colnames(metadata_df)[[1]]
  sample_ids <- as.character(metadata_df[[sample_col]])
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  sample_ids <- intersect(sample_ids, colnames(vst_mat))
  if (length(sample_ids) < 2) {
    return(NULL)
  }

  mat <- vst_mat[, sample_ids, drop = FALSE]
  vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  ord <- order(vars, decreasing = TRUE)
  top_n <- min(as.integer(top_n), nrow(mat))
  mat <- mat[ord[seq_len(max(2, top_n))], , drop = FALSE]
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)

  scores <- as.data.frame(pca$x, stringsAsFactors = FALSE)
  scores$sample_id <- rownames(scores)
  md_idx <- match(scores$sample_id, as.character(metadata_df[[sample_col]]))
  md_sub <- metadata_df[md_idx, , drop = FALSE]
  if (sample_col %in% colnames(md_sub)) {
    md_sub[[sample_col]] <- NULL
  }
  out <- cbind(scores, md_sub)
  colnames(out) <- make.unique(colnames(out))
  rownames(out) <- NULL
  out
}

plot_pca_scores <- function(scores_df, pc_x = "PC1", pc_y = "PC2", color_by = NULL, shape_by = NULL) {
  if (is.null(scores_df) || !(pc_x %in% colnames(scores_df)) || !(pc_y %in% colnames(scores_df))) {
    return(NULL)
  }
  map_color <- !is.null(color_by) && nzchar(color_by) && color_by %in% colnames(scores_df) &&
    length(unique(as.character(scores_df[[color_by]]))) > 1
  map_shape <- !is.null(shape_by) && nzchar(shape_by) && shape_by %in% colnames(scores_df) &&
    length(unique(as.character(scores_df[[shape_by]]))) > 1

  color_map <- if (map_color) color_by else NULL
  shape_map <- if (map_shape) shape_by else NULL

  plot_df <- scores_df
  if (!is.null(color_map)) {
    plot_df[[color_map]] <- as.factor(as.character(plot_df[[color_map]]))
  }
  if (!is.null(shape_map)) {
    plot_df[[shape_map]] <- as.factor(as.character(plot_df[[shape_map]]))
  }

  p <- tryCatch({
    if (!is.null(color_map) && !is.null(shape_map)) {
      ggplot2::ggplot(plot_df, ggplot2::aes_string(x = pc_x, y = pc_y, colour = color_map, shape = shape_map))
    } else if (!is.null(color_map)) {
      ggplot2::ggplot(plot_df, ggplot2::aes_string(x = pc_x, y = pc_y, colour = color_map))
    } else if (!is.null(shape_map)) {
      ggplot2::ggplot(plot_df, ggplot2::aes_string(x = pc_x, y = pc_y, shape = shape_map))
    } else {
      ggplot2::ggplot(plot_df, ggplot2::aes_string(x = pc_x, y = pc_y))
    }
  }, error = function(e) {
    ggplot2::ggplot(plot_df, ggplot2::aes_string(x = pc_x, y = pc_y))
  })

  p +
    ggplot2::geom_point(size = 2.5, alpha = 0.9) +
    ggplot2::theme_light() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::xlab(pc_x) +
    ggplot2::ylab(pc_y)
}

build_distance_matrix <- function(vst_mat, metadata_df, top_n = 1000) {
  if (is.null(vst_mat) || nrow(vst_mat) < 2 || ncol(vst_mat) < 2 || is.null(metadata_df) || !nrow(metadata_df)) {
    return(NULL)
  }
  sample_col <- if ("sample_id" %in% colnames(metadata_df)) "sample_id" else colnames(metadata_df)[[1]]
  sample_ids <- as.character(metadata_df[[sample_col]])
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  sample_ids <- intersect(sample_ids, colnames(vst_mat))
  if (length(sample_ids) < 2) {
    return(NULL)
  }

  mat <- vst_mat[, sample_ids, drop = FALSE]
  vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  ord <- order(vars, decreasing = TRUE)
  top_n <- min(as.integer(top_n), nrow(mat))
  mat <- mat[ord[seq_len(max(2, top_n))], , drop = FALSE]
  as.matrix(stats::dist(t(mat)))
}

build_topvar_heatmap_matrix <- function(vst_mat, metadata_df, top_n = 1000) {
  if (is.null(vst_mat) || nrow(vst_mat) < 2 || ncol(vst_mat) < 2 || is.null(metadata_df) || !nrow(metadata_df)) {
    return(NULL)
  }
  sample_col <- if ("sample_id" %in% colnames(metadata_df)) "sample_id" else colnames(metadata_df)[[1]]
  sample_ids <- as.character(metadata_df[[sample_col]])
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  sample_ids <- intersect(sample_ids, colnames(vst_mat))
  if (length(sample_ids) < 2) {
    return(NULL)
  }
  mat <- vst_mat[, sample_ids, drop = FALSE]
  vars <- apply(mat, 1, stats::var, na.rm = TRUE)
  vars[is.na(vars)] <- 0
  ord <- order(vars, decreasing = TRUE)
  top_n <- min(as.integer(top_n), nrow(mat))
  mat[ord[seq_len(max(2, top_n))], , drop = FALSE]
}

read_qc_table_if_exists <- function(qc_dir, file_name) {
  path <- file.path(qc_dir, file_name)
  if (!file.exists(path)) {
    return(NULL)
  }
  utils::read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
