read_enrichment_tables <- function(contrast_dir) {
  if (is.null(contrast_dir) || !nzchar(contrast_dir) || !dir.exists(contrast_dir)) {
    return(list(clusterprofiler = NULL, decoupler = NULL))
  }

  rowbind_fill <- function(df_list) {
    df_list <- df_list[!vapply(df_list, is.null, logical(1))]
    if (length(df_list) == 0) {
      return(NULL)
    }
    all_cols <- unique(unlist(lapply(df_list, colnames), use.names = FALSE))
    aligned <- lapply(df_list, function(one_df) {
      missing_cols <- setdiff(all_cols, colnames(one_df))
      if (length(missing_cols) > 0) {
        for (mc in missing_cols) {
          one_df[[mc]] <- NA
        }
      }
      one_df <- one_df[, all_cols, drop = FALSE]
      rownames(one_df) <- NULL
      one_df
    })
    out <- do.call(rbind, aligned)
    rownames(out) <- NULL
    out
  }

  cp_dir <- file.path(contrast_dir, "clusterProfiler")
  cp_tbl <- NULL
  if (dir.exists(cp_dir)) {
    cp_files <- list.files(
      cp_dir,
      pattern = "\\.(ORA|GSEA)\\.tsv$",
      full.names = TRUE,
      recursive = TRUE
    )
    cp_rows <- lapply(cp_files, function(one_file) {
      tab <- tryCatch(utils::read.delim(one_file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      if (is.null(tab) || !nrow(tab)) {
        return(NULL)
      }
      tab$source_file <- substring(one_file, nchar(cp_dir) + 2L)
      tab
    })
    cp_tbl <- rowbind_fill(cp_rows)
  }

  dc_dir <- file.path(contrast_dir, "decoupler")
  dc_tbl <- NULL
  if (dir.exists(dc_dir)) {
    dc_files <- list.files(dc_dir, pattern = "\\.tsv$", full.names = TRUE)
    dc_rows <- lapply(dc_files, function(one_file) {
      tab <- tryCatch(utils::read.delim(one_file, stringsAsFactors = FALSE, check.names = FALSE), error = function(e) NULL)
      if (is.null(tab) || !nrow(tab)) {
        return(NULL)
      }
      if (!all(c("source", "score") %in% colnames(tab))) {
        return(NULL)
      }
      tab$source_file <- basename(one_file)
      tab
    })
    dc_tbl <- rowbind_fill(dc_rows)
  }

  list(clusterprofiler = cp_tbl, decoupler = dc_tbl)
}

plot_cp_top_terms <- function(cp_tbl, top_n = 15, colors = c(positive = "#14796F", negative = "#D77A63")) {
  if (is.null(cp_tbl) || !nrow(cp_tbl)) {
    return(NULL)
  }
  desc_col <- if ("Description" %in% colnames(cp_tbl)) "Description" else if ("ID" %in% colnames(cp_tbl)) "ID" else NULL
  p_col <- if ("p.adjust" %in% colnames(cp_tbl)) "p.adjust" else if ("qvalue" %in% colnames(cp_tbl)) "qvalue" else NULL
  if (is.null(desc_col) || is.null(p_col)) {
    return(NULL)
  }
  df <- cp_tbl
  df$term <- as.character(df[[desc_col]])
  df$padj <- suppressWarnings(as.numeric(df[[p_col]]))
  method <- if ("method" %in% colnames(df)) unique(as.character(df$method)) else "ORA"
  if (length(method) != 1) stop("Plot ORA and GSEA results separately.")
  is_gsea <- identical(method, "GSEA")
  if (is_gsea && !("NES" %in% colnames(df))) return(NULL)
  df$score <- if (is_gsea) suppressWarnings(as.numeric(df$NES)) else -log10(pmax(df$padj, 1e-300))
  df <- df[is.finite(df$score) & is.finite(df$padj) & df$padj >= 0 & df$padj <= 1 & !is.na(df$term) & nzchar(df$term), , drop = FALSE]
  if (!nrow(df)) {
    return(NULL)
  }
  df <- df[order(df$padj, -abs(df$score)), , drop = FALSE]
  df <- df[!duplicated(df$term), , drop = FALSE]
  df <- utils::head(df, top_n)
  df <- df[order(df$score), , drop = FALSE]
  df$term <- factor(df$term, levels = unique(df$term))
  if (is_gsea) {
    df$direction <- ifelse(df$score >= 0, "Positive NES", "Negative NES")
    p <- ggplot2::ggplot(df, ggplot2::aes(x = score, y = term, fill = direction)) +
      ggplot2::geom_vline(xintercept = 0, color = "#B9C6CF", linewidth = 0.4) +
      ggplot2::geom_col(width = 0.64) +
      ggplot2::scale_fill_manual(name = NULL, values = c("Positive NES" = unname(colors[["positive"]]), "Negative NES" = unname(colors[["negative"]]))) +
      ggplot2::xlab("Normalized enrichment score (NES)")
  } else {
    direction <- if ("direction" %in% colnames(df)) df$direction[[1]] else "up"
    color <- if (identical(direction, "down")) colors[["negative"]] else colors[["positive"]]
    p <- ggplot2::ggplot(df, ggplot2::aes(x = score, y = term)) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = score, yend = term), color = "#DCE5E8", linewidth = 1.3)
    if ("Count" %in% colnames(df)) {
      df$Count <- suppressWarnings(as.numeric(df$Count))
      p <- p + ggplot2::geom_point(data = df, ggplot2::aes(size = Count), color = color) +
        ggplot2::scale_size_continuous(name = "Genes", range = c(3.5, 7))
    } else {
      p <- p + ggplot2::geom_point(size = 4.5, color = color)
    }
    p <- p + ggplot2::xlab("-log10 adjusted P-value")
  }
  p + ggplot2::ylab(NULL) +
    ggplot2::scale_y_discrete(labels = function(labels) vapply(labels, function(label) paste(strwrap(gsub("_", " ", label), width = 45), collapse = "\n"), character(1))) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "bottom", legend.title = ggplot2::element_text(size = 10),
                   axis.text.y = ggplot2::element_text(color = "#304957", size = 12),
                   axis.text.x = ggplot2::element_text(size = 11), plot.margin = ggplot2::margin(12, 18, 12, 12))
}

plot_decoupler_top <- function(dc_tbl, top_n = 20) {
  if (is.null(dc_tbl) || !nrow(dc_tbl) || !all(c("source", "score") %in% colnames(dc_tbl))) {
    return(NULL)
  }
  df <- dc_tbl
  df$score <- suppressWarnings(as.numeric(df$score))
  df <- df[is.finite(df$score), , drop = FALSE]
  if (!nrow(df)) {
    return(NULL)
  }
  df <- df[order(abs(df$score), decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, top_n)
  df$source <- as.character(df$source)
  df$source <- factor(df$source, levels = rev(unique(df$source)))
  ggplot2::ggplot(df, ggplot2::aes(x = source, y = score, fill = source_file)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_light() +
    ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::xlab(NULL) +
    ggplot2::ylab("Activity score")
}
