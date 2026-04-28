read_de_table <- function(table_path) {
  if (is.null(table_path) || !nzchar(table_path) || !file.exists(table_path)) {
    stop("DE table not found: ", table_path)
  }
  utils::read.table(
    table_path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

detect_gene_column <- function(tbl) {
  candidates <- c("symbol", "gene_symbol", "SYMBOL", "external_gene_name", "gene", "Geneid", "gene_id", "Geneid_noVersion")
  hit <- candidates[candidates %in% colnames(tbl)]
  if (length(hit) > 0) {
    return(hit[[1]])
  }
  colnames(tbl)[[1]]
}

filter_de_table <- function(tbl, padj_max = 1, lfc_min = 0, base_mean_min = 0, gene_query = "") {
  if (!("padj" %in% colnames(tbl)) || !("log2FoldChange" %in% colnames(tbl))) {
    return(tbl)
  }
  keep <- rep(TRUE, nrow(tbl))
  keep <- keep & !is.na(tbl$padj) & tbl$padj <= as.numeric(padj_max)
  keep <- keep & !is.na(tbl$log2FoldChange) & abs(tbl$log2FoldChange) >= as.numeric(lfc_min)
  if ("baseMean" %in% colnames(tbl)) {
    keep <- keep & !is.na(tbl$baseMean) & tbl$baseMean >= as.numeric(base_mean_min)
  }

  query <- trimws(as.character(gene_query))
  if (nzchar(query)) {
    gene_col <- detect_gene_column(tbl)
    query_idx <- grepl(query, as.character(tbl[[gene_col]]), ignore.case = TRUE, fixed = FALSE)
    keep <- keep & query_idx
  }
  tbl[keep, , drop = FALSE]
}

summarize_contrast <- function(tbl, alpha = 0.05, lfc_cutoff = 1) {
  if (nrow(tbl) == 0 || !("padj" %in% colnames(tbl)) || !("log2FoldChange" %in% colnames(tbl))) {
    return(data.frame(
      n_rows = nrow(tbl),
      sig_total = 0L,
      sig_up = 0L,
      sig_down = 0L,
      stringsAsFactors = FALSE
    ))
  }
  sig <- !is.na(tbl$padj) & tbl$padj <= as.numeric(alpha) &
    !is.na(tbl$log2FoldChange) & abs(tbl$log2FoldChange) >= as.numeric(lfc_cutoff)
  sig_up <- sum(sig & tbl$log2FoldChange > 0, na.rm = TRUE)
  sig_down <- sum(sig & tbl$log2FoldChange < 0, na.rm = TRUE)
  data.frame(
    n_rows = nrow(tbl),
    sig_total = as.integer(sig_up + sig_down),
    sig_up = as.integer(sig_up),
    sig_down = as.integer(sig_down),
    stringsAsFactors = FALSE
  )
}
