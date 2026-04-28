resolve_gene_label_column <- function(tbl) {
  for (candidate in c("gene_symbol", "symbol", "SYMBOL", "gene_id", "Geneid", "gene")) {
    if (candidate %in% colnames(tbl)) {
      return(candidate)
    }
  }
  colnames(tbl)[[1]]
}

top_genes_for_boxplot <- function(tbl, n_top = 10, alpha = 0.05, lfc_cutoff = 1) {
  if (!all(c("padj", "log2FoldChange") %in% colnames(tbl))) {
    return(character(0))
  }
  sub <- tbl[!is.na(tbl$padj) & tbl$padj <= alpha & abs(tbl$log2FoldChange) >= lfc_cutoff, , drop = FALSE]
  if (!nrow(sub)) {
    return(character(0))
  }
  sub <- sub[order(sub$padj, -abs(sub$log2FoldChange)), , drop = FALSE]
  label_col <- resolve_gene_label_column(sub)
  unique(as.character(utils::head(sub[[label_col]], n_top)))
}
