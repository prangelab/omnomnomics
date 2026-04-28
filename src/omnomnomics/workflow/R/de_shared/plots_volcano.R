build_volcano_data <- function(tbl, alpha = 0.05, lfc_cutoff = 1) {
  required_cols <- c("log2FoldChange", "padj")
  if (!all(required_cols %in% colnames(tbl))) {
    return(NULL)
  }
  out <- tbl
  out$padj[is.na(out$padj)] <- 1
  out$log2FoldChange[is.na(out$log2FoldChange)] <- 0
  out$neglog10padj <- -log10(pmax(out$padj, .Machine$double.xmin))
  out$is_sig <- out$padj <= alpha & abs(out$log2FoldChange) >= lfc_cutoff
  out
}

plot_volcano <- function(vdf, alpha = 0.05, lfc_cutoff = 1) {
  if (is.null(vdf) || !nrow(vdf)) {
    return(NULL)
  }
  ggplot2::ggplot(vdf, ggplot2::aes(x = log2FoldChange, y = neglog10padj, color = is_sig)) +
    ggplot2::geom_point(alpha = 0.7, size = 1.2) +
    ggplot2::scale_color_manual(values = c("FALSE" = "grey65", "TRUE" = "firebrick")) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, 0, lfc_cutoff), linetype = c("dashed", "solid", "dashed")) +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
    ggplot2::theme_light() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none") +
    ggplot2::xlab("log2 Fold Change") +
    ggplot2::ylab("-log10 adjusted P-value")
}
