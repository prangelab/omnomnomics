explorer_palette_presets <- list(
  coastal = c("#14796F", "#D77A63", "#5479A2", "#B49B55", "#866694", "#68A99E", "#A56676", "#667680"),
  okabe_ito = c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#333333"),
  muted = c("#4878A8", "#D08A52", "#6C9A77", "#B26879", "#8A78A8", "#8F8877", "#83ADB5", "#B5A45A")
)

explorer_palette <- function(n, preset = "coastal", start = "#175A81", end = "#DA845C") {
  if (n < 1) {
    return(character(0))
  }
  if (identical(preset, "gradient")) {
    if (requireNamespace("coloRspacer", quietly = TRUE)) {
      return(coloRspacer::coloRspacer(startcolor = start, endcolor = end, steps = n))
    }
    return(grDevices::colorRampPalette(c(start, end))(n))
  }
  colors <- explorer_palette_presets[[preset]]
  if (is.null(colors)) {
    colors <- explorer_palette_presets$coastal
  }
  if (n <= length(colors)) colors[seq_len(n)] else grDevices::colorRampPalette(colors)(n)
}

explorer_color_id <- function(column, group) {
  identity <- charToRaw(enc2utf8(paste(column, group, sep = "\u001f")))
  paste0("group_", paste(sprintf("%02x", as.integer(identity)), collapse = ""))
}

explorer_color_input <- function(id, label, value) {
  shiny::tags$div(
    class = "explorer-color-row",
    shiny::tags$label(`for` = id, label),
    shiny::tags$div(
      class = "explorer-color-control",
      shiny::tags$input(id = id, type = "color", value = value, class = "explorer-color-input"),
      shiny::tags$span(class = "explorer-color-hex", value)
    )
  )
}

explorer_enrichment_title <- function(tbl) {
  category <- tools::toTitleCase(gsub("_", " ", sub("^msig_", "", tbl$category[[1]])))
  direction <- c(up = "Upregulated", down = "Downregulated", all = "All features", ranked = "Ranked genes")
  paste(category, tbl$method[[1]], unname(direction[tbl$direction[[1]]]), sep = " · ")
}
