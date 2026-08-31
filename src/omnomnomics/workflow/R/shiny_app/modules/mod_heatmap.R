mod_heatmap_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Heatmap"),
    shiny::fluidRow(
      shiny::column(3, shiny::numericInput(ns("top_n"), "Top features", value = 50, min = 2, step = 1)),
      shiny::column(3, shiny::checkboxInput(ns("zscore_rows"), "Row z-score", value = TRUE)),
      shiny::column(6, shiny::downloadButton(ns("download_heatmap"), "Download Heatmap PDF"))
    ),
    shiny::selectInput(ns("annotation_cols"), "Annotation columns", choices = character(0), multiple = TRUE),
    shiny::verbatimTextOutput(ns("status")),
    shiny::plotOutput(ns("heatmap_plot"), height = "700px")
  )
}

mod_heatmap_server <- function(id, project_index, selected_contrast, filtered_table, metadata) {
  shiny::moduleServer(id, function(input, output, session) {
    observe_annotation_choices <- shiny::observe({
      md <- metadata()
      if (is.null(md) || !ncol(md)) {
        shiny::updateSelectInput(session, "annotation_cols", choices = character(0), selected = character(0))
        return(invisible(NULL))
      }
      sample_col <- if ("sample_id" %in% colnames(md)) "sample_id" else colnames(md)[[1]]
      candidates <- setdiff(colnames(md), sample_col)
      candidates <- candidates[vapply(candidates, function(one_col) {
        vals <- as.character(md[[one_col]])
        vals <- vals[!is.na(vals) & nzchar(vals)]
        length(unique(vals)) >= 2
      }, logical(1))]
      selected <- input$annotation_cols
      selected <- selected[selected %in% candidates]
      shiny::updateSelectInput(session, "annotation_cols", choices = candidates, selected = selected)
      invisible(NULL)
    })

    heatmap_data <- shiny::reactive({
      contrast_obj <- selected_contrast()
      project_obj <- project_index()
      tbl <- filtered_table()
      md <- metadata()
      if (is.null(contrast_obj) || nrow(contrast_obj) == 0 || is.null(project_obj) || is.null(tbl) || !nrow(tbl)) {
        return(NULL)
      }

      vst_mat <- load_vst_matrix(project_obj)
      if (is.null(vst_mat) || !nrow(vst_mat)) {
        return(NULL)
      }

      label_col <- resolve_gene_label_column(tbl)
      labels <- resolve_heatmap_genes(tbl, vst_mat, top_n = input$top_n)
      if (length(labels) < 2) {
        return(list(mat = NULL, n_genes = 0L, n_samples = 0L, reason = paste0("No overlap between selected DE rows (first key column: ", label_col, ") and qc/vst_matrix.tsv rownames.")))
      }

      sub <- vst_mat[labels, , drop = FALSE]
      label_map <- resolve_gene_display_labels(tbl, labels)
      rownames(sub) <- make.unique(as.character(label_map[labels]))
      if (isTRUE(input$zscore_rows)) {
        sub <- row_zscore(sub)
      }

      annotation_df <- NULL
      selected_ann_cols <- input$annotation_cols
      if (!is.null(md) && nrow(md) > 0 && length(selected_ann_cols) > 0) {
        sample_col <- if ("sample_id" %in% colnames(md)) "sample_id" else colnames(md)[[1]]
        md_samples <- as.character(md[[sample_col]])
        md_idx <- match(colnames(sub), md_samples)
        valid <- !is.na(md_idx)
        if (any(valid)) {
          sub <- sub[, valid, drop = FALSE]
          md_idx <- md_idx[valid]
          ann_cols <- selected_ann_cols[selected_ann_cols %in% colnames(md)]
          if (length(ann_cols) > 0) {
            annotation_df <- md[md_idx, ann_cols, drop = FALSE]
            rownames(annotation_df) <- colnames(sub)
          }
        }
      }

      if (ncol(sub) < 2) {
        return(list(mat = NULL, n_genes = 0L, n_samples = 0L, reason = "Not enough samples after metadata annotation matching."))
      }

      list(mat = sub, annotation_df = annotation_df, n_genes = nrow(sub), n_samples = ncol(sub), reason = "")
    })

    output$status <- shiny::renderText({
      hmd <- heatmap_data()
      if (is.null(hmd)) {
        return("Heatmap unavailable. Requires qc/vst_matrix.tsv and overlapping gene labels in selected table.")
      }
      if (is.null(hmd$mat) || !nrow(hmd$mat)) {
        return(hmd$reason)
      }
      paste0("Rendering ", hmd$n_genes, " features across ", hmd$n_samples, " samples.")
    })

    output$heatmap_plot <- shiny::renderPlot({
      hmd <- heatmap_data()
      shiny::req(hmd)
      shiny::req(hmd$mat)
      pheatmap::pheatmap(
        hmd$mat,
        scale = "none",
        annotation_col = hmd$annotation_df,
        show_colnames = FALSE,
        show_rownames = TRUE
      )
    })

    output$download_heatmap <- shiny::downloadHandler(
      filename = function() {
        contrast_obj <- selected_contrast()
        label <- if (is.null(contrast_obj) || nrow(contrast_obj) == 0) "contrast" else as.character(contrast_obj$contrast_id[[1]])
        paste0(label, ".sig_heatmap.shiny.pdf")
      },
      content = function(file) {
        hmd <- heatmap_data()
        if (is.null(hmd) || is.null(hmd$mat) || !nrow(hmd$mat)) {
          grDevices::pdf(file, width = 8, height = 6)
          graphics::plot.new()
          graphics::text(0.5, 0.5, "No heatmap data available.")
          grDevices::dev.off()
          return(invisible(NULL))
        }
        grDevices::pdf(file, width = 10, height = 10)
            pheatmap::pheatmap(
              hmd$mat,
              scale = "none",
              annotation_col = hmd$annotation_df,
              show_colnames = FALSE,
              show_rownames = TRUE
            )
        grDevices::dev.off()
        invisible(NULL)
      }
    )

    session$onSessionEnded(function() {
      observe_annotation_choices$destroy()
    })
  })
}
