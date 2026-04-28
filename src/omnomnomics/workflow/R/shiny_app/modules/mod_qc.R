mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Shared QC"),
    shiny::fluidRow(
      shiny::column(3, shiny::numericInput(ns("top_n"), "Top variable genes", value = 1000, min = 100, step = 100)),
      shiny::column(3, shiny::selectInput(ns("pc_x"), "X axis", choices = c("PC1", "PC2", "PC3"), selected = "PC1")),
      shiny::column(3, shiny::selectInput(ns("pc_y"), "Y axis", choices = c("PC1", "PC2", "PC3"), selected = "PC2")),
      shiny::column(3, shiny::downloadButton(ns("download_pca"), "Download PCA PDF"))
    ),
    shiny::fluidRow(
      shiny::column(6, shiny::selectInput(ns("color_by"), "Color by", choices = character(0))),
      shiny::column(6, shiny::selectInput(ns("shape_by"), "Shape by", choices = character(0)))
    ),
    shiny::selectInput(ns("topvar_annotation_cols"), "TopVar heatmap annotation columns", choices = character(0), multiple = TRUE),
    shiny::verbatimTextOutput(ns("status")),
    shiny::plotOutput(ns("pca_plot"), height = "560px"),
    shiny::h4("Sample Distance Heatmap"),
    shiny::plotOutput(ns("distance_plot"), height = "620px"),
    shiny::h4("Top Variable Genes Heatmap"),
    shiny::plotOutput(ns("topvar_heatmap"), height = "700px"),
    shiny::h4("Filtering Summary"),
    shiny::tableOutput(ns("filtering_tbl")),
    shiny::h4("Cook's Global Summary"),
    shiny::tableOutput(ns("cooks_global_tbl")),
    shiny::h4("Cook's Sample Summary"),
    shiny::tableOutput(ns("cooks_sample_tbl")),
    shiny::h4("Cook's Outlier Counts"),
    shiny::plotOutput(ns("cooks_counts_plot"), height = "420px"),
    shiny::h4("Cook's q99 per Sample"),
    shiny::plotOutput(ns("cooks_q99_plot"), height = "420px")
  )
}

mod_qc_server <- function(id, project_index, metadata) {
  shiny::moduleServer(id, function(input, output, session) {
    observe_metadata_cols <- shiny::observe({
      md <- metadata()
      if (is.null(md) || !ncol(md)) {
        shiny::updateSelectInput(session, "color_by", choices = character(0), selected = character(0))
        shiny::updateSelectInput(session, "shape_by", choices = character(0), selected = character(0))
        return(invisible(NULL))
      }
      sample_col <- if ("sample_id" %in% colnames(md)) "sample_id" else colnames(md)[[1]]
      choices <- setdiff(colnames(md), sample_col)
      choices <- choices[vapply(choices, function(one_col) {
        vals <- as.character(md[[one_col]])
        vals <- vals[!is.na(vals) & nzchar(vals)]
        length(unique(vals)) >= 2
      }, logical(1))]
      selected_color <- if (length(choices) > 0) choices[[1]] else character(0)
            selected_shape <- if (length(choices) > 1) choices[[2]] else character(0)
            empty_choice <- setNames("", "")
            shiny::updateSelectInput(session, "color_by", choices = c(empty_choice, choices), selected = selected_color)
            shiny::updateSelectInput(session, "shape_by", choices = c(empty_choice, choices), selected = selected_shape)
            shiny::updateSelectInput(session, "topvar_annotation_cols", choices = choices, selected = choices)
            invisible(NULL)
        })

    pca_payload <- shiny::reactive({
      project_obj <- project_index()
      md <- metadata()
      if (is.null(project_obj) || is.null(md)) {
        return(NULL)
      }
      vst_mat <- load_vst_matrix(project_obj)
      if (is.null(vst_mat)) {
        return(NULL)
      }
      scores <- build_pca_scores(vst_mat, md, top_n = input$top_n)
      if (is.null(scores)) {
        return(NULL)
      }
            p <- tryCatch(
              plot_pca_scores(
                scores_df = scores,
                pc_x = input$pc_x,
                pc_y = input$pc_y,
                color_by = input$color_by,
                shape_by = input$shape_by
              ),
              error = function(e) e
            )
            if (inherits(p, "error")) {
              return(list(scores = scores, plot = NULL, error = conditionMessage(p)))
            }
            list(scores = scores, plot = p, error = "")
        })

    observe_pc_choices <- shiny::observe({
      payload <- pca_payload()
      if (is.null(payload) || is.null(payload$scores)) {
        return(invisible(NULL))
      }
      pc_choices <- grep("^PC[0-9]+$", colnames(payload$scores), value = TRUE)
      if (length(pc_choices) < 2) {
        return(invisible(NULL))
      }
      selected_x <- input$pc_x
      selected_y <- input$pc_y
      if (!nzchar(selected_x) || !(selected_x %in% pc_choices)) {
        selected_x <- pc_choices[[1]]
      }
      if (!nzchar(selected_y) || !(selected_y %in% pc_choices) || identical(selected_y, selected_x)) {
        selected_y <- if (length(pc_choices) > 1) pc_choices[[2]] else pc_choices[[1]]
      }
      shiny::updateSelectInput(session, "pc_x", choices = pc_choices, selected = selected_x)
      shiny::updateSelectInput(session, "pc_y", choices = pc_choices, selected = selected_y)
      invisible(NULL)
    })

        output$status <- shiny::renderText({
            payload <- pca_payload()
            if (is.null(payload)) {
                return("PCA unavailable. Requires qc/vst_matrix.tsv and metadata sample overlap.")
            }
            if (is.null(payload$plot)) {
              return(paste0("PCA render error: ", payload$error))
            }
            paste0("PCA rendered for ", nrow(payload$scores), " samples.")
        })

        output$pca_plot <- shiny::renderPlot({
            payload <- pca_payload()
            shiny::req(payload)
            shiny::req(payload$plot)
            tryCatch(
              print(payload$plot),
              error = function(e) {
                graphics::plot.new()
                graphics::text(0.5, 0.5, paste0("PCA render error: ", conditionMessage(e)))
              }
            )
        })

    output$download_pca <- shiny::downloadHandler(
      filename = function() {
        paste0("qc.pca.", input$pc_x, "_", input$pc_y, ".pdf")
      },
      content = function(file) {
        payload <- pca_payload()
        if (is.null(payload) || is.null(payload$plot)) {
          grDevices::pdf(file, width = 8, height = 6)
          graphics::plot.new()
          graphics::text(0.5, 0.5, "No PCA available.")
          grDevices::dev.off()
          return(invisible(NULL))
        }
        ggplot2::ggsave(file, payload$plot, width = 8, height = 6)
        invisible(NULL)
      }
    )

    output$distance_plot <- shiny::renderPlot({
      project_obj <- project_index()
      md <- metadata()
      shiny::req(project_obj)
      shiny::req(md)
      vst_mat <- load_vst_matrix(project_obj)
      shiny::req(vst_mat)
      dmat <- build_distance_matrix(vst_mat, md, top_n = input$top_n)
      shiny::req(dmat)
      pheatmap::pheatmap(
        dmat,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = TRUE,
        show_rownames = TRUE
      )
    })

    output$topvar_heatmap <- shiny::renderPlot({
      project_obj <- project_index()
      md <- metadata()
      shiny::req(project_obj)
      shiny::req(md)
      vst_mat <- load_vst_matrix(project_obj)
      shiny::req(vst_mat)
      top_mat <- build_topvar_heatmap_matrix(vst_mat, md, top_n = input$top_n)
      shiny::req(top_mat)
      annotation_df <- NULL
      selected_ann <- input$topvar_annotation_cols
      if (!is.null(selected_ann) && length(selected_ann) > 0) {
        sample_col <- if ("sample_id" %in% colnames(md)) "sample_id" else colnames(md)[[1]]
        md_idx <- match(colnames(top_mat), as.character(md[[sample_col]]))
        valid <- !is.na(md_idx)
        if (any(valid)) {
          top_mat <- top_mat[, valid, drop = FALSE]
          md_idx <- md_idx[valid]
          ann_cols <- selected_ann[selected_ann %in% colnames(md)]
          if (length(ann_cols) > 0) {
            annotation_df <- md[md_idx, ann_cols, drop = FALSE]
            rownames(annotation_df) <- colnames(top_mat)
          }
        }
      }
      pheatmap::pheatmap(
        row_zscore(top_mat),
        annotation_col = annotation_df,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_colnames = FALSE,
        show_rownames = FALSE
      )
    })

    output$filtering_tbl <- shiny::renderTable({
      project_obj <- project_index()
      shiny::req(project_obj)
      tbl <- read_qc_table_if_exists(project_obj$qc_dir, "filtering_summary.tsv")
      shiny::req(tbl)
      tbl
    })

    output$cooks_global_tbl <- shiny::renderTable({
      project_obj <- project_index()
      shiny::req(project_obj)
      tbl <- read_qc_table_if_exists(project_obj$qc_dir, "cooks_distance_global_summary.tsv")
      shiny::req(tbl)
      tbl
    })

        output$cooks_sample_tbl <- shiny::renderTable({
            project_obj <- project_index()
            shiny::req(project_obj)
            tbl <- read_qc_table_if_exists(project_obj$qc_dir, "cooks_distance_sample_summary.tsv")
            shiny::req(tbl)
            utils::head(tbl, 40)
        })

        output$cooks_counts_plot <- shiny::renderPlot({
          project_obj <- project_index()
          shiny::req(project_obj)
          tbl <- read_qc_table_if_exists(project_obj$qc_dir, "cooks_distance_sample_summary.tsv")
          shiny::req(tbl)
          shiny::req(all(c("sample", "genes_cooks_gt_10") %in% colnames(tbl)))
          plot_tbl <- tbl[order(tbl$genes_cooks_gt_10, decreasing = TRUE), , drop = FALSE]
          plot_tbl$sample <- factor(plot_tbl$sample, levels = plot_tbl$sample)
          ggplot2::ggplot(plot_tbl, ggplot2::aes(x = sample, y = genes_cooks_gt_10)) +
            ggplot2::geom_col(fill = "steelblue") +
            ggplot2::coord_flip() +
            ggplot2::theme_light() +
            ggplot2::theme(panel.grid = ggplot2::element_blank()) +
            ggplot2::xlab("Sample") +
            ggplot2::ylab("Genes with Cook's distance > 10")
        })

        output$cooks_q99_plot <- shiny::renderPlot({
          project_obj <- project_index()
          shiny::req(project_obj)
          tbl <- read_qc_table_if_exists(project_obj$qc_dir, "cooks_distance_sample_summary.tsv")
          shiny::req(tbl)
          shiny::req(all(c("sample", "q99_cooks") %in% colnames(tbl)))
          plot_tbl <- tbl[order(tbl$q99_cooks, decreasing = TRUE), , drop = FALSE]
          plot_tbl$sample <- factor(plot_tbl$sample, levels = plot_tbl$sample)
          ggplot2::ggplot(plot_tbl, ggplot2::aes(x = sample, y = q99_cooks)) +
            ggplot2::geom_col(fill = "darkorange3") +
            ggplot2::coord_flip() +
            ggplot2::theme_light() +
            ggplot2::theme(panel.grid = ggplot2::element_blank()) +
            ggplot2::xlab("Sample") +
            ggplot2::ylab("99th percentile Cook's distance")
        })

        session$onSessionEnded(function() {
            observe_metadata_cols$destroy()
            observe_pc_choices$destroy()
        })
  })
}
