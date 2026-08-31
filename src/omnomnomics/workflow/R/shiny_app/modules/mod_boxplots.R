mod_boxplots_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Feature Boxplots"),
    shiny::fluidRow(
      shiny::column(3, shiny::numericInput(ns("top_n"), "Top features", value = 10, min = 1, step = 1)),
      shiny::column(3, shiny::radioButtons(ns("gene_mode"), "Feature source", choices = c("Top N" = "top_n", "Custom feature" = "custom"), selected = "top_n", inline = TRUE)),
      shiny::column(3, shiny::selectInput(ns("group_col"), "Group column", choices = character(0))),
      shiny::column(3, shiny::downloadButton(ns("download_boxplot"), "Download Boxplot PDF"))
    ),
    shiny::verbatimTextOutput(ns("status")),
    shiny::selectizeInput(ns("gene_pick"), "Feature", choices = character(0), options = list(placeholder = "Select or search a feature")),
    shiny::plotOutput(ns("boxplot"), height = "520px")
  )
}

mod_boxplots_server <- function(id, project_index, selected_contrast, filtered_table, metadata) {
  shiny::moduleServer(id, function(input, output, session) {
    gene_choices_state <- shiny::reactiveVal(list(key = NULL, selected = NULL))
    observe_group_cols <- shiny::observe({
      md <- metadata()
      if (is.null(md) || !ncol(md)) {
        shiny::updateSelectInput(session, "group_col", choices = character(0), selected = character(0))
        return(invisible(NULL))
      }
      categorical <- colnames(md)[vapply(md, function(col) {
        !is.numeric(col) && length(unique(as.character(col[!is.na(col)]))) >= 2
      }, logical(1))]
      if (length(categorical) == 0) {
        categorical <- colnames(md)
      }
      selected <- input$group_col
      if (!nzchar(selected) || !(selected %in% categorical)) {
        selected <- categorical[[1]]
      }
      shiny::updateSelectInput(session, "group_col", choices = categorical, selected = selected)
      invisible(NULL)
    })

    all_gene_ids <- shiny::reactive({
      project_obj <- project_index()
      tbl <- filtered_table()
      if (is.null(project_obj) || is.null(tbl) || !nrow(tbl)) {
        return(character(0))
      }
      vst_mat <- load_vst_matrix(project_obj)
      if (is.null(vst_mat) || !nrow(vst_mat)) {
        return(character(0))
      }
      id_col <- if ("gene_id" %in% colnames(tbl)) "gene_id" else colnames(tbl)[[1]]
      ids <- as.character(tbl[[id_col]])
      ids <- ids[!is.na(ids) & nzchar(ids)]
      ids <- unique(ids)
      intersect(ids, rownames(vst_mat))
    })

    top_gene_ids <- shiny::reactive({
      project_obj <- project_index()
      tbl <- filtered_table()
      if (is.null(tbl) || !nrow(tbl) || is.null(project_obj)) {
        return(character(0))
      }
      vst_mat <- load_vst_matrix(project_obj)
      if (is.null(vst_mat) || !nrow(vst_mat)) {
        return(character(0))
      }
      resolved <- resolve_heatmap_genes(tbl, vst_mat, top_n = input$top_n)
      as.character(resolved)
    })

    gene_choice_map <- shiny::reactive({
      project_obj <- project_index()
      tbl <- filtered_table()
      mode <- if (nzchar(input$gene_mode)) input$gene_mode else "top_n"
      ids <- if (identical(mode, "custom")) all_gene_ids() else top_gene_ids()
      if (is.null(project_obj) || is.null(tbl) || !length(ids)) {
        return(list(choices = character(0), selected = character(0), id_by_label = character(0)))
      }
      vst_mat <- load_vst_matrix(project_obj)
      if (is.null(vst_mat)) {
        return(list(choices = ids, selected = ids[[1]], id_by_label = stats::setNames(ids, ids)))
      }
      labels <- as.character(resolve_gene_display_labels(tbl, ids)[ids])
      labels[is.na(labels) | !nzchar(labels)] <- ids[is.na(labels) | !nzchar(labels)]
      display <- make.unique(labels)
      id_by_label <- stats::setNames(ids, display)
      list(
        choices = display,
        selected = if (length(display)) display[[1]] else character(0),
        id_by_label = id_by_label
      )
    })

    observe_gene_choices <- shiny::observeEvent(gene_choice_map(), {
      choice_map <- gene_choice_map()
      choices <- choice_map$choices
      selected <- isolate(input$gene_pick)
      state <- gene_choices_state()
      choices_key <- paste(choices, collapse = "\u001f")
      if (length(choices) == 0) {
        if (!is.null(state$key) || (length(selected) > 0 && nzchar(selected))) {
          shiny::freezeReactiveValue(input, "gene_pick")
          shiny::updateSelectizeInput(session, "gene_pick", choices = character(0), selected = character(0), server = TRUE)
          gene_choices_state(list(key = NULL, selected = NULL))
        }
        return(invisible(NULL))
      }

      if (identical(state$key, choices_key) && length(selected) > 0 && nzchar(selected) && (selected %in% choices)) {
        return(invisible(NULL))
      }

      if (!(length(selected) > 0 && nzchar(selected) && (selected %in% choices))) {
        if (!is.null(state$selected) && nzchar(state$selected) && (state$selected %in% choices)) {
          selected <- state$selected
        } else {
          selected <- choices[[1]]
        }
      }
      shiny::freezeReactiveValue(input, "gene_pick")
      shiny::updateSelectizeInput(session, "gene_pick", choices = choices, selected = selected, server = TRUE)
      gene_choices_state(list(key = choices_key, selected = selected))
      invisible(NULL)
    }, ignoreInit = FALSE)

    observe_gene_pick <- shiny::observeEvent(input$gene_pick, {
      picked <- input$gene_pick
      if (!length(picked) || !nzchar(picked)) {
        return(invisible(NULL))
      }
      state <- gene_choices_state()
      gene_choices_state(list(key = state$key, selected = picked))
      invisible(NULL)
    }, ignoreInit = TRUE)

    boxplot_df <- shiny::reactive({
      project_obj <- project_index()
      md <- metadata()
      choice_map <- gene_choice_map()
      gene_label <- input$gene_pick
      if (is.null(choice_map$id_by_label) || !length(choice_map$id_by_label)) {
        return(NULL)
      }
      if (!nzchar(gene_label) || !(gene_label %in% names(choice_map$id_by_label))) {
        gene_label <- as.character(choice_map$selected)
      }
      if (!nzchar(gene_label) || !(gene_label %in% names(choice_map$id_by_label))) {
        return(NULL)
      }
      gene <- as.character(unname(choice_map$id_by_label[[gene_label]]))
      if (!length(gene) || is.na(gene[[1]]) || !nzchar(gene[[1]])) {
        return(NULL)
      }
      gene <- gene[[1]]
      group_col <- input$group_col
      if (is.null(project_obj) || is.null(md) || !nzchar(gene) || !nzchar(group_col) || !(group_col %in% colnames(md))) {
        return(NULL)
      }

      vst_mat <- load_vst_matrix(project_obj)
      if (is.null(vst_mat)) {
        return(NULL)
      }
      gene_idx <- match(gene, rownames(vst_mat))
      if (is.na(gene_idx)) {
        return(NULL)
      }

      sample_col <- if ("sample_id" %in% colnames(md)) "sample_id" else colnames(md)[[1]]
      sample_ids <- as.character(md[[sample_col]])
      sample_ids <- sample_ids[nzchar(sample_ids) & !is.na(sample_ids)]
      sample_ids <- intersect(sample_ids, colnames(vst_mat))
      if (length(sample_ids) < 2) {
        return(NULL)
      }

      md_idx <- match(sample_ids, as.character(md[[sample_col]]))
      expr <- tryCatch(
        as.numeric(vst_mat[gene_idx, sample_ids, drop = TRUE]),
        error = function(e) NULL
      )
      if (is.null(expr) || length(expr) != length(sample_ids)) {
        return(NULL)
      }
      plot_df <- data.frame(
        sample_id = sample_ids,
        expression = expr,
        group = as.character(md[[group_col]][md_idx]),
        gene_label = if (nzchar(gene_label)) gene_label else gene,
        stringsAsFactors = FALSE
      )
      plot_df <- plot_df[!is.na(plot_df$group) & nzchar(plot_df$group), , drop = FALSE]
      if (nrow(plot_df) < 2) {
        return(NULL)
      }
      plot_df
    })

    make_boxplot <- function(df, gene_label, group_col) {
      ggplot2::ggplot(df, ggplot2::aes(x = group, y = expression, color = group)) +
        ggplot2::geom_boxplot(outlier.shape = NA) +
        ggplot2::geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
        ggplot2::xlab(group_col) +
        ggplot2::ylab("VST abundance") +
        ggplot2::ggtitle(gene_label) +
        ggplot2::theme_light() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none")
    }

    output$status <- shiny::renderText({
      df <- boxplot_df()
      if (is.null(df)) {
        return("Boxplot unavailable for the current selection. Check that qc/vst_matrix.tsv exists and metadata sample IDs overlap VST column names.")
      }
      paste0("Rendering ", nrow(df), " samples across ", length(unique(df$group)), " groups.")
    })

    output$boxplot <- shiny::renderPlot({
      df <- boxplot_df()
      shiny::req(df)
      p <- make_boxplot(df, unique(df$gene_label)[[1]], input$group_col)
      p
    })

    output$download_boxplot <- shiny::downloadHandler(
      filename = function() {
        contrast_obj <- selected_contrast()
        contrast_label <- if (is.null(contrast_obj) || nrow(contrast_obj) == 0) "contrast" else as.character(contrast_obj$contrast_id[[1]])
        gene <- if (nzchar(input$gene_pick)) input$gene_pick else "feature"
        paste0(contrast_label, ".", gene, ".boxplot.shiny.pdf")
        },
      content = function(file) {
        df <- boxplot_df()
        if (is.null(df)) {
          grDevices::pdf(file, width = 7, height = 5)
          graphics::plot.new()
          graphics::text(0.5, 0.5, "No boxplot data available.")
          grDevices::dev.off()
          return(invisible(NULL))
        }
        p <- make_boxplot(df, unique(df$gene_label)[[1]], input$group_col)
        ggplot2::ggsave(file, p, width = 7, height = 5)
        invisible(NULL)
      }
    )

    session$onSessionEnded(function() {
      observe_group_cols$destroy()
      observe_gene_choices$destroy()
      observe_gene_pick$destroy()
    })
  })
}
