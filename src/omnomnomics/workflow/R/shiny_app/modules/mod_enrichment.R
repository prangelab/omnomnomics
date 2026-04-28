mod_enrichment_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Enrichment"),
    shiny::fluidRow(
      shiny::column(4, shiny::numericInput(ns("cp_top_n"), "Top CP terms", value = 15, min = 5, step = 1)),
      shiny::column(4, shiny::numericInput(ns("dc_top_n"), "Top decoupleR regulators", value = 20, min = 5, step = 1)),
      shiny::column(4, shiny::downloadButton(ns("download_cp_plot"), "Download CP Plot PDF"))
    ),
    shiny::checkboxInput(ns("dc_include_all_genes"), "Include decoupleR all-genes outputs", value = FALSE),
    shiny::verbatimTextOutput(ns("status")),
    shiny::h4("clusterProfiler (Hallmark/Reactome/etc.)"),
    shiny::uiOutput(ns("cp_plots_ui")),
    shiny::tableOutput(ns("cp_table")),
    shiny::h4("decoupleR"),
    shiny::uiOutput(ns("dc_plots_ui")),
    shiny::tableOutput(ns("dc_table"))
  )
}

mod_enrichment_server <- function(id, selected_contrast) {
  shiny::moduleServer(id, function(input, output, session) {
    sanitize_id <- function(x) {
      gsub("[^A-Za-z0-9_]+", "_", tolower(as.character(x)))
    }

    annotate_categories <- function(tbl, contrast_id, kind = c("cp", "dc")) {
      kind <- match.arg(kind)
      if (is.null(tbl) || !nrow(tbl) || !("source_file" %in% colnames(tbl))) {
        return(tbl)
      }
      out <- tbl
      base <- sub("\\.tsv$", "", as.character(out$source_file))
      base <- sub(paste0("^", contrast_id, "\\."), "", base)
      if (kind == "cp") {
        category <- ifelse(
          grepl("\\.(up|down)\\.ORA$", base),
          sub("\\.(up|down)\\.ORA$", "", base),
          ifelse(grepl("\\.GSEA$", base), sub("\\.GSEA$", "", base), "other")
        )
      } else {
        scope <- ifelse(grepl("^all_genes\\.", base), "all_genes", "DE_genes")
        category <- sub("^(all_genes|DE_genes)\\.", "", base)
        out$scope <- scope
      }
      out$category <- category
      out
    }

    enrichment_payload <- shiny::reactive({
      contrast_obj <- selected_contrast()
      if (is.null(contrast_obj) || nrow(contrast_obj) == 0) {
        return(NULL)
      }
      contrast_id <- as.character(contrast_obj$contrast_id[[1]])
      contrast_dir <- as.character(contrast_obj$contrast_dir[[1]])
      raw <- read_enrichment_tables(contrast_dir)
      raw$clusterprofiler <- annotate_categories(raw$clusterprofiler, contrast_id, kind = "cp")
      raw$decoupler <- annotate_categories(raw$decoupler, contrast_id, kind = "dc")
      raw
    })

    cp_split <- shiny::reactive({
      payload <- enrichment_payload()
      if (is.null(payload) || is.null(payload$clusterprofiler) || !nrow(payload$clusterprofiler)) {
        return(list())
      }
      split(payload$clusterprofiler, payload$clusterprofiler$category)
    })

    dc_split <- shiny::reactive({
      payload <- decoupler_payload()
      if (is.null(payload) || !nrow(payload)) {
        return(list())
      }
      split(payload, payload$category)
    })

    decoupler_payload <- shiny::reactive({
      payload <- enrichment_payload()
      if (is.null(payload) || is.null(payload$decoupler) || !nrow(payload$decoupler)) {
        return(NULL)
      }
      dc_tbl <- payload$decoupler
      include_all <- isTRUE(input$dc_include_all_genes)
      if (!include_all && ("scope" %in% colnames(dc_tbl))) {
        dc_tbl <- dc_tbl[dc_tbl$scope == "DE_genes", , drop = FALSE]
      }
      dc_tbl
    })

    output$status <- shiny::renderText({
      payload <- enrichment_payload()
      if (is.null(payload)) {
        return("Select a contrast to load enrichment outputs.")
      }
      cp_n <- if (is.null(payload$clusterprofiler)) 0 else nrow(payload$clusterprofiler)
      dc_tbl <- decoupler_payload()
      dc_n <- if (is.null(dc_tbl)) 0 else nrow(dc_tbl)
      dc_mode <- if (isTRUE(input$dc_include_all_genes)) "all_genes+DE_genes" else "DE_genes_only"
      paste0("Loaded enrichment rows: clusterProfiler=", cp_n, ", decoupleR=", dc_n, " (", dc_mode, ")")
    })

    output$cp_plots_ui <- shiny::renderUI({
      split_list <- cp_split()
      if (length(split_list) == 0) {
        return(shiny::tags$em("No clusterProfiler enrichment tables found for this contrast."))
      }
      ui_items <- lapply(names(split_list), function(cat_name) {
        plot_id <- paste0("cp_plot_", sanitize_id(cat_name))
        shiny::tagList(
          shiny::h5(cat_name),
          shiny::plotOutput(session$ns(plot_id), height = "420px")
        )
      })
      do.call(shiny::tagList, ui_items)
    })

    observe_cp_plots <- shiny::observe({
      split_list <- cp_split()
      if (length(split_list) == 0) {
        return(invisible(NULL))
      }
      for (cat_name in names(split_list)) {
        local({
          one_cat <- cat_name
          one_df <- split_list[[one_cat]]
          plot_id <- paste0("cp_plot_", sanitize_id(one_cat))
          output[[plot_id]] <- shiny::renderPlot({
            p <- plot_cp_top_terms(one_df, top_n = input$cp_top_n)
            shiny::req(p)
            p + ggplot2::ggtitle(one_cat)
          })
        })
      }
      invisible(NULL)
    })

    output$cp_table <- shiny::renderTable({
      payload <- enrichment_payload()
      shiny::req(payload)
      shiny::req(payload$clusterprofiler)
      utils::head(payload$clusterprofiler, 40)
    })

    output$dc_plots_ui <- shiny::renderUI({
      split_list <- dc_split()
      if (length(split_list) == 0) {
        return(shiny::tags$em("No decoupleR enrichment tables found for this contrast."))
      }
      ui_items <- lapply(names(split_list), function(cat_name) {
        plot_id <- paste0("dc_plot_", sanitize_id(cat_name))
        shiny::tagList(
          shiny::h5(cat_name),
          shiny::plotOutput(session$ns(plot_id), height = "420px")
        )
      })
      do.call(shiny::tagList, ui_items)
    })

    observe_dc_plots <- shiny::observe({
      split_list <- dc_split()
      if (length(split_list) == 0) {
        return(invisible(NULL))
      }
      for (cat_name in names(split_list)) {
        local({
          one_cat <- cat_name
          one_df <- split_list[[one_cat]]
          plot_id <- paste0("dc_plot_", sanitize_id(one_cat))
          output[[plot_id]] <- shiny::renderPlot({
            p <- plot_decoupler_top(one_df, top_n = input$dc_top_n)
            shiny::req(p)
            p + ggplot2::ggtitle(one_cat)
          })
        })
      }
      invisible(NULL)
    })

    output$dc_table <- shiny::renderTable({
      dc_tbl <- decoupler_payload()
      shiny::req(dc_tbl)
      utils::head(dc_tbl, 40)
    })

    output$download_cp_plot <- shiny::downloadHandler(
      filename = function() {
        contrast_obj <- selected_contrast()
        label <- if (is.null(contrast_obj) || nrow(contrast_obj) == 0) "contrast" else as.character(contrast_obj$contrast_id[[1]])
        paste0(label, ".clusterprofiler.shiny.pdf")
      },
      content = function(file) {
        split_list <- cp_split()
        if (length(split_list) == 0) {
          grDevices::pdf(file, width = 8, height = 6)
          graphics::plot.new()
          graphics::text(0.5, 0.5, "No clusterProfiler plot available.")
          grDevices::dev.off()
          return(invisible(NULL))
        }
        first_cat <- names(split_list)[[1]]
        p <- plot_cp_top_terms(split_list[[first_cat]], top_n = input$cp_top_n)
        if (is.null(p)) {
          grDevices::pdf(file, width = 8, height = 6)
          graphics::plot.new()
          graphics::text(0.5, 0.5, "No clusterProfiler plot available.")
          grDevices::dev.off()
          return(invisible(NULL))
        }
        ggplot2::ggsave(file, p + ggplot2::ggtitle(first_cat), width = 9, height = 6)
        invisible(NULL)
      }
    )

    session$onSessionEnded(function() {
      observe_cp_plots$destroy()
      observe_dc_plots$destroy()
    })
  })
}
