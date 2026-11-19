source("global.R", local = TRUE)

app_ui <- shiny::fluidPage(
  title = "omnomnomics Differential Explorer",
  shiny::tags$head(
    shiny::tags$link(rel = "stylesheet", type = "text/css", href = "explorer.css")
  ),
  shiny::tags$div(
    class = "explorer-shell",
    shiny::tags$header(
      class = "explorer-header",
      shiny::tags$div(
        shiny::tags$div(class = "explorer-brand", "omnomnomics"),
        shiny::tags$h1("Differential Explorer")
      )
    ),
    shiny::tags$div(
      class = "explorer-workspace",
      shiny::tags$aside(
        class = "explorer-sidebar",
        shiny::tags$section(
          class = "explorer-sidebar-section",
          mod_project_loader_ui("project_loader")
        ),
        shiny::tags$section(
          class = "explorer-sidebar-section",
          mod_analysis_selector_ui("analysis_selector")
        )
      ),
      shiny::tags$main(
        class = "explorer-main",
        shiny::uiOutput("contrast_heading", class = "explorer-contrast"),
        shiny::uiOutput("project_metrics", class = "explorer-metrics"),
        shiny::tags$section(
          class = "explorer-content",
          shiny::tabsetPanel(
            id = "explorer_view",
            shiny::tabPanel(
              "Overview",
              shiny::h3("Sample metadata"),
              shiny::tableOutput("metadata_preview"),
              shiny::tags$details(
                class = "explorer-details",
                shiny::tags$summary("Project and contrast details"),
                shiny::verbatimTextOutput("overview_text"),
                shiny::h4("Selected analysis"),
                shiny::tableOutput("selected_analysis_tbl"),
                shiny::h4("Selected contrast"),
                shiny::tableOutput("selected_contrast_tbl")
              ),
              shiny::tags$details(
                class = "explorer-details",
                shiny::tags$summary("QC files"),
                shiny::tableOutput("qc_preview")
              )
            ),
            shiny::tabPanel("QC", mod_qc_ui("qc")),
            shiny::tabPanel("DE Table", mod_de_table_ui("de_table")),
            shiny::tabPanel("Volcano", mod_volcano_ui("volcano")),
            shiny::tabPanel("Heatmap", mod_heatmap_ui("heatmap")),
            shiny::tabPanel("Boxplots", mod_boxplots_ui("boxplots")),
            shiny::tabPanel("Enrichment", mod_enrichment_ui("enrichment"))
          )
        )
      )
    )
  )
)

app_server <- function(input, output, session) {
  loader <- mod_project_loader_server("project_loader")
  selector <- mod_analysis_selector_server("analysis_selector", loader$project_index)
  de_table <- mod_de_table_server("de_table", selector$selected_contrast)
  volcano <- mod_volcano_server("volcano", selector$selected_contrast, de_table$filtered_table)
  heatmap <- mod_heatmap_server("heatmap", loader$project_index, selector$selected_contrast, de_table$filtered_table, loader$metadata)
  boxplots <- mod_boxplots_server("boxplots", loader$project_index, selector$selected_contrast, de_table$filtered_table, loader$metadata)
  qc <- mod_qc_server("qc", loader$project_index, loader$metadata)
  enrichment <- mod_enrichment_server("enrichment", selector$selected_contrast)

  output$contrast_heading <- shiny::renderUI({
    contrast_obj <- selector$selected_contrast()
    if (is.null(contrast_obj) || nrow(contrast_obj) == 0) {
      return(shiny::tags$h2("Select a project and contrast"))
    }
    shiny::tags$h2(as.character(contrast_obj$contrast_id[[1]]))
  })

  output$project_metrics <- shiny::renderUI({
    metadata_tbl <- loader$metadata()
    result_tbl <- de_table$de_table()
    filtered_tbl <- de_table$filtered_table()
    shiny::req(!is.null(metadata_tbl), !is.null(result_tbl), !is.null(filtered_tbl))
    sample_column <- if ("sample_id" %in% colnames(metadata_tbl)) "sample_id" else colnames(metadata_tbl)[[1]]
    sample_ids <- unique(as.character(metadata_tbl[[sample_column]]))
    sample_count <- sum(!is.na(sample_ids) & nzchar(sample_ids))
    shiny::tagList(
      shiny::tags$div(
        class = "explorer-metric",
        shiny::tags$span("Samples"),
        shiny::tags$strong(format(sample_count, big.mark = ",", trim = TRUE))
      ),
      shiny::tags$div(
        class = "explorer-metric",
        shiny::tags$span("Features in table"),
        shiny::tags$strong(format(nrow(result_tbl), big.mark = ",", trim = TRUE))
      ),
      shiny::tags$div(
        class = "explorer-metric explorer-metric-accent",
        shiny::tags$span("After table filters"),
        shiny::tags$strong(format(nrow(filtered_tbl), big.mark = ",", trim = TRUE))
      )
    )
  })

  output$overview_text <- shiny::renderText({
    loaded <- loader$project_index()
    if (is.null(loaded)) {
      return("Load a DE_calling directory to start exploring outputs.")
    }
    analysis_obj <- selector$selected_analysis()
    contrast_obj <- selector$selected_contrast()
    analysis_id <- if (is.null(analysis_obj)) "none" else analysis_obj$analysis_id
    contrast_id <- if (is.null(contrast_obj)) "none" else as.character(contrast_obj$contrast_id[[1]])
    paste0(
      "Project: ", loaded$de_root, "\n",
      "Metadata rows: ", nrow(loaded$metadata), "\n",
      "Shared QC files: ", nrow(loaded$qc_index), "\n",
      "Analyses detected: ", length(loaded$analyses), "\n",
      "Selected analysis: ", analysis_id, "\n",
      "Selected contrast: ", contrast_id
    )
  })

  output$selected_analysis_tbl <- shiny::renderTable({
    analysis_obj <- selector$selected_analysis()
    shiny::req(analysis_obj)
    data.frame(
      analysis_id = analysis_obj$analysis_id,
      analysis_dir = analysis_obj$analysis_dir,
      n_contrasts = nrow(analysis_obj$contrast_index),
      stringsAsFactors = FALSE
    )
  })

  output$selected_contrast_tbl <- shiny::renderTable({
    contrast_obj <- selector$selected_contrast()
    shiny::req(contrast_obj)
    out <- contrast_obj
    ctype <- if ("contrast_type" %in% colnames(out)) as.character(out$contrast_type[[1]]) else "factor"
    numerator <- if ("numerator" %in% colnames(out)) as.character(out$numerator[[1]]) else ""
    denominator <- if ("denominator" %in% colnames(out)) as.character(out$denominator[[1]]) else ""
    coef_name <- if ("coefficient_name" %in% colnames(out)) as.character(out$coefficient_name[[1]]) else ""
    if (identical(ctype, "factor") && nzchar(numerator) && nzchar(denominator)) {
      out$positive_direction <- paste0("higher in ", numerator, " vs ", denominator)
      out$negative_direction <- paste0("higher in ", denominator, " vs ", numerator)
    } else if (identical(ctype, "coefficient") && nzchar(coef_name)) {
      out$positive_direction <- paste0("positive coefficient direction: ", coef_name)
      out$negative_direction <- "opposite coefficient direction"
    } else {
      out$positive_direction <- "modeled numerator direction"
      out$negative_direction <- "opposite direction"
    }
    out
  })

  output$metadata_preview <- shiny::renderTable({
    metadata_tbl <- loader$metadata()
    shiny::req(metadata_tbl)
    utils::head(metadata_tbl, 12)
  })

  output$qc_preview <- shiny::renderTable({
    qc_tbl <- loader$qc_index()
    shiny::req(qc_tbl)
    utils::head(qc_tbl, 20)
  })
}

shiny::shinyApp(ui = app_ui, server = app_server)
