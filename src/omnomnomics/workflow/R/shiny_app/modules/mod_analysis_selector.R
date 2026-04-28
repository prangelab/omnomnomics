mod_analysis_selector_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Selection"),
    shiny::selectInput(ns("analysis_id"), "Analysis", choices = character(0)),
    shiny::selectInput(ns("contrast_id"), "Contrast", choices = character(0)),
    shiny::verbatimTextOutput(ns("selection_status")),
    shiny::h4("Direction Legend"),
    shiny::verbatimTextOutput(ns("direction_legend"))
  )
}

mod_analysis_selector_server <- function(id, project_index) {
  shiny::moduleServer(id, function(input, output, session) {
    observe_analysis_choices <- shiny::observe({
      loaded <- project_index()
      if (is.null(loaded) || length(loaded$analyses) == 0) {
        shiny::updateSelectInput(session, "analysis_id", choices = character(0), selected = character(0))
        shiny::updateSelectInput(session, "contrast_id", choices = character(0), selected = character(0))
        return(invisible(NULL))
      }

      analysis_ids <- names(loaded$analyses)
      selected_analysis <- input$analysis_id
      if (!nzchar(selected_analysis) || !(selected_analysis %in% analysis_ids)) {
        selected_analysis <- analysis_ids[[1]]
      }
      shiny::updateSelectInput(session, "analysis_id", choices = analysis_ids, selected = selected_analysis)
      invisible(NULL)
    })

    observe_contrast_choices <- shiny::observe({
      loaded <- project_index()
      selected_analysis <- input$analysis_id
      if (is.null(loaded) || !nzchar(selected_analysis) || !(selected_analysis %in% names(loaded$analyses))) {
        shiny::updateSelectInput(session, "contrast_id", choices = character(0), selected = character(0))
        return(invisible(NULL))
      }

      contrast_tbl <- loaded$analyses[[selected_analysis]]$contrast_index
      contrast_ids <- if (nrow(contrast_tbl) > 0) contrast_tbl$contrast_id else character(0)
      contrast_ids <- as.character(contrast_ids)

      selected_contrast <- input$contrast_id
      if (!nzchar(selected_contrast) || !(selected_contrast %in% contrast_ids)) {
        selected_contrast <- if (length(contrast_ids) > 0) contrast_ids[[1]] else character(0)
      }
      shiny::updateSelectInput(session, "contrast_id", choices = contrast_ids, selected = selected_contrast)
      invisible(NULL)
    })

    output$selection_status <- shiny::renderText({
      loaded <- project_index()
      if (is.null(loaded) || length(loaded$analyses) == 0) {
        return("Load a project to enable analysis and contrast selection.")
      }

      analysis_id <- input$analysis_id
      if (!nzchar(analysis_id) || !(analysis_id %in% names(loaded$analyses))) {
        return("Select an analysis.")
      }

      contrast_tbl <- loaded$analyses[[analysis_id]]$contrast_index
      contrast_id <- input$contrast_id
      if (nrow(contrast_tbl) == 0) {
        return(paste0("Analysis selected: ", analysis_id, "\nNo contrasts detected."))
      }
      if (!nzchar(contrast_id) || !(contrast_id %in% contrast_tbl$contrast_id)) {
        return(paste0("Analysis selected: ", analysis_id, "\nSelect a contrast."))
      }
      paste0(
        "Analysis: ", analysis_id, "\n",
        "Contrast: ", contrast_id
      )
    })

    output$direction_legend <- shiny::renderText({
      contrast_obj <- selected_contrast()
      if (is.null(contrast_obj) || nrow(contrast_obj) == 0) {
        return("Select a contrast to show log2FC direction interpretation.")
      }
      ctype <- if ("contrast_type" %in% colnames(contrast_obj)) as.character(contrast_obj$contrast_type[[1]]) else "factor"
      numerator <- if ("numerator" %in% colnames(contrast_obj)) as.character(contrast_obj$numerator[[1]]) else ""
      denominator <- if ("denominator" %in% colnames(contrast_obj)) as.character(contrast_obj$denominator[[1]]) else ""
      coef_name <- if ("coefficient_name" %in% colnames(contrast_obj)) as.character(contrast_obj$coefficient_name[[1]]) else ""
      if (identical(ctype, "factor") && nzchar(numerator) && nzchar(denominator)) {
        return(
          paste0(
            "Positive log2FC/score: higher in ", numerator, " vs ", denominator, ".\n",
            "Negative log2FC/score: higher in ", denominator, " vs ", numerator, "."
          )
        )
      }
      if (identical(ctype, "coefficient") && nzchar(coef_name)) {
        return(
          paste0(
            "Coefficient contrast: ", coef_name, ".\n",
            "Positive values follow positive coefficient direction; negative values follow the opposite direction."
          )
        )
      }
      "Positive values indicate enrichment in the modeled numerator direction; negative values indicate the opposite direction."
    })

    session$onSessionEnded(function() {
      observe_analysis_choices$destroy()
      observe_contrast_choices$destroy()
    })

    selected_analysis <- shiny::reactive({
      loaded <- project_index()
      if (is.null(loaded) || length(loaded$analyses) == 0) {
        return(NULL)
      }
      analysis_id <- input$analysis_id
      if (!nzchar(analysis_id) || !(analysis_id %in% names(loaded$analyses))) {
        return(NULL)
      }
      loaded$analyses[[analysis_id]]
    })

    selected_contrast <- shiny::reactive({
      analysis_obj <- selected_analysis()
      if (is.null(analysis_obj)) {
        return(NULL)
      }
      contrast_tbl <- analysis_obj$contrast_index
      if (nrow(contrast_tbl) == 0) {
        return(NULL)
      }
      contrast_id <- input$contrast_id
      if (!nzchar(contrast_id) || !(contrast_id %in% contrast_tbl$contrast_id)) {
        return(NULL)
      }
      contrast_tbl[contrast_tbl$contrast_id == contrast_id, , drop = FALSE]
    })

    list(
      selected_analysis = selected_analysis,
      selected_contrast = selected_contrast
    )
  })
}
