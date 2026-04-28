mod_project_loader_ui <- function(id) {
  ns <- shiny::NS(id)
  default_root <- Sys.getenv("OMNOMNOMICS_DE_APP_PROJECT", unset = "")
  if (!nzchar(default_root)) {
    default_root <- getwd()
  }
  shiny::tagList(
    shiny::h3("Project"),
    shiny::textInput(ns("de_root"), "Project directory or DE_calling directory", value = default_root),
    shiny::actionButton(ns("load"), "Load Project"),
    shiny::verbatimTextOutput(ns("status")),
    shiny::h4("Analyses"),
    shiny::tableOutput(ns("analysis_table"))
  )
}

mod_project_loader_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {
    project_index <- shiny::reactiveVal(NULL)
    status_text <- shiny::reactiveVal("No project loaded.")

    observe_load <- shiny::observeEvent(input$load, {
      selected_root <- trimws(as.character(input$de_root))
      if (!nzchar(selected_root)) {
        status_text("Please enter a DE_calling directory.")
        return(invisible(NULL))
      }

      loaded <- tryCatch(
        discover_de_calling(selected_root),
        error = function(e) e
      )
      if (inherits(loaded, "error")) {
        project_index(NULL)
        status_text(paste("Load failed:", conditionMessage(loaded)))
        return(invisible(NULL))
      }

      project_index(loaded)
      n_analyses <- length(loaded$analyses)
      n_contrasts <- sum(vapply(loaded$analyses, function(one_analysis) {
        nrow(one_analysis$contrast_index)
      }, integer(1)))
      status_text(
        paste0(
          "Loaded ", loaded$de_root,
          "\nAnalyses: ", n_analyses,
          "\nContrasts: ", n_contrasts
        )
      )
      invisible(NULL)
    })

    output$status <- shiny::renderText({
      status_text()
    })

    output$analysis_table <- shiny::renderTable({
      loaded <- project_index()
      shiny::req(loaded)
      if (length(loaded$analyses) == 0) {
        return(data.frame(
          analysis = character(0),
          contrasts = integer(0),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        analysis = vapply(loaded$analyses, function(one_analysis) one_analysis$analysis_id, character(1)),
        contrasts = vapply(loaded$analyses, function(one_analysis) nrow(one_analysis$contrast_index), integer(1)),
        stringsAsFactors = FALSE
      )
    })

    session$onSessionEnded(function() {
      observe_load$destroy()
    })

    auto_load_once <- shiny::observe({
      selected_root <- trimws(as.character(input$de_root))
      if (!nzchar(selected_root)) {
        return(invisible(NULL))
      }
      if (!is.null(project_index())) {
        return(invisible(NULL))
      }
      loaded <- tryCatch(
        discover_de_calling(selected_root),
        error = function(e) e
      )
      if (inherits(loaded, "error")) {
        status_text(paste("Auto-load skipped:", conditionMessage(loaded)))
        return(invisible(NULL))
      }
      project_index(loaded)
      n_analyses <- length(loaded$analyses)
      n_contrasts <- sum(vapply(loaded$analyses, function(one_analysis) {
        nrow(one_analysis$contrast_index)
      }, integer(1)))
      status_text(
        paste0(
          "Loaded ", loaded$de_root,
          "\nAnalyses: ", n_analyses,
          "\nContrasts: ", n_contrasts
        )
      )
      auto_load_once$destroy()
      invisible(NULL)
    })

    list(
      project_index = shiny::reactive(project_index()),
      metadata = shiny::reactive({
        loaded <- project_index()
        if (is.null(loaded)) {
          return(NULL)
        }
        loaded$metadata
      }),
      qc_index = shiny::reactive({
        loaded <- project_index()
        if (is.null(loaded)) {
          return(NULL)
        }
        loaded$qc_index
      })
    )
  })
}
