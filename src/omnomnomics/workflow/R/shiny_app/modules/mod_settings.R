mod_settings_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Plot settings"),
    shiny::p("Group colours are shared by PCA, boxplots and heatmap annotations."),
    shiny::fluidRow(
      shiny::column(6, shiny::selectInput(ns("group_by"), "Sample grouping", choices = character(0))),
      shiny::column(6, shiny::selectInput(ns("preset"), "Palette preset", choices = c(
        "Teal & coral" = "coastal", "Okabe–Ito" = "okabe_ito", "Muted" = "muted", "Custom gradient" = "gradient"
      )))
    ),
    shiny::conditionalPanel(
      condition = "input.preset === 'gradient'", ns = ns,
      shiny::tags$div(
        class = "explorer-gradient-controls",
        explorer_color_input(ns("gradient_start"), "Gradient start", "#175A81"),
        explorer_color_input(ns("gradient_end"), "Gradient end", "#DA845C")
      ),
      shiny::textOutput(ns("generator"))
    ),
    shiny::h4("Colours per group"),
    shiny::uiOutput(ns("group_colors"), class = "explorer-color-grid"),
    shiny::actionButton(ns("reset_palette"), "Reset to preset"),
    shiny::h4("PCA points"),
    shiny::numericInput(ns("pca_size"), "Point size", value = 6, min = 2, max = 10, step = 0.5)
  )
}

mod_settings_server <- function(id, metadata) {
  shiny::moduleServer(id, function(input, output, session) {
    observe_group_choices <- shiny::observeEvent(metadata(), {
      md <- metadata()
      shiny::req(md)
      choices <- colnames(md)[vapply(md, function(values) {
        any(!is.na(values) & nzchar(as.character(values)))
      }, logical(1))]
      preferred <- c("condition", "sample_type", "sample_color", choices)
      selected <- preferred[preferred %in% choices][1]
      if (!length(selected) || is.na(selected)) selected <- character(0)
      shiny::updateSelectInput(session, "group_by", choices = choices, selected = selected)
    }, ignoreNULL = TRUE)

    group_definition <- shiny::reactive({
      md <- metadata()
      column <- input$group_by
      if (is.null(md) || !length(column) || !(column %in% colnames(md))) {
        return(list(column = "", groups = character(0)))
      }
      values <- as.character(md[[column]])
      list(column = column, groups = sort(unique(values[!is.na(values) & nzchar(values)])))
    })

    palette_for <- function(n) {
      preset <- input$preset
      if (!length(preset)) preset <- "coastal"
      start <- input$gradient_start
      end <- input$gradient_end
      if (!length(start) || !grepl("^#[0-9A-Fa-f]{6}$", start)) start <- "#175A81"
      if (!length(end) || !grepl("^#[0-9A-Fa-f]{6}$", end)) end <- "#DA845C"
      explorer_palette(n, preset, start, end)
    }

    colors_for <- function(values, column) {
      md <- metadata()
      groups <- as.character(values)
      if (!is.null(md) && column %in% colnames(md)) {
        groups <- as.character(md[[column]])
      }
      groups <- sort(unique(groups[!is.na(groups) & nzchar(groups)]))
      colors <- stats::setNames(palette_for(length(groups)), groups)
      for (group in groups) {
        picked <- input[[explorer_color_id(column, group)]]
        if (length(picked) && grepl("^#[0-9A-Fa-f]{6}$", picked)) colors[[group]] <- picked
      }
      colors
    }

    output$group_colors <- shiny::renderUI({
      definition <- group_definition()
      if (!length(definition$groups)) return(shiny::tags$em("Load a project to choose group colours."))
      colors <- shiny::isolate(colors_for(definition$groups, definition$column))
      shiny::tagList(lapply(definition$groups, function(group) {
        explorer_color_input(session$ns(explorer_color_id(definition$column, group)), group, colors[[group]])
      }))
    })

    observe_preset <- shiny::observeEvent(list(input$preset, input$gradient_start, input$gradient_end, input$reset_palette), {
      definition <- group_definition()
      colors <- palette_for(length(definition$groups))
      for (index in seq_along(definition$groups)) {
        session$sendInputMessage(explorer_color_id(definition$column, definition$groups[[index]]), list(value = colors[[index]]))
      }
    }, ignoreInit = TRUE)

    output$generator <- shiny::renderText({
      if (requireNamespace("coloRspacer", quietly = TRUE)) {
        paste("Gradient generator: coloRspacer", utils::packageVersion("coloRspacer"))
      } else {
        "Gradient generator: R colorRampPalette. Uses coloRspacer when installed."
      }
    })

    annotation_colors <- function(annotation) {
      if (is.null(annotation) || !ncol(annotation)) return(NA)
      categorical <- names(annotation)[!vapply(annotation, is.numeric, logical(1))]
      colors <- lapply(categorical, function(column) colors_for(annotation[[column]], column))
      stats::setNames(colors, categorical)
    }

    session$onSessionEnded(function() {
      observe_group_choices$destroy()
      observe_preset$destroy()
    })

    list(
      group_by = shiny::reactive(input$group_by),
      colors_for = colors_for,
      annotation_colors = annotation_colors,
      pca_size = shiny::reactive({
        value <- suppressWarnings(as.numeric(input$pca_size))
        if (length(value) != 1 || !is.finite(value)) 6 else max(2, min(10, value))
      })
    )
  })
}
