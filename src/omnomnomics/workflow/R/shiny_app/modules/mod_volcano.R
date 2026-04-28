mod_volcano_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("Volcano"),
    shiny::fluidRow(
      shiny::column(3, shiny::numericInput(ns("alpha"), "alpha", value = 0.05, min = 0, max = 1, step = 0.001)),
      shiny::column(3, shiny::numericInput(ns("lfc_cutoff"), "|log2FC| cutoff", value = 1, min = 0, step = 0.1)),
      shiny::column(6, shiny::downloadButton(ns("download_plot"), "Download Volcano PDF"))
    ),
    shiny::verbatimTextOutput(ns("status")),
    shiny::plotOutput(ns("volcano_plot"), height = "560px")
  )
}

mod_volcano_server <- function(id, selected_contrast, filtered_table) {
  shiny::moduleServer(id, function(input, output, session) {
    volcano_plot_obj <- shiny::reactive({
      contrast_obj <- selected_contrast()
      if (is.null(contrast_obj) || nrow(contrast_obj) == 0) {
        return(NULL)
      }

      contrast_id <- as.character(contrast_obj$contrast_id[[1]])
      contrast_dir <- as.character(contrast_obj$contrast_dir[[1]])
      existing_pdf <- file.path(contrast_dir, paste0(contrast_id, ".volcano_plot_labels.pdf"))
      if (!file.exists(existing_pdf)) {
        existing_pdf <- file.path(contrast_dir, paste0(contrast_id, ".volcano_plot.pdf"))
      }

      tbl <- filtered_table()
      if (is.null(tbl) || !nrow(tbl)) {
        return(list(plot = NULL, source = if (file.exists(existing_pdf)) existing_pdf else NA_character_))
      }
      vdf <- build_volcano_data(tbl, alpha = input$alpha, lfc_cutoff = input$lfc_cutoff)
      p <- plot_volcano(vdf, alpha = input$alpha, lfc_cutoff = input$lfc_cutoff)
      list(plot = p, source = if (file.exists(existing_pdf)) existing_pdf else NA_character_)
    })

    output$status <- shiny::renderText({
      vp <- volcano_plot_obj()
      if (is.null(vp)) {
        return("Select a contrast to view volcano outputs.")
      }
      if (!is.na(vp$source) && nzchar(vp$source)) {
        paste0("Rendered from table filters. Existing exported PDF: ", vp$source)
      } else {
        "Rendered from table filters. No existing volcano PDF detected for this contrast."
      }
    })

    output$volcano_plot <- shiny::renderPlot({
      vp <- volcano_plot_obj()
      shiny::req(vp)
      shiny::req(vp$plot)
      vp$plot
    })

    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        contrast_obj <- selected_contrast()
        label <- if (is.null(contrast_obj) || nrow(contrast_obj) == 0) "contrast" else as.character(contrast_obj$contrast_id[[1]])
        paste0(label, ".volcano.shiny.pdf")
      },
      content = function(file) {
        vp <- volcano_plot_obj()
        if (is.null(vp) || is.null(vp$plot)) {
          grDevices::pdf(file, width = 8, height = 6)
          graphics::plot.new()
          graphics::text(0.5, 0.5, "No volcano plot available.")
          grDevices::dev.off()
          return(invisible(NULL))
        }
        ggplot2::ggsave(file, vp$plot, width = 8, height = 6)
        invisible(NULL)
      }
    )
  })
}
