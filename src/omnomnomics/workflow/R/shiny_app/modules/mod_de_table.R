mod_de_table_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::h3("DE Table"),
    shiny::fluidRow(
      shiny::column(3, shiny::numericInput(ns("padj_max"), "padj <=", value = 0.05, min = 0, max = 1, step = 0.001)),
      shiny::column(3, shiny::numericInput(ns("lfc_min"), "|log2FC| >=", value = 1, min = 0, step = 0.1)),
      shiny::column(3, shiny::numericInput(ns("base_mean_min"), "baseMean >=", value = 0, min = 0, step = 1)),
      shiny::column(3, shiny::textInput(ns("gene_query"), "Gene search", value = ""))
    ),
    shiny::fluidRow(
      shiny::column(6, shiny::downloadButton(ns("download_filtered"), "Download Filtered TSV")),
      shiny::column(6, shiny::verbatimTextOutput(ns("table_status")))
    ),
    shiny::h4("Contrast Summary"),
    shiny::tableOutput(ns("summary_tbl")),
    shiny::h4("Filtered Table (Preview)"),
    shiny::tableOutput(ns("filtered_tbl"))
  )
}

mod_de_table_server <- function(id, selected_contrast) {
  shiny::moduleServer(id, function(input, output, session) {
    current_de_table <- shiny::reactive({
      contrast_obj <- selected_contrast()
      if (is.null(contrast_obj) || nrow(contrast_obj) == 0) {
        return(NULL)
      }
      de_path <- as.character(contrast_obj$de_table[[1]])
      if (!nzchar(de_path) || is.na(de_path) || !file.exists(de_path)) {
        return(NULL)
      }
      read_de_table(de_path)
    })

    filtered_table <- shiny::reactive({
      tbl <- current_de_table()
      if (is.null(tbl)) {
        return(NULL)
      }
      filter_de_table(
        tbl = tbl,
        padj_max = input$padj_max,
        lfc_min = input$lfc_min,
        base_mean_min = input$base_mean_min,
        gene_query = input$gene_query
      )
    })

    output$table_status <- shiny::renderText({
      tbl <- current_de_table()
      filtered <- filtered_table()
      if (is.null(tbl)) {
        return("No DE table available for the selected contrast.")
      }
      paste0("Rows: ", nrow(tbl), " total, ", nrow(filtered), " after filters.")
    })

    output$summary_tbl <- shiny::renderTable({
      tbl <- filtered_table()
      shiny::req(tbl)
      summarize_contrast(tbl, alpha = input$padj_max, lfc_cutoff = input$lfc_min)
    })

    output$filtered_tbl <- shiny::renderTable({
      tbl <- filtered_table()
      shiny::req(tbl)
      utils::head(tbl, 50)
    })

    output$download_filtered <- shiny::downloadHandler(
      filename = function() {
        contrast_obj <- selected_contrast()
        label <- if (is.null(contrast_obj) || nrow(contrast_obj) == 0) "contrast" else as.character(contrast_obj$contrast_id[[1]])
        paste0(label, ".filtered.tsv")
      },
      content = function(file) {
        tbl <- filtered_table()
        if (is.null(tbl)) {
          utils::write.table(
            data.frame(message = "No filtered DE table available.", stringsAsFactors = FALSE),
            file = file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
          )
          return(invisible(NULL))
        }
        utils::write.table(tbl, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
        invisible(NULL)
      }
    )

    list(
      de_table = current_de_table,
      filtered_table = filtered_table
    )
  })
}
