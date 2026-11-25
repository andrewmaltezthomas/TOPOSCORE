# TOPOSCORE Calculator: Shiny app
library(shiny)
library(bslib)
library(DT)
library(ggplot2)

source("R/toposcore_calc.R")
source("R/plots.R")
options(shiny.maxRequestSize = 100 * 1024^2)

APP_VERSION <- "2.0.0"
XW <- load_toposcore_crosswalk("data/toposcore_crosswalk.tsv")
REF_SCORES <- read.csv("data/scores_disc.csv", check.names = FALSE)
REF_CLIN <- read.csv("data/clin_disc.csv", check.names = FALSE)
DB_CHOICES <- setNames(toposcore_databases$id, toposcore_databases$label)
FORMAT_LABEL <- c(metaphlan = "MetaPhlAn (k__|p__|...)", gtdb = "GTDB-style (d__;p__;...)", meteor = "Meteor2 (msp)")
db_family <- function(fmt) sub("_(species|sgb)$", "", fmt)          # metaphlan_sgb -> metaphlan

theme_app <- bs_theme(version = 5, primary = BRAND[["coral"]], secondary = BRAND[["grey"]],
                      success = "#009E73", danger = "#C0392B",
                      base_font = font_collection("Helvetica Neue", "Arial", "sans-serif"))

logo <- function(src, alt, href, height = "40px")
  tags$a(href = href, target = "_blank", rel = "noopener", class = "me-3",
         tags$img(src = src, alt = alt, style = paste0("height:", height, ";")))

# ---- UI ------------------------------------------------------------------------------------------
ui <- page_navbar(
  title = tags$span(tags$strong("TOPOSCORE"), " Calculator", class = "fs-4"),
  theme = theme_app, window_title = "TOPOSCORE Calculator", fillable = FALSE, bg = "#FFFFFF", inverse = FALSE,
  header = tags$head(tags$style(HTML("
    .navbar { border-bottom: 3px solid #F24E57; }
    .navbar .nav-link.active { color: #F24E57 !important; font-weight: 600; }
    .bslib-value-box .value-box-value { font-size: 1.9rem; }
    .value-box-title { font-size: .95rem; }
    .dt-taxa { max-width: 340px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }
    .not-assessable li { font-style: italic; }"))),

  nav_panel("Calculator", icon = icon("calculator"),
    layout_sidebar(fillable = FALSE,
      sidebar = sidebar(width = 360,
        fileInput("file", "1. Upload a merged profile table (tab-separated)",
                  accept = c(".tsv", ".txt", ".tab", "text/tab-separated-values")),
        selectInput("database", "2. Profiler and database", choices = DB_CHOICES, selected = "SGB_Jan21"),
        uiOutput("format_hint"),
        actionButton("calculate", "3. Calculate TOPOSCORE", class = "btn-primary w-100", icon = icon("play")),
        hr(),
        accordion(open = FALSE,
          accordion_panel("Input format", icon = icon("circle-info"),
            HTML("<p>Taxa in rows, samples in columns, relative abundances.</p>
                  <ul><li><b>MetaPhlAn 4</b>: output of <code>merge_metaphlan_tables.py</code>
                  (column <code>clade_name</code>). SGB options need <code>t__SGB</code> rows.</li>
                  <li><b>GTDB r207/r220, UHGG v2.0.2</b>: column <code>clade_name</code> with
                  <code>d__;p__;c__;o__;f__;g__;s__</code> lineages (e.g. sylph, mOTUs4, Kraken2/Bracken).
                  Fractions or percentages are both accepted.</li>
                  <li><b>Meteor2</b>: columns <code>msp</code> and <code>msp_taxonomy</code>.</li></ul>
                  <p>For Kraken2/Bracken, correct read fractions for genome length first: the
                  <i>Akkermansia</i> cut-off is defined on taxonomic (cell) abundance.</p>")),
          accordion_panel("Example files", icon = icon("download"),
            downloadButton("dl_example_mpa", "MetaPhlAn 4 vJan21", class = "btn-sm btn-outline-secondary w-100 mb-2"),
            downloadButton("dl_example_gtdb", "sylph GTDB r220", class = "btn-sm btn-outline-secondary w-100"))
        )
      ),
      uiOutput("before_run"),
      uiOutput("summary_boxes"),
      uiOutput("run_details"),
      conditionalPanel("output.has_results",
        navset_card_tab(id = "result_tabs", full_screen = TRUE,
          nav_panel("Results table",
            div(class = "d-flex justify-content-end mb-2",
                downloadButton("dl_results", "Download results (CSV)", class = "btn-sm btn-outline-secondary")),
            DTOutput("results_table", height = "auto", fill = FALSE)),
          nav_panel("Score distribution",
            layout_columns(col_widths = c(8, 4),
              selectizeInput("highlight", "Highlight samples (or select rows in the table)", choices = NULL,
                             multiple = TRUE, options = list(placeholder = "Type sample IDs", plugins = list("remove_button"))),
              div(class = "pt-4", actionLink("clear_highlight", "Clear highlighted samples"))),
            plotOutput("score_plot", height = "460px"),
            p(class = "text-muted small", "Reference: discovery cohort of Derosa et al. (2024). Vertical lines: S-score cut-offs ",
              TOPOSCORE_CUTOFFS[["SIG1"]], " and ", TOPOSCORE_CUTOFFS[["SIG2"]], ". Ticks at the bottom: your samples."),
            div(downloadButton("dl_score_plot", "Download plot (PDF)", class = "btn-sm btn-outline-secondary"))),
          nav_panel("Signature species",
            checkboxInput("heatmap_highlight_only", "Show only highlighted samples", FALSE),
            plotOutput("heatmap", height = "auto"),
            div(downloadButton("dl_heatmap", "Download heatmap (PDF)", class = "btn-sm btn-outline-secondary")))
        )
      )
    )
  ),

  nav_panel("Signature species", icon = icon("list"),
    card(card_header("SIG1 and SIG2 species and their names in each profiler"),
         p("TOPOSCORE uses 37 SIG1 (resistance-associated) and 45 SIG2 (response-associated) species. ",
           "The crosswalk lists the identifier(s) each profiler and database uses for every species; a species ",
           "counts as present if any of its identifiers has a non-zero abundance."),
         DTOutput("crosswalk_table", height = "auto", fill = FALSE),
         div(class = "mt-2", downloadButton("dl_crosswalk", "Download crosswalk (TSV)", class = "btn-sm")))
  ),

  nav_panel("About", icon = icon("book-open"),
    layout_columns(col_widths = c(7, 5),
      card(card_header("About TOPOSCORE"),
        HTML(sprintf("
        <p>TOPOSCORE classifies gut microbiome profiles of cancer patients by the detection of 37 SIG1 and
        45 SIG2 species (Derosa et al., <i>Cell</i> 2024).</p>
        <ol>
          <li><b>S-score</b> = 0.5 &times; (SIG2 detected / 45 &minus; SIG1 detected / 37 + 1)</li>
          <li><b>SIG class</b>: SIG1 if S &le; %.4f; SIG2 if S &ge; %.4f; otherwise grey zone.</li>
          <li><b><i>Akkermansia muciniphila</i></b>: Zero, Normal, or High (&ge; 4.799%%; &ge; 4.52%% for Meteor2).</li>
          <li><b>TOPOSCORE</b>: SIG1 &rarr; SIG1+; SIG2 &rarr; SIG2+; grey zone &rarr; SIG1+ if
          <i>Akkermansia</i> is Zero or High, SIG2+ if Normal.</li>
        </ol>
        <p>Profiles from other profilers are mapped to the original MetaPhlAn 4 (vJan21) species through a
        genome-based crosswalk (GTDB r207/r220, UHGG v2.0.2, Meteor2 MSPs). Species a database cannot represent
        are listed after each run.</p>
        <p><b>Citation</b>: Derosa L, Iebba V, Silva CAC, et al. Custom scoring based on ecological topology of gut
        microbiota associated with cancer immunotherapy outcome. <i>Cell</i> 2024.
        <a href='https://doi.org/10.1016/j.cell.2024.05.029' target='_blank'>doi:10.1016/j.cell.2024.05.029</a></p>",
        TOPOSCORE_CUTOFFS[["SIG1"]], TOPOSCORE_CUTOFFS[["SIG2"]]))),
      card(card_header("Data and privacy"),
        p("Uploaded files are processed in memory for the duration of your session and are not stored."),
        p("This tool is for research use only and is not a medical device."),
        hr(),
        p(class = "text-muted small", "App version ", APP_VERSION, ". Developed with the Seerave Foundation and Gustave Roussy."),
        div(class = "d-flex align-items-center gap-3 mt-2",
            tags$img(src = "seerave.png", alt = "Seerave Foundation", style = "height:44px;"),
            tags$img(src = "gustave_roussy.png", alt = "Gustave Roussy", style = "height:56px;")))
    )
  ),

  nav_spacer(),
  nav_item(logo("seerave.png", "Seerave Foundation", "https://www.seerave.org", "34px")),
  nav_item(logo("gustave_roussy.png", "Gustave Roussy", "https://www.gustaveroussy.fr", "46px"))
)

# ---- server --------------------------------------------------------------------------------------
server <- function(input, output, session) {
  results <- reactiveVal(NULL)
  output$has_results <- reactive(!is.null(results()))
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)

  detected_format <- reactive({ req(input$file); detect_profile_format(input$file$datapath) })
  selected_db <- reactive(toposcore_databases[toposcore_databases$id == input$database, ])

  output$format_hint <- renderUI({
    req(input$file)
    f <- detected_format()
    if (is.na(f)) return(div(class = "alert alert-warning py-2 small", "Could not recognise the table format."))
    ok <- db_family(selected_db()$format) == f
    div(class = paste("alert py-2 small", if (ok) "alert-success" else "alert-danger"),
        "File looks like ", tags$b(FORMAT_LABEL[[f]]), ".",
        if (!ok) " This does not match the selected profiler: choose a matching option.")
  })

  observeEvent(input$calculate, {
    if (is.null(input$file)) {
      showNotification("Please upload a profile table first.", type = "warning"); return()
    }
    f <- detected_format()
    if (!is.na(f) && db_family(selected_db()$format) != f) {
      showNotification(sprintf("The file looks like %s, which does not match '%s'.", FORMAT_LABEL[[f]],
                               selected_db()$label), type = "error", duration = 10)
      return()
    }
    res <- withProgress(message = "Calculating TOPOSCORE", value = 0.5, tryCatch(
      calculate_toposcore(input$file$datapath, input$database, crosswalk = XW),
      error = function(e) { showNotification(paste("Could not calculate:", conditionMessage(e)), type = "error", duration = NULL); NULL }))
    results(res)
    if (!is.null(res)) {
      updateSelectizeInput(session, "highlight", choices = res$Sample_id, selected = character(0), server = TRUE)
      showNotification(sprintf("TOPOSCORE calculated for %d samples.", nrow(res)), type = "message")
    }
  })

  # Changing the database or file invalidates previous results
  observeEvent(list(input$database, input$file), results(NULL), ignoreInit = TRUE)

  output$before_run <- renderUI({
    if (!is.null(results())) return(NULL)
    card(card_body(
      h4("Calculate TOPOSCORE from a taxonomic profile"),
      p("Upload a merged profile table, choose the profiler and database it was produced with, and press ",
        tags$b("Calculate"), ". Results, score distribution and a heatmap of detected signature species will appear here."),
      p(class = "text-muted small mb-0", "Supported: MetaPhlAn 4 (vJan21, vJun23, vJan25), GTDB r207/r220 profilers (sylph, mOTUs4, ",
        "Kraken2/Bracken), Kraken2/Bracken with UHGG v2.0.2, and Meteor2.")))
  })

  output$summary_boxes <- renderUI({
    res <- req(results())
    n <- function(col, v) sum(res[[col]] == v)
    layout_columns(fill = FALSE,
      value_box("Samples", nrow(res), theme = value_box_theme(bg = "#FFFFFF", fg = BRAND[["grey"]])),
      value_box("TOPOSCORE SIG1+", n("Toposcore", "SIG1+"), p(sprintf("%.0f%% of samples", 100 * mean(res$Toposcore == "SIG1+"))), theme = value_box_theme(bg = TOPO_COLOURS[["SIG1"]], fg = "#FFFFFF")),
      value_box("TOPOSCORE SIG2+", n("Toposcore", "SIG2+"), p(sprintf("%.0f%% of samples", 100 * mean(res$Toposcore == "SIG2+"))), theme = value_box_theme(bg = TOPO_COLOURS[["SIG2"]], fg = "#FFFFFF")),
      value_box("SIG class", sprintf("%d / %d / %d", n("SIG_class", "SIG1"), n("SIG_class", "Gray"), n("SIG_class", "SIG2")),
                p("SIG1 / grey zone / SIG2"),
                theme = value_box_theme(bg = "#FFFFFF", fg = BRAND[["grey"]]))
    )
  })

  output$run_details <- renderUI({
    res <- req(results()); d <- attr(res, "diagnostics")
    miss <- c(d$sig1_not_in_table, d$sig2_not_in_table)
    lst <- function(x) if (length(x)) tags$ul(class = "not-assessable small mb-1", lapply(gsub("_", " ", x), tags$li)) else span(class = "small", "none")
    title <- if (length(miss)) sprintf("Run details: %d SIG1 and %d SIG2 species not found in this table",
                                       length(d$sig1_not_in_table), length(d$sig2_not_in_table)) else
      "Run details: all signature species found"
    accordion(class = "mb-3", open = FALSE,
      accordion_panel(title, value = "details", icon = icon(if (length(miss)) "triangle-exclamation" else "circle-check"),
        layout_columns(col_widths = c(5, 7),
          div(tags$dl(class = "small mb-0",
            tags$dt("Profiler / database"), tags$dd(d$database),
            tags$dt("Akkermansia"), tags$dd(HTML(sprintf("%s; High at &ge; %.3g%%", d$akk_source, d$akk_cutoff))),
            tags$dt("Abundance units"), tags$dd(if (isTRUE(d$scaled_to_percent)) "fractions, converted to %" else "percent"),
            tags$dt("Assessable species"), tags$dd(sprintf("SIG1 %d/37, SIG2 %d/45",
                                                     37 - length(d$sig1_not_in_table), 45 - length(d$sig2_not_in_table))))),
          div(p(class = "small mb-1", tags$b("Signature species not found in this table"),
                " (absent from every sample, or not represented in the database). They count as not detected,",
                " so S-scores can be biased when many are missing."),
              layout_columns(div(tags$b(class = "small", "SIG1"), lst(d$sig1_not_in_table)),
                             div(tags$b(class = "small", "SIG2"), lst(d$sig2_not_in_table))))))
    )
  })

  # Results table: compact taxa columns with full text on hover
  output$results_table <- renderDT({
    res <- req(results())
    show <- res[, c("Sample_id", "S_score", "SIG_class", "Toposcore", "Akk_status", "Akk_abundance",
                    "SIG1_count", "SIG2_count", "SIG1_taxa", "SIG2_taxa")]
    show$Akk_abundance <- signif(show$Akk_abundance, 3)
    show$SIG1_taxa <- gsub("_", " ", gsub(";", "; ", show$SIG1_taxa))
    show$SIG2_taxa <- gsub("_", " ", gsub(";", "; ", show$SIG2_taxa))
    datatable(show, rownames = FALSE, selection = "multiple", fillContainer = FALSE,
              colnames = c("Sample", "S-score", "SIG class", "TOPOSCORE", "Akkermansia", "Akkermansia (%)",
                           "SIG1 detected", "SIG2 detected", "SIG1 species", "SIG2 species"),
              options = list(pageLength = 15, scrollX = TRUE, order = list(list(1, "asc")),
                             columnDefs = list(list(targets = 8:9, render = JS(
                               "function(d, t) { return t === 'display' ? '<div class=\"dt-taxa\" title=\"' + d + '\">' + d + '</div>' : d; }"))))) |>
      formatStyle("Toposcore", color = "white", fontWeight = "bold",
                  backgroundColor = styleEqual(c("SIG1+", "SIG2+"), unname(TOPO_COLOURS[c("SIG1+", "SIG2+")]))) |>
      formatStyle("SIG_class", color = styleEqual(c("SIG1", "Gray", "SIG2"), unname(TOPO_COLOURS[c("SIG1", "Gray", "SIG2")])),
                  fontWeight = "bold")
  })

  # Table row selection adds samples to the highlight list
  observeEvent(input$results_table_rows_selected, {
    res <- req(results())
    ids <- union(input$highlight, res$Sample_id[input$results_table_rows_selected])
    updateSelectizeInput(session, "highlight", choices = res$Sample_id, selected = ids, server = TRUE)
  })
  observeEvent(input$clear_highlight, {
    res <- req(results())
    updateSelectizeInput(session, "highlight", choices = res$Sample_id, selected = character(0), server = TRUE)
    selectRows(dataTableProxy("results_table"), NULL)
  })

  score_plot <- reactive(plot_reference_density(REF_SCORES, REF_CLIN, results(), input$highlight))
  output$score_plot <- renderPlot(score_plot(), res = 96)

  heatmap_samples <- reactive(if (isTRUE(input$heatmap_highlight_only)) input$highlight else NULL)
  heatmap_plot <- reactive({
    res <- req(results()); d <- attr(res, "diagnostics")
    plot_sig_heatmap(res, attr(res, "presence"), XW, c(d$sig1_not_in_table, d$sig2_not_in_table), heatmap_samples())
  })
  heatmap_height <- function() 220 + 13 * (nrow(XW) + 1)
  output$heatmap <- renderPlot(heatmap_plot(), height = heatmap_height, res = 96)

  # Downloads
  stamp <- function() format(Sys.time(), "%Y%m%d_%H%M")
  output$dl_results <- downloadHandler(
    filename = function() sprintf("toposcore_%s_%s.csv", input$database, stamp()),
    content = function(file) {
      res <- results(); d <- attr(res, "diagnostics")
      out <- as.data.frame(res); attr(out, "presence") <- NULL; attr(out, "diagnostics") <- NULL
      out$database <- d$database; out$app_version <- APP_VERSION
      write.csv(out, file, row.names = FALSE)
    })
  output$dl_score_plot <- downloadHandler(
    filename = function() sprintf("toposcore_distribution_%s.pdf", stamp()),
    content = function(file) ggsave(file, score_plot(), width = 10, height = 5.5))
  output$dl_heatmap <- downloadHandler(
    filename = function() sprintf("toposcore_signature_heatmap_%s.pdf", stamp()),
    content = function(file) ggsave(file, heatmap_plot(), width = max(8, 4 + 0.12 * nrow(req(results()))),
                                    height = heatmap_height() / 96, limitsize = FALSE))
  output$dl_example_mpa <- downloadHandler(function() "example_metaphlan4_vJan21.tsv",
                                           function(file) file.copy("data/mpa4_vJan21_test_data.tsv", file))
  output$dl_example_gtdb <- downloadHandler(function() "example_sylph_GTDBr220.tsv",
                                            function(file) file.copy("data/sylph_GTDBr220_test_data.tsv", file))

  # Crosswalk tab: identifiers for the selected database
  output$crosswalk_table <- renderDT({
    db <- input$database
    tab <- data.frame(Species = gsub("_", " ", XW$species_Jan21), Signature = XW$SIG,
                      Identifiers = gsub(",", ", ", XW[[db]]), check.names = FALSE)
    names(tab)[3] <- paste("Identifier(s) in", names(DB_CHOICES)[DB_CHOICES == db])
    datatable(tab, rownames = FALSE, filter = "top", fillContainer = FALSE, options = list(pageLength = 25, scrollX = TRUE)) |>
      formatStyle("Signature", color = styleEqual(c("SIG1", "SIG2"), unname(TOPO_COLOURS[c("SIG1", "SIG2")])), fontWeight = "bold")
  })
  output$dl_crosswalk <- downloadHandler(function() "toposcore_crosswalk.tsv",
                                         function(file) file.copy("data/toposcore_crosswalk.tsv", file))
}

shinyApp(ui, server)
