library(shiny)
library(shinyWidgets)
library(shinyjs)
library(bslib)
library(bsicons)
library(shinybusy)
library(digest)
library(commonmark)

source('global.R')

script_path <- Sys.getenv('SCRIPT_PATH', unset = '/usr/local/bin/generate-ont-report.py')

app_tmp_dir <- file.path(getwd(), "tmp")
dir.create(app_tmp_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
sidebar <- sidebar(
  title = 'Controls', width = 300, gap = '0.5rem',
  div(id = 'controls',

    selectInput('num_runs',
      label = tags$span(
        "Number of flow cells", class = "small fw-semibold",
        tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                "How many ONT flow cells to include in the report.", placement = "right")
      ),
      choices = 1:8, selected = 2),

    uiOutput('file_inputs_ui'),

    layout_columns(
      col_widths = c(8, 4), gap = '0.5rem', class = 'mt-1',
      textInput('report_title',
        label = tags$span(
          "Report title", class = "small fw-semibold",
          tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                  "Title shown at the top of the generated HTML report.", placement = "right")
        ),
        value = 'ONT Sequencing Report'),

      numericInput('sample_hz',
        label = tags$span(
          "Sampling (min)", class = "small fw-semibold",
          tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                  "Sampling interval for the Active-Pores chart in minutes.", placement = "right")
        ),
        value = 5, min = 1, max = 60, step = 1)
    ),

    div(class = "d-flex gap-2 mt-1",
      shinyjs::disabled(
        actionButton('start', 'Generate Report', class = 'btn-success flex-fill')
      ),
      actionButton('reset', 'Reset', class = 'flex-shrink-0')
    )
  )
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- page_navbar(
  useShinyjs(),
  fillable = TRUE,
  title    = 'ONT Run Report',
  theme    = bs_theme(bootswatch = 'yeti', primary = '#2E4053',
                       font_scale = 0.82, spacer = '0.7rem'),
  sidebar  = sidebar,
  nav_panel(
    use_busy_spinner(spin = "double-bounce", position = 'top-right', color = '#E67E22'),
    title = 'Generate Report',
    layout_column_wrap(
      width = NULL, fill = TRUE, gap = '0.5rem',
      style = htmltools::css(grid_template_columns = "1fr"),
      card(
        full_screen = TRUE,
        card_header(
          class = 'py-1 d-flex align-items-center gap-1 small fw-semibold',
          'Terminal output',
          tooltip(bs_icon("question-circle", class = "text-muted small"),
                  "Live execution log or command preview.", placement = "right")
        ),
        card_body(class = 'p-2', verbatimTextOutput('stdout'))
      ),
      shinyjs::hidden(
        div(id = 'download_panel',
          card(
            card_header(class = 'py-1 small fw-semibold', 'Report ready'),
            card_body(class = 'p-2',
              downloadButton('download_report', 'Download HTML Report',
                             class = 'btn-primary', style = 'width:100%')
            )
          )
        )
      )
    )
  ),
  nav_panel(
    title = 'README',
    card(
      style = "height:100%;overflow-y:auto;padding:12px;",
      tags$style(HTML("
        .readme-container table { width:100%; margin:15px 0; border-collapse:collapse; }
        .readme-container th, .readme-container td { padding:8px 12px; border:1px solid #dee2e6; }
        .readme-container th { background:#f8f9fa; font-weight:bold; }
        .readme-container blockquote { padding:10px 20px; margin:15px 0; border-left:5px solid #E67E22; background:#fcf8f2; }
        .readme-container code { background:#f8f9fa; padding:2px 4px; border-radius:4px; }
        .readme-container pre code { background:transparent; padding:0; }
      ")),
      div(class = "readme-container",
        htmltools::HTML(commonmark::markdown_html(
          paste(readLines('README.md', warn = FALSE), collapse = '\n'),
          extensions = TRUE
        ))
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  rv <- reactiveValues(
    log_file    = NULL,
    status_file = NULL,
    report_file = NULL,
    is_running  = FALSE,
    show_log    = FALSE
  )

  # Check script on startup
  observe({
    if (!script_exists(script_path)) {
      showNotification('generate-ont-report.py not found or not executable!', type = 'error', duration = NULL)
    }
  })

  # ---------------------------------------------------------------------------
  # Dynamic file input UI
  # ---------------------------------------------------------------------------
  output$file_inputs_ui <- renderUI({
    n <- req(as.integer(input$num_runs))
    panels <- lapply(seq_len(n), function(i) {
      accordion_panel(
        title = paste0("Flow cell ", i), icon = bs_icon("hdd-stack"),
        textInput(paste0("label_", i), label = "Label", value = paste0("Run ", i)),
        layout_columns(
          col_widths = c(6, 6), gap = '0.4rem',
          fileInput(paste0("pore_activity_", i),
            label = tags$span("pore_activity", class = "small",
              tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:3px;cursor:pointer;"),
                      "pore_activity_*.csv exported by MinKNOW.", placement = "right")),
            accept = ".csv", buttonLabel = "Browse…",
            placeholder = "CSV"),
          fileInput(paste0("throughput_", i),
            label = tags$span("throughput", class = "small",
              tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:3px;cursor:pointer;"),
                      "throughput_*.csv exported by MinKNOW.", placement = "right")),
            accept = ".csv", buttonLabel = "Browse…",
            placeholder = "CSV")
        )
      )
    })
    do.call(accordion, c(
      list(id = 'flowcell_accordion', open = 'Flow cell 1', multiple = TRUE, class = 'mb-1'),
      panels
    ))
  })

  # ---------------------------------------------------------------------------
  # Enable Start button when all files uploaded
  # ---------------------------------------------------------------------------
  observe({
    n <- req(as.integer(input$num_runs))
    all_ready <- all(vapply(seq_len(n), function(i) {
      pa <- input[[paste0("pore_activity_", i)]]
      tp <- input[[paste0("throughput_", i)]]
      !is.null(pa) && nrow(pa) > 0 && !is.null(tp) && nrow(tp) > 0
    }, logical(1)))

    if (all_ready && !rv$is_running) shinyjs::enable('start')
    else                             shinyjs::disable('start')
  })

  # ---------------------------------------------------------------------------
  # Poll the log file
  # ---------------------------------------------------------------------------
  poll_log <- reactivePoll(
    500, session,
    checkFunc = function() {
      lf <- rv$log_file
      if (is.null(lf) || !file.exists(lf)) return(0)
      file.info(lf)$size
    },
    valueFunc = function() {
      lf <- rv$log_file
      if (is.null(lf) || !file.exists(lf)) return("")
      paste(readLines(lf, warn = FALSE), collapse = "\n")
    }
  )

  # ---------------------------------------------------------------------------
  # Terminal output
  # ---------------------------------------------------------------------------
  output$stdout <- renderPrint({
    if (rv$show_log) {
      cat(poll_log())
    } else {
      n <- req(as.integer(input$num_runs))
      cat("Command Preview:\n\n")
      cat(sprintf("generate-ont-report.py \\\n"))
      cat(sprintf("  --title %s \\\n", shQuote(input$report_title %||% "ONT Sequencing Report")))
      cat(sprintf("  --sample-hz %s \\\n", input$sample_hz %||% 5))
      cat("  --out report.html \\\n")
      pa_names  <- character(n); tp_names  <- character(n); lbls <- character(n)
      for (i in seq_len(n)) {
        pa  <- input[[paste0("pore_activity_", i)]]
        tp  <- input[[paste0("throughput_", i)]]
        lbl <- input[[paste0("label_", i)]]
        pa_names[i] <- if (!is.null(pa) && nrow(pa) > 0) pa$name[1] else paste0("[pore_activity_", i, ".csv]")
        tp_names[i] <- if (!is.null(tp) && nrow(tp) > 0) tp$name[1] else paste0("[throughput_",    i, ".csv]")
        lbls[i]     <- if (nzchar(lbl %||% "")) lbl else paste0("Run ", i)
      }
      cat("  --pore-activity", paste(shQuote(pa_names), collapse = " "), "\\\n")
      cat("  --throughput",    paste(shQuote(tp_names), collapse = " "), "\\\n")
      cat("  --labels",        paste(shQuote(lbls),     collapse = " "), "\n")
    }
  })

  # Auto-scroll terminal
  observe({
    req(rv$is_running); poll_log()
    runjs("var el=document.getElementById('stdout'); if(el) el.parentElement.scrollTo({top:1e9,behavior:'smooth'});")
  })

  # ---------------------------------------------------------------------------
  # Poll status file for completion
  # ---------------------------------------------------------------------------
  observe({
    req(rv$is_running)
    invalidateLater(1000, session)
    sf <- rv$status_file
    if (is.null(sf) || !file.exists(sf)) return()
    status <- trimws(readLines(sf, warn = FALSE)[1])
    isolate({
      rv$is_running <- FALSE
      shinyjs::enable('controls')
      shinyjs::html(id = 'start', 'Generate Report')
      hide_spinner()
      if (status == "0") {
        sendSweetAlert(title = NULL, text = "Report generated successfully!",
                       type = "success", btn_labels = NA, closeOnClickOutside = TRUE)
        shinyjs::show('download_panel')
      } else {
        sendSweetAlert(title = "Error", text = "Report generation failed! Check the terminal output.",
                       type = "error", btn_labels = "OK")
      }
    })
  })

  # ---------------------------------------------------------------------------
  # Start
  # ---------------------------------------------------------------------------
  observeEvent(input$start, {
    n <- as.integer(input$num_runs)

    run_id       <- digest::digest(Sys.time(), algo = 'crc32')
    work_dir     <- file.path(app_tmp_dir, paste0("run_", run_id))
    dir.create(work_dir, showWarnings = FALSE)

    log_file     <- file.path(app_tmp_dir, paste0("run_", run_id, ".log"))
    status_file  <- file.path(app_tmp_dir, paste0("run_", run_id, ".status"))
    report_file  <- file.path(work_dir, "report.html")
    inner_script <- file.path(app_tmp_dir, paste0("run_", run_id, ".sh"))

    pa_paths <- character(n); tp_paths <- character(n); lbl_args <- character(n)

    for (i in seq_len(n)) {
      pa  <- input[[paste0("pore_activity_", i)]]
      tp  <- input[[paste0("throughput_",    i)]]
      lbl <- input[[paste0("label_",         i)]]

      pa_dest <- file.path(work_dir, pa$name[1])
      tp_dest <- file.path(work_dir, tp$name[1])
      file.copy(pa$datapath[1], pa_dest, overwrite = TRUE)
      file.copy(tp$datapath[1], tp_dest, overwrite = TRUE)

      pa_paths[i] <- pa_dest
      tp_paths[i] <- tp_dest
      lbl_args[i] <- if (nzchar(lbl %||% "")) lbl else paste0("Run ", i)
    }

    title_arg <- if (nzchar(input$report_title %||% "")) input$report_title else "ONT Sequencing Report"

    py_cmd <- paste(
      "python3", shQuote(script_path),
      "--pore-activity", paste(shQuote(pa_paths), collapse = " "),
      "--throughput",    paste(shQuote(tp_paths), collapse = " "),
      "--labels",        paste(shQuote(lbl_args), collapse = " "),
      "--title",         shQuote(title_arg),
      "--sample-hz",     as.integer(input$sample_hz),
      "--out",           shQuote(report_file)
    )

    writeLines(c(
      "#!/bin/bash",
      py_cmd,
      paste0("echo $? > ", shQuote(status_file))
    ), inner_script)
    system2("chmod", c("+x", inner_script))

    system(paste0(
      "nohup bash ", shQuote(inner_script),
      " >> ", shQuote(log_file),
      " 2>&1 </dev/null &"
    ))

    rv$log_file    <- log_file
    rv$status_file <- status_file
    rv$report_file <- report_file
    rv$is_running  <- TRUE
    rv$show_log    <- TRUE

    shinyjs::disable('controls')
    shinyjs::html(id = 'start', 'Generating...')
    shinyjs::hide('download_panel')
    show_spinner()
  })

  # ---------------------------------------------------------------------------
  # Download
  # ---------------------------------------------------------------------------
  output$download_report <- downloadHandler(
    filename = function() paste0("ont-report-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".html"),
    content  = function(file) {
      rf <- rv$report_file
      if (!is.null(rf) && file.exists(rf)) file.copy(rf, file)
    },
    contentType = "text/html"
  )

  # ---------------------------------------------------------------------------
  # Reset
  # ---------------------------------------------------------------------------
  observeEvent(input$reset, {
    rv$log_file    <- NULL
    rv$status_file <- NULL
    rv$report_file <- NULL
    rv$is_running  <- FALSE
    rv$show_log    <- FALSE

    # Clean up tmp files
    tmp_files <- list.files(app_tmp_dir, pattern = "^run_", full.names = TRUE, recursive = TRUE)
    unlink(tmp_files, recursive = TRUE)

    updateSelectInput(session,  "num_runs",      selected = 2)
    updateTextInput(session,    "report_title",  value = "ONT Sequencing Report")
    updateNumericInput(session, "sample_hz",     value = 5)
    shinyjs::hide('download_panel')
  })
}

shinyApp(ui, server)
