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
addResourcePath("reports", app_tmp_dir)

# Sweep stale run artifacts left behind by crashed/killed sessions (e.g. a
# session that never fired onSessionEnded). Age-gated so it never deletes
# files from a run that's still active in another session on app restart.
stale_run_files <- list.files(app_tmp_dir, pattern = "^run_", full.names = TRUE)
if (length(stale_run_files)) {
  age_mins <- as.numeric(difftime(Sys.time(), file.info(stale_run_files)$mtime, units = "mins"))
  unlink(stale_run_files[age_mins > 60], recursive = TRUE)
}

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
sidebar <- sidebar(
  title = 'Controls', width = 400, gap = '0.5rem',
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
      col_widths = c(7, 5), gap = '0.5rem', class = 'mt-1',
      textInput('report_title',
        label = tags$span(
          "Report title", class = "small fw-semibold",
          tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                  "Title shown at the top of the generated HTML report.", placement = "right")
        ),
        value = 'ONT Run Report'),

      numericInput('sample_hz',
        label = tags$span(
          "Sampling (min)", class = "small fw-semibold",
          tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                  "Sampling interval for the Active-Pores chart in minutes.", placement = "right")
        ),
        value = 5, min = 1, max = 60, step = 1)
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
                       font_scale = 1.0, spacer = '0.7rem'),
  header   = tags$style(HTML("
    .status-dot {
      display: inline-block; width: 12px; height: 12px; border-radius: 50%;
      margin-right: 6px; vertical-align: middle;
    }
    .status-red    { background-color: #e74c3c; }
    .status-yellow { background-color: #f1c40f; }
    .status-green  { background-color: #2ecc71; }
  ")),
  sidebar  = sidebar,
  nav_panel(
    use_busy_spinner(spin = "double-bounce", position = 'top-right', color = '#E67E22'),
    title = 'Generate Report',
    div(class = "d-flex gap-2 mb-2",
      shinyjs::disabled(
        actionButton('start', 'Generate Report', class = 'btn-success flex-fill')
      ),
      actionButton('reset', 'Reset', class = 'flex-shrink-0')
    ),
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
            card_body(class = 'p-2', uiOutput('report_link_ui'))
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
    log_file     = NULL,
    status_file  = NULL,
    report_file  = NULL,
    report_url   = NULL,
    work_dir     = NULL,
    inner_script = NULL,
    is_running   = FALSE,
    show_log     = FALSE
  )

  # Remove only the files/dirs belonging to THIS session. Never glob the
  # shared app_tmp_dir by pattern - other sessions may have runs in flight
  # there, and deleting by pattern is what corrupts them.
  cleanup_session_files <- function() {
    paths <- c(rv$log_file, rv$status_file, rv$inner_script, rv$work_dir)
    paths <- paths[!vapply(paths, is.null, logical(1))]
    if (length(paths)) unlink(paths, recursive = TRUE)
  }

  # Catch-all: if the user closes the tab/browser without hitting Reset,
  # clean up this session's run artifacts automatically.
  session$onSessionEnded(function() {
    isolate(cleanup_session_files())
  })

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
        value = paste0("Run ", i),
        title = tagList(
          tags$span(id = paste0("status_dot_", i), class = "status-dot status-red"),
          paste0("Run ", i)
        ),
        icon = bs_icon("hdd-stack"),
        textInput(paste0("label_", i), label = "Label", value = paste0("Run ", i)),
        layout_columns(
          col_widths = c(6, 6), gap = '0.4rem',
          fileInput(paste0("pore_activity_", i),
            label = tags$span("pore_activity", class = "small",
              tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:3px;cursor:pointer;"),
                      "pore_activity_*.csv generated by MinKNOW during run.", placement = "right")),
            accept = ".csv", buttonLabel = "Browse…",
            placeholder = "CSV"),
          fileInput(paste0("throughput_", i),
            label = tags$span("throughput", class = "small",
              tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:3px;cursor:pointer;"),
                      "throughput_*.csv generated by MinKNOW during run.", placement = "right")),
            accept = ".csv", buttonLabel = "Browse…",
            placeholder = "CSV")
        )
      )
    })
    do.call(accordion, c(
      list(id = 'flowcell_accordion', open = FALSE, multiple = TRUE, class = 'mb-1'),
      panels
    ))
  })

  # ---------------------------------------------------------------------------
  # Per-run status dot (yellow = files needed, green = ready)
  # Created once for the max possible number of runs so they persist across
  # changes to num_runs (the accordion panels get re-rendered each time).
  # ---------------------------------------------------------------------------
  lapply(1:8, function(i) {
    observe({
      pa <- input[[paste0("pore_activity_", i)]]
      tp <- input[[paste0("throughput_", i)]]
      pa_ok <- !is.null(pa) && nrow(pa) > 0
      tp_ok <- !is.null(tp) && nrow(tp) > 0
      n_ok  <- sum(pa_ok, tp_ok)

      mismatch <- FALSE
      if (pa_ok && tp_ok) {
        pa_suffix <- sub("^pore_activity_", "", pa$name[1])
        tp_suffix <- sub("^throughput_",    "", tp$name[1])
        mismatch  <- !identical(pa_suffix, tp_suffix)
      }

      state <- if (mismatch)     "status-red"
               else if (n_ok == 2) "status-green"
               else if (n_ok == 1) "status-yellow"
               else                "status-red"
      all_states <- c("status-red", "status-yellow", "status-green")

      dot_id <- paste0("status_dot_", i)
      shinyjs::removeClass(id = dot_id, class = paste(setdiff(all_states, state), collapse = " "))
      shinyjs::addClass(id = dot_id, class = state)

      # Warn if the pore_activity/throughput filenames' hash suffixes don't match
      if (mismatch) {
        sendSweetAlert(
          title = "File mismatch",
          text  = sprintf(
            "The pore_activity and throughput files for Run %d don't look like they're from the same run:\n%s\n%s",
            i, pa$name[1], tp$name[1]
          ),
          type = "warning", btn_labels = "OK", showCloseButton = TRUE
        )
      } else if (n_ok == 2) {
        accordion_panel_close("flowcell_accordion", paste0("Run ", i))
      }
    })
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
      cat(sprintf("  --title %s \\\n", shQuote(input$report_title %||% "ONT Run Report")))
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
      shinyjs::enable('reset')
      shinyjs::html(id = 'start', 'Generate Report')
      hide_spinner()
      if (status == "0") {
        showNotification("Report generated successfully!", type = "message")
        # sendSweetAlert(title = NULL, text = "Report generated successfully!",
        #                type = "success", btn_labels = NA, closeOnClickOutside = TRUE, showCloseButton = TRUE)
        shinyjs::show('download_panel')
      } else {
        sendSweetAlert(title = "Error", text = "Report generation failed! Check the terminal output.",
                       type = "error", btn_labels = "OK", showCloseButton = TRUE)
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
    report_url   <- file.path("reports", basename(work_dir), "report.html")
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

    title_arg <- if (nzchar(input$report_title %||% "")) input$report_title else "ONT Report"

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

    rv$log_file     <- log_file
    rv$status_file  <- status_file
    rv$report_file  <- report_file
    rv$report_url   <- report_url
    rv$work_dir     <- work_dir
    rv$inner_script <- inner_script
    rv$is_running   <- TRUE
    rv$show_log     <- TRUE

    shinyjs::disable('controls')
    shinyjs::disable('reset')
    shinyjs::html(id = 'start', 'Generating...')
    shinyjs::hide('download_panel')
    show_spinner()
  })

  # ---------------------------------------------------------------------------
  # Report link (opens the generated HTML report in a new tab)
  # ---------------------------------------------------------------------------
  output$report_link_ui <- renderUI({
    req(rv$report_url)
    div(class = "d-flex gap-2",
      tags$a(
        href = rv$report_url, target = "_blank", rel = "noopener noreferrer",
        class = "btn btn-primary flex-fill",
        bs_icon("box-arrow-up-right", class = "me-1"), "Open HTML Report"
      ),
      downloadButton('download_report', 'Download',
        class = "btn btn-outline-primary flex-fill")
    )
  })

  output$download_report <- downloadHandler(
    filename = function() {
      title_arg <- if (nzchar(input$report_title %||% "")) input$report_title else "ONT Run Report"
      paste0(gsub("[^A-Za-z0-9_-]+", "_", title_arg), ".html")
    },
    content = function(file) {
      req(rv$report_file)
      file.copy(rv$report_file, file)
    }
  )

  # ---------------------------------------------------------------------------
  # Reset
  # ---------------------------------------------------------------------------
  observeEvent(input$reset, {
    cleanup_session_files()

    rv$log_file     <- NULL
    rv$status_file  <- NULL
    rv$report_file  <- NULL
    rv$report_url   <- NULL
    rv$work_dir     <- NULL
    rv$inner_script <- NULL
    rv$is_running   <- FALSE
    rv$show_log     <- FALSE

    updateSelectInput(session,  "num_runs",      selected = 2)
    updateTextInput(session,    "report_title",  value = "ONT Run Report")
    updateNumericInput(session, "sample_hz",     value = 5)
    for (i in 1:8) {
      shinyjs::reset(paste0("pore_activity_", i))
      shinyjs::reset(paste0("throughput_", i))
      updateTextInput(session, paste0("label_", i), value = paste0("Run ", i))
      dot_id <- paste0("status_dot_", i)
      shinyjs::removeClass(id = dot_id, class = "status-yellow status-green")
      shinyjs::addClass(id = dot_id, class = "status-red")
    }
    shinyjs::hide('download_panel')
    shinyjs::enable('reset')
  })
}

shinyApp(ui, server)
