library(shiny)
library(shinyWidgets)
library(shinyjs)
library(bslib)
library(bsicons)
library(shinybusy)
library(digest)
library(commonmark)

source('global.R')

# Local bin/ scripts take priority; env vars are a fallback for deployment.
app_dir          <- normalizePath(getwd())
local_ont  <- file.path(app_dir, 'bin', 'generate-ont-report.py')
local_hifi <- file.path(app_dir, 'bin', 'generate-hifi-report.py')
script_path      <- if (file.exists(local_ont))  local_ont  else Sys.getenv('SCRIPT_PATH',      unset = '/usr/local/bin/generate-ont-report.py')
hifi_script_path <- if (file.exists(local_hifi)) local_hifi else Sys.getenv('HIFI_SCRIPT_PATH', unset = '/usr/local/bin/generate-hifi-report.py')

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
  title = NULL, width = 390, gap = '0.1rem',
  div(id = 'controls',
    # ---- Sequencer type toggle ----
    radioButtons('seq_type',
      label = tags$span(
        "Sequencer", class = "small fw-semibold",
        tooltip(
          bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
          tags$span("ONT: pore_activity_*.csv + throughput_*.csv, with optional browser-processed sequencing summaries.", tags$br(),
                    "PacBio HiFi: *run-qc-export*.csv (SMRT Link), with optional *.barcodes.summary.csv demultiplexing files."),
          placement = "left"
        )
      ),
      choices = c('ONT' = 'ont', 'PacBio HiFi' = 'hifi'),
      selected = 'ont', inline = TRUE
    ),

    tags$hr(class = 'my-1'),

    # ---- ONT: all MinKNOW CSVs at once, paired by run ID ----
    conditionalPanel(
      condition = "input.seq_type == 'ont'",
      fileInput('ont_files',
        label = tags$span(
          "MinKNOW CSV files", class = "small fw-semibold",
          tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                  "Select the pore_activity_*.csv and throughput_*.csv files for all flow cells at once.",
                  "They are paired by the run ID in their file names.", placement = "right")
        ),
        accept = ".csv", buttonLabel = "Browse…",
        placeholder = "pore_activity + throughput CSVs", multiple = TRUE
      ),
      uiOutput('ont_pairs_ui')
    ),

    # ---- HiFi: run-qc-export CSVs ----
    conditionalPanel(
      condition = "input.seq_type == 'hifi'",
      fileInput('hifi_files',
        label = tags$span(
          "run-qc-export and barcodes CSV file(s)", class = "small fw-semibold",
          tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                  "One or more *run-qc-export*.csv files exported from SMRT Link.",
                  "Optionally add *.barcodes.summary.csv files to get a Demultiplexing section;",
                  "each is matched to its SMRT Cell by the instrument, movie date and sample name in its file name.",
                  placement = "right")
        ),
        accept = ".csv", buttonLabel = "Browse…",
        placeholder = "Select CSV file(s)", multiple = TRUE
      )
    ),

    layout_columns(
      col_widths = c(7, 5), gap = '0.5rem', class = 'mt-1',
      textInput('report_title',
        label = tags$span(
          "Report title", class = "small fw-semibold",
          tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                  "Title shown at the top of the generated HTML report.", placement = "left")
        ),
        value = 'ONT Run Report'),

      # Sampling interval is ONT-only
      conditionalPanel(
        condition = "input.seq_type == 'ont'",
        numericInput('sample_hz',
          label = tags$span(
            "Sampling (min)", class = "small fw-semibold",
            tooltip(bs_icon("question-circle", class = "text-muted small", style = "margin-left:4px;cursor:pointer;"),
                    "Sampling interval for the Active-Pores chart in minutes.", placement = "left")
          ),
          value = 5, min = 1, max = 60, step = 1)
      )
    )
  )
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- page_navbar(
  useShinyjs(),
  fillable = TRUE,
  title    = 'TGS Run Report',
  theme    = bs_theme(
    version = 5,
    bg = '#f4f6fb', fg = '#172238',
    primary = '#465cf5', secondary = '#66738a',
    success = '#129b78'
  ),
  header   = tags$head(
    tags$link(rel = 'stylesheet', type = 'text/css', href = 'app.css'),
    tags$script(src = "nx-summary.js")
  ),
  sidebar  = sidebar,
  nav_panel(
    use_busy_spinner(spin = "double-bounce", position = 'top-right', color = '#E67E22'),
    title = 'Generate Report',
    div(class = "workspace-shell",
      div(class = "action-toolbar",
        shinyjs::disabled(
          actionButton('start', 'Generate Report', class = 'btn-success flex-fill')
        ),
        actionButton('reset', 'Reset', class = 'btn-outline-secondary flex-shrink-0')
      ),
      card(
        class = 'console-card',
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
            class = 'ready-card',
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
    show_log     = FALSE,
    # fileInput values survive shinyjs::reset(); this hides them until the next upload
    ont_cleared  = FALSE,
    nx_slots     = 0
  )

  observeEvent(input$ont_files, rv$ont_cleared <- FALSE)

  ont_pairs <- reactive({
    f <- input$ont_files
    if (rv$ont_cleared || is.null(f) || nrow(f) == 0) {
      return(list(pairs = list(), unmatched = character(0)))
    }
    pair_ont_files(f)
  })

  # Labels typed in the flow-cell panels, falling back to the flow cell ID.
  ont_labels <- function(pairs) {
    vapply(seq_along(pairs), function(i) {
      lbl <- input[[paste0("label_", i)]]
      if (nzchar(lbl %||% "")) lbl else pairs[[i]]$label
    }, character(1))
  }

  # Clear browser-computed Nx values for every flow-cell slot used so far.
  reset_nx_summaries <- function() {
    for (i in seq_len(max(rv$nx_slots, 1))) session$sendCustomMessage("resetNxSummary", i)
  }

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

  # When switching sequencer type, reset the report title default
  observeEvent(input$seq_type, {
    reset_nx_summaries()
    current_title <- input$report_title %||% ""
    if (input$seq_type == 'hifi' && current_title %in% c('', 'ONT Run Report')) {
      updateTextInput(session, 'report_title', value = 'HiFi Sequencing Report')
    } else if (input$seq_type == 'ont' && current_title %in% c('', 'HiFi Sequencing Report')) {
      updateTextInput(session, 'report_title', value = 'ONT Run Report')
    }
  }, ignoreInit = TRUE)

  # Flow-cell panels are re-rendered when the pairing changes, so Nx values
  # computed for the previous panels no longer belong to the same flow cells.
  observeEvent(ont_pairs(), {
    isolate(reset_nx_summaries())
    rv$nx_slots <- max(isolate(rv$nx_slots), length(ont_pairs()$pairs))
  }, ignoreInit = TRUE)

  # Check scripts on startup
  observe({
    if (!script_exists(script_path)) {
      showNotification('generate-ont-report.py not found or not executable!', type = 'error', duration = NULL)
    }
  })

  observe({
    req(input$seq_type == 'hifi')
    if (!script_exists(hifi_script_path)) {
      showNotification('generate-hifi-report.py not found or not executable!', type = 'warning', duration = NULL)
    }
  })

  # ---------------------------------------------------------------------------
  # ONT flow cells: one panel per paired run (label + optional summary)
  # ---------------------------------------------------------------------------
  output$ont_pairs_ui <- renderUI({
    p <- ont_pairs()
    if (!length(p$pairs) && !length(p$unmatched)) {
      return(tags$small(class = "form-text text-muted d-block mb-2",
        "Select every pore_activity and throughput CSV together; files from the same run are paired automatically."))
    }

    panels <- lapply(seq_along(p$pairs), function(i) {
      pair <- p$pairs[[i]]
      accordion_panel(
        value = paste0("Run ", i),
        title = tagList(tags$span(class = "status-dot status-green"), pair$flow_cell),
        icon = bs_icon("hdd-stack"),
        textInput(paste0("label_", i), label = "Label", value = pair$label),
        tags$small(class = "form-text text-muted d-block mb-2 text-break",
          style = "font-size:10.5px;line-height:1.35;",
          pair$pa_name, tags$br(), pair$tp_name),
        div(class = "mb-2 nx-summary-control",
          tags$label(
            "sequencing_summary (optional)",
            `for` = paste0("sequencing_summary_file_", i),
            class = "form-label small"
          ),
          tags$input(
            type = "file",
            id = paste0("sequencing_summary_file_", i),
            accept = ".txt,.tsv",
            class = "form-control form-control-sm nx-summary-file"
          ),
          tags$small(
            id = paste0("sequencing_summary_status_", i),
            class = "form-text text-muted",
            ""
          )
        )
      )
    })

    tagList(
      if (length(p$unmatched)) {
        div(class = "alert alert-warning py-1 px-2 small mb-2",
          bs_icon("exclamation-triangle", class = "me-1"),
          tags$strong("Not paired, will be ignored: "),
          paste(p$unmatched, collapse = ", "))
      },
      if (length(panels)) {
        tagList(
          tags$small(class = "form-text text-muted d-block mb-1",
            sprintf("%d flow cell%s paired", length(panels), if (length(panels) == 1) "" else "s")),
          do.call(accordion, c(
            list(id = 'flowcell_accordion', open = FALSE, multiple = TRUE, class = 'mb-1'),
            panels
          ))
        )
      }
    )
  })

  # ---------------------------------------------------------------------------
  # Enable Start button when all files uploaded
  # ---------------------------------------------------------------------------
  observe({
    seq_type <- input$seq_type %||% 'ont'
    all_ready <- if (seq_type == 'hifi') {
      hf <- input$hifi_files
      !is.null(hf) && nrow(hf) > 0
    } else {
      length(ont_pairs()$pairs) > 0
    }

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
      seq_type <- input$seq_type %||% 'ont'
      cat("Command Preview:\n\n")
      if (seq_type == 'hifi') {
        hf <- input$hifi_files
        file_args <- if (!is.null(hf) && nrow(hf) > 0) hf$name else "[run-qc-export.csv ...]"
        cat("generate-hifi-report.py \\\n")
        cat(sprintf("  --title %s \\\n", shQuote(input$report_title %||% "HiFi Sequencing Report")))
        cat("  --out hifi-report.html \\\n")
        cat("  ", paste(shQuote(file_args), collapse = " "), "\n")
      } else {
        pairs <- ont_pairs()$pairs
        n <- length(pairs)
        cat(sprintf("generate-ont-report.py \\\n"))
        cat(sprintf("  --title %s \\\n", shQuote(input$report_title %||% "ONT Run Report")))
        cat(sprintf("  --sample-hz %s \\\n", input$sample_hz %||% 5))
        cat("  --out ont-report.html \\\n")
        if (n == 0) {
          cat("  --pore-activity [pore_activity_*.csv ...] \\\n")
          cat("  --throughput [throughput_*.csv ...]\n")
          return(invisible())
        }
        pa_names  <- vapply(pairs, `[[`, character(1), "pa_name")
        tp_names  <- vapply(pairs, `[[`, character(1), "tp_name")
        lbls      <- ont_labels(pairs)
        nx_values <- lapply(seq_len(n), function(i) input[[paste0("sequencing_summary_nx_", i)]])
        cat("  --pore-activity", paste(shQuote(pa_names), collapse = " "), "\\\n")
        cat("  --throughput",    paste(shQuote(tp_names), collapse = " "), "\\\n")
        has_nx <- !vapply(nx_values, is.null, logical(1))
        if (any(has_nx) && all(has_nx)) {
          nx_json <- paste0(
            "[", paste(vapply(nx_values, as.character, character(1)), collapse = ","), "]"
          )
          cat("  --nx-values", shQuote(nx_json), "\\\n")
        } else if (any(has_nx)) {
          cat("  # Process a sequencing summary for every run to enable the Nx plot\n")
        }
        cat("  --labels",        paste(shQuote(lbls),     collapse = " "), "\n")
      }
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
    seq_type <- input$seq_type %||% 'ont'
    if (seq_type != 'hifi') {
      n <- length(ont_pairs()$pairs)
      nx_values <- lapply(seq_len(n), function(i) {
        input[[paste0("sequencing_summary_nx_", i)]]
      })
      has_nx <- !vapply(nx_values, is.null, logical(1))
      if (any(has_nx) && !all(has_nx)) {
        showNotification(
          "Process a sequencing summary for every run, or remove the partial selections.",
          type = "error"
        )
        return()
      }
    }

    run_id       <- digest::digest(Sys.time(), algo = 'crc32')
    work_dir     <- file.path(app_tmp_dir, paste0("run_", run_id))
    dir.create(work_dir, showWarnings = FALSE)

    log_file     <- file.path(app_tmp_dir, paste0("run_", run_id, ".log"))
    status_file  <- file.path(app_tmp_dir, paste0("run_", run_id, ".status"))
    report_name  <- report_file_name(seq_type)
    report_file  <- file.path(work_dir, report_name)
    report_url   <- file.path("reports", basename(work_dir), report_name)
    inner_script <- file.path(app_tmp_dir, paste0("run_", run_id, ".sh"))

    title_arg <- if (nzchar(input$report_title %||% "")) input$report_title else
      if (seq_type == 'hifi') "HiFi Sequencing Report" else "ONT Run Report"

    if (seq_type == 'hifi') {
      # ---- HiFi mode ----
      hf <- input$hifi_files
      hifi_paths <- vapply(seq_len(nrow(hf)), function(i) {
        dest <- file.path(work_dir, hf$name[i])
        file.copy(hf$datapath[i], dest, overwrite = TRUE)
        dest
      }, character(1))

      py_cmd <- paste(
        "python3", shQuote(hifi_script_path),
        paste(shQuote(hifi_paths), collapse = " "),
        "--title", shQuote(title_arg),
        "--out",   shQuote(report_file)
      )
    } else {
      # ---- ONT mode ----
      pairs <- ont_pairs()$pairs
      n <- length(pairs)
      pa_paths <- character(n); tp_paths <- character(n)
      lbl_args <- ont_labels(pairs)
      nx_values <- lapply(seq_len(n), function(i) {
        input[[paste0("sequencing_summary_nx_", i)]]
      })
      use_nx <- !is.null(nx_values[[1]])

      for (i in seq_len(n)) {
        pa_dest <- file.path(work_dir, pairs[[i]]$pa_name)
        tp_dest <- file.path(work_dir, pairs[[i]]$tp_name)
        file.copy(pairs[[i]]$pa_path, pa_dest, overwrite = TRUE)
        file.copy(pairs[[i]]$tp_path, tp_dest, overwrite = TRUE)
        pa_paths[i] <- pa_dest
        tp_paths[i] <- tp_dest
      }

      nx_arg <- if (use_nx) {
        nx_json <- paste0(
          "[", paste(vapply(nx_values, as.character, character(1)), collapse = ","), "]"
        )
        paste("--nx-values", shQuote(nx_json))
      } else {
        ""
      }
      py_cmd <- paste(
        "python3", shQuote(script_path),
        "--pore-activity", paste(shQuote(pa_paths), collapse = " "),
        "--throughput",    paste(shQuote(tp_paths), collapse = " "),
        nx_arg,
        "--labels",        paste(shQuote(lbl_args), collapse = " "),
        "--title",         shQuote(title_arg),
        "--sample-hz",     as.integer(input$sample_hz),
        "--out",           shQuote(report_file)
      )
    }

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
    filename = function() basename(rv$report_file %||% report_file_name(input$seq_type %||% 'ont')),
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

    updateRadioButtons(session, "seq_type",      selected = 'ont')
    updateTextInput(session,    "report_title",  value = "ONT Run Report")
    updateNumericInput(session, "sample_hz",     value = 5)
    shinyjs::reset('hifi_files')
    shinyjs::reset('ont_files')
    rv$ont_cleared <- TRUE
    reset_nx_summaries()
    shinyjs::hide('download_panel')
    shinyjs::enable('reset')
  })
}

shinyApp(ui, server)
