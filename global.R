# helper functions for ont-run-report app

script_exists <- function(path) {
  file.exists(path) && file.access(path, 1) == 0
}

# Pair MinKNOW uploads pore_activity_<run>.csv / throughput_<run>.csv by their
# shared <run> suffix (flow cell ID + run hashes). `files` is a fileInput value.
# Returns the pairs (sorted by run) and the names of files that couldn't be paired.
pair_ont_files <- function(files) {
  kind <- ifelse(startsWith(files$name, "pore_activity_"), "pa",
          ifelse(startsWith(files$name, "throughput_"), "tp", NA))
  run  <- sub("\\.csv$", "", sub("^(pore_activity|throughput)_", "", files$name),
              ignore.case = TRUE)

  pairs <- list()
  unmatched <- files$name[is.na(kind)]
  for (r in sort(unique(run[!is.na(kind)]))) {
    pa <- which(run == r & kind %in% "pa")
    tp <- which(run == r & kind %in% "tp")
    if (length(pa) == 1 && length(tp) == 1) {
      pairs[[length(pairs) + 1]] <- list(
        run = r, flow_cell = sub("_.*$", "", r),
        pa_name = files$name[pa], pa_path = files$datapath[pa],
        tp_name = files$name[tp], tp_path = files$datapath[tp]
      )
    } else {
      unmatched <- c(unmatched, files$name[c(pa, tp)])
    }
  }

  # Default labels: the flow cell ID, numbered if a flow cell was run twice.
  fc <- vapply(pairs, `[[`, character(1), "flow_cell")
  for (i in seq_along(pairs)) {
    dup <- sum(fc[seq_len(i)] == fc[i])
    pairs[[i]]$label <- if (sum(fc == fc[i]) > 1) sprintf("%s (%d)", fc[i], dup) else fc[i]
  }
  list(pairs = pairs, unmatched = unmatched)
}

# Name of the generated report file for a sequencer type ('ont' or 'hifi').
report_file_name <- function(seq_type) {
  if (identical(seq_type, 'hifi')) 'hifi-report.html' else 'ont-report.html'
}
