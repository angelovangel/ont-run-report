# helper functions for ont-run-report app

script_exists <- function(path) {
  file.exists(path) && file.access(path, 1) == 0
}
