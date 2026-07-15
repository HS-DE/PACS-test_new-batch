#!/usr/bin/env Rscript

parts <- c(
  "script/parts/01_setup_and_functions.R",
  "script/parts/02_load_and_audit.R",
  "script/parts/03_help_sets_and_factors.R",
  "script/parts/04_strategies.R",
  "script/parts/05_evaluate_strategies.R",
  "script/parts/06_plots_and_report.R"
)

missing_parts <- parts[!file.exists(parts)]
if (length(missing_parts)) {
  stop("Missing analysis script parts: ", paste(missing_parts, collapse = ", "))
}

for (part in parts) {
  message("Running ", part)
  sys.source(part, envir = globalenv())
}
