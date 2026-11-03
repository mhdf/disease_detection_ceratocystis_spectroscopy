# Tables 1, A1 and A2: final Random Forest and SVM models with 75 features, for the validation plants.
# With --leaf and --method only that part runs, which takes about 15 minutes for each combination.
source("R/setup.R")

args <- parse_args(commandArgs(TRUE))
final_sets <- read_result(FINAL_SETS_FILE)
combinations <- final_sets %>% filter(target_size == ARTICLE_SIZE) %>% distinct(method, leaf)
if (!is.null(args$leaf)) combinations <- combinations %>% filter(leaf == args$leaf)
if (!is.null(args$method)) combinations <- combinations %>% filter(method == args$method)

metrics <- list(); validation <- list(); calibration <- list()
for (i in seq_len(nrow(combinations))) for (pretreatment in PRETREATMENTS) {
  run <- list(method = combinations$method[i], leaf = combinations$leaf[i], pretreatment = pretreatment)
  bands <- final_sets %>%
    filter(method == run$method, leaf == run$leaf, pretreatment == !!pretreatment, target_size == ARTICLE_SIZE) %>%
    pull(bands) %>%
    split_bands()
  result <- run_final_models(run, bands)
  metrics[[length(metrics) + 1]] <- result$metrics
  validation[[length(validation) + 1]] <- result$validation
  calibration[[length(calibration) + 1]] <- result$calibration
  message(sprintf("%s %s %s finished", run$method, run$leaf, pretreatment))
}
tag <- paste(c(args$method, args$leaf), collapse = "_")
suffix <- if (nzchar(tag)) paste0("_", tag) else ""
write_table(bind_rows(metrics), file.path(OUTPUTS_DIR, paste0("validation_metrics", suffix, ".csv")))
write_table(bind_rows(validation), file.path(OUTPUTS_DIR, paste0("validation_predictions", suffix, ".csv")))
write_table(bind_rows(calibration), file.path(OUTPUTS_DIR, paste0("calibration_predictions", suffix, ".csv")))
cat("Tables written to", OUTPUTS_DIR, "\n")
