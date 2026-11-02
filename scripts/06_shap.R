# Figure 4: importance of the wavelengths with SHAP values.
# The simulations run in parallel without a fixed random stream, so the values change slightly between runs.
#
#   Rscript scripts/06_shap.R --cores 8
source("R/setup.R")

args <- parse_args(commandArgs(TRUE))
cores <- as.integer(if (is.null(args$cores)) max(1, parallel::detectCores(logical = FALSE) - 1) else args$cores)
final_sets <- read_result(FINAL_SETS_FILE)
metrics <- read_result(VALIDATION_METRICS_FILE)

rows <- list()
for (method in METHODS) for (pretreatment in PRETREATMENTS) {
  selected <- final_sets %>%
    filter(method == !!method, leaf == "apical", pretreatment == !!pretreatment, target_size == ARTICLE_SIZE)
  parameters <- metrics %>% filter(method == !!method, leaf == "apical", pretreatment == !!pretreatment, learner == "ranger")
  rows[[length(rows) + 1]] <- shap_rf(method, pretreatment, split_bands(selected$bands), parameters, cores)
  message(method, " ", pretreatment, " finished")
}
for (pretreatment in PRETREATMENTS) {
  rows[[length(rows) + 1]] <- shap_plsda(pretreatment, cores)
  message("PLS-DA ", pretreatment, " finished")
}
write_table(bind_rows(rows), file.path(OUTPUTS_DIR, basename(SHAP_FILE)))
cat("Table written to", OUTPUTS_DIR, "\n")
