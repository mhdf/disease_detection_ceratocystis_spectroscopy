# Table 1: PLS-DA with all wavelengths. Number of components in the 25 folds, then the final models.
source("R/setup.R")

components <- list(); metrics <- list(); validation <- list(); calibration <- list()
for (pretreatment in PRETREATMENTS) {
  for (repetition in seq_len(N_REPETITIONS)) {
    components[[length(components) + 1]] <- run_plsda_components(pretreatment, repetition)
  }
  final <- run_plsda_final(pretreatment)
  metrics[[pretreatment]] <- final$metrics
  validation[[pretreatment]] <- final$validation
  calibration[[pretreatment]] <- final$calibration
  message(pretreatment, " finished")
}
write_table(bind_rows(components), file.path(OUTPUTS_DIR, basename(PLSDA_CV_FILE)))
write_table(bind_rows(metrics), file.path(OUTPUTS_DIR, "validation_metrics_plsda.csv"))
write_table(bind_rows(validation), file.path(OUTPUTS_DIR, "validation_predictions_plsda.csv"))
write_table(bind_rows(calibration), file.path(OUTPUTS_DIR, "calibration_predictions_plsda.csv"))
cat("Tables written to", OUTPUTS_DIR, "\n")
