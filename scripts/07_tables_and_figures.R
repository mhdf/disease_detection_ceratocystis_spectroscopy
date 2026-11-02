# Tables 1, A1 and A2 and Figures 1, 3, 4, 5 and A1, from the tables in results/ and data/. Takes about a minute.
source("R/setup.R")
source("R/summaries.R")
source("R/figures.R")

metrics <- read_result(VALIDATION_METRICS_FILE)
refits <- read_result(REFITS_FILE)
plsda_cv <- read_result(PLSDA_CV_FILE)
final_sets <- read_result(FINAL_SETS_FILE)
shap <- read_result(SHAP_FILE)
calibration <- read_result(CALIBRATION_PREDICTIONS_FILE)
validation <- read_result(VALIDATION_PREDICTIONS_FILE)
severity <- read_result(SEVERITY_FILE)

performance <- performance_table(metrics, refits, plsda_cv)
write_table(table_1(performance), file.path(RESULTS_DIR, "table_1.csv"))
write_table(table_a(performance, "ranger"), file.path(RESULTS_DIR, "table_a1.csv"))
write_table(table_a(performance, "svm"), file.path(RESULTS_DIR, "table_a2.csv"))

published <- figure_5_values(calibration, validation, as_published = TRUE)
intended <- figure_5_values(calibration, validation, as_published = FALSE)
write_table(published, file.path(RESULTS_DIR, "figure_5_values.csv"))
write_table(intended, file.path(RESULTS_DIR, "figure_5_values_intended.csv"))

save_figure <- function(plot, name, width, height) {
  dir.create(FIGURES_DIR, showWarnings = FALSE)
  ggsave(file.path(FIGURES_DIR, name), plot = plot, width = width, height = height, dpi = 200)
}
data_apical <- prepare_fit_data("reflect", "apical")
split_table <- read_result(SPLIT_FILE)
meta <- read_spectra("reflect") %>%
  filter(leaf == "apical") %>%
  mutate(dataset = split_table$split[match(plant, split_table$plant)]) %>%
  select(plant, clone, inoculation, time, dataset)
save_figure(figure_1(meta, make_outer_folds(data_apical$strata)), "figure_1.png", 11, 8)
save_figure(figure_3(read_spectra("reflect"), read_spectra("CR"), final_sets, shap), "figure_3.png", 9, 9)
save_figure(figure_4(shap), "figure_4.png", 14, 5)
save_figure(figure_5(published), "figure_5.png", 9, 4)
save_figure(figure_5(intended), "figure_5_intended.png", 9, 4)
save_figure(figure_a1(severity), "figure_a1.png", 10, 6)
cat("Tables written to", RESULTS_DIR, "and figures to", FIGURES_DIR, "\n")
