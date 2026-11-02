# Export the data and result tables of this repository from the original project folder.
# It is kept to document where these files came from: it needs the original files, which are not shared.
#
#   Rscript tools/export_original_results.R --source <project folder>

suppressPackageStartupMessages(library(dplyr))
source("R/config.R")
args <- parse_args(commandArgs(TRUE))
source_dir <- args$source
if (is.null(source_dir)) stop("Indicate the original project folder with --source")
raw_dir <- file.path(source_dir, "2_data")
out_dir <- file.path(source_dir, "5_outputs")
dir.create(RESULTS_DIR, showWarnings = FALSE)

write_csv_plain <- function(x, path) write.table(x, path, sep = ",", row.names = FALSE, na = "")
as_labels <- function(x) {
  if (is.factor(x)) return(as.character(x))
  if (is.numeric(x)) return(c("no", "yes")[x])
  x
}

# Runs of the article: folder of the original project, feature selection method and canopy layer
RUNS <- data.frame(
  folder = c("model_fitting_v17", "model_fitting_v18", "model_fitting_v19", "model_fitting_v21"),
  method = c("rfe", "rfe", "rfe", "genetic_algorithm"),
  leaf = c("apical", "middle", "basal", "apical"),
  stringsAsFactors = FALSE
)
PLSDA_FOLDER <- "model_fitting_v22"

# ---- Spectra: Experiment A, 12 to 60 hours after inoculation, plants in the analysis ----
SOURCE_TABLES <- c(reflect = "data_to_model.csv", SNV = "data_to_model_SNV.csv", CR = "data_to_model_CR.csv")
for (pretreatment in PRETREATMENTS) {
  table_5nm <- read.csv(file.path(out_dir, SOURCE_TABLES[[pretreatment]]), stringsAsFactors = FALSE) %>%
    filter(experiment == "A", time %in% PRE_SYMPTOM_TIMES) %>%
    select(-experiment)
  write_csv_plain(table_5nm, spectra_path(pretreatment))
}

raw <- read.csv(file.path(raw_dir, "Dados_experimento_v2.csv"), stringsAsFactors = FALSE)
names(raw)[1:6] <- c("plant", "clone", "inoculation", "leaf", "time", "sample")
raw_bands <- paste0("C", RAW_WAVELENGTHS)
names(raw)[7:ncol(raw)] <- raw_bands
kept <- read.csv(spectra_path("reflect"), stringsAsFactors = FALSE)
averages <- raw %>%
  filter(time %in% PRE_SYMPTOM_TIMES) %>%
  mutate(
    plant = as.integer(plant),
    clone = as.integer(clone),
    inoculation = c(nao = "no", sim = "yes")[inoculation],
    leaf = c(base = "basal", media = "middle", apical = "apical")[leaf],
    across(all_of(raw_bands), as.numeric)
  ) %>%
  semi_join(kept, by = c("plant", "leaf", "time")) %>%
  group_by(plant, clone, inoculation, leaf, time) %>%
  summarise(n_readings = n(), across(all_of(raw_bands), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
averages <- averages[match(paste(kept$plant, kept$leaf), paste(averages$plant, averages$leaf)), ]
averages_csv <- tempfile(fileext = ".csv")
write_csv_plain(averages, averages_csv)
connection <- xzfile(RAW_MEANS_FILE, "wb", compression = 9)
writeBin(readBin(averages_csv, "raw", file.size(averages_csv)), connection)
close(connection)

# ---- Severity: times reported in the article, without the three excluded plants ----
severity <- read.csv(file.path(out_dir, "severity_experiment_A_to_model.csv"), stringsAsFactors = FALSE) %>%
  filter(time %in% c(PRE_SYMPTOM_TIMES, POST_SYMPTOM_TIMES))
write_csv_plain(severity, SEVERITY_FILE)

# ---- Plants used for calibration and validation, as drawn in the original run ----
validation_plants <- readRDS(file.path(out_dir, PLSDA_FOLDER, "validation_results_all_nbr_vars_all.rds"))$predictions[[1]]$plant
plants <- kept %>% filter(leaf == "apical") %>% pull(plant)
write_csv_plain(
  data.frame(plant = plants, split = ifelse(plants %in% validation_plants, "validation", "calibration")),
  SPLIT_FILE
)

# ---- Candidate feature sets and their cross-validation results ----
candidates <- list()
for (r in seq_len(nrow(RUNS))) for (learner in LEARNERS) for (pretreatment in PRETREATMENTS) {
  for (size in TARGET_SIZES) for (repetition in seq_len(N_REPETITIONS)) for (fold in seq_len(N_FOLDS)) {
    file_name <- sprintf(
      "outer_results_test_kfold_seq_%02d_%s_%s_%s_%s_1_rem_pre_A_all_outer_loop_%d_nbr_vars_%d.rds",
      repetition, pretreatment, RUNS$leaf[r], RUNS$method[r], learner, fold, size
    )
    run <- readRDS(file.path(out_dir, RUNS$folder[r], file_name))
    parameters <- as.list(run$best_parameters)
    candidates[[length(candidates) + 1]] <- data.frame(
      method = RUNS$method[r], leaf = RUNS$leaf[r], pretreatment = pretreatment, target_size = size, learner = learner,
      repetition = repetition, fold = fold,
      precision = run$performance[1], recall = run$performance[2], bacc = run$performance[3], fbeta = run$performance[4],
      num_trees = if (learner == "ranger") parameters$num.trees else NA,
      mtry = if (learner == "ranger") parameters$mtry else NA,
      min_node_size = if (learner == "ranger") parameters$min.node.size else NA,
      cost = if (learner == "svm") parameters$cost else NA,
      gamma = if (learner == "svm") parameters$gamma else NA,
      bands = if (learner == "ranger") paste(run$selected_vars, collapse = ";") else NA,
      stringsAsFactors = FALSE
    )
  }
}
write_csv_plain(bind_rows(candidates), CANDIDATES_FILE)

# ---- Final feature sets and the 25 refits that give the standard deviations of Table 1 ----
final_sets <- list()
refits <- list()
for (r in seq_len(nrow(RUNS))) for (pretreatment in PRETREATMENTS) for (size in TARGET_SIZES) {
  rf <- readRDS(file.path(out_dir, RUNS$folder[r], sprintf("features_selected_final_ranger_%s_nbr_vars_%d.rds", pretreatment, size)))
  svm <- readRDS(file.path(out_dir, RUNS$folder[r], sprintf("features_selected_final_svm_%s_nbr_vars_%d.rds", pretreatment, size)))
  stopifnot(setequal(rf$selected_vars, svm$selected_vars))
  final_sets[[length(final_sets) + 1]] <- data.frame(
    method = RUNS$method[r], leaf = RUNS$leaf[r], pretreatment = pretreatment, target_size = size,
    n_bands = length(rf$selected_vars), bands = paste(rf$selected_vars, collapse = ";"), stringsAsFactors = FALSE
  )
  for (learner in LEARNERS) {
    performance <- if (learner == "ranger") rf$performance else svm$performance
    refits[[length(refits) + 1]] <- data.frame(
      method = RUNS$method[r], leaf = RUNS$leaf[r], pretreatment = pretreatment, target_size = size, learner = learner,
      resample = seq_len(nrow(performance)), as.data.frame(performance), stringsAsFactors = FALSE
    )
  }
}
write_csv_plain(bind_rows(final_sets), FINAL_SETS_FILE)
write_csv_plain(bind_rows(refits), REFITS_FILE)

# ---- Models of the article: validation metrics and predictions, out-of-fold predictions ----
metrics <- list()
validation <- list()
calibration <- list()
key <- function(method, leaf, pretreatment, learner) {
  data.frame(method = method, leaf = leaf, pretreatment = pretreatment, learner = learner, stringsAsFactors = FALSE)
}
clean_validation <- function(table, id) {
  cbind(id, as.data.frame(table)[, c("plant", "clone", "inoculation", "time")], response = as_labels(table$response))
}
clean_calibration <- function(table, id) {
  table <- as.data.frame(table)
  responses <- grep("^response_", names(table), value = TRUE)
  table[responses] <- lapply(table[responses], as_labels)
  cbind(id, table[, c("plant", "clone", "inoculation", "time", responses)])
}
for (r in seq_len(nrow(RUNS))) {
  table <- readRDS(file.path(out_dir, RUNS$folder[r], sprintf("results_validation_all_nbr_vars_%d.rds", ARTICLE_SIZE)))
  for (learner in LEARNERS) {
    # The runs of the middle and basal layers kept only the metrics
    has_predictions <- RUNS$leaf[r] == "apical"
    if (has_predictions) {
      predictions_validation <- readRDS(file.path(out_dir, RUNS$folder[r], sprintf("predictions_validation_%s_all_nbr_vars_%d.rds", learner, ARTICLE_SIZE)))
      predictions_calibration <- readRDS(file.path(out_dir, RUNS$folder[r], sprintf("predictions_calibration_%s_all_nbr_vars_%d.rds", learner, ARTICLE_SIZE)))
    }
    for (k in seq_along(PRETREATMENTS)) {
      id <- key(RUNS$method[r], RUNS$leaf[r], PRETREATMENTS[k], learner)
      row <- table %>% filter(dataset == PRETREATMENTS[k], .data$learner == !!learner)
      metrics[[length(metrics) + 1]] <- cbind(id, data.frame(
        target_size = ARTICLE_SIZE,
        num_trees = as.numeric(row$num.trees), mtry = as.numeric(row$mtry), min_node_size = as.numeric(row$min.node.size),
        cost = as.numeric(row$cost), gamma = as.numeric(row$gamma), ncomp = NA,
        precision = as.numeric(row$precision), recall = as.numeric(row$recall), bacc = as.numeric(row$bacc), fbeta = as.numeric(row$fbeta)
      ))
      if (has_predictions) {
        validation[[length(validation) + 1]] <- clean_validation(predictions_validation[[k]], id)
        calibration[[length(calibration) + 1]] <- clean_calibration(predictions_calibration[[k]], id)
      }
    }
  }
}
plsda_validation <- readRDS(file.path(out_dir, PLSDA_FOLDER, "validation_results_all_nbr_vars_all.rds"))
plsda_calibration <- readRDS(file.path(out_dir, PLSDA_FOLDER, "calibration_results_all_nbr_vars_all.rds"))
for (k in seq_along(PRETREATMENTS)) {
  id <- key("none", "apical", PRETREATMENTS[k], "plsda")
  row <- plsda_validation$performance %>% filter(dataset == PRETREATMENTS[k])
  metrics[[length(metrics) + 1]] <- cbind(id, data.frame(
    target_size = NA, num_trees = NA, mtry = NA, min_node_size = NA, cost = NA, gamma = NA, ncomp = PLSDA_COMPONENTS[[PRETREATMENTS[k]]],
    precision = row$precision, recall = row$recall, bacc = row$bacc, fbeta = row$fbeta
  ))
  validation[[length(validation) + 1]] <- clean_validation(plsda_validation$predictions[[k]], id)
  calibration[[length(calibration) + 1]] <- clean_calibration(plsda_calibration[[k]], id)
}
write_csv_plain(bind_rows(metrics), VALIDATION_METRICS_FILE)
write_csv_plain(bind_rows(validation), VALIDATION_PREDICTIONS_FILE)
write_csv_plain(bind_rows(calibration), CALIBRATION_PREDICTIONS_FILE)

# ---- PLS-DA: cross-validation of the number of components ----
plsda_cv <- list()
for (pretreatment in PRETREATMENTS) for (repetition in seq_len(N_REPETITIONS)) for (fold in seq_len(N_FOLDS)) {
  file_name <- sprintf("results_plsda_test_kfold_seq_%02d_%s_apical_vip_1_rem_pre_A_all_outer_loop_%d_nbr_vars_all.rds", repetition, pretreatment, fold)
  run <- readRDS(file.path(out_dir, PLSDA_FOLDER, file_name))
  plsda_cv[[length(plsda_cv) + 1]] <- data.frame(
    pretreatment = pretreatment, repetition = repetition, fold = fold, as.data.frame(run$performance), stringsAsFactors = FALSE
  )
}
write_csv_plain(bind_rows(plsda_cv), PLSDA_CV_FILE)

# ---- SHAP importance, in percent of the total of each model ----
shap <- list()
for (method in c("rfe", "genetic_algorithm")) for (pretreatment in PRETREATMENTS) {
  file_name <- sprintf("relative_shap_nbr_vars_%d_apical_%s_%s.rds", ARTICLE_SIZE, method, pretreatment)
  importance <- unlist(as.data.frame(readRDS(file.path(out_dir, file_name))))
  shap[[length(shap) + 1]] <- data.frame(
    method = method, pretreatment = pretreatment, learner = "ranger",
    wavelength = as.numeric(sub("C", "", names(importance))), importance = unname(importance)
  )
}
plsda_shap <- readRDS(file.path(out_dir, PLSDA_FOLDER, "relative_shap_all_nbr_vars_all.rds"))
for (pretreatment in PRETREATMENTS) {
  importance <- plsda_shap %>% filter(dataset == pretreatment) %>% select(-dataset) %>% unlist()
  shap[[length(shap) + 1]] <- data.frame(
    method = "none", pretreatment = pretreatment, learner = "plsda",
    wavelength = as.numeric(sub("C", "", names(importance))), importance = unname(importance)
  )
}
write_csv_plain(bind_rows(shap), SHAP_FILE)

cat("Exported to", DATA_DIR, "and", RESULTS_DIR, "\n")
