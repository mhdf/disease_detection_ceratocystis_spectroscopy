# Checks of this repository against the article and the original run. From the repository folder:
#
#   Rscript tests/test_reproduction.R            data, folds, tables and Figure 5 (about a minute)
#   Rscript tests/test_reproduction.R --models   also refits models and compares with the original run (about an hour)

source("R/setup.R")
source("R/summaries.R")

run_models <- "--models" %in% commandArgs(TRUE)
failures <- 0
check <- function(name, ok, detail = "") {
  cat(sprintf("%-4s %s %s\n", if (isTRUE(ok)) "ok" else "FAIL", name, detail))
  if (!isTRUE(ok)) failures <<- failures + 1
}
TOLERANCE <- 1e-8

# ---- Data ----
split_table <- read_result(SPLIT_FILE)
check("77 plants: 61 for calibration and 16 for validation",
      nrow(split_table) == 77 && sum(split_table$split == "calibration") == 61 && sum(split_table$split == "validation") == 16)
for (pretreatment in PRETREATMENTS) {
  spectra <- read_spectra(pretreatment)
  check(sprintf("%s table has 77 plants in 3 canopy layers and 411 bands", pretreatment),
        nrow(spectra) == 231 && all(BANDS %in% names(spectra)) && length(BANDS) == 411 &&
          all(table(spectra$leaf) == 77) && all(spectra$time %in% PRE_SYMPTOM_TIMES))
}
apical <- read_spectra("reflect") %>% filter(leaf == "apical")
counts <- apical %>% count(clone, inoculation) %>% pull(n)
check("Apical samples per clone and treatment are 20, 20, 20 and 17", identical(as.integer(counts), c(20L, 20L, 20L, 17L)))
check("Excluded plants are not in the data", !any(c(39, 57, 97) %in% apical$plant) && !any(c(39, 57, 97) %in% read_result(SEVERITY_FILE)$plant))

prepared <- prepare_spectra()
for (pretreatment in PRETREATMENTS) {
  difference <- max(abs(as.matrix(prepared[[pretreatment]][, BANDS]) - as.matrix(as.data.frame(read_spectra(pretreatment))[, BANDS])))
  check(sprintf("%s table is rebuilt from the mean spectra", pretreatment), difference < 1e-12, sprintf("(largest difference %.1e)", difference))
}

severity <- read_result(SEVERITY_FILE)
range_of <- function(times, clone_value, inoculated) {
  range(severity$severity[severity$time %in% times & severity$clone == clone_value & severity$inoculation == inoculated])
}
check("Severity of the susceptible clone before symptoms is 0.70-5.63%", all(round(range_of(PRE_SYMPTOM_TIMES, 223, "yes"), 2) == c(0.70, 5.63)))
check("Severity of the resistant clone before symptoms is at most 2.36%", round(range_of(PRE_SYMPTOM_TIMES, 236, "yes")[2], 2) == 2.36)
check("Severity after symptoms reaches 84.42% and 3.57%",
      round(range_of(POST_SYMPTOM_TIMES, 223, "yes")[2], 2) == 84.42 && round(range_of(POST_SYMPTOM_TIMES, 236, "yes")[2], 2) == 3.57)

candidates <- read_result(CANDIDATES_FILE)
final_sets <- read_result(FINAL_SETS_FILE)
refits <- read_result(REFITS_FILE)

# ---- Tables and figures ----
performance <- performance_table(read_result(VALIDATION_METRICS_FILE), refits, read_result(PLSDA_CV_FILE))
published <- read_result(file.path("tests", "expected", "published_tables.csv"))
compared <- published %>% inner_join(performance, by = c("method", "leaf", "pretreatment", "learner"), suffix = c("_article", "_here"))
names_compared <- c("precision", "precision_sd", "recall", "recall_sd", "bacc", "bacc_sd", "fbeta", "fbeta_sd")
differences <- sapply(names_compared, function(name) abs(compared[[paste0(name, "_article")]] - compared[[paste0(name, "_here")]]))
check(sprintf("The %d values of Tables 1, A1 and A2 match the article to three decimals", length(differences)),
      nrow(compared) == nrow(published) && max(differences) < 0.00051, sprintf("(largest difference %.4f)", max(differences)))

calibration <- read_result(CALIBRATION_PREDICTIONS_FILE)
validation <- read_result(VALIDATION_PREDICTIONS_FILE)
figure_values <- figure_5_values(calibration, validation, as_published = TRUE)
expected_figure <- read_result(file.path("tests", "expected", "figure_5_published.csv")) %>%
  inner_join(figure_values, by = c("model", "clone", "time"), suffix = c("_article", "_here"))
check("The 24 bars of Figure 5 match the article", nrow(expected_figure) == 24 && max(abs(expected_figure$bacc_article - expected_figure$bacc_here)) < 1e-8)

pre_classes <- severity_classes(severity %>% filter(time %in% PRE_SYMPTOM_TIMES))
post_classes <- severity_classes(severity %>% filter(time %in% POST_SYMPTOM_TIMES))
check("The legends of Figure A1 match the article",
      identical(pre_classes$labels, c("58.44%", "79.22%", "92.21%", "97.40%", "100.00%")) &&
        identical(post_classes$labels, c("77.50%", "82.50%", "87.50%", "95.00%", "100.00%")))

plsda_cv <- read_result(PLSDA_CV_FILE)
best_components <- sapply(PRETREATMENTS, function(pretreatment) {
  bounds <- plsda_cv %>%
    filter(pretreatment == !!pretreatment) %>%
    group_by(ncomp) %>%
    summarise(bound = mean(fbeta) - sd(fbeta))
  bounds$ncomp[which.max(bounds$bound)]
})
check("The PLS-DA components of Table 1 have the highest lower bound of the F1-score", identical(best_components, PLSDA_COMPONENTS))

gsets <- final_sets %>% filter(method == "genetic_algorithm", leaf == "apical", target_size == 75)
check("Final GA sets have 74, 75 and 65 bands, as described in the README",
      identical(as.integer(gsets$n_bands[match(PRETREATMENTS, gsets$pretreatment)]), c(74L, 75L, 65L)))

if (run_models) {
  compare_values <- function(new, old) max(abs(as.numeric(unlist(new)) - as.numeric(unlist(old))))

  # One fold of RFE at all target sizes
  run <- list(method = "rfe", leaf = "apical", pretreatment = "SNV")
  rf_rows <- run_feature_selection(run, repetition = 1, folds = 1, sizes = RUN_SIZES)
  for (size in TARGET_SIZES) {
    new <- rf_rows %>% filter(target_size == size)
    old <- candidates %>% filter(method == "rfe", leaf == "apical", pretreatment == "SNV", learner == "ranger", target_size == size, repetition == 1, fold == 1)
    check(sprintf("RFE, %d features: same bands and metrics as the original run", size),
          identical(new$bands, old$bands) && compare_values(new[, METRICS], old[, METRICS]) < TOLERANCE)
  }

  # The same fold with SVM
  svm_rows <- run_svm_candidates(run, repetition = 1, candidates = candidates, folds = 1, sizes = 75)
  old <- candidates %>% filter(method == "rfe", leaf == "apical", pretreatment == "SNV", learner == "svm", target_size == 75, repetition == 1, fold == 1)
  check("SVM on the RFE set: same metrics as the original run", compare_values(svm_rows[, METRICS], old[, METRICS]) < TOLERANCE)

  # The 25 refits of the final set, with the seeds of its position among the candidates
  for (spec in list(list(method = "rfe", pretreatment = "SNV"), list(method = "genetic_algorithm", pretreatment = "reflect"))) {
    run <- list(method = spec$method, leaf = "apical", pretreatment = spec$pretreatment)
    bands <- final_sets %>% filter(method == spec$method, leaf == "apical", pretreatment == spec$pretreatment, target_size == 75) %>% pull(bands) %>% split_bands()
    position <- which(sapply(candidate_sets(candidates, run, 75), function(x) identical(x, bands)))
    best <- choose_final_set(run, 75, candidates, only_candidate = position)
    old <- refits %>% filter(method == spec$method, leaf == "apical", pretreatment == spec$pretreatment, target_size == 75, learner == "ranger")
    check(sprintf("%s %s: 25 refits of the final set as in the original run", spec$method, spec$pretreatment),
          compare_values(best$performance[, METRICS], old[, METRICS]) < TOLERANCE)
  }

  # Final models for the validation plants
  metrics <- read_result(VALIDATION_METRICS_FILE)
  for (method in METHODS) {
    bands <- final_sets %>% filter(method == !!method, leaf == "apical", pretreatment == "SNV", target_size == 75) %>% pull(bands) %>% split_bands()
    result <- run_final_models(list(method = method, leaf = "apical", pretreatment = "SNV"), bands)
    old <- metrics %>% filter(method == !!method, leaf == "apical", pretreatment == "SNV")
    check(sprintf("%s SNV: final RF and SVM models as in the original run", method),
          compare_values(result$metrics[, METRICS], old[match(result$metrics$learner, old$learner), METRICS]) < TOLERANCE)
    old_predictions <- validation %>% filter(method == !!method, leaf == "apical", pretreatment == "SNV")
    check(sprintf("%s SNV: predictions for the validation plants as in the original run", method),
          identical(result$validation$response, old_predictions$response))
  }

  # PLS-DA
  components <- bind_rows(lapply(seq_len(N_REPETITIONS), function(r) run_plsda_components("SNV", r)))
  old <- read_result(PLSDA_CV_FILE) %>% filter(pretreatment == "SNV")
  check("PLS-DA SNV: components as in the original run", compare_values(components[, METRICS], old[, METRICS]) < TOLERANCE)
  final <- run_plsda_final("SNV")
  check("PLS-DA SNV: final model as in the original run",
        compare_values(final$metrics[, METRICS], metrics %>% filter(learner == "plsda", pretreatment == "SNV") %>% select(all_of(METRICS))) < TOLERANCE)
}

cat(if (failures == 0) "\nAll checks passed\n" else sprintf("\n%d checks failed\n", failures))
quit(status = if (failures == 0) 0 else 1, save = "no")
