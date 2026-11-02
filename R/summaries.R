# Tables and values of the figures, from the result tables.

METRICS <- c("precision", "recall", "bacc", "fbeta")

# Validation value and standard deviation across the 25 refits of the final feature set
performance_table <- function(metrics, refits, plsda_cv) {
  sd_rf_svm <- refits %>%
    filter(target_size == ARTICLE_SIZE) %>%
    group_by(method, leaf, pretreatment, learner) %>%
    summarise(across(all_of(METRICS), ~ sd(.x, na.rm = TRUE), .names = "{.col}_sd"), .groups = "drop")
  components <- data.frame(pretreatment = names(PLSDA_COMPONENTS), ncomp = unname(PLSDA_COMPONENTS), stringsAsFactors = FALSE)
  sd_plsda <- plsda_cv %>%
    inner_join(components, by = c("pretreatment", "ncomp")) %>%
    group_by(pretreatment) %>%
    summarise(across(all_of(METRICS), ~ sd(.x, na.rm = TRUE), .names = "{.col}_sd"), .groups = "drop") %>%
    mutate(method = "none", leaf = "apical", learner = "plsda")
  metrics %>%
    select(method, leaf, pretreatment, learner, all_of(METRICS)) %>%
    inner_join(bind_rows(sd_rf_svm, sd_plsda), by = c("method", "leaf", "pretreatment", "learner")) %>%
    select(method, leaf, pretreatment, learner, precision, precision_sd, recall, recall_sd, bacc, bacc_sd, fbeta, fbeta_sd)
}

order_rows <- function(x) {
  x %>%
    mutate(pretreatment = factor(pretreatment, levels = PRETREATMENTS),
           leaf = factor(leaf, levels = LEAVES),
           method = factor(method, levels = c("none", METHODS)),
           learner = factor(learner, levels = c("plsda", LEARNERS))) %>%
    arrange(method, learner, leaf, pretreatment) %>%
    mutate(across(c(pretreatment, leaf, method, learner), as.character))
}

table_1 <- function(performance) order_rows(filter(performance, leaf == "apical"))
table_a <- function(performance, selected_learner) order_rows(filter(performance, method == "rfe", learner == selected_learner))

# Mode of the five out-of-fold predictions of each plant in a repetition
modal_response <- function(values) {
  values <- values[!is.na(values)]
  unique_values <- unique(values)
  unique_values[which.max(tabulate(match(values, unique_values)))]
}

balanced_accuracy <- function(truth, response) {
  mlr3measures::bacc(truth = as_outcome(truth), response = as_outcome(response), positive = "yes")
}

calibration_modes <- function(calibration, responses) {
  values <- as.matrix(calibration[, responses, drop = FALSE])
  apply(values, 1, modal_response)
}

# Standard deviation across the five repetitions, using the calibration plants only
repetition_sd <- function(calibration) {
  responses <- grep("^response_", names(calibration), value = TRUE)
  repetitions <- split(responses, ceiling(seq_along(responses) / N_FOLDS))
  per_repetition <- lapply(repetitions, function(columns) {
    calibration %>%
      mutate(response = calibration_modes(calibration, columns)) %>%
      group_by(clone, time) %>%
      summarise(bacc = balanced_accuracy(inoculation, response), .groups = "drop")
  })
  bind_rows(per_repetition) %>%
    group_by(clone, time) %>%
    summarise(bacc_sd = sd(bacc, na.rm = TRUE), .groups = "drop")
}

# Balanced accuracy by clone and time: mode of all out-of-fold predictions for the calibration
# plants, together with the predictions for the validation plants
clone_time_accuracy <- function(calibration, validation) {
  responses <- grep("^response_", names(calibration), value = TRUE)
  all_plants <- bind_rows(
    calibration %>% mutate(response = calibration_modes(calibration, responses)) %>% select(plant, clone, inoculation, time, response),
    validation %>% select(plant, clone, inoculation, time, response)
  )
  all_plants %>%
    group_by(clone, time) %>%
    summarise(bacc = balanced_accuracy(inoculation, response), .groups = "drop")
}

# Values of Figure 5. In the article, the bars of a model use the pre-treatment in the position of the
# model in the list (reflectance for PLS-DA, SNV for RF+RFE and CR for SVM+GA), while the error bars use
# the intended pre-treatments (SNV, SNV and reflectance). Both versions are available here.
figure_5_values <- function(calibration, validation, as_published = TRUE) {
  rows <- list()
  for (i in seq_len(nrow(FIGURE_5_MODELS))) {
    chosen <- FIGURE_5_MODELS[i, ]
    bars_pretreatment <- if (as_published) PRETREATMENTS[i] else chosen$pretreatment
    select_model <- function(table, pretreatment_used) {
      table %>% filter(learner == chosen$learner, method == chosen$method, leaf == "apical", pretreatment == pretreatment_used)
    }
    accuracy <- clone_time_accuracy(select_model(calibration, bars_pretreatment), select_model(validation, bars_pretreatment))
    deviation <- repetition_sd(select_model(calibration, chosen$pretreatment))
    rows[[i]] <- accuracy %>%
      inner_join(deviation, by = c("clone", "time")) %>%
      mutate(model = chosen$label, bars_pretreatment = bars_pretreatment, error_bars_pretreatment = chosen$pretreatment)
  }
  bind_rows(rows)
}

# Severity classes of Figure A1: five intervals of equal width and the share of plants up to each one
severity_classes <- function(severity) {
  intervals <- cut(severity$severity, breaks = 5)
  limits <- unique(as.numeric(unlist(strsplit(gsub("[][()]", "", levels(intervals)), ","))))
  shares <- cumsum(table(intervals)) / length(intervals)
  list(limits = limits, class = as.integer(intervals), labels = sprintf("%.2f%%", 100 * shares))
}
