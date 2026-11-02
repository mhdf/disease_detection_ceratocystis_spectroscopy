# PLS-DA with all wavelengths: number of components and final models.

plsda_fit <- function(data, ncomp, scaler) {
  train(
    x = predict(scaler, as.data.frame(select(data, -inoculation))),
    y = as_outcome(data$inoculation),
    method = "pls", metric = "F1", trControl = train_control(), tuneGrid = data.frame(ncomp = ncomp)
  )
}

plsda_scaler <- function(data) preProcess(as.data.frame(select(data, -inoculation)), method = c("center", "scale"))

plsda_predict <- function(model, data, scaler) {
  predict_labels(model, predict(scaler, as.data.frame(select(data, -inoculation))))
}

# Performance of each number of components in the folds of one repetition
run_plsda_components <- function(pretreatment, repetition, folds = seq_len(N_FOLDS)) {
  data <- prepare_fit_data(pretreatment, "apical")
  ids <- repetition_ids(make_outer_folds(data$strata), repetition)
  set.seed(24886)
  seed_sequence <- sample(0:1000000, N_REPETITIONS, replace = FALSE)
  set.seed(seed_sequence[repetition])
  seeds <- lapply(seq_len(N_REPETITIONS), function(ii) {
    replicate(N_FOLDS, sample(0:1000000, 2, replace = FALSE), simplify = FALSE)
  })[[repetition]]

  rows <- list()
  for (fold in folds) {
    train_data <- data$data_calib[ids$train[[fold]], ]
    test_data <- data$data_calib[ids$test[[fold]], ]
    scaler <- plsda_scaler(train_data)
    set.seed(seeds[[fold]][1])
    for (ncomp in seq_len(PLSDA_MAX_COMPONENTS)) {
      model <- plsda_fit(train_data, ncomp, scaler)
      values <- performance_values(as_outcome(test_data$inoculation), plsda_predict(model, test_data, scaler))
      rows[[length(rows) + 1]] <- data.frame(
        pretreatment = pretreatment, repetition = repetition, fold = fold, ncomp = ncomp,
        precision = values[1], recall = values[2], bacc = values[3], fbeta = values[4], stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(rows)
}

plsda_seeds <- function() {
  set.seed(2358)
  final <- sample(0:1000000, length(PRETREATMENTS), replace = FALSE)
  set.seed(2458)
  shap <- sample(0:1000000, length(PRETREATMENTS), replace = FALSE)
  set.seed(5698)
  calibration <- sample(0:1000000, length(PRETREATMENTS), replace = FALSE)
  list(final = final, shap = shap, calibration = calibration)
}

# Final model of a pre-treatment, with the number of components given in the article
run_plsda_final <- function(pretreatment) {
  data <- prepare_fit_data(pretreatment, "apical")
  ids <- all_fold_ids(make_outer_folds(data$strata))
  ncomp <- PLSDA_COMPONENTS[[pretreatment]]
  index <- match(pretreatment, PRETREATMENTS)
  seeds <- plsda_seeds()
  set.seed(seeds$calibration[index])
  fold_seeds <- lapply(seq_along(ids$train), function(jj) sample(0:1000000, 1, replace = FALSE))

  calibration_predictions <- as.data.frame(data$meta_calibration[, c("plant", "clone", "inoculation", "time")])
  for (ii in seq_along(ids$train)) {
    train_data <- data$data_calib[ids$train[[ii]], ]
    test_data <- data$data_calib[ids$test[[ii]], ]
    scaler <- plsda_scaler(train_data)
    set.seed(fold_seeds[[ii]][1])
    model <- plsda_fit(train_data, ncomp, scaler)
    rows <- match(data$meta_calibration$plant[ids$test[[ii]]], calibration_predictions$plant)
    calibration_predictions[[response_name(ii)]] <- NA_character_
    calibration_predictions[[response_name(ii)]][rows] <- as.character(plsda_predict(model, test_data, scaler))
  }

  scaler <- plsda_scaler(data$data_calib)
  set.seed(seeds$final[index])
  model <- plsda_fit(data$data_calib, ncomp, scaler)
  response <- plsda_predict(model, data$data_valid, scaler)
  values <- performance_values(as_outcome(data$data_valid$inoculation), response)
  id <- data.frame(method = "none", leaf = "apical", pretreatment = pretreatment, learner = "plsda", stringsAsFactors = FALSE)
  list(
    model = model, scaler = scaler, data = data,
    metrics = cbind(id, data.frame(target_size = NA, num_trees = NA, mtry = NA, min_node_size = NA, cost = NA, gamma = NA,
                                   ncomp = ncomp, precision = values[1], recall = values[2], bacc = values[3], fbeta = values[4])),
    validation = cbind(id, as.data.frame(data$meta_validation[, c("plant", "clone", "inoculation", "time")]), response = as.character(response)),
    calibration = cbind(id, calibration_predictions)
  )
}
