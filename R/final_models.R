# Final Random Forest and SVM models: hyperparameters from the folds of the calibration set,
# out-of-fold predictions for Figure 5, and predictions for the validation plants.

final_model_seeds <- function(pretreatment, n_folds_total) {
  index <- match(pretreatment, PRETREATMENTS)
  set.seed(25473)
  seed_all <- sample(0:1000000, length(PRETREATMENTS), replace = FALSE)
  set.seed(seed_all[index])
  seed_model <- lapply(seq_len(n_folds_total), function(jj) sample(0:1000000, 2, replace = FALSE))
  seed_final <- sample(0:1000000, 2, replace = FALSE)
  set.seed(6536)
  seed_all <- sample(0:1000000, length(PRETREATMENTS), replace = FALSE)
  set.seed(seed_all[index])
  seed_calibration <- lapply(seq_len(n_folds_total), function(jj) sample(0:1000000, 2, replace = FALSE))
  list(model = seed_model, final = seed_final, calibration = seed_calibration)
}

response_name <- function(ii) sprintf("response_%02d", ii)

fit_final_learner <- function(data, ids, bands, seeds, learner) {
  position <- match(learner, LEARNERS)
  calibration_predictions <- data$meta_calibration[, c("plant", "clone", "inoculation", "time")]
  best_measure <- 0
  selected <- NULL
  scaler_all <- if (learner == "svm") scaler_for(data$data_calib, bands) else NULL

  for (ii in seq_along(ids$train)) {
    train_data <- data$data_calib[ids$train[[ii]], ]
    test_data <- data$data_calib[ids$test[[ii]], ]
    test_meta <- data$meta_calibration[ids$test[[ii]], ]
    scaler <- if (learner == "svm") scaler_for(train_data, bands) else NULL

    set.seed(seeds$model[[ii]][position])
    best <- if (learner == "ranger") {
      best_rf(tune_rf(train_data, test_data, bands))
    } else {
      best_svm(tune_svm(train_data, test_data, bands, scaler))
    }
    if (best$measure > best_measure) {
      selected <- best
      best_measure <- selected$measure
    }

    set.seed(seeds$calibration[[ii]][position])
    response <- if (learner == "ranger") {
      predict_labels(fit_rf(train_data, bands, selected), as_predictors(test_data, bands))
    } else {
      predict_labels(fit_svm(train_data, bands, selected, scaler), predict(scaler, as_predictors(test_data, bands)))
    }
    rows <- match(test_meta$plant, calibration_predictions$plant)
    calibration_predictions[[response_name(ii)]] <- NA_character_
    calibration_predictions[[response_name(ii)]][rows] <- as.character(response)
  }

  set.seed(seeds$final[position])
  response <- if (learner == "ranger") {
    predict_labels(fit_rf(data$data_calib, bands, selected), as_predictors(data$data_valid, bands))
  } else {
    predict_labels(fit_svm(data$data_calib, bands, selected, scaler_all), predict(scaler_all, as_predictors(data$data_valid, bands)))
  }
  list(selected = selected, response = response, calibration_predictions = calibration_predictions,
       values = performance_values(as_outcome(data$data_valid$inoculation), response))
}

run_final_models <- function(run, bands) {
  data <- prepare_fit_data(run$pretreatment, run$leaf)
  ids <- all_fold_ids(make_outer_folds(data$strata))
  seeds <- final_model_seeds(run$pretreatment, length(ids$train))
  metrics <- list(); validation <- list(); calibration <- list()
  for (learner in LEARNERS) {
    result <- fit_final_learner(data, ids, bands, seeds, learner)
    id <- data.frame(method = run$method, leaf = run$leaf, pretreatment = run$pretreatment, learner = learner, stringsAsFactors = FALSE)
    metrics[[learner]] <- cbind(id, data.frame(
      target_size = ARTICLE_SIZE,
      num_trees = if (learner == "ranger") result$selected$num.trees else NA,
      mtry = if (learner == "ranger") result$selected$mtry else NA,
      min_node_size = if (learner == "ranger") result$selected$min.node.size else NA,
      cost = if (learner == "svm") result$selected$cost else NA,
      gamma = if (learner == "svm") result$selected$gamma else NA,
      ncomp = NA,
      precision = result$values[1], recall = result$values[2], bacc = result$values[3], fbeta = result$values[4]
    ))
    validation[[learner]] <- cbind(id, as.data.frame(data$meta_validation[, c("plant", "clone", "inoculation", "time")]),
                                   response = as.character(result$response))
    calibration[[learner]] <- cbind(id, as.data.frame(result$calibration_predictions))
  }
  list(metrics = bind_rows(metrics), validation = bind_rows(validation), calibration = bind_rows(calibration))
}
