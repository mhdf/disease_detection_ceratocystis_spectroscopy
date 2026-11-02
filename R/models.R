# Model fitting, prediction and evaluation shared by the steps of the analysis.

f1_summary <- function(data, lev = NULL, model = NULL) {
  precision <- caret::posPredValue(data$pred, data$obs, positive = "yes")
  recall <- caret::sensitivity(data$pred, data$obs, positive = "yes")
  if (is.na(precision) | is.na(recall) | (precision + recall) == 0) {
    f1 <- 0
  } else {
    f1 <- 2 * (precision * recall) / (precision + recall)
  }
  c(F1 = f1, Accuracy = mean(data$pred == data$obs))
}

train_control <- function() {
  trainControl(method = "none", summaryFunction = f1_summary, classProbs = TRUE,
               savePredictions = "final", allowParallel = FALSE)
}

# Radial SVM with more iterations than the caret default
CUSTOM_SVM <- getModelInfo("svmRadial", regex = FALSE)[[1]]
CUSTOM_SVM$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
  kernlab::ksvm(
    x = as.matrix(x), y = y, type = "C-svc", kernel = "rbfdot", kpar = list(sigma = param$sigma),
    C = param$C, prob.model = classProbs, scaled = FALSE, maxiter = 10000, tol = 0.01, ...
  )
}

as_outcome <- function(x) factor(x, levels = c("no", "yes"))
as_predictors <- function(data, bands) as.data.frame(data[, bands, drop = FALSE])

predict_labels <- function(model, newdata) {
  probabilities <- predict(model, newdata = newdata, type = "prob")
  factor(ifelse(probabilities$yes > 0.5, "yes", "no"), levels = c("no", "yes"))
}

# Precision, recall, balanced accuracy and F1-score, in this order
performance_values <- function(truth, response) {
  c(
    mlr3measures::precision(truth = truth, response = response, positive = "yes"),
    mlr3measures::recall(truth = truth, response = response, positive = "yes"),
    mlr3measures::bacc(truth = truth, response = response, positive = "yes"),
    mlr3measures::fbeta(truth = truth, response = response, beta = 1, na_value = 0, positive = "yes")
  )
}

f1_or_zero <- function(truth, response) {
  score <- F1_Score(y_true = truth, y_pred = response, positive = "yes")
  if (is.na(score)) 0 else score
}

# ---- Random Forest ----

fit_rf <- function(data, bands, parameters) {
  train(
    x = as_predictors(data, bands),
    y = as_outcome(data$inoculation),
    method = "ranger", metric = "F1", trControl = train_control(),
    tuneGrid = data.frame(mtry = parameters$mtry, splitrule = "gini", min.node.size = parameters$min.node.size),
    num.trees = parameters$num.trees
  )
}

rf_grid <- function(n_bands) {
  mtry_grid <- floor(n_bands * MTRY_RATIO_GRID) %>% unique(.) %>% .[. > 0]
  if (n_bands <= 20) {
    mtry_grid <- c(mtry_grid, floor(sqrt(n_bands)), n_bands) %>% unique(.) %>% .[. > 0 & . <= n_bands]
  }
  grid <- as.data.frame(expand.grid(mtry = mtry_grid, num.trees = NUM_TREES_GRID,
                                    min.node.size = MIN_NODE_SIZE_GRID, stringsAsFactors = FALSE))
  grid$measure <- NA
  grid
}

# F1-score of every combination of the grid; the random seed is set before calling it
tune_rf <- function(fit_data, eval_data, bands) {
  grid <- rf_grid(length(bands))
  for (kk in seq_len(nrow(grid))) {
    model <- fit_rf(fit_data, bands, grid[kk, ])
    grid[kk, "measure"] <- f1_or_zero(eval_data$inoculation, predict_labels(model, as_predictors(eval_data, bands)))
  }
  grid
}

# Ties go to the most complex model: more trees, then larger mtry, then larger node size
best_rf <- function(grid) {
  best <- grid[which(grid$measure == max(grid$measure)), ]
  element <- 1
  while (nrow(best) > 1) {
    best <- best[best[, RF_PARAMETERS[element]] == max(best[, RF_PARAMETERS[element]]), , drop = FALSE]
    element <- element + 1
  }
  best
}

# ---- Support Vector Machine ----

fit_svm <- function(data, bands, parameters, scaler) {
  train(
    x = predict(scaler, as_predictors(data, bands)),
    y = as_outcome(data$inoculation),
    method = CUSTOM_SVM, metric = "F1", trControl = train_control(),
    tuneGrid = data.frame(C = parameters$cost, sigma = parameters$gamma)
  )
}

scaler_for <- function(data, bands) preProcess(as_predictors(data, bands), method = c("center", "scale"))

tune_svm <- function(fit_data, eval_data, bands, scaler) {
  grid <- as.data.frame(expand.grid(kernel = "svmRadial", cost = COST_GRID, gamma = GAMMA_GRID, stringsAsFactors = FALSE))
  grid$measure <- NA
  for (kk in seq_len(nrow(grid))) {
    model <- fit_svm(fit_data, bands, grid[kk, ], scaler)
    response <- predict_labels(model, predict(scaler, as_predictors(eval_data, bands)))
    grid[kk, "measure"] <- f1_or_zero(eval_data$inoculation, response)
  }
  grid
}

# Ties go to the largest cost, then the largest gamma
best_svm <- function(grid) {
  best <- grid[which(grid$measure == max(grid$measure)), ]
  element <- 1
  parameters <- c("cost", "gamma")
  while (nrow(best) > 1) {
    best <- best[best[, parameters[element]] == max(best[, parameters[element]]), , drop = FALSE]
    element <- element + 1
  }
  best
}
