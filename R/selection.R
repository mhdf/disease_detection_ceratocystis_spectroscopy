# Feature selection with Random Forest (RFE and genetic algorithm), and choice of the final feature set.

# ---- Random seeds of a repetition ----
run_seeds <- function(repetition) {
  set.seed(24886)
  seed_sequence <- sample(0:1000000, N_REPETITIONS, replace = FALSE)
  set.seed(seed_sequence[repetition])
  seed_outer <- lapply(seq_len(N_REPETITIONS), function(ii) {
    replicate(N_FOLDS, sample(0:1000000, 3, replace = FALSE), simplify = FALSE)
  })
  seed_ga <- lapply(seq_len(N_REPETITIONS), function(ii) {
    lapply(seq_len(N_FOLDS), function(jj) lapply(seq_along(RUN_SIZES), function(kk) sample(0:1000000, 6, replace = FALSE)))
  })
  seed_rfe <- lapply(seq_len(N_REPETITIONS), function(ii) {
    lapply(seq_len(N_FOLDS), function(jj) lapply(seq_along(RUN_SIZES), function(kk) sample(0:1000000, 4, replace = FALSE)))
  })
  list(outer = seed_outer[[repetition]], ga = seed_ga[[repetition]], rfe = seed_rfe[[repetition]])
}

# Seeds of the grid search and of the model: the third and fourth of the GA, the first and second of RFE
step_seeds <- function(seeds, method, fold, size_index) {
  if (method == "genetic_algorithm") seeds$ga[[fold]][[size_index]][3:4] else seeds$rfe[[fold]][[size_index]][1:2]
}

# ---- Genetic algorithm ----
initialize_population_k <- function(popSize, num_features, target_k) {
  population <- matrix(0, nrow = popSize, ncol = num_features)
  for (i in 1:popSize) {
    population[i, sample(num_features, target_k)] <- 1
  }
  population
}

# Fitness: F1-score and closeness to the target number of features, with equal weights
fitness_fn <- function(chromosome, target_k, fit_data, eval_data, train_ctrl, mtry_ratio_grid,
                       num_trees_grid, min_node_size_grid, ...) {
  nbr_selected_features <- sum(chromosome)
  nbr_features <- length(chromosome)
  if (nbr_selected_features == 0) {
    return(0)
  }
  selected_features <- names(fit_data)[which(chromosome == 1) + 1] # the first column is the outcome

  inner_X_fit <- fit_data[, selected_features, drop = FALSE] %>% as.data.frame()
  inner_y_fit <- fit_data$inoculation %>% factor(., levels = c("no", "yes"))
  inner_X_eval <- eval_data[, selected_features, drop = FALSE] %>% as.data.frame()
  inner_y_eval <- eval_data$inoculation %>% factor(., levels = c("no", "yes"))

  mtry_grid <- floor(nbr_selected_features * mtry_ratio_grid) %>% unique(.) %>% .[. > 0]
  if (nbr_selected_features <= 20) {
    mtry_grid <- c(mtry_grid, floor(sqrt(nbr_selected_features)), nbr_selected_features) %>%
      unique(.) %>% .[. > 0 & . <= nbr_selected_features]
  }
  results_grid <- expand.grid(mtry = mtry_grid, num.trees = num_trees_grid,
                              min.node.size = min_node_size_grid, stringsAsFactors = FALSE) %>% as.data.frame()
  results_grid$measure <- NA

  for (elmnt in 1:nrow(results_grid)) {
    model <- train(
      x = inner_X_fit, y = inner_y_fit, method = "ranger", metric = "F1", trControl = train_ctrl,
      tuneGrid = data.frame(mtry = results_grid$mtry[elmnt], splitrule = "gini",
                            min.node.size = results_grid$min.node.size[elmnt]),
      num.trees = results_grid$num.trees[elmnt]
    )
    pred_probs <- predict(model, newdata = inner_X_eval, type = "prob")
    pred_labels <- factor(ifelse(pred_probs$yes > 0.5, "yes", "no"), levels = c("no", "yes"))
    f1_score <- F1_Score(y_true = inner_y_eval, y_pred = pred_labels, positive = "yes")
    if (is.na(f1_score)) {
      f1_score <- 0
    }
    results_grid[elmnt, "measure"] <- max(f1_score)
  }

  max_dev <- max((target_k - 1), (nbr_features - target_k))
  feature_score <- 1 - (abs(nbr_selected_features - target_k) / max_dev)
  0.5 * max(results_grid$measure, na.rm = TRUE) + 0.5 * feature_score
}

select_with_ga <- function(fit_data, eval_data, target_k, size_seeds, cores) {
  if (target_k == length(BANDS)) {
    return(BANDS)
  }
  set.seed(size_seeds[1])
  suggestions <- initialize_population_k(GA_POPULATION, ncol(fit_data) - 1, target_k)
  result <- ga(
    type = "binary", fitness = fitness_fn, nBits = ncol(fit_data) - 1, popSize = GA_POPULATION,
    pcrossover = 0.8, pmutation = 0.1, elitism = GA_POPULATION * 0.1, maxiter = GA_ITERATIONS,
    run = floor(0.2 * GA_ITERATIONS), maxFitness = 1, suggestions = suggestions, keepBest = TRUE,
    parallel = cores, seed = size_seeds[2],
    target_k = target_k, fit_data = fit_data, eval_data = eval_data, train_ctrl = train_control(),
    mtry_ratio_grid = MTRY_RATIO_GRID, num_trees_grid = NUM_TREES_GRID, min_node_size_grid = MIN_NODE_SIZE_GRID
  )
  BANDS[unname(which(result@solution[1, ] == 1))]
}

# ---- Recursive feature elimination, with the hyperparameters searched inside each model ----
custom_rfe_function <- caretFuncs
custom_rfe_function$fit <- function(x, y, first, last, mtry_ratio_grid, num_trees_grid, min_node_size_grid,
                                    inner_folds_nbr_hyper, seeds_rfe_inner_train, ...) {
  n_features <- ncol(x)
  mtry_grid <- floor(n_features * mtry_ratio_grid) %>% unique(.) %>% .[. > 0]
  if (n_features <= 20) {
    mtry_grid <- c(mtry_grid, floor(sqrt(n_features)), n_features) %>% unique(.) %>% .[. > 0 & . <= n_features]
  }
  results_grid <- expand.grid(mtry = mtry_grid, num.trees = num_trees_grid,
                              min.node.size = min_node_size_grid, stringsAsFactors = FALSE)
  results_grid$measure <- NA
  inner_train_ctrl <- trainControl(method = "cv", number = inner_folds_nbr_hyper, seeds = seeds_rfe_inner_train,
                                   summaryFunction = f1_summary, classProbs = TRUE,
                                   savePredictions = "all", allowParallel = FALSE)
  best_perf <- -Inf
  best_model <- NULL
  for (elmnt in seq(1, dim(results_grid)[1], 1)) {
    model <- train(
      x = x, y = y, method = "ranger", metric = "F1", trControl = inner_train_ctrl,
      tuneGrid = data.frame(mtry = results_grid$mtry[elmnt], splitrule = "gini",
                            min.node.size = results_grid$min.node.size[elmnt]),
      num.trees = results_grid$num.trees[elmnt], importance = "permutation"
    )
    perf <- max(model$results$F1, na.rm = TRUE)
    results_grid$measure[elmnt] <- perf
    if (best_perf < perf) {
      best_perf <- perf
      best_model <- model
    }
  }
  best_model
}
custom_rfe_function$summary <- f1_summary

# All the sizes come from one call: the bands are ranked once and the top ones are kept at each size
run_rfe <- function(train_data, fit_idx, seeds_outer_fold) {
  rfe_ctrl <- rfeControl(functions = custom_rfe_function, method = "LGOCV",
                         index = list(HouldoutTrain = as.vector(fit_idx)), verbose = FALSE, returnResamp = "all")
  n_combinations <- nrow(expand.grid(mtry = MTRY_RATIO_GRID, num.trees = NUM_TREES_GRID,
                                     min.node.size = MIN_NODE_SIZE_GRID, stringsAsFactors = FALSE))
  inner_seeds <- vector("list", INNER_FOLDS + 1)
  set.seed(seeds_outer_fold[2])
  for (element in seq_len(INNER_FOLDS)) {
    inner_seeds[[element]] <- sample(0:1000000, n_combinations, replace = FALSE)
  }
  inner_seeds[[INNER_FOLDS + 1]] <- sample(0:1000000, 1, replace = FALSE)

  set.seed(seeds_outer_fold[3])
  result <- caret::rfe(
    x = as.data.frame(select(train_data, -inoculation)),
    y = as_outcome(train_data$inoculation),
    sizes = RUN_SIZES, rfeControl = rfe_ctrl, metric = "F1",
    mtry_ratio_grid = MTRY_RATIO_GRID, num_trees_grid = NUM_TREES_GRID, min_node_size_grid = MIN_NODE_SIZE_GRID,
    inner_folds_nbr_hyper = INNER_FOLDS, seeds_rfe_inner_train = inner_seeds
  )
  result$variables
}

rfe_bands <- function(rfe_variables, size) {
  selected <- rfe_variables %>% filter(Variables == size)
  selected[order(as.numeric(sub("C", "", selected$var))), "var"]
}

candidate_row <- function(run, size, learner, repetition, fold, performance, best, bands) {
  data.frame(
    method = run$method, leaf = run$leaf, pretreatment = run$pretreatment, target_size = size, learner = learner,
    repetition = repetition, fold = fold,
    precision = performance[1], recall = performance[2], bacc = performance[3], fbeta = performance[4],
    num_trees = if (learner == "ranger") best$num.trees else NA,
    mtry = if (learner == "ranger") best$mtry else NA,
    min_node_size = if (learner == "ranger") best$min.node.size else NA,
    cost = if (learner == "svm") best$cost else NA,
    gamma = if (learner == "svm") best$gamma else NA,
    bands = if (learner == "ranger") paste(bands, collapse = ";") else NA,
    stringsAsFactors = FALSE
  )
}

# Train and test rows of a fold, with the 60/40 split used to compare hyperparameters
fold_data <- function(data, ids, fold, partition_seed) {
  train_data <- data$data_calib[ids$train[[fold]], ]
  set.seed(partition_seed)
  fit_idx <- createDataPartition(data$strata[ids$train[[fold]]], p = 0.60, list = FALSE)
  eval_idx <- setdiff(seq_len(nrow(train_data)), fit_idx)
  list(
    train_data = train_data,
    test_data = data$data_calib[ids$test[[fold]], ],
    fit_idx = fit_idx,
    fit_data = train_data[fit_idx, ],
    eval_data = train_data[eval_idx, ]
  )
}

# One repetition of the feature selection with Random Forest
run_feature_selection <- function(run, repetition, folds = seq_len(N_FOLDS), sizes = RUN_SIZES, cores = 12) {
  data <- prepare_fit_data(run$pretreatment, run$leaf)
  ids <- repetition_ids(make_outer_folds(data$strata), repetition)
  seeds <- run_seeds(repetition)
  rows <- list()
  for (fold in folds) {
    parts <- fold_data(data, ids, fold, seeds$outer[[fold]][1])
    rfe_variables <- if (run$method == "rfe") run_rfe(parts$train_data, parts$fit_idx, seeds$outer[[fold]]) else NULL
    for (size in sizes) {
      size_index <- match(size, RUN_SIZES)
      bands <- if (run$method == "rfe") {
        rfe_bands(rfe_variables, size)
      } else {
        select_with_ga(parts$fit_data, parts$eval_data, size, seeds$ga[[fold]][[size_index]], cores)
      }
      model_seeds <- step_seeds(seeds, run$method, fold, size_index)
      set.seed(model_seeds[1])
      best <- best_rf(tune_rf(parts$fit_data, parts$eval_data, bands))
      set.seed(model_seeds[2])
      model <- fit_rf(parts$train_data, bands, best)
      response <- predict_labels(model, as_predictors(parts$test_data, bands))
      rows[[length(rows) + 1]] <- candidate_row(
        run, size, "ranger", repetition, fold,
        performance_values(as_outcome(parts$test_data$inoculation), response), best, bands
      )
    }
    message(sprintf("%s %s %s: repetition %d, fold %d finished", run$method, run$leaf, run$pretreatment, repetition, fold))
  }
  bind_rows(rows)
}

# The same candidate feature sets, with SVM models
run_svm_candidates <- function(run, repetition, candidates, folds = seq_len(N_FOLDS), sizes = TARGET_SIZES) {
  data <- prepare_fit_data(run$pretreatment, run$leaf)
  ids <- repetition_ids(make_outer_folds(data$strata), repetition)
  seeds <- run_seeds(repetition)
  rows <- list()
  for (fold in folds) {
    parts <- fold_data(data, ids, fold, seeds$outer[[fold]][1])
    for (size in sizes) {
      size_index <- match(size, RUN_SIZES)
      bands <- candidate_bands(candidates, run, size, repetition, fold)
      scaler_fit <- scaler_for(parts$fit_data, bands)
      scaler_train <- scaler_for(parts$train_data, bands)
      model_seeds <- step_seeds(seeds, run$method, fold, size_index)
      set.seed(model_seeds[1])
      best <- best_svm(tune_svm(parts$fit_data, parts$eval_data, bands, scaler_fit))
      set.seed(model_seeds[2])
      model <- fit_svm(parts$train_data, bands, best, scaler_train)
      response <- predict_labels(model, predict(scaler_train, as_predictors(parts$test_data, bands)))
      rows[[length(rows) + 1]] <- candidate_row(
        run, size, "svm", repetition, fold,
        performance_values(as_outcome(parts$test_data$inoculation), response), best, bands
      )
    }
    message(sprintf("SVM %s %s %s: repetition %d, fold %d finished", run$method, run$leaf, run$pretreatment, repetition, fold))
  }
  bind_rows(rows)
}

split_bands <- function(text) strsplit(text, ";", fixed = TRUE)[[1]]

candidate_bands <- function(candidates, run, size, repetition, fold) {
  row <- candidates %>% filter(
    method == run$method, leaf == run$leaf, pretreatment == run$pretreatment, target_size == size,
    learner == "ranger", repetition == !!repetition, fold == !!fold
  )
  split_bands(row$bands)
}

# Candidate feature sets of a target size, from repetition 1 fold 1 to repetition 5 fold 5
candidate_sets <- function(candidates, run, size) {
  rows <- candidates %>%
    filter(method == run$method, leaf == run$leaf, pretreatment == run$pretreatment, target_size == size, learner == "ranger") %>%
    arrange(repetition, fold)
  lapply(rows$bands, split_bands)
}

# ---- Choice of the final feature set ----
choice_seeds <- function(pretreatment, n_candidates, n_folds_total) {
  set.seed(25473)
  seed_all <- sample(0:1000000, length(PRETREATMENTS), replace = FALSE)[match(pretreatment, PRETREATMENTS)]
  set.seed(seed_all)
  seed_fit_data <- lapply(seq_len(n_folds_total), function(ii) sample(0:1000000, 1, replace = FALSE))
  seed_model <- lapply(seq_len(n_candidates), function(ii) {
    lapply(seq_len(n_folds_total), function(jj) sample(0:1000000, 2, replace = FALSE))
  })
  list(fit_data = seed_fit_data, model = seed_model)
}

# Each candidate is refitted on the 25 folds of the calibration set, and the final set is the one
# with the highest lower bound (mean - SD) of the F1-score
refit_candidate <- function(data, ids, bands, candidate_seeds, fit_seeds, learner) {
  rows <- list()
  for (jj in seq_along(ids$train)) {
    train_data <- data$data_calib[ids$train[[jj]], ]
    test_data <- data$data_calib[ids$test[[jj]], ]
    set.seed(fit_seeds[[jj]][1])
    fit_idx <- createDataPartition(data$strata[ids$train[[jj]]], p = 0.60, list = FALSE)
    fit_data <- train_data[fit_idx, ]
    eval_data <- train_data[setdiff(seq_len(nrow(train_data)), fit_idx), ]
    if (learner == "ranger") {
      set.seed(candidate_seeds[[jj]][1])
      best <- best_rf(tune_rf(fit_data, eval_data, bands))
      set.seed(candidate_seeds[[jj]][2])
      model <- fit_rf(train_data, bands, best)
      response <- predict_labels(model, as_predictors(test_data, bands))
    } else {
      scaler_fit <- scaler_for(fit_data, bands)
      scaler_train <- scaler_for(train_data, bands)
      set.seed(candidate_seeds[[jj]][1])
      best <- best_svm(tune_svm(fit_data, eval_data, bands, scaler_fit))
      set.seed(candidate_seeds[[jj]][2])
      model <- fit_svm(train_data, bands, best, scaler_train)
      response <- predict_labels(model, predict(scaler_train, as_predictors(test_data, bands)))
    }
    values <- performance_values(as_outcome(test_data$inoculation), response)
    rows[[jj]] <- data.frame(resample = jj, precision = values[1], recall = values[2], bacc = values[3], fbeta = values[4])
  }
  bind_rows(rows)
}

choose_final_set <- function(run, size, candidates, only_candidate = NULL) {
  data <- prepare_fit_data(run$pretreatment, run$leaf)
  ids <- all_fold_ids(make_outer_folds(data$strata))
  sets <- candidate_sets(candidates, run, size)
  seeds <- choice_seeds(run$pretreatment, length(sets), length(ids$train))
  best_value <- 0
  best <- NULL
  for (ii in if (is.null(only_candidate)) seq_along(sets) else only_candidate) {
    performance <- refit_candidate(data, ids, sets[[ii]], seeds$model[[ii]], seeds$fit_data, "ranger")
    value <- mean(performance$fbeta, na.rm = TRUE) - sd(performance$fbeta, na.rm = TRUE)
    if (best_value < value) {
      best_value <- value
      best <- list(bands = sets[[ii]], performance = performance)
    }
    message(sprintf("%s %s %s, %d features: candidate %d of %d", run$method, run$leaf, run$pretreatment, size, ii, length(sets)))
  }
  best
}

# Performance of the chosen feature set with SVM
refit_final_svm <- function(run, size, candidates, bands) {
  data <- prepare_fit_data(run$pretreatment, run$leaf)
  ids <- all_fold_ids(make_outer_folds(data$strata))
  sets <- candidate_sets(candidates, run, size)
  seeds <- choice_seeds(run$pretreatment, length(sets), length(ids$train))
  ii <- which(sapply(sets, function(x) identical(x, bands)))
  refit_candidate(data, ids, bands, seeds$model[[ii]], seeds$fit_data, "svm")
}
