# Importance of the wavelengths with SHAP values, in percent of the total of each model.
# The values are estimated with random simulations in parallel, so they change slightly between runs.

relative_importance <- function(explanation) {
  mean_absolute <- apply(abs(explanation), 2, mean)
  100 * mean_absolute / sum(mean_absolute)
}

explain_in_parallel <- function(object, X, newdata, predict_function, cores, packages) {
  cluster <- makeCluster(cores)
  registerDoParallel(cluster)
  on.exit(stopCluster(cluster))
  clusterExport(cluster, "packages", envir = environment())
  clusterEvalQ(cluster, for (p in packages) library(p, character.only = TRUE))
  fastshap::explain(object = object, X = X, newdata = newdata, pred_wrapper = predict_function,
                    nsim = 100, parallel = TRUE)
}

shap_rf <- function(method, pretreatment, bands, parameters, cores) {
  data <- prepare_fit_data(pretreatment, "apical")
  training <- as.data.frame(data$data_calib[, c(bands, "inoculation")])
  training$inoculation <- as_outcome(training$inoculation)
  # The model explained continues the random numbers of the final model, as in the original run
  set.seed(final_model_seeds(pretreatment, N_FOLDS * N_REPETITIONS)$final[1])
  fit_rf(data$data_calib, bands, list(mtry = parameters$mtry, min.node.size = parameters$min_node_size, num.trees = parameters$num_trees))
  model <- ranger::ranger(
    formula = inoculation ~ ., data = training, num.trees = parameters$num_trees, mtry = parameters$mtry,
    min.node.size = parameters$min_node_size, splitrule = "gini", classification = TRUE, probability = TRUE,
    importance = "impurity"
  )
  # Seeds for reflectance, SNV and CR, alternating RFE and GA as in the original run
  set.seed(2458)
  seeds <- sample(0:1000000, length(METHODS) * length(PRETREATMENTS), replace = FALSE)
  labels <- as.vector(outer(METHODS, PRETREATMENTS, paste, sep = "_"))
  set.seed(seeds[match(paste0(method, "_", pretreatment), labels)])
  explanation <- explain_in_parallel(
    model, as_predictors(data$data_calib, bands), as_predictors(data$data_valid, bands),
    function(object, newdata) predict(object, data = newdata)$predictions[, "yes"], cores, "ranger"
  )
  importance <- relative_importance(explanation)
  data.frame(method = method, pretreatment = pretreatment, learner = "ranger",
             wavelength = as.numeric(sub("C", "", names(importance))), importance = unname(importance))
}

shap_plsda <- function(pretreatment, cores) {
  final <- run_plsda_final(pretreatment)
  seeds <- plsda_seeds()
  set.seed(seeds$shap[match(pretreatment, PRETREATMENTS)])
  explanation <- explain_in_parallel(
    final$model,
    predict(final$scaler, as.data.frame(select(final$data$data_calib, -inoculation))),
    predict(final$scaler, as.data.frame(select(final$data$data_valid, -inoculation))),
    function(object, newdata) predict(object, newdata = as.data.frame(newdata), type = "prob")[, "yes"],
    cores, "caret"
  )
  importance <- relative_importance(explanation)
  data.frame(method = "none", pretreatment = pretreatment, learner = "plsda",
             wavelength = as.numeric(sub("C", "", names(importance))), importance = unname(importance))
}
