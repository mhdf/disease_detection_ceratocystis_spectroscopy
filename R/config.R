# Paths, settings of the article and fixed random seeds.

DATA_DIR <- "data"
RESULTS_DIR <- "results"
OUTPUTS_DIR <- "outputs"
FIGURES_DIR <- "figures"

RAW_MEANS_FILE <- file.path(DATA_DIR, "spectra_mean_1nm.csv.xz")
SEVERITY_FILE <- file.path(DATA_DIR, "severity.csv")
SPLIT_FILE <- file.path(DATA_DIR, "sample_split.csv")
spectra_path <- function(pretreatment, folder = DATA_DIR) {
  file.path(folder, sprintf("spectra_%s_5nm.csv", tolower(pretreatment)))
}

CANDIDATES_FILE <- file.path(RESULTS_DIR, "cv_candidate_feature_sets.csv")
FINAL_SETS_FILE <- file.path(RESULTS_DIR, "final_feature_sets.csv")
REFITS_FILE <- file.path(RESULTS_DIR, "cv_final_feature_sets.csv")
PLSDA_CV_FILE <- file.path(RESULTS_DIR, "cv_plsda_components.csv")
VALIDATION_METRICS_FILE <- file.path(RESULTS_DIR, "validation_metrics.csv")
VALIDATION_PREDICTIONS_FILE <- file.path(RESULTS_DIR, "validation_predictions.csv")
CALIBRATION_PREDICTIONS_FILE <- file.path(RESULTS_DIR, "calibration_predictions.csv")
SHAP_FILE <- file.path(RESULTS_DIR, "shap_importance.csv")

# Evaluation times: t1 to t4 are 12, 24, 48 and 60 hours after inoculation, t5 and t6 are 30 and 40 days
PRE_SYMPTOM_TIMES <- c("t1", "t2", "t3", "t4")
POST_SYMPTOM_TIMES <- c("t5", "t6")
TIME_LABELS <- c(t1 = "12 hours", t2 = "24 hours", t3 = "48 hours", t4 = "60 hours")
CLONE_LABELS <- c("223" = "Susceptible", "236" = "Resistant")

PRETREATMENTS <- c("reflect", "SNV", "CR")
PRETREATMENT_LABELS <- c(reflect = "Reflectance", SNV = "SNV", CR = "CR")
LEAVES <- c("apical", "middle", "basal")
METHODS <- c("rfe", "genetic_algorithm")
METHOD_LABELS <- c(rfe = "RFE", genetic_algorithm = "GA")
LEARNERS <- c("ranger", "svm")

# Spectra
RAW_WAVELENGTHS <- seq(350, 2500, 1)
WAVELENGTHS <- seq(400, 2450, 5)
BANDS <- paste0("C", WAVELENGTHS)
RESAMPLING_FWHM <- c(rep(3, length(seq(400, 700, 5))), rep(10, length(seq(705, 2450, 5))))

# Cross-validation
N_FOLDS <- 5
N_REPETITIONS <- 5
INNER_FOLDS <- 5
FOLD_SEEDS <- c(2458, 5768, 44726, 235688, 745)

# Numbers of features evaluated in each run. Only TARGET_SIZES were considered in the article, but the
# whole sequence is kept because it defines the RFE steps and the random seeds of the original run.
RUN_SIZES <- c(411, 200, 100, 75, 50, 25, 10, 5)
TARGET_SIZES <- c(25, 50, 75, 100)
ARTICLE_SIZE <- 75

# Hyperparameter grids
NUM_TREES_GRID <- c(20, 50, 100, 200)
MTRY_RATIO_GRID <- c(0.05, 0.10, round(seq(0.2, 0.5, length.out = 3), 2))
MIN_NODE_SIZE_GRID <- c(1, 2, 5)
RF_PARAMETERS <- c("num.trees", "mtry", "min.node.size")
COST_GRID <- round(sort(c(10^seq(-4, 3, length.out = 7), 200, 500)), 4)
GAMMA_GRID <- round(sort(c(10^seq(-4, -1, by = 0.25), 0.023, 0.042)), 4)
PLSDA_MAX_COMPONENTS <- 15
PLSDA_COMPONENTS <- c(reflect = 7, SNV = 12, CR = 6)

# Genetic algorithm
GA_POPULATION <- 50
GA_ITERATIONS <- 100

# Models shown in Figure 5
FIGURE_5_MODELS <- data.frame(
  learner = c("plsda", "ranger", "svm"),
  method = c("none", "rfe", "genetic_algorithm"),
  pretreatment = c("SNV", "SNV", "reflect"),
  label = c("PLS-DA", "RF+RFE", "SVM+GA"),
  stringsAsFactors = FALSE
)

# Simple reader of arguments given as --name value
parse_args <- function(arguments) {
  names_index <- grep("^--", arguments)
  values <- as.list(arguments[names_index + 1])
  names(values) <- sub("^--", "", arguments[names_index])
  values
}
