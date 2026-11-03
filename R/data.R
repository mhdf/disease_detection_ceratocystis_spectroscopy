# Spectra, calibration and validation sets, and cross-validation folds.

# From the mean spectrum of each plant and canopy layer at 1 nm to the three tables used by the models
prepare_spectra <- function() {
  means <- read.csv(RAW_MEANS_FILE, stringsAsFactors = FALSE)
  metadata <- means[, c("plant", "clone", "inoculation", "leaf", "time")]
  raw <- as.matrix(means[, paste0("C", RAW_WAVELENGTHS)])

  reflectance <- prospectr::resample2(X = raw, wav = RAW_WAVELENGTHS, new.wav = WAVELENGTHS, fwhm = RESAMPLING_FWHM)
  colnames(reflectance) <- BANDS
  # Standard Normal Variate: each spectrum minus its own mean, divided by its own standard deviation
  snv <- t(scale(t(reflectance), center = TRUE, scale = TRUE))
  continuum <- prospectr::continuumRemoval(
    X = reflectance, wav = WAVELENGTHS, type = "R", interpol = "linear", method = "substraction"
  )
  colnames(continuum) <- BANDS

  list(
    reflect = cbind(metadata, reflectance),
    SNV = cbind(metadata, snv),
    CR = cbind(metadata, continuum)
  )
}

read_spectra <- function(pretreatment) {
  as_tibble(read.csv(spectra_path(pretreatment), stringsAsFactors = FALSE))
}

# Data used by every model of a given pre-treatment and canopy layer
prepare_fit_data <- function(pretreatment, leaf) {
  split_table <- read.csv(SPLIT_FILE, stringsAsFactors = FALSE)
  selected_leaf <- leaf
  data_to_fit <- read_spectra(pretreatment) %>%
    filter(leaf == selected_leaf, time %in% PRE_SYMPTOM_TIMES) %>%
    mutate(dataset = split_table$split[match(plant, split_table$plant)])
  calibration <- data_to_fit %>% filter(dataset == "calibration")
  validation <- data_to_fit %>% filter(dataset == "validation")
  list(
    strata = interaction(calibration$clone, calibration$time, calibration$inoculation, drop = TRUE),
    data_calib = calibration[, c("inoculation", BANDS)],
    data_valid = validation[, c("inoculation", BANDS)],
    meta_calibration = calibration %>% select(-all_of(BANDS)),
    meta_validation = validation %>% select(-all_of(BANDS))
  )
}

# Five repetitions of a five-fold split, stratified by clone, time and inoculation, with no repeated folds
make_outer_folds <- function(strata) {
  is_vector_in_list <- function(vector, list) any(sapply(list, function(x) identical(vector, x)))
  all_vectors_different <- function(list1, list2) all(sapply(list1, function(vector1) !is_vector_in_list(vector1, list2)))

  for (ii in seq_along(FOLD_SEEDS)) {
    if (ii == 1) {
      set.seed(FOLD_SEEDS[ii])
      folds <- caret::createFolds(strata, k = N_FOLDS, list = TRUE, returnTrain = FALSE)
    } else {
      seed_to_update <- FOLD_SEEDS[ii]
      current <- folds
      while (!all_vectors_different(current, folds)) {
        set.seed(seed_to_update)
        current <- caret::createFolds(strata, k = N_FOLDS, list = TRUE, returnTrain = FALSE)
        seed_to_update <- seed_to_update + 1
      }
      names(current) <- paste0("Fold", seq(length(folds) + 1, length(folds) + N_FOLDS, 1))
      folds <- c(folds, current)
    }
  }
  folds
}

# Training and test rows of each fold of a repetition. The folds are used from last to first, as in the original run
repetition_ids <- function(folds, repetition) {
  last <- repetition * N_FOLDS
  selected <- folds[paste0("Fold", seq(last - (N_FOLDS - 1), last, 1))]
  list(
    train = lapply(rev(seq_len(N_FOLDS)), function(k) sort(unname(unlist(selected[-k])))),
    test = lapply(rev(seq_len(N_FOLDS)), function(k) sort(unname(unlist(selected[k]))))
  )
}

all_fold_ids <- function(folds) {
  ids <- lapply(seq_len(N_REPETITIONS), function(repetition) repetition_ids(folds, repetition))
  list(train = unlist(lapply(ids, `[[`, "train"), recursive = FALSE),
       test = unlist(lapply(ids, `[[`, "test"), recursive = FALSE))
}
