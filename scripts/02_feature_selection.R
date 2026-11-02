# Section 2.4: feature selection with Random Forest, and the same feature sets with SVM.
# One run takes one combination:
#
#   Rscript scripts/02_feature_selection.R --method rfe --leaf apical --pretreatment SNV --repetition 1
#
# Optional: --folds 1,2 --sizes 75,50 --cores 12 (genetic algorithm) --learner ranger|svm (both by default)
source("R/setup.R")

args <- parse_args(commandArgs(TRUE))
run <- list(method = args$method, leaf = args$leaf, pretreatment = args$pretreatment)
repetition <- as.integer(args$repetition)
numbers <- function(text, default) if (is.null(text)) default else as.numeric(strsplit(text, ",")[[1]])
folds <- numbers(args$folds, seq_len(N_FOLDS))
learners <- if (is.null(args$learner)) LEARNERS else args$learner
output <- file.path(OUTPUTS_DIR, "candidates", sprintf("%s_%s_%s_%d", run$method, run$leaf, run$pretreatment, repetition))

candidates <- NULL
if ("ranger" %in% learners) {
  sizes <- numbers(args$sizes, RUN_SIZES)
  candidates <- run_feature_selection(run, repetition, folds, sizes, cores = as.integer(if (is.null(args$cores)) 12 else args$cores))
  write_table(candidates, paste0(output, "_ranger.csv"))
}
if ("svm" %in% learners) {
  source_table <- if (is.null(candidates)) read_result(CANDIDATES_FILE) else candidates
  sizes <- numbers(args$sizes, TARGET_SIZES)
  write_table(run_svm_candidates(run, repetition, source_table, folds, sizes), paste0(output, "_svm.csv"))
}
cat("Tables written to", dirname(output), "\n")
