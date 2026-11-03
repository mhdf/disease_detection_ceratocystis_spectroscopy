# Section 2.4: choice of the final feature set of each target size, from the candidate sets.
#
#   Rscript scripts/03_final_feature_sets.R --method rfe --leaf apical --pretreatment SNV --size 75
#
# The candidates come from results/, or from outputs/candidates/ with --use-outputs yes.
source("R/setup.R")

args <- parse_args(commandArgs(TRUE))
run <- list(method = args$method, leaf = args$leaf, pretreatment = args$pretreatment)
size <- as.numeric(args$size)
candidates <- if (identical(args$`use-outputs`, "yes")) {
  bind_rows(lapply(list.files(file.path(OUTPUTS_DIR, "candidates"), full.names = TRUE), read_result))
} else {
  read_result(CANDIDATES_FILE)
}

best <- choose_final_set(run, size, candidates)
svm_performance <- refit_final_svm(run, size, candidates, best$bands)

id <- data.frame(method = run$method, leaf = run$leaf, pretreatment = run$pretreatment, target_size = size, stringsAsFactors = FALSE)
tag <- sprintf("%s_%s_%s_%d", run$method, run$leaf, run$pretreatment, size)
write_table(cbind(id, n_bands = length(best$bands), bands = paste(best$bands, collapse = ";")),
            file.path(OUTPUTS_DIR, "final_sets", paste0(tag, ".csv")))
write_table(bind_rows(cbind(id, learner = "ranger", best$performance), cbind(id, learner = "svm", svm_performance)),
            file.path(OUTPUTS_DIR, "refits", paste0(tag, ".csv")))
cat("Tables written to", OUTPUTS_DIR, "\n")
