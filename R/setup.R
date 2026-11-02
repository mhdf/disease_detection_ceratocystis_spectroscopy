# Load the packages and the code of the analysis. Scripts are run from the repository folder.

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(GA)
  library(doParallel)
  library(MLmetrics)
  library(prospectr)
})
options(scipen = 999)

source("R/config.R")
source("R/data.R")
source("R/models.R")
source("R/selection.R")
source("R/final_models.R")
source("R/plsda.R")
source("R/shap.R")

dir.create(OUTPUTS_DIR, showWarnings = FALSE)

write_table <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  write.table(x, path, sep = ",", row.names = FALSE, na = "")
}

# Empty cells of the result tables are missing values
read_result <- function(path) read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
