# Section 2.3: from the mean spectra at 1 nm to the tables with reflectance, SNV and continuum removal at 5 nm.
source("R/setup.R")

tables <- prepare_spectra()
for (pretreatment in PRETREATMENTS) {
  write_table(tables[[pretreatment]], spectra_path(pretreatment, OUTPUTS_DIR))
}
cat("Tables written to", OUTPUTS_DIR, "\n")
