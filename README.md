<div align="justify">

# Disease detection of Ceratocystis wilt with leaf spectroscopy

Code and data to reproduce the analysis of the article:

> Santiago, M.I.L., Franceschini, M.H.D., Del Ponte, E.M., do Amaral, C.H. and Alfenas, R.F., 2025. Leaf Diffuse Reflectance Spectroscopy for Early Detection of Ceratocystis Wilt in Eucalyptus Cuttings. *Forest Pathology*, 55(6), e70052. https://doi.org/10.1111/efp.70052

The study evaluates whether leaf reflectance from 400 to 2450 nm can detect Ceratocystis wilt, caused by *Ceratocystis fimbriata*, in *Eucalyptus* cuttings before visible symptoms appear. Two clones, one susceptible and one resistant, were inoculated in a greenhouse and assessed 12 to 60 hours after inoculation.

All values in Tables 1, A1 and A2 of the article are reproduced by the results in this repository. Some descriptions in the article and the bars of its Figure 5 differ from what the code was intended to do, as explained at the end of this page.

## How the method works

1. **Spectra.** The readings of each plant and canopy layer are averaged, resampled every 5 nm and used as reflectance, after Standard Normal Variate (SNV, which corrects differences in overall intensity) and after continuum removal (CR, which compares absorption features against a common baseline).
2. **Data split.** The 77 plants are split into calibration (61) and validation (16). The calibration plants are divided again into five folds, five times, to compare models.
3. **Feature selection.** Recursive Feature Elimination (RFE) and a Genetic Algorithm (GA), both with Random Forest (RF), select the wavelengths in each of the 25 folds. The final set for a target number of features is the candidate with the highest F1-score minus its standard deviation (the lower bound) when refitted on the 25 folds.
4. **Number of features.** The authors chose a target of 75 features from the cross-validation results, judging model stability, and used the same target for all pre-treatments so they could be compared.
5. **Models.** RF and Support Vector Machine (SVM) models use the selected wavelengths, and Partial Least Squares Discriminant Analysis (PLS-DA) uses all wavelengths. The final models are tested on the validation plants.
6. **Importance.** SHAP values estimate how much each wavelength contributes to the predictions.

## Main results

Best model of each approach for the apical canopy layer. The values are for the 16 validation plants, and the value after ± is the standard deviation in the cross-validation of the calibration plants.

| Model | Pre-treatment | Balanced accuracy | F1-score |
|---|---|---|---|
| PLS-DA, all wavelengths | SNV | 0.625 ± 0.113 | 0.625 ± 0.119 |
| RF with RFE | SNV | 0.750 ± 0.107 | 0.750 ± 0.131 |
| SVM with GA | Reflectance | 0.750 ± 0.125 | 0.714 ± 0.203 |

The complete values are in `results/` and the figures in `figures/`.

## Reproducing the results

Install [R](https://www.r-project.org/) 4.4 and the packages with the versions in `renv.lock`, using [renv](https://rstudio.github.io/renv/) (`renv::restore()`). Run the scripts from the repository folder.

The tables and figures take about a minute, because they use the saved results:

```text
Rscript scripts/07_tables_and_figures.R
Rscript tests/test_reproduction.R
```

Each step can also be run again. Its tables are written to `outputs/`, so the saved results are kept for comparison.

| Script | Article | Time |
|---|---|---|
| `01_prepare_spectra.R` | Section 2.3 | Seconds |
| `02_feature_selection.R` | Section 2.4, candidate feature sets | About 3 minutes per RFE fold, and about 8 minutes per GA fold and target size |
| `03_final_feature_sets.R` | Section 2.4, final feature sets | About 30 minutes per combination |
| `04_final_models.R` | Tables 1, A1 and A2, RF and SVM | About 15 minutes per feature selection method and canopy layer |
| `05_plsda.R` | Table 1, PLS-DA | About 3 minutes |
| `06_shap.R` | Figure 4 | About an hour with 4 processor cores |
| `07_tables_and_figures.R` | Tables and Figures 1, 3, 4, 5 and A1 | About a minute |

The options of each script are described at its beginning. The figures have a simpler layout than in the article, where some were assembled by hand. Figure 2 is a schematic and is not included.

The random seeds of the original run are kept, so the same results are obtained, with two exceptions:

- **Genetic algorithm.** The search runs in parallel and is not repeatable: in two tests, a new search shared only 16 and 17 of its wavelengths with the original run. The candidate sets of the original run are in `results/`, and all the following steps reproduce the article from them.
- **SHAP values.** They come from parallel simulations, so they change slightly between runs. In a test, the correlation with the original values was at least 0.99 for every model, and the largest difference was 0.4 percentage points.

`Rscript tests/test_reproduction.R --models` refits parts of the analysis and compares them with the original run, which takes about 15 minutes.

## Repository contents

| Folder | Content |
|---|---|
| `data/` | Spectra, severity and the data split, explained in `data/README.md` |
| `R/` | The analysis code |
| `scripts/` | The steps to run |
| `results/` | The tables of the article and the results of the original run |
| `figures/` | The figures obtained with this code |
| `tests/` | Checks, with the values printed in the article in `tests/expected/` |
| `tools/` | The script that exported the data and results, which needs the original files |

| Table in `results/` | Content |
|---|---|
| `table_1.csv`, `table_a1.csv`, `table_a2.csv` | Tables of the article |
| `validation_metrics.csv`, `validation_predictions.csv` | Final models on the validation plants, with their hyperparameters. Predictions were saved only for the apical layer |
| `calibration_predictions.csv` | Out-of-fold predictions of the calibration plants, used in Figure 5 |
| `final_feature_sets.csv` | Final wavelengths of each target number of features |
| `cv_final_feature_sets.csv` | The 25 refits of each final feature set, which give the standard deviations of the tables |
| `cv_candidate_feature_sets.csv` | The 25 candidate feature sets of each combination and their results |
| `cv_plsda_components.csv` | PLS-DA with 1 to 15 components in the 25 folds |
| `shap_importance.csv` | Importance of the wavelengths, in percent of the total of each model |
| `figure_5_values.csv`, `figure_5_values_intended.csv` | Values of Figure 5, explained below |

The original scripts also tested each candidate model on the validation plants. Those results were not used to select models or reported in the article, so they are not included.

## Citation and license

Please cite the article above when using this code or these data. Citation details are in `CITATION.cff`.

The code is released under the [MIT license](LICENSE). The data are released under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license, as described in `data/README.md`.

Parts of the code, tests and documentation of this repository were reorganized with the help of an AI coding assistant (Claude, by Anthropic). The author directed the work, checked the results and is responsible for the content.

## Differences from the published article

Preparing this repository showed that some descriptions in the article differ from the code that produced its results. The values in Tables 1, A1 and A2 are not affected, except where noted.

| Topic | Article | Code and results |
|---|---|---|
| Resampling (Section 2.3) | Spline function | `prospectr::resample2`, based on the band width (full width at half maximum): 3 nm up to 700 nm and 10 nm above |
| Raw bands (Section 2.3) | 2150 | 2151, from 350 to 2500 nm |
| Readings (Section 2.3) | Three per leaf, two leaves per canopy layer | From 2 to 10 per plant and canopy layer, 6 in 67 of 231 cases, shown in `n_readings` |
| RFE (Section 2.4) | One feature removed at a time | Wavelengths ranked once and reduced to fixed sizes (411, 200, 100, 75, 50, 25, 10 and 5), of which 25 to 100 were considered |
| GA mutation (Section 2.4) | 0.10 for each gene | 0.10 for each candidate (`pmutation` of the GA package) |
| Target number of features (Section 2.4) | Highest lower bound of the F1-score | Chosen by the authors, as described in *How the method works* |
| Number of GA features (Table 1, Figure 3) | 75 in each case | 74, 75 and 65 for reflectance, SNV and CR, because the GA treats the target as a penalty, not a limit |
| Standard deviations (Table 1) | 25 complete feature selection and model fitting procedures | 25 refits of the single final feature set, which gives smaller values: 0.107 instead of 0.128 for the balanced accuracy of RF with RFE and SNV |
| Bars of Figure 5 | PLS-DA, RF+RFE and SVM+GA | See below |

**Figure 5.** The code selected the predictions of each model by its position in a list, not by its pre-treatment. As a result, the bars show PLS-DA with reflectance and SVM+GA with CR, while the error bars show the intended pre-treatments (SNV and reflectance). RF+RFE is correct. `figures/figure_5.png` reproduces the published figure, and `figures/figure_5_intended.png` shows the intended pre-treatments in the bars. With the intended pre-treatments, the accuracy of the susceptible clone still generally increases at later times, and that of the resistant clone is still lower overall and stable or declining, as the article states.

</div>
