<div align="justify">

# Data

These are the data of the greenhouse experiment used in the article, limited to what the article analyses. They are shared under the [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) license. The MIT license of the code does not cover them. Please cite the article when using them.

| File | Rows | Content |
|---|---|---|
| `spectra_mean_1nm.csv.xz` | 231 | Mean reflectance of each plant and canopy layer, every 1 nm from 350 to 2500 nm |
| `spectra_reflect_5nm.csv` | 231 | Reflectance every 5 nm from 400 to 2450 nm, used by the models |
| `spectra_snv_5nm.csv` | 231 | The same spectra after Standard Normal Variate (SNV): each spectrum minus its own mean, divided by its own standard deviation |
| `spectra_cr_5nm.csv` | 231 | The same spectra after continuum removal (CR), which compares absorption features against a common baseline |
| `severity.csv` | 117 | Plant height, lesion length and disease severity |
| `sample_split.csv` | 77 | Plants used for calibration (61) and validation (16) |

R reads the compressed `.xz` file directly, for example with `read.csv()`.

## Columns

| Column | Meaning |
|---|---|
| `plant` | Plant number. Each plant was assessed once, because the assessment is destructive |
| `clone` | Clone: 223 is susceptible and 236 is resistant |
| `inoculation` | `yes` for plants inoculated with *Ceratocystis fimbriata*, `no` for non-inoculated plants |
| `leaf` | Canopy layer: `apical`, `middle` or `basal` |
| `time` | Time after inoculation: `t1` 12 hours, `t2` 24 hours, `t3` 48 hours, `t4` 60 hours, `t5` 30 days, `t6` 40 days |
| `n_readings` | Number of spectral readings averaged for that plant and canopy layer |
| `C350` to `C2500`, `C400` to `C2450` | Value at that wavelength, in nm |
| `height`, `lesion` | Plant height and lesion length, in cm |
| `severity` | Lesion length divided by plant height, in % |
| `split` | `calibration` or `validation` |

## Before using the tables

- **Plants.** The spectra cover the 77 plants assessed before visual symptoms (12 to 60 hours). `severity.csv` also includes the 40 plants assessed at 30 and 40 days, which are used only in Figure A1.
- **Excluded plants.** Three inoculated plants of the resistant clone (39, 57 and 97) are not included, as described in the article: their disease development was abnormally high, which is suspected to result from genetic material leakage during propagation.
- **Calibration and validation.** The split was drawn from the complete experiment, so it cannot be drawn again from these tables. It is therefore stored in `sample_split.csv`.
- **Averages.** The individual spectral readings are not included. The averages give the three 5-nm tables again, with differences below 1e-14, as checked by the tests.

</div>
