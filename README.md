# ccpl_data_nlFit2026a

Data and analysis code for *Relations Among nlFit Parameters, Gender, and Math Anxiety in Three Samples of U.S. Adults* (Thompson, Cohen, and colleagues; Journal of Educational Psychology, revised submission).

## Layout

| Path | Contents |
|---|---|
| `runAnalysis.r` | Driver. Runs every dataset, with and without the outlier filters. |
| `nlFitRun.r` | Wrapper that fits the nlFit model to one participant's estimates. |
| `samClarissa*.r` | Study 1 pipeline: data prep, nlFit fits, traditional analyses, parameter analyses. |
| `Exp 1 Kent/` | Study 1 University sample: raw workbook `SSZpilotdata_for_dale_Sep17.xlsx` and the demographics table the prep step writes. |
| `Exp 2 Prolific/` | Study 1 Prolific sample: raw workbook `SSZ_SM1data_fordale_Sep17.xlsx` and its demographics table. |
| `mielicki 2022/` | Study 2 pipeline (`mielicki*.r`), its demographics table, and `OSF/`, the original data and codebook from Mielicki, Fitzsimmons, et al. (2023). |

## Running

Open R in this directory and run:

```r
source("runAnalysis.r")
```

The driver creates, inside each dataset folder, a `withoutOutliers/` and a `withOutliers/` folder holding the prepared data, the per-participant fits (in a folder named for the installed nlFit version, e.g. `nlFit_v0.3.0/`), and the analysis output. The flags at the top of `runAnalysis.r` (`runDataPrep`, `runNlFit`, `runTradAnalysis`, `runParamAnalysis`) select which steps run; `datasets` selects the datasets. The nlFit step is the slow one.

The demographics tables (`*Demographics.txt`, `mielickiDemo.txt`) are written by the prep step and read by the parameter-analysis step. They are included so the later steps can be run without re-running prep.

## Requirements

R packages: `nlFit` (0.3.0 or later; the scripts use the `pConceptual` parameter name introduced in 0.3.0), `smartGridSearch2`, `chutils`, `readxl`, `dplyr`, `tidyr`, `emmeans`, `lmerTest`.

`nlFit`, `smartGridSearch2` and `chutils` are at https://github.com/ccpluncw. The fits reported in the manuscript were produced with nlFit 0.2.0; 0.3.0 renamed one parameter and added a convenience function without changing the model, so fits from 0.3.0 reproduce the reported values.

## Citation

Thompson, C. A., & Cohen, D. J. (2025). A new model explaining number-line estimation performance across development: An adversarial collaboration. *American Psychologist, 80*(9), 1473. (The nlFit model.)

Mielicki, M. K., Fitzsimmons, C. J., Schiller, L. K., Scheibe, D., Taber, J. M., Sidney, P. G., Matthews, P. G., Waters, E. A., Coifman, K. G., & Thompson, C. A. (2023). Number lines can be more effective at facilitating adults' performance on health-related ratio problems than risk ladders and icon arrays. *Journal of Experimental Psychology: Applied, 29*(3), 529–543. https://doi.org/10.1037/xap0000456 (Study 2 data, in `mielicki 2022/OSF/`.)
