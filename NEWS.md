# Icare 1.0.1

### Dependency Refactoring
*   **Streamlined Dependencies**: Moved heavy/optional packages from `Imports` to `Suggests`, including the `mlr3` ecosystem (`mlr3`, `mlr3proba`, `mlr3tuning`, `mlr3pipelines`, `mlr3tuningspaces`, `paradox`), `shiny`, `shinydashboard`, `DT`, `NMF`, `Rtsne`, `umap`, `mclust`, `NbClust`, `factoextra`, `mice`, `rms`, `timeROC`, `nricens`, `survminer`, `maxstat`, `survcomp`, `risksetROC`, `shapviz`, `survex`, `rBayesianOptimization`, `autoReg`, `officer`, `flextable`, `randomForestSRC`, `Hmisc`, and others. Only core data-handling and plotting packages remain hard dependencies.
*   **Lazy Dependency Checks**: Added internal `.require_pkgs()` helper that validates optional packages at function entry and raises actionable `install.packages()` hints; applied across the model, prognosis, subtyping, stat, and viz modules.
*   **Namespaced Calls**: Removed blanket `@import` / `@importFrom` directives for optional packages and switched to explicit `pkg::fun()` calls (e.g. `shiny::`, `shinydashboard::`, `mlr3::`, `mlr3proba::`, `mlr3pipelines::`, `mlr3tuning::`, `paradox::`) in deployment apps and mlr3-based functions.
*   **Roxygen & Docs**: Upgraded `roxygen2` to 8.1.0 and regenerated `NAMESPACE` and man pages; wrapped runnable examples in `\dontrun{}`; declared `utils::globalVariables(c("BIC", "time"))`; dropped `mdbrown/rmda` from `Remotes`.
*   **`DMwR` via GitHub Remotes**: `DMwR` is archived on CRAN, so it is installed from the official CRAN read-only mirror `cran/DMwR` (added to `Remotes`, kept in `Suggests`). This fixes the GitHub Actions dependency step, which failed because `setup-r-dependencies` could not resolve `DMwR` from CRAN. `R/model_05_train.R`'s SMOTE path now prints `remotes::install_github("cran/DMwR")` when the package is missing.
*   **Housekeeping**: Added `.reasonix` / `reasonix.toml` to both `.Rbuildignore` and `.gitignore`.

# Icare 1.0.0

### New Features & Improvements
*   **Version 1.0 Release**: Completed version 1.0 of the Icare R package — functions and tests refined and officially updated.
*   **Non-binary Downstream Analysis**: Added `R/model_05b_nonbinary_downstream.R` for multiclass and regression downstream tasks (multi-class ROC / confusion matrices, regression fit diagnostics, preprocessing helpers).
*   **Prognosis Extensions**: Added `R/prognosis_04_Extensions.R` covering time-dependent ROC comparison, repeated cross-validation evaluation, DeLong-style model comparison, and Bayesian hyperparameter tuning.
*   **Bundled Reference Data**: Added built-in `allmodel.rda` and `inst/extdata/all_model.csv` with model metadata.
*   **`sub_dir()` Utility**: New helper creating subdirectories under the global output root, with a session-temp fallback when no root is set.

### Bug Fixes
*   **Encoding Cleanup**: Replaced Unicode dashes and non-ASCII characters in code, docs, and messages with ASCII equivalents across deploy and prognosis modules.
*   **Ensemble Prediction**: Hardened error handling in ensemble prediction dispatch and non-binary probability conversion in `ModelDeployment`.
*   **Dependencies**: Added `maxstat` and `survcomp` to `Imports`, and `bbotk` / `mlr3mbo` to `Suggests`.

### Packaging & Housekeeping 
*   **GitHub Packaging**: Added `mlr-org/mlr3proba` to `Remotes` so a current version is installed from GitHub.
*   **Ignore Temporary Files**: Removed stray `.Rhistory` and added an ignore rule.
*   **Merge PR #7**: Resolved merge conflicts, keeping both sides' complementary changes.

# Icare 0.1.2

### Refactoring & Improvements
*   **Code Consolidation**: Merged `run_univariate_cox_analysis` from `module_3` and `module_4` into a single, enhanced implementation in `module_4_pn5_0_cox_univariate.R`.
    *   Added support for **covariate adjustment** and custom **P-value thresholds** (migrated from module 3).
    *   Retained **robust status column handling** and data validation (from module 4).
    *   Deprecated the redundant function in `module_3` to prevent namespace conflicts and installation warnings.

# Icare 0.1.1

### Bug Fixes
*   **Scope Issues**: Fixed variable scope errors in `survfit` and `survdiff` calls across multiple models (Lasso, PLS, RSF, CoxPH, SuperPC) by implementing formula injection.
*   **CoxBoost**: Resolved row mismatch errors in hazard ratio calculation caused by missing feature names.
*   **RSF**: Fixed execution errors in Random Survival Forest model when variable importance is empty; added robust handling for empty result data frames.
*   **SuperPC**: Fixed "logical subscript too long" error by ensuring proper numeric conversion of transposed feature matrices.
*   **Theme Conflicts**: Resolved conflicts between `ggprism` and `ggplot2` themes by standardizing on `theme_classic` and removing unstable dependencies.

### New Features & Improvements
*   **Dependencies**: Added `randomForestSRC`, `plsRcox`, `superpc`, and `rmda` to package dependencies to support advanced modeling and DCA.
*   **DCA Support**: Enhanced Decision Curve Analysis (DCA) integration, ensuring robust extraction of best model results.


# Icare 0.1.0

Icare is a comprehensive R package designed for survival analysis, clinical prediction modeling, and bioinformatics data analysis. It provides a unified framework for:

*   **Data Preprocessing**: Missing value handling, outlier detection, normalization, and descriptive statistics.
*   **Survival Modeling**: Implementation of various survival models including Lasso-Cox, CoxPH, Random Survival Forests (RSF), PLS-Cox, CoxBoost, and SuperPC.
*   **Model Evaluation**: Time-dependent ROC analysis, Kaplan-Meier survival curves, Decision Curve Analysis (DCA), and calibration assessment.
*   **Subtyping Analysis**: Molecular subtyping using K-means, NMF, and consensus clustering with t-SNE/UMAP visualization.
*   **Visualization**: High-quality, publication-ready plots for all analysis steps.


