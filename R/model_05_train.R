# ---------------------------------------------------------------------------
# Defensive summaryFunction for metric = "ROC" (binary classification).
#
# Root cause this guards against: caret::twoClassSummary() calls
# pROC::roc() once per resample (fold). If, within a SINGLE fold, the
# predicted probabilities for the positive class come out with (near-)zero
# variance -- e.g. a perfectly-separating glm on a small/high-dimensional
# fold, or an rpart fold whose chosen tree collapses to a single leaf --
# pROC::roc() returns NA for that fold's ROC/Sens/Spec. caret's internal
# check `"Something is wrong; all the ROC metric values are missing"` is
# only triggered when EVERY fold for that model came back NA, but a single
# resample producing NA is enough to poison bestTune selection (NA
# propagates through mean/SD aggregation for that tuning row), and with
# small tuning grids (e.g. glm/rpart with 1-3 candidate rows) it is common
# for this to cascade into the entire model being dropped, as observed with
# glm/rpart in practice.
#
# Fix: fall back to a neutral, "no-information" score (ROC = Sens = Spec =
# 0.5, i.e. equivalent to random guessing) for the offending fold ONLY,
# instead of propagating NA. This keeps that one resample's contribution to
# the model's CV mean numerically well-defined (dragging the mean toward
# 0.5, which is the honest, conservative interpretation of "this fold gave
# no usable signal") rather than causing caret to discard the model
# outright. It does not fabricate an optimistic result, and it never
# affects folds that computed a valid AUC.
# ---------------------------------------------------------------------------
.safeTwoClassSummary <- function(data, lev = NULL, model = NULL) {
  res <- tryCatch({
    caret::twoClassSummary(data = data, lev = lev, model = model)
  }, error = function(e) NULL)
  
  if (!is.null(res) && !any(is.na(res))) {
    return(res)
  }
  
  # Defensive fallback: happens when a fold's predicted probabilities have
  # (near-)zero variance, or pROC::roc() errored/returned NA for another
  # reason. Score this fold as no-better-than-random rather than NA.
  c(ROC = 0.5, Sens = 0.5, Spec = 0.5)
}


# ---------------------------------------------------------------------------
# Defensive summaryFunction for metric = "Accuracy" / "Kappa"
# (classification, no class probabilities required).
#
# Root cause this guards against: caret::defaultSummary() calls
# confusionMatrix() internally, which errors/returns NA for the WHOLE fold
# if even a single row's predicted class is NA. In practice this most often
# happens with `glm` (a categorical predictor has a level present in the
# held-out fold's rows but absent from the fold actually used to fit the
# model, so the design matrix can't build a contrast for it and predict()
# returns NA for those rows) or `rpart` (similar unseen-level / degenerate
# split issue). This is exactly the failure mode observed with glm/rpart:
# "Something is wrong; all the Accuracy metric values are missing".
#
# Fix, in order of preference:
#   1. Try caret::defaultSummary() as normal.
#   2. If it errored or returned NA, drop only the rows with NA pred/obs
#      and recompute on the remaining valid rows (so one bad row doesn't
#      void an otherwise-fine fold).
#   3. If literally every row in the fold is unusable, fall back to a
#      neutral chance-level score (Accuracy = 1 / number of classes,
#      Kappa = 0) instead of NA, so this fold drags the model's CV mean
#      toward "no better than guessing" rather than poisoning it outright.
# ---------------------------------------------------------------------------
.safeDefaultSummary <- function(data, lev = NULL, model = NULL) {
  res <- tryCatch({
    caret::defaultSummary(data = data, lev = lev, model = model)
  }, error = function(e) NULL)
  
  if (!is.null(res) && !any(is.na(res))) {
    return(res)
  }
  
  # Retry after dropping rows with NA prediction/observation.
  ok <- stats::complete.cases(data[, c("pred", "obs")])
  if (any(ok) && !all(ok)) {
    res2 <- tryCatch({
      caret::defaultSummary(data = data[ok, , drop = FALSE], lev = lev, model = model)
    }, error = function(e) NULL)
    if (!is.null(res2) && !any(is.na(res2))) return(res2)
  }
  
  # Every row unusable (or partial-row recovery still failed): chance-level
  # fallback rather than NA.
  n_lev <- length(lev)
  if (is.null(n_lev) || n_lev < 1) n_lev <- tryCatch(nlevels(data$obs), error = function(e) 2)
  c(Accuracy = if (n_lev > 0) 1 / n_lev else 0.5, Kappa = 0)
}


# ---------------------------------------------------------------------------
# Defensive summaryFunction for regression metrics (RMSE / Rsquared / MAE).
#
# Same failure mode as .safeDefaultSummary(), applied to numeric outcomes:
# if any row's prediction is NA, caret::defaultSummary() returns NA for
# the whole fold. Here we drop unusable rows and recompute; if literally
# nothing is usable we return a deliberately BAD (not NA) score so this
# fold penalizes rather than crashes the model -- RMSE/MAE = Inf and
# Rsquared = 0 both correctly sort as "worst" during model comparison,
# whereas NA would have been silently dropped or crashed ranking.
# ---------------------------------------------------------------------------
.safeRegressionSummary <- function(data, lev = NULL, model = NULL) {
  res <- tryCatch({
    caret::defaultSummary(data = data, lev = lev, model = model)
  }, error = function(e) NULL)
  
  if (!is.null(res) && !any(is.na(res))) {
    return(res)
  }
  
  ok <- stats::complete.cases(data[, c("pred", "obs")])
  if (any(ok) && !all(ok)) {
    res2 <- tryCatch({
      caret::defaultSummary(data = data[ok, , drop = FALSE], lev = lev, model = model)
    }, error = function(e) NULL)
    if (!is.null(res2) && !any(is.na(res2))) return(res2)
  }
  
  c(RMSE = Inf, Rsquared = 0, MAE = Inf)
}


# ---------------------------------------------------------------------------
# Fallback list of known two-class-only caret methods (package-level
# constant). The PRIMARY check (in train_and_evaluate_models() and
# ModelTrainAnalysis()) reads caret's own `getModelInfo(method)$tags` for
# "Two Class Only" -- that is the authoritative, self-updating source (it's
# literally what generates https://topepo.github.io/caret/train-models-by-tag.html#Two_Class_Only)
# and covers all 200+ caret methods without manual upkeep. This list only
# matters as a fallback if a model's `tags` field is missing/empty.
# ---------------------------------------------------------------------------
BINARY_ONLY_METHODS <- c("glm", "glmStepAIC", "glm.nb", "plr", "LMT", "bayesglm")

#' Extract Filtered Training and Testing Data from a Model Data Object
#'
#' This function extracts the filtered training and testing datasets stored in the
#' 'filtered.set' slot of a 'Train_Model' object. If the 'filtered.set' slot is
#' not available or there is an error in accessing it, the function returns NULL for
#' both training and testing datasets.
#'
#' @param object An object of class 'Train_Model' containing the 'filtered.set' slot,
#'   which holds the filtered training and testing data.
#'
#' @returns A list with two elements:
#'   - `training`: The filtered training dataset (if available), otherwise NULL.
#'   - `testing`: The filtered testing dataset (if available), otherwise NULL.
#'
#' @export
#' @examples
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' # Normally split later; for demo, we just show extraction (empty)
#' extracted <- Extract_filtered.set(model)
#' str(extracted)
Extract_filtered.set <- function(object) {
  train <- tryCatch(object@filtered.set$training,
                    error = function(e) NULL)
  test <- tryCatch(object@filtered.set$testing,
                   error = function(e) NULL)
  return(list(training = train, testing = test))
}

#' Train Multiple caret Models with Automatic Tuning and Imbalance Handling
#'
#' Supports binary classification, multiclass classification, and regression.
#' Automatically selects the appropriate summary function. User can supply
#' custom tune grids or use `tuneLength` for automatic grid generation.
#'
#' All methods requested in a single call are evaluated on the SAME
#' cross-validation folds (built once via [caret::createFolds()] /
#' [caret::createMultiFolds()] and passed to `trainControl(index = ...)`),
#' so that CV-based comparisons across models are fair rather than each
#' model drawing its own random resamples. This is only done automatically
#' for `control$method %in% c("cv", "repeatedcv")`; for other resampling
#' schemes (e.g. `"boot"`), folds are not pre-shared and a message is
#' printed to make that explicit.
#'
#' NOTE on parallel processing: this function runs strictly sequentially.
#' Parallel backends (`doParallel`/`foreach`) were removed because, combined
#' with a shared RNG seed and per-fold sampling (`smote`/`rose`), parallel
#' workers do not reliably reproduce the same random state across runs/
#' platforms, and because failures inside a parallel worker are much harder
#' to diagnose/report per-model than a sequential `tryCatch`. If you need
#' speed, parallelize across `methods` yourself at a higher level (e.g. one
#' \code{train_and_evaluate_models()} call per method, dispatched with your
#' own scheduler) with `set.seed()` fixed inside each worker.
#'
#' NOTE on robustness: each method is checked BEFORE training for (a)
#' availability in \code{caret::modelLookup()}, (b) installed dependency
#' packages (via \code{caret::getModelInfo()$library}), (c) task-type
#' compatibility (Classification vs. Regression), and (d) known
#' binary-classification-only methods (e.g. `glm`, `glmStepAIC`) being used
#' on outcomes with more than 2 classes. Methods failing any check are
#' skipped with an informative message rather than aborting the whole run.
#'
#' @param data Training data frame.
#' @param methods Character vector of caret model methods.
#' @param control List with elements `method`, `number`, `repeats` for trainControl.
#' @param tune_grids Named list of data frames of tuning parameters (optional).
#' @param tuneLength Integer, number of automatic tuning levels (default 3).
#' @param preProcess Character vector of preprocessing methods (passed to caret::train).
#'   Set to `NULL` if data already preprocessed externally.
#' @param loocv_threshold Numeric, sample size below which LOOCV is forced.
#' @param classProbs Logical, compute class probabilities.
#' @param allowParallel Deprecated and ignored (kept only for backward
#'   compatibility with existing call sites). This function always runs
#'   sequentially now; see Details.
#' @param group_col Name of outcome column.
#' @param metric Evaluation metric (e.g., "ROC", "Accuracy", "RMSE").
#' @param imbalance_handling Character: "none" (default), "auto", "up", "down",
#'   "smote", "rose", or "weights". Defines how to handle class imbalance.
#'   "auto" applies "smote" if minority class proportion < imbalance_threshold
#'   AND the task is binary; for multiclass tasks "auto" is not supported by
#'   caret's `trainControl(sampling = ...)` and is downgraded to "none" with
#'   a warning.
#'   NOTE: "smote" requires the (CRAN-archived) 'DMwR' package and "rose"
#'   requires the 'ROSE' package; both are checked explicitly and a clear
#'   error (with install instructions) is raised if missing, rather than
#'   letting caret fail with an opaque internal error.
#' @param imbalance_threshold Numeric (0-1), proportion below which auto
#'   handling is triggered. Default 0.2.
#' @param seed Optional integer. If supplied, set ONCE at the start of this
#'   function (before the CV fold indices are built), so that (a) fold
#'   construction is reproducible and (b) every model in `methods` starts
#'   from the same RNG state for anything not already fixed by the shared
#'   folds (e.g. `smote`/`rose` resampling inside a fold, or stochastic
#'   algorithms like `rf`/`gbm`/`nnet`). The seed is intentionally NOT reset
#'   again inside the per-model loop: re-seeding before every single model
#'   would make later models in `methods` replay the exact same random draws
#'   as earlier ones for anything caret does internally between folds (tree
#'   bootstraps, random restarts, etc.), which is not what "reproducible" is
#'   meant to guarantee -- it should reproduce the whole run end-to-end, not
#'   force every model onto an identical internal random trajectory.
#' @return List of caret train objects.
#' @importFrom caret train trainControl getModelInfo modelLookup defaultSummary
#'   twoClassSummary createFolds createMultiFolds
#' @importFrom stats as.formula
#' @export
#'
#' @examples
#' \dontrun{
#' mtcars$am <- as.factor(mtcars$am)
#' control <- list(method = "cv", number = 3)
#' results <- train_and_evaluate_models(data = mtcars,
#'                                      methods = c("glm", "rf"),
#'                                      control = control,
#'                                      group_col = "am",
#'                                      metric = "Accuracy",
#'                                      seed = 123)
#' print(names(results))
#' }
train_and_evaluate_models <- function(data,
                                      methods,
                                      control = list(method = "repeatedcv", number = 10, repeats = 3),
                                      tune_grids = NULL,
                                      tuneLength = 3,
                                      preProcess = NULL,
                                      loocv_threshold = 100,
                                      classProbs = TRUE,
                                      allowParallel = FALSE,
                                      group_col = "group",
                                      metric = "ROC",
                                      imbalance_handling = c("none", "auto", "up", "down", "smote", "rose", "weights"),
                                      imbalance_threshold = 0.2,
                                      seed = NULL) {
  
  # ---- Dependency checks (always required) ----
  .check_required_pkgs <- function() {
    required <- c("caret")
    missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
    if (length(missing) > 0) {
      stop("Missing packages: ", paste(missing, collapse = ", "),
           ". Please install them before using this function.")
    }
  }
  .check_required_pkgs()
  
  # ---- Dependency check for imbalance-handling sampling methods (conditional) ----
  # caret::trainControl(sampling = ...) silently delegates "smote" to the
  # 'DMwR' package and "rose" to the 'ROSE' package. Neither is a hard
  # dependency of this function, so we check for them only when actually
  # requested, and fail with an actionable message instead of letting
  # caret::train() error out deep inside its internals.
  .check_sampling_pkg <- function(sampling_method) {
    if (is.null(sampling_method)) return(invisible(TRUE))
    if (sampling_method == "smote" && !requireNamespace("DMwR", quietly = TRUE)) {
      stop(
        "imbalance_handling = 'smote' requires the 'DMwR' package, which is not installed.\n",
        "'DMwR' was archived on CRAN. Install an archived source build, e.g.:\n",
        "  install.packages(\n",
        "    'https://cran.r-project.org/src/contrib/Archive/DMwR/DMwR_0.4.1.tar.gz',\n",
        "    repos = NULL, type = 'source')\n",
        "Or avoid the dependency entirely with imbalance_handling = 'up', 'down', or 'weights'."
      )
    }
    if (sampling_method == "rose" && !requireNamespace("ROSE", quietly = TRUE)) {
      stop(
        "imbalance_handling = 'rose' requires the 'ROSE' package, which is not installed.\n",
        "Install it with install.packages('ROSE'), or use imbalance_handling = ",
        "'up', 'down', or 'weights' instead."
      )
    }
    invisible(TRUE)
  }
  
  imbalance_handling <- match.arg(imbalance_handling)
  
  if (isTRUE(allowParallel)) {
    message("allowParallel = TRUE was requested but parallel processing has ",
            "been removed from this function; running sequentially instead.")
  }
  allowParallel <- FALSE
  
  # ---- Determine task type ----
  is_regression <- metric %in% c("RMSE", "Rsquared", "MAE")
  
  if (is_regression) {
    data[[group_col]] <- as.numeric(data[[group_col]])
    classProbs <- FALSE
    summary_fn <- .safeRegressionSummary
    imbalance_handling <- "none"
  } else {
    data[[group_col]] <- factor(data[[group_col]])
    levels(data[[group_col]]) <- make.names(levels(data[[group_col]]))
    n_lv <- length(levels(data[[group_col]]))
    if (metric == "ROC" && n_lv != 2) {
      stop("metric = 'ROC' is only valid for binary classification (found ",
           n_lv, " levels). Use metric = 'Accuracy' or 'Kappa' for multiclass.")
    }
    summary_fn <- if (metric == "ROC") .safeTwoClassSummary else .safeDefaultSummary
  }
  
  # ---- Auto imbalance handling ----
  if (!is_regression && imbalance_handling == "auto") {
    if (n_lv != 2) {
      # caret::trainControl(sampling = "smote"/"rose") only implements the
      # binary case; silently leaving imbalance_handling == "auto" here would
      # later fall through with sampling_method left NULL while the variable
      # itself still reads "auto" (not "none"), which is a latent bug users
      # would only notice by the *absence* of any imbalance correction.
      warning("imbalance_handling = 'auto' is only supported for binary ",
              "classification (found ", n_lv, " classes). Downgrading to 'none'. ",
              "Consider up/down-sampling manually, or class weights via ",
              "imbalance_handling = 'weights', for multiclass imbalance.")
      imbalance_handling <- "none"
    } else {
      tbl <- table(data[[group_col]])
      min_prop <- min(tbl) / sum(tbl)
      if (min_prop < imbalance_threshold) {
        message("Auto imbalance handling triggered: minority class proportion = ",
                round(min_prop, 3), " < ", imbalance_threshold, ". Applying SMOTE.")
        imbalance_handling <- "smote"
      } else {
        message("Auto imbalance handling not triggered: minority proportion = ",
                round(min_prop, 3), " >= ", imbalance_threshold, ". No resampling.")
        imbalance_handling <- "none"
      }
    }
  }
  
  # ---- Class weights ----
  use_weights <- FALSE
  weights <- NULL
  if (!is_regression && imbalance_handling == "weights") {
    tbl <- table(data[[group_col]])
    wts <- 1 / tbl[as.character(data[[group_col]])]
    wts <- as.numeric(wts) * (length(wts) / sum(wts))  # normalize: mean weight = 1
    weights <- wts
    use_weights <- TRUE
  }
  
  # ---- Build trainControl ----
  n_samples <- nrow(data)
  formula <- stats::as.formula(paste(group_col, "~ ."))
  
  if (!is.null(seed)) set.seed(seed)
  
  if (control$method == "LOOCV" || n_samples <= loocv_threshold) {
    fitControl <- caret::trainControl(
      method = "LOOCV",
      classProbs = classProbs,
      summaryFunction = summary_fn,
      savePredictions = "final",
      allowParallel = allowParallel
    )
    cat("Using LOOCV (n =", n_samples, "samples)\n")
  } else {
    sampling_method <- NULL
    if (!is_regression && imbalance_handling %in% c("up", "down", "smote", "rose")) {
      sampling_method <- imbalance_handling
    }
    .check_sampling_pkg(sampling_method)
    
    repeats_val <- if (!is.null(control$repeats)) control$repeats else 1
    
    if (control$method %in% c("cv", "repeatedcv")) {
      # Pre-generate resampling indices ONCE so every method in `methods` is
      # evaluated on identical folds. Without this, each caret::train() call
      # draws its own random folds and CV metrics across models are not
      # directly comparable (a "better" mean CV score could just reflect an
      # easier random split rather than a better model).
      resample_index <- if (control$method == "repeatedcv") {
        caret::createMultiFolds(data[[group_col]], k = control$number, times = repeats_val)
      } else {
        caret::createFolds(data[[group_col]], k = control$number, returnTrain = TRUE)
      }
      fitControl <- caret::trainControl(
        method = control$method,
        number = control$number,
        repeats = repeats_val,
        classProbs = classProbs,
        summaryFunction = summary_fn,
        savePredictions = "final",
        allowParallel = allowParallel,
        sampling = sampling_method,
        index = resample_index
      )
      fold_note <- "| folds shared across all methods for fair comparison"
    } else {
      # e.g. "boot", "adaptive_cv", etc. -- caret builds its own resamples
      # internally for these and does not accept a simple pre-built `index`
      # in the same way; we fall back to per-model resampling and are
      # explicit about the resulting caveat.
      fitControl <- caret::trainControl(
        method = control$method,
        number = control$number,
        repeats = repeats_val,
        classProbs = classProbs,
        summaryFunction = summary_fn,
        savePredictions = "final",
        allowParallel = allowParallel,
        sampling = sampling_method
      )
      fold_note <- paste0(
        "| NOTE: control$method = '", control$method, "' does not use ",
        "pre-shared CV folds here -- CV metrics across models may not be ",
        "perfectly comparable. Use 'cv' or 'repeatedcv' for guaranteed fold-sharing."
      )
    }
    
    cat("Using", control$method, "with", control$number, "folds",
        if (repeats_val > 1) paste("and", repeats_val, "repeats"),
        "(n =", n_samples, ")",
        if (!is.null(sampling_method)) paste("| sampling:", sampling_method),
        fold_note, "\n")
  }
  
  results <- list()
  needed_type <- if (is_regression) "Regression" else "Classification"
  
  # Small helper: attempt caret::train() with a given argument list, retrying
  # once with `verbose` dropped if the underlying fit function doesn't accept
  # it (many caret methods forward `...` straight to a fitting function with
  # no `verbose` formal, e.g. some via do.call-style dispatch). Different
  # caret/model versions phrase this error differently -- base R's usual
  # "unused argument (verbose = ...)" for some methods, but e.g. rpart can
  # instead raise "Argument verbose not matched" -- so we match on any error
  # message that mentions `verbose` rather than requiring one exact phrasing.
  # Not every caret method supports `verbose`, so we no longer pass it
  # unconditionally.
  .safe_train <- function(train_args) {
    tryCatch(
      do.call(caret::train, train_args),
      error = function(e) {
        if (grepl("verbose", e$message, ignore.case = TRUE) &&
            "verbose" %in% names(train_args)) {
          train_args$verbose <- NULL
          tryCatch(do.call(caret::train, train_args),
                   error = function(e2) stop(e2))
        } else {
          stop(e)
        }
      }
    )
  }
  
  # ---- Train each method ----
  for (method in methods) {
    
    # -- 1. Method known to caret at all? --
    if (!method %in% caret::modelLookup()$model) {
      message("  Skipping '", method, "': not in caret::modelLookup().")
      next
    }
    
    minfo <- tryCatch(caret::getModelInfo(method, regex = FALSE)[[1]],
                      error = function(e) NULL)
    if (is.null(minfo)) {
      message("  Skipping '", method, "': caret::getModelInfo() failed to return details.")
      next
    }
    
    # -- 2. Required package(s) installed? --
    pkgs <- minfo$library
    if (!is.null(pkgs) && length(pkgs) > 0) {
      missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
      if (length(missing_pkgs) > 0) {
        message("  Skipping '", method, "': requires package(s) not installed: ",
                paste(missing_pkgs, collapse = ", "), ".")
        next
      }
    }
    
    # -- 3. Task-type compatible (Classification vs Regression)? --
    if (!is.null(minfo$type) && !(needed_type %in% minfo$type)) {
      message("  Skipping '", method, "': does not support ", needed_type,
              " (supports: ", paste(minfo$type, collapse = ", "), ").")
      next
    }
    
    # -- 4. Two-class-only method used on >2 classes? --
    # caret tags every model definition with a `tags` character vector (this
    # is exactly what generates the official "train Models By Tag" /
    # "Two Class Only" page: https://topepo.github.io/caret/train-models-by-tag.html#Two_Class_Only).
    # Reading it programmatically is authoritative and self-updating -- it
    # covers every current AND future caret model without us having to hand-
    # maintain a name list (which is necessarily incomplete: caret ships 200+
    # methods, and new ones are added over time). BINARY_ONLY_METHODS is kept
    # only as a fallback for the rare case a model's `tags` element is
    # missing/empty (e.g. some custom/older model definitions).
    is_two_class_only <- (!is.null(minfo$tags) && "Two Class Only" %in% minfo$tags) ||
      method %in% BINARY_ONLY_METHODS
    if (!is_regression && exists("n_lv") && n_lv > 2 && is_two_class_only) {
      message("  Skipping '", method, "': tagged 'Two Class Only' by caret ",
              "(fits a two-class model internally), but the outcome has ",
              n_lv, " classes. Use a multiclass-capable method instead ",
              "(e.g. 'multinom', 'rf', 'nnet', 'xgbTree').")
      next
    }
    
    # -- 5. Class-probability capability, if the run needs it --
    needs_probs <- isTRUE(classProbs) && (metric == "ROC" || identical(summary_fn, .safeTwoClassSummary))
    if (needs_probs && is.null(minfo$prob)) {
      message("  Skipping '", method, "': cannot output class probabilities, ",
              "but metric = '", metric, "' requires them. Use an Accuracy/Kappa-",
              "based metric for probability-incapable methods.")
      next
    }
    
    cat("  Training", method, "... ")
    
    if (!is.null(tune_grids) && method %in% names(tune_grids) &&
        !is.null(tune_grids[[method]])) {
      tuneGrid <- tune_grids[[method]]
      useTuneLength <- FALSE
    } else {
      tuneGrid <- NULL
      useTuneLength <- TRUE
    }
    
    base_args <- list(
      form = formula,
      data = data,
      method = method,
      trControl = fitControl,
      metric = metric,
      preProcess = preProcess,
      tuneGrid = tuneGrid,
      tuneLength = if (useTuneLength) tuneLength else NULL,
      weights = weights,
      verbose = FALSE
    )
    
    model <- tryCatch({
      .safe_train(base_args)
    }, error = function(e) {
      if (use_weights && grepl("weight", e$message, ignore.case = TRUE)) {
        message("  Model '", method, "' does not support weights. Retrying without weights.")
        no_wt_args <- base_args
        no_wt_args$weights <- NULL
        return(tryCatch(.safe_train(no_wt_args), error = function(e2) {
          message("  Model '", method, "' failed even without weights: ", e2$message)
          NULL
        }))
      } else {
        message("  Model '", method, "' failed: ", e$message)
        NULL
      }
    })
    
    if (!is.null(model)) {
      results[[method]] <- model
      cat("Done.\n")
    } else {
      cat("Skipped.\n")
    }
  }
  
  if (length(results) == 0) {
    stop("No models could be trained. Check installed packages, method/task-type ",
         "compatibility, and data format.")
  }
  
  return(results)
}



#' Evaluate Model Performance Metrics (Caret-based)
#'
#' Uses \code{\link[caret]{confusionMatrix}} to compute all classification
#' metrics and \code{\link[pROC]{roc}} for AUC. Automatically detects group
#' levels - no hard-coded \code{"0"} / \code{"1"}. Features robust defensive 
#' probability extraction and type coercion to ensure downstream stability.
#'
#' @param data          A data frame containing test features and the group column.
#' @param model_result  A single \code{train} object or a named list of
#'   \code{train} objects.
#' @param group_col     Name of the group column (default \code{"group"}).
#' @param custom_cutoff Optional numeric probability cutoff; default = 0.5.
#'
#' @return A data frame with one row per model, or \code{NULL} if all models fail.
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc auc
#' @importFrom stats predict
#' @export
#' @examples
#' \dontrun{
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' set.seed(123)
#' idx <- sample(1:nrow(mtcars), 20)
#' model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#' trained <- ModelTrainAnalysis(model, methods = c("rf"), 
#' control = list(method = "cv", number = 3), save_plots = FALSE)
#' perf <- evaluate_model_performance(trained@filtered.set$testing,
#'                                    trained@train.models$rf,
#'                                    group_col = "am")
#' print(perf)
#' }
evaluate_model_performance <- function(data,
                                       model_result,
                                       group_col   = "group",
                                       custom_cutoff = NULL) {
  
  # ---- Helper that processes a single caret model with robust probability extraction ----
  .eval_one <- function(model, model_name, data, group_col, custom_cutoff) {
    
    # 1. Ensure response variable factor levels are well-defined
    ref_levels <- levels(data[[group_col]])
    pos_class_col <- ref_levels[2] # Identify the positive class (second level)
    
    # 2. Predict probabilities from the model with error handling
    prob_matrix <- tryCatch(
      predict(model, newdata = data, type = "prob"),
      error = function(e) NULL
    )
    
    if (is.null(prob_matrix)) {
      warning("Could not predict probabilities for model '", model_name, "'")
      return(NULL)
    }
    
    # 3. Robust probability extraction with defensive validation and as.numeric coercion
    prob <- tryCatch({
      if (is.data.frame(prob_matrix) || is.matrix(prob_matrix)) {
        if (pos_class_col %in% colnames(prob_matrix)) {
          # Preferred method: match by explicit column name and ensure numeric
          as.numeric(prob_matrix[[pos_class_col]])
        } else if (ncol(prob_matrix) >= 2) {
          # Fallback method with a warning if column names are misaligned
          warning("Positive class '", pos_class_col, "' not found in prediction column names for model '", 
                  model_name, "'. Falling back to column index 2.")
          as.numeric(prob_matrix[, 2])
        } else {
          stop("Prediction probability matrix has fewer than 2 columns.")
        }
      } else {
        # Handle edge cases where models return a numeric vector directly
        as.numeric(prob_matrix)
      }
    }, error = function(e) {
      warning("Failed to extract probability vector for model '", model_name, "': ", e$message)
      NULL
    })
    
    if (is.null(prob)) {
      return(NULL)
    }
    
    cutoff <- if (!is.null(custom_cutoff)) custom_cutoff else 0.5
    pred_class <- factor(
      ifelse(prob > cutoff, levels(data[[group_col]])[2], levels(data[[group_col]])[1]),
      levels = levels(data[[group_col]])
    )
    true_class <- factor(data[[group_col]], levels = levels(data[[group_col]]))
    
    cm <- tryCatch(
      caret::confusionMatrix(pred_class, true_class, positive = levels(true_class)[2]),
      error = function(e) NULL
    )
    if (is.null(cm)) {
      warning("confusionMatrix failed for model '", model_name, "'")
      return(NULL)
    }
    
    roc_obj <- tryCatch(
      pROC::roc(true_class, prob, levels = levels(true_class), direction = "auto", quiet = TRUE),
      error = function(e) NULL
    )
    auc_val <- if (!is.null(roc_obj)) as.numeric(pROC::auc(roc_obj)) else NA_real_
    
    # Coerce all metrics to numeric explicitly to strip names and prevent rbind warnings
    data.frame(
      Model                     = model_name,
      Cutoff_Used               = cutoff,
      Sensitivity               = as.numeric(cm$byClass["Sensitivity"]),
      Specificity               = as.numeric(cm$byClass["Specificity"]),
      Positive_predictive_value = as.numeric(cm$byClass["Pos Pred Value"]),
      Negative_predictive_value = as.numeric(cm$byClass["Neg Pred Value"]),
      accuracy_score            = as.numeric(cm$overall["Accuracy"]),
      Precision                 = as.numeric(cm$byClass["Precision"]),
      recall_score              = as.numeric(cm$byClass["Recall"]),
      f1_score                  = as.numeric(cm$byClass["F1"]),
      auc                       = as.numeric(auc_val),
      row.names                 = NULL
    )
  }
  
  if (inherits(model_result, "train")) {
    return(.eval_one(model_result, model_result$method, data, group_col, custom_cutoff))
  } else if (is.list(model_result) && !inherits(model_result, "train")) {
    res_list <- lapply(names(model_result), function(nm) {
      .eval_one(model_result[[nm]], nm, data, group_col, custom_cutoff)
    })
    res_list <- res_list[!sapply(res_list, is.null)]
    if (length(res_list) == 0) return(NULL)
    return(do.call(rbind, res_list))
  } else {
    stop("model_result must be a caret train object or a named list of such objects.")
  }
}
#' Plot ROC Curves for Multiple Models
#' 
#' This function generates ROC curves with AUC and confidence intervals for a list of
#' binary classification models evaluated on validation data. It supports customizable
#' color palettes, plot saving, and data export.
#'
#' @param model_list A named list of model objects.
#' @param validation_data A data frame containing the validation dataset.
#' @param group_col A character string specifying the outcome column. Default `"group"`.
#' @param palette_name Wes Anderson palette name. Default `"AsteroidCity1"`.
#' @param base_size Base font size. Default `14`.
#' @param save_plots Logical. Save the ROC plot as PDF. Default `FALSE`.
#' @param save_dir Directory to save plots/data.
#' @param plot_width Width in inches. Default `5`.
#' @param plot_height Height in inches. Default `5`.
#' @param alpha Line transparency. Default `1`.
#' @param save_data Logical. Save ROC curve data as CSV. Default `FALSE`.
#' @return A list with `roc_objects`, `plot_data`, `auc_results`.
#' @importFrom pROC roc auc ci.auc coords
#' @importFrom ggplot2 ggplot geom_line geom_abline scale_color_manual labs scale_x_continuous scale_y_continuous theme element_rect element_text
#' @importFrom ggprism theme_prism
#' @importFrom wesanderson wes_palette
#' @importFrom viridis viridis
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("pROC", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   set.seed(123)
#'   idx <- sample(1:nrow(mtcars), 20)
#'   model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#'   trained <- ModelTrainAnalysis(model, methods = c("glm", "rf"),
#'                                 control = list(method = "cv", number = 3),
#'                                 save_plots = FALSE)
#'   if (interactive()) {
#'     plot_roc_curve(trained@train.models, 
#'     validation_data = trained@filtered.set$testing, 
#'     group_col = "am", save_plots = FALSE)
#'   }
#' }
#' }
plot_roc_curve <- function(model_list,
                           validation_data,
                           group_col = "group",
                           palette_name = "AsteroidCity1",
                           base_size = 14,
                           save_plots = FALSE,
                           save_dir = NULL,
                           plot_width = 5,
                           plot_height = 5,
                           alpha = 1,
                           save_data = FALSE) {

  roc_list <- list()
  plot_list <- list()
  auc_results <- numeric()

  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    predictions <- predict(model, validation_data, type = "prob")[, 2]

    actual_levels <- levels(factor(validation_data[[group_col]]))
    if (length(actual_levels) != 2) {
      stop("Group column must have exactly 2 levels. Found: ",
           paste(actual_levels, collapse = ", "))
    }

    roc_obj <- roc(validation_data[[group_col]], predictions,
                   levels = actual_levels,
                   direction = "auto")
    auc_value <- auc(roc_obj)
    roc_list[[model_name]] <- roc_obj
    auc_results[model_name] <- auc_value

    plot_data <- data.frame(
      Specificity = 1 - roc_obj$specificities,
      Sensitivity = roc_obj$sensitivities,
      Dataset = model_name
    )
    plot_list[[model_name]] <- plot_data
  }

  auc_results <- sort(auc_results, decreasing = TRUE)
  plot_list <- plot_list[names(auc_results)]

  combined_plot_data <- do.call(rbind, plot_list)

  dataset_levels <- paste0(names(auc_results),
                           " (AUC = ", round(auc_results, 3),
                           ", CI = [",
                           sapply(names(auc_results), function(nm) round(ci.auc(roc_list[[nm]])[1], 3)),
                           ", ",
                           sapply(names(auc_results), function(nm) round(ci.auc(roc_list[[nm]])[3], 3)),
                           "])")

  combined_plot_data$Dataset <- factor(combined_plot_data$Dataset, levels = names(auc_results))
  levels(combined_plot_data$Dataset) <- dataset_levels

  n_colors <- length(dataset_levels)
  palette_colors <- tryCatch({
    cols <- wesanderson::wes_palette(palette_name, type = "discrete")
    if (length(cols) < n_colors) {
      rep(cols, length.out = n_colors)
    } else {
      cols[1:n_colors]
    }
  }, error = function(e) {
    cat("Failed to use specified palette, falling back to viridis default colors.\n")
    viridis::viridis(n_colors)
  })

  p <- ggplot(combined_plot_data, aes(x = Specificity, y = Sensitivity, color = Dataset)) +
    geom_line(size = 1.25, alpha = alpha) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    scale_color_manual(values = palette_colors) +
    labs(title = "ROC Curves for Best Model on Training Data",
         subtitle = "Including AUC and 95% Confidence Intervals",
         x = "1 - Specificity",
         y = "Sensitivity",
         color = "Dataset (AUC and CI)") +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    ggprism::theme_prism(base_size = base_size) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", alpha = 0.8),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  print(p)

  if (save_plots) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggsave(filename = file.path(save_dir, "roc_curves.pdf"),
           plot = p, width = plot_width, height = plot_height, device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "roc_curves.pdf"), "\n")
  }

  if (save_data) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    csv_path <- file.path(save_dir, "roc_curve_data.csv")
    write.csv(combined_plot_data, csv_path, row.names = FALSE)
    cat("Plot data saved to:", csv_path, "\n")
  }

  print(auc_results)

  return(list(
    roc_objects = roc_list,
    plot_data = combined_plot_data,
    auc_results = auc_results
  ))
}


#' Check the Number of Levels in the Factor Variable
#'
#' @param data A data frame or tibble containing the dataset.
#' @param group_col The name of the column to check.
#' @returns TRUE if the column contains at least two levels, otherwise stops.
#' @export
check_factor_level <- function(data, group_col) {
  levels_present <- levels(as.factor(data[[group_col]]))

  if (length(levels_present) < 2) {
    stop("Validation data must contain at least two class levels.")
  }

  return(TRUE)
}

# ---------------------------------------------------------------------------
# Internal helper (not exported): a small, reusable CV-comparison bar chart
# used by ModelTrainAnalysis() when save_plots = TRUE, for ANY task_type
# (binary / multiclass / regression). Kept intentionally simple and in the
# same visual style (wesanderson palette + ggprism theme) as plot_roc_curve()
# so figures look consistent across the package.
# ---------------------------------------------------------------------------
#' Plot Cross-Validation Metric Comparison (internal)
#'
#' Creates a bar chart comparing a chosen cross-validation metric across
#' models, with error bars representing the standard deviation. Optionally
#' saves the plot (as PDF) and underlying data (as CSV) to disk.
#'
#' @param cv_metrics A data.frame containing at least the columns \code{Model},
#'   \code{CV_Mean}, and \code{CV_SD}.
#' @param metric_selection Character string naming the metric being plotted
#'   (used only for axis/title labeling).
#' @param palette_name Character string giving a \pkg{wesanderson} palette
#'   name. Falls back to \pkg{viridis} if the palette is unavailable or has
#'   fewer colors than needed.
#' @param base_size Base font size passed to \code{ggprism::theme_prism()}.
#' @param save_plots Logical. If \code{TRUE}, saves the plot and data to
#'   \code{save_dir}.
#' @param save_dir Character string. Directory to save outputs to. Required if
#'   \code{save_plots = TRUE}.
#' @param plot_width Numeric. Width (inches) for the saved PDF.
#' @param plot_height Numeric. Height (inches) for the saved PDF.
#' @param filename Character string. Filename for the saved PDF (a matching
#'   \code{.csv} is derived from it).
#'
#' @return Invisibly returns the \code{ggplot} object.
#'
#' @keywords internal
#' @noRd
.plot_cv_bar <- function(cv_metrics, metric_selection, palette_name = "AsteroidCity1",
                         base_size = 14, save_plots = FALSE, save_dir = NULL,
                         plot_width = 6, plot_height = 5, filename = "cv_comparison.pdf") {
  if (is.null(cv_metrics) || nrow(cv_metrics) == 0) {
    message("No CV metrics available to plot.")
    return(invisible(NULL))
  }

  cv_metrics$Model <- factor(cv_metrics$Model, levels = cv_metrics$Model)

  n_colors <- nrow(cv_metrics)
  palette_colors <- tryCatch({
    cols <- wesanderson::wes_palette(palette_name, type = "discrete")
    if (length(cols) < n_colors) rep(cols, length.out = n_colors) else cols[1:n_colors]
  }, error = function(e) viridis::viridis(n_colors))

  p <- ggplot2::ggplot(cv_metrics, ggplot2::aes(x = .data$Model, y = .data$CV_Mean, fill = .data$Model)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$CV_Mean - ifelse(is.na(.data$CV_SD), 0, .data$CV_SD),
                   ymax = .data$CV_Mean + ifelse(is.na(.data$CV_SD), 0, .data$CV_SD)),
      width = 0.2
    ) +
    ggplot2::scale_fill_manual(values = palette_colors) +
    ggplot2::labs(title = paste("Cross-Validation Comparison (", metric_selection, ")"),
                  x = NULL, y = metric_selection) +
    ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

  print(p)

  if (isTRUE(save_plots)) {
    if (is.null(save_dir)) stop("save_dir must be provided when save_plots = TRUE.")
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(filename = file.path(save_dir, filename),
                     plot = p, width = plot_width, height = plot_height, device = "pdf")
    write.csv(cv_metrics, file.path(save_dir, sub("\\.pdf$", ".csv", filename)), row.names = FALSE)
    cat("CV comparison plot + data saved to:", save_dir, "\n")
  }

  invisible(p)
}

#' Comprehensive Model Training and Analysis (Binary-Compatible Mode)
#'
#' Supports binary classification, multiclass classification, and regression.
#' Data is assumed to already be clean and analysis-ready.
#'
#' For binary classification, this function stores test-set performance as a
#' data frame in `@all.results`, identical to the original `ModelTrainAnalysis`,
#' so all downstream functions (SelectBestModel, plots, ensembles) work unchanged.
#' For other task types, `@all.results` is stored as a list; see
#' `model_05b_nonbinary_downstream.R` for a matching set of downstream
#' plotting/summary functions built specifically for that list shape (this
#' function's *output* is left exactly as-is -- multiclass/regression results
#' are still stored in the same list layout as before).
#'
#' @param object Train_Model object (must have filtered.set with training/testing)
#' @param methods Character vector of caret model methods
#' @param control List of trainControl parameters (passed to train_and_evaluate_models)
#' @param tune_grids List of tuning grids per model (if NULL, uses tuneLength)
#' @param tuneLength Integer, number of automatic tuning levels
#' @param preProcess Character vector of preprocessing methods. Default c("center","scale").
#'   Set to NULL if your data truly needs no transformation at all (e.g. it was
#'   already centered/scaled upstream -- passing a non-NULL preProcess here on
#'   top of already-scaled data will re-transform it a second time).
#' @param loocv_threshold Threshold for LOOCV
#' @param classProbs Logical. Ignored (forced FALSE) for regression.
#' @param allowParallel Logical
#' @param group_col Target column name (will be taken from object if possible)
#' @param task_type One of "auto" (default), "binary", "multiclass", "regression".
#' @param metric_selection Metric to select best model on CV. If NULL, a sensible
#'   default is chosen per task_type: "ROC" (binary), "Accuracy" (multiclass),
#'   "RMSE" (regression).
#' @param imbalance_handling Character: "none", "auto", "up", "down", "smote",
#'   "rose", or "weights". Passed to train_and_evaluate_models.
#' @param imbalance_threshold Numeric (0-1), passed to train_and_evaluate_models.
#' @param seed Random seed (set once at the start, and re-used per model
#'   inside train_and_evaluate_models() for reproducible, fold-shared CV).
#' @param save_plots Logical. If `TRUE`, generates and saves diagnostic plots
#'   to `save_dir`: for binary tasks, ROC curves for every trained model via
#'   [plot_roc_curve()] on the held-out test set; for multiclass/regression,
#'   a cross-validation comparison bar chart. CSVs of the underlying metrics
#'   are saved alongside the plots.
#' @param save_dir Directory to save plots/CSVs to. Required if `save_plots = TRUE`.
#' @param ... Additional arguments passed to train_and_evaluate_models
#' @return Updated Train_Model object (or list if input was list). When
#'   preprocessing is learned internally (`preProcess` non-NULL), the fitted
#'   `caret::preProcess` object is now also stored -- in
#'   `object@process.info$preprocess_model` for `Train_Model` input, or as
#'   `$preProcess_model` in the returned list otherwise -- so it can be
#'   reapplied later (e.g. to new/external validation data at deployment
#'   time) instead of being learned-and-discarded.
#' @importFrom caret preProcess train postResample confusionMatrix
#' @importFrom stats predict as.formula na.omit
#' @export
#' @examples
#' \dontrun{
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' # Split (simulated)
#' set.seed(123)
#' idx <- sample(1:nrow(mtcars), 20)
#' model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#' trained <- ModelTrainAnalysis(model,
#'                               methods = c("glm", "rf"),
#'                               control = list(method = "cv", number = 3),
#'                               save_plots = FALSE)
#'                               
#' }
ModelTrainAnalysis <- function(object,
                               methods = c("glmnet", "rf", "knn", "svmRadial", "rpart"),
                               control = list(method = "repeatedcv", number = 10, repeats = 3),
                               tune_grids = NULL,
                               tuneLength = 3,
                               preProcess = c("center", "scale"),
                               loocv_threshold = 100,
                               classProbs = TRUE,
                               allowParallel = FALSE,
                               group_col = NULL,
                               task_type = c("auto", "binary", "multiclass", "regression"),
                               metric_selection = NULL,
                               imbalance_handling = c("none", "auto", "up", "down", "smote", "rose", "weights"),
                               imbalance_threshold = 0.2,
                               seed = 123,
                               save_plots = FALSE,
                               save_dir = NULL,
                               ...) {
  
  task_type <- match.arg(task_type)
  imbalance_handling <- match.arg(imbalance_handling)
  set.seed(seed)
  
  if (isTRUE(save_plots) && is.null(save_dir)) {
    stop("save_dir must be provided when save_plots = TRUE.")
  }
  
  # Metrics direction
  higher_is_better <- c("ROC", "Accuracy", "Kappa", "Sensitivity", "Specificity",
                        "Precision", "F1", "AUC", "Rsquared")
  lower_is_better  <- c("RMSE", "MAE", "logLoss")
  
  # ---------- 1. Extract train/test sets ----------
  if (inherits(object, "Train_Model")) {
    data_sets <- Extract_filtered.set(object)
    train_data <- data_sets$training
    test_data  <- data_sets$testing
    if (is.null(group_col)) group_col <- object@group_col
  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    train_data <- object$train
    test_data  <- object$test
  } else {
    stop("Input must be Train_Model or list with 'train' and 'test' elements")
  }
  
  if (is.null(test_data) || nrow(test_data) == 0) {
    stop("Test data is missing. Cannot evaluate generalization.")
  }
  if (is.null(group_col) || !(group_col %in% names(train_data))) {
    stop("group_col ('", group_col, "') not found in training data.")
  }
  
  cat("Train samples:", nrow(train_data), " | Test samples:", nrow(test_data), "\n")
  
  # ---------- 2. Determine task_type ----------
  outcome_raw <- train_data[[group_col]]
  
  if (task_type == "auto") {
    n_unique <- length(unique(stats::na.omit(outcome_raw)))
    if (is.numeric(outcome_raw) && n_unique > 10) {
      task_type <- "regression"
    } else if (n_unique == 2) {
      task_type <- "binary"
    } else if (n_unique > 2) {
      task_type <- "multiclass"
    } else {
      warning("Could not confidently detect task type. Defaulting to 'binary'.")
      task_type <- "binary"
    }
  }
  cat(">>> Task type:", task_type, "\n")
  
  # ---------- 3. Prepare outcome ----------
  if (task_type %in% c("binary", "multiclass")) {
    train_data[[group_col]] <- factor(train_data[[group_col]])
    test_data[[group_col]]  <- factor(test_data[[group_col]], levels = levels(train_data[[group_col]]))
    
    n_lv <- nlevels(train_data[[group_col]])
    if (task_type == "binary" && n_lv != 2) {
      stop("task_type = 'binary' but outcome has ", n_lv, " levels.")
    }
    if (task_type == "multiclass" && n_lv < 3) {
      stop("task_type = 'multiclass' but outcome has only ", n_lv, " level(s).")
    }
  } else {
    train_data[[group_col]] <- as.numeric(train_data[[group_col]])
    test_data[[group_col]]  <- as.numeric(test_data[[group_col]])
    classProbs <- FALSE
  }
  
  # ---------- 4. Default metric_selection ----------
  if (is.null(metric_selection)) {
    metric_selection <- switch(task_type,
                               binary     = "Accuracy",  # more robust default; pass metric_selection = "ROC" explicitly to rank by AUC
                               multiclass = "Accuracy",
                               regression = "RMSE")
  }
  if (task_type == "multiclass" && metric_selection == "ROC") {
    warning("ROC not supported for multiclass; switching to 'Accuracy'.")
    metric_selection <- "Accuracy"
  }
  if (task_type == "regression" && metric_selection %in% c("ROC", "Accuracy", "Kappa")) {
    warning("metric_selection '", metric_selection, "' invalid for regression; using 'RMSE'.")
    metric_selection <- "RMSE"
  }
  cat(">>> Selection metric:", metric_selection, "\n")
  rank_decreasing <- !(metric_selection %in% lower_is_better)
  
  # ---------- 5. Preprocessing (learned on training; reusable core) ----------
  # Delegates to FitApplyPreprocess() (see that function's docs for why
  # imbalance/resampling is intentionally NOT done here) so the same logic
  # is available standalone via PreprocessModelData() for use outside this
  # function too.
  pp_result <- FitApplyPreprocess(
    training = train_data,
    testing  = test_data,
    group_col = group_col,
    method = preProcess
  )
  train_data_pp <- pp_result$training
  test_data_pp  <- pp_result$testing
  pp            <- pp_result$preProc   # <- fitted caret::preProcess object, if any
  
  # Also mirror the result into split.scale.data (kept alongside filtered.set
  # below), so the scaled data + fitted preprocessing model are available for
  # reuse from a single, predictable slot -- not just recomputed/discarded.
  if (inherits(object, "Train_Model")) {
    object@split.scale.data <- list(training = train_data_pp, testing = test_data_pp)
  }
  
  # ---------- 6. Train models ----------
  if (task_type == "multiclass") {
    # Same authoritative check as inside train_and_evaluate_models(): prefer
    # caret's own "Two Class Only" tag over the hand-curated fallback list,
    # so this preview message matches what will actually be skipped.
    is_two_class_only_method <- function(m) {
      minfo <- tryCatch(caret::getModelInfo(m, regex = FALSE)[[1]], error = function(e) NULL)
      (!is.null(minfo) && !is.null(minfo$tags) && "Two Class Only" %in% minfo$tags) ||
        m %in% BINARY_ONLY_METHODS
    }
    binary_only_requested <- Filter(is_two_class_only_method, methods)
    if (length(binary_only_requested) > 0) {
      message(">>> Note: ", paste(binary_only_requested, collapse = ", "),
              " are tagged 'Two Class Only' by caret and will be skipped ",
              "(outcome has ", n_lv, " classes). See train_and_evaluate_models().")
    }
  }
  
  cat("Training models using CV...\n")
  model_list <- train_and_evaluate_models(
    data = train_data_pp,
    methods = methods,
    control = control,
    tune_grids = tune_grids,
    tuneLength = tuneLength,
    preProcess = NULL,   # already manually preprocessed above
    loocv_threshold = loocv_threshold,
    classProbs = classProbs,
    allowParallel = allowParallel,
    group_col = group_col,
    metric = metric_selection,
    imbalance_handling = imbalance_handling,
    imbalance_threshold = imbalance_threshold,
    seed = seed          # honoured inside train_and_evaluate_models (shared folds, seeded once before they're built)
  )
  
  if (length(model_list) == 0) stop("No models trained successfully.")
  
  # ---------- 7. Extract CV metrics (robust matching) ----------
  cv_metrics <- data.frame()
  for (nm in names(model_list)) {
    model <- model_list[[nm]]
    if (!is.null(model$bestTune) && !is.null(model$results)) {
      best_row <- NULL
      for (i in seq_len(nrow(model$results))) {
        row_vals  <- model$results[i, colnames(model$bestTune), drop = FALSE]
        tune_vals <- model$bestTune[1, , drop = FALSE]
        # Use all.equal rather than == to avoid floating-point mismatches
        # between numeric tuning-parameter columns (e.g. 0.1 vs 0.10000001).
        if (isTRUE(all.equal(as.list(row_vals), as.list(tune_vals), check.attributes = FALSE))) {
          best_row <- model$results[i, ]
          break
        }
      }
      if (!is.null(best_row) && metric_selection %in% colnames(best_row)) {
        mean_val <- best_row[[metric_selection]]
        sd_col <- paste0(metric_selection, "SD")
        sd_val <- if (sd_col %in% colnames(best_row)) best_row[[sd_col]] else NA_real_
        cv_metrics <- rbind(cv_metrics, data.frame(Model = nm, CV_Mean = mean_val, CV_SD = sd_val))
      }
    } else if (!is.null(model$resample) && metric_selection %in% colnames(model$resample)) {
      mean_val <- mean(model$resample[[metric_selection]], na.rm = TRUE)
      sd_val   <- sd(model$resample[[metric_selection]], na.rm = TRUE)
      cv_metrics <- rbind(cv_metrics, data.frame(Model = nm, CV_Mean = mean_val, CV_SD = sd_val))
    }
  }
  
  if (nrow(cv_metrics) == 0) stop("No CV metrics available for selection.")
  cv_metrics <- cv_metrics[order(cv_metrics$CV_Mean, decreasing = rank_decreasing), ]
  cat("\n========== Cross-Validation Performance ==========\n")
  print(cv_metrics)
  
  best_model_name <- cv_metrics$Model[1]
  best_model <- model_list[[best_model_name]]
  
  # ---------- 8. Final evaluation (binary uses data.frame, others list) ----------
  .evaluate_final <- function(data, model_result, group_col, task_type) {
    if (task_type == "binary") {
      return(evaluate_model_performance(data = data,
                                        model_result = model_result,
                                        group_col = group_col))
    } else if (task_type == "regression") {
      preds <- as.numeric(stats::predict(model_result, newdata = data))
      truth <- data[[group_col]]
      perf <- caret::postResample(pred = preds, obs = truth)
      return(list(RMSE = unname(perf["RMSE"]),
                  Rsquared = unname(perf["Rsquared"]),
                  MAE = unname(perf["MAE"])))
    } else { # multiclass
      pred_class <- stats::predict(model_result, newdata = data)
      truth <- data[[group_col]]
      cm <- caret::confusionMatrix(pred_class, truth)
      result <- list(accuracy = unname(cm$overall["Accuracy"]),
                     kappa = unname(cm$overall["Kappa"]))
      probs <- tryCatch(stats::predict(model_result, newdata = data, type = "prob"),
                        error = function(e) NULL)
      if (!is.null(probs) && requireNamespace("pROC", quietly = TRUE)) {
        roc_multi <- tryCatch(pROC::multiclass.roc(truth, probs),
                              error = function(e) NULL)
        if (!is.null(roc_multi)) result$auc_macro <- as.numeric(roc_multi$auc)
      }
      return(result)
    }
  }
  
  # ---- For binary: evaluate all models on test set to produce data frame ----
  if (task_type == "binary") {
    test_perf_df <- evaluate_model_performance(
      data = test_data_pp,
      model_result = model_list,
      group_col = group_col
    )
    if (is.null(test_perf_df)) stop("Test performance evaluation failed.")
    test_perf_df <- test_perf_df[order(test_perf_df$auc, decreasing = TRUE), ]
    
    train_perf_df <- evaluate_model_performance(
      data = train_data_pp,
      model_result = model_list,
      group_col = group_col
    )
    
    best_test_perf <- test_perf_df[test_perf_df$Model == best_model_name, ]
    best_train_perf <- if (!is.null(train_perf_df)) {
      train_perf_df[train_perf_df$Model == best_model_name, ]
    } else NULL
    
    # ---------- 9. Print diagnostics ----------
    cat("\n========== Final Evaluation on Holdout Test Set ==========\n")
    if (!is.null(best_test_perf) && nrow(best_test_perf) > 0) {
      cat("Test ROC     :", round(best_test_perf$auc, 4), "\n")
      cat("Test Accuracy:", round(best_test_perf$accuracy_score, 4), "\n")
    }
    if (!is.null(best_train_perf) && nrow(best_train_perf) > 0) {
      gap <- best_train_perf$auc - best_test_perf$auc
      cat("Train AUC    :", round(best_train_perf$auc, 4), "\n")
      cat("Gap          :", round(gap, 4))
      if (gap > 0.1) cat(" (High! Consider simpler models)\n") else cat(" (Acceptable)\n")
    }
    
    # ---------- 9b. save_plots (now fully wired) ----------
    if (isTRUE(save_plots)) {
      cat("\nSaving diagnostic plots to:", save_dir, "\n")
      plot_roc_curve(
        model_list = model_list,
        validation_data = test_data_pp,
        group_col = group_col,
        save_plots = TRUE,
        save_dir = save_dir,
        save_data = TRUE
      )
      .plot_cv_bar(cv_metrics, metric_selection, save_plots = TRUE, save_dir = save_dir,
                   filename = "cv_comparison.pdf")
      write.csv(test_perf_df, file.path(save_dir, "test_performance.csv"), row.names = FALSE)
      write.csv(train_perf_df, file.path(save_dir, "train_performance.csv"), row.names = FALSE)
    }
    
    # ---------- 10. Store in Train_Model (binary: data.frame) ----------
    if (inherits(object, "Train_Model")) {
      object@filtered.set <- list(training = train_data_pp, testing = test_data_pp)
      object@train.models <- model_list
      object@all.results  <- test_perf_df   # data frame
      object@process.info$cv_metrics <- cv_metrics
      object@process.info$task_type   <- task_type
      object@process.info$preProcess  <- preProcess
      object@process.info$preprocess_model <- pp   # <- fitted preProcess object (NULL if none learned)
      object@process.info$imbalance_handling <- imbalance_handling
      cat("\n>>> Train_Model object updated (binary mode).\n")
      return(object)
    } else {
      return(list(
        task_type = task_type,
        models = model_list,
        cv_performance = cv_metrics,
        test_performance = test_perf_df,
        train_performance = train_perf_df,
        preProcess = preProcess,
        preProcess_model = pp,
        imbalance_handling = imbalance_handling
      ))
    }
  } else {
    # ---- Non-binary: evaluate only the best model (list) ----
    final_test_perf <- .evaluate_final(test_data_pp, best_model, group_col, task_type)
    final_train_perf <- .evaluate_final(train_data_pp, best_model, group_col, task_type)
    
    cat("\n========== Final Evaluation on Holdout Test Set ==========\n")
    print(final_test_perf)
    
    # ---------- save_plots (now fully wired) for multiclass/regression ----------
    if (isTRUE(save_plots)) {
      cat("\nSaving diagnostic plots to:", save_dir, "\n")
      .plot_cv_bar(cv_metrics, metric_selection, save_plots = TRUE, save_dir = save_dir,
                   filename = "cv_comparison.pdf")
      # A lightweight, task-agnostic performance dump; richer, task-specific
      # plots (confusion-matrix heatmaps, per-class ROC, predicted-vs-actual)
      # live in model_05b_nonbinary_downstream.R and should be called
      # separately on the returned/updated object.
      perf_df <- data.frame(
        Set = c("train", "test"),
        t(rbind(unlist(final_train_perf), unlist(final_test_perf)))
      )
      write.csv(perf_df, file.path(save_dir, paste0(task_type, "_best_model_performance.csv")),
                row.names = FALSE)
    }
    
    if (inherits(object, "Train_Model")) {
      object@filtered.set <- list(training = train_data_pp, testing = test_data_pp)
      object@train.models <- model_list
      object@all.results <- list(
        task_type = task_type,
        cv_selection = cv_metrics,
        train_final = final_train_perf,
        test_final  = final_test_perf,
        preProcess  = preProcess,
        imbalance_handling = imbalance_handling
      )
      object@best.model.result <- list(
        model = best_model,
        model_type = best_model_name,
        task_type = task_type,
        cv_metric = cv_metrics$CV_Mean[1],
        test_performance = final_test_perf,
        train_performance = final_train_perf,
        selection_basis = paste0("Cross-Validation (", metric_selection, ")")
      )
      object@process.info$preprocess_model <- pp   # <- fitted preProcess object (NULL if none learned)
      cat("\n>>> Train_Model object updated (non-binary mode).\n")
      return(object)
    } else {
      return(list(
        task_type = task_type,
        models = model_list,
        cv_performance = cv_metrics,
        best_model = best_model,
        test_performance = final_test_perf,
        train_performance = final_train_perf,
        preProcess = preProcess,
        preProcess_model = pp,
        imbalance_handling = imbalance_handling
      ))
    }
  }
}

#' Plot ROC Curves for Best Model Across Multiple Datasets
#'
#' Generates and visualizes Receiver Operating Characteristic (ROC) curves 
#' with Area Under the Curve (AUC) and confidence intervals for a trained 
#' classification model evaluated across multiple datasets (training, testing, 
#' validation, and external validation).
#'
#' @param best_model A trained classification model object with a valid `predict` method.
#' @param train_data Data frame containing training data with the response variable 
#'   (required to establish authoritative factor levels). Default `NULL`.
#' @param test_data Data frame containing testing data with the response variable. Default `NULL`.
#' @param validation_data Data frame containing validation data with the response variable. Default `NULL`.
#' @param external_validation Data frame containing external validation data with the response variable. Default `NULL`.
#' @param group_col Character string specifying the response variable column name. Default `"group"`.
#' @param palette_name Character string specifying a Wes Anderson palette name. Default `"AsteroidCity1"`.
#' @param base_size Numeric base font size for the plot theme. Default `14`.
#' @param save_plots Logical indicating whether to save the plot as a PDF. Default `FALSE`.
#' @param save_dir Character string specifying the directory to save outputs. Default `NULL`.
#' @param plot_width Numeric plot width in inches. Default `5`.
#' @param plot_height Numeric plot height in inches. Default `5`.
#' @param alpha Numeric transparency value for ROC curve lines. Default `1`.
#' @param subtitle Character string for the plot subtitle. Default `"Training and Testing Datasets"`.
#' @param save_data Logical indicating whether to save ROC curve data as a CSV. Default `FALSE`.
#' @param csv_filename Character string for the saved data filename. Default `"best_model_roc_data.csv"`.
#' 
#' @return A named list of data frames containing ROC coordinates and performance metrics per dataset.
#' 
#' @importFrom pROC roc auc ci.auc
#' @importFrom ggplot2 ggplot geom_line geom_abline scale_color_manual labs scale_x_continuous scale_y_continuous theme element_rect element_text expansion ggsave aes
#' @importFrom ggprism theme_prism
#' @importFrom wesanderson wes_palette
#' @importFrom viridis viridis
#' @export
#' @examples
#' \dontrun{
#' if (requireNamespace("pROC", quietly = TRUE)) {
#'   mtcars$am <- as.factor(mtcars$am)
#'   model <- CreateModelObject(data = mtcars, group_col = "am")
#'   set.seed(123)
#'   idx <- sample(1:nrow(mtcars), 20)
#'   model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#'   trained <- ModelTrainAnalysis(model, methods = c("glm"), 
#'   control = list(method = "cv", number = 3), save_plots = FALSE)
#'   best_model <- trained@train.models$glm
#'   if (interactive()) {
#'     plot_best_model_roc(best_model,
#'     train_data = trained@filtered.set$training,
#'     test_data = trained@filtered.set$testing,
#'     group_col = "am",
#'     save_plots = FALSE)
#'   }
#' }
#' }
plot_best_model_roc <- function(best_model,
                                train_data = NULL,
                                test_data = NULL,
                                validation_data = NULL,
                                external_validation = NULL,
                                group_col = "group",
                                palette_name = "AsteroidCity1",
                                base_size = 14,
                                save_plots = FALSE,
                                save_dir = NULL,
                                plot_width = 5,
                                plot_height = 5,
                                alpha = 1,
                                subtitle = "Training and Testing Datasets",
                                save_data = FALSE,
                                csv_filename = "best_model_roc_data.csv") {
  
  # 1. Validation: train_data must be provided to set a strict factor level benchmark
  if (is.null(train_data)) {
    stop("`train_data` must be provided to establish authoritative factor levels.")
  }
  
  # Establish authoritative factor levels from training data
  train_data[[group_col]] <- factor(train_data[[group_col]])
  ref_levels <- levels(train_data[[group_col]])
  
  if (length(ref_levels) != 2) {
    stop("Group column must have exactly 2 levels for ROC analysis. Found: ", 
         length(ref_levels), " (levels: ", paste(ref_levels, collapse = ", "), ")")
  }
  
  plot_data_list <- list()
  roc_objects <- list()
  
  # 2. Internal helper function to process and evaluate each dataset robustly
  .process_dataset <- function(df, dataset_name) {
    if (is.null(df)) return(NULL)
    
    # Remove missing values and align factor levels with the training baseline
    df <- stats::na.omit(df)
    df[[group_col]] <- factor(df[[group_col]], levels = ref_levels)
    
    if (any(is.na(df[[group_col]]))) {
      stop("Dataset '", dataset_name, "' contains values not present in training data levels: ", 
           paste(ref_levels, collapse = ", "))
    }
    
    # Predict probabilities for the positive class (second level)
    prob_matrix <- tryCatch(
      predict(best_model, newdata = df, type = "prob"),
      error = function(e) stop("Failed to predict probabilities for ", dataset_name, ": ", e$message)
    )
    
    pos_class_col <- ref_levels[2]
    preds <- if (pos_class_col %in% colnames(prob_matrix)) {
      prob_matrix[[pos_class_col]]
    } else {
      prob_matrix[, 2] # Fallback to second column if column name matching fails
    }
    
    truth <- df[[group_col]]
    
    # Compute ROC curve using direction = "auto" to prevent manual direction inversion errors
    roc_obj <- pROC::roc(truth, preds, levels = ref_levels, direction = "auto", quiet = TRUE)
    
    auc_val <- pROC::auc(roc_obj)
    auc_ci <- pROC::ci.auc(roc_obj)
    ci_half_width <- (auc_ci[3] - auc_ci[1]) / 2
    
    # Format label including AUC and 95% confidence interval half-width
    dataset_label <- paste0(dataset_name, " (AUC = ", sprintf("%.3f", auc_val),
                            "+/-", sprintf("%.3f", ci_half_width), ")")
    
    plot_df <- data.frame(
      Specificity = 1 - roc_obj$specificities,
      Sensitivity = roc_obj$sensitivities,
      Dataset = dataset_label
    )
    
    list(plot_data = plot_df, roc_obj = roc_obj)
  }
  
  # 3. Process all available datasets sequentially
  datasets_to_process <- list(
    "Training Set" = train_data,
    "Testing Set" = test_data,
    "Validation Set" = validation_data,
    "External Validation Set" = external_validation
  )

for (name in names(datasets_to_process)) {
  res <- .process_dataset(datasets_to_process[[name]], name)
  if (!is.null(res)) {
    plot_data_list[[name]] <- res$plot_data
    roc_objects[[name]] <- res$roc_obj
  }
}

if (length(plot_data_list) == 0) {
  stop("No datasets provided or successfully processed for ROC plotting.")
}

# Combine plot data and enforce correct factor level ordering
combined_plot_data <- do.call(rbind, plot_data_list)
dataset_order_labels <- unlist(lapply(plot_data_list, function(x) unique(x$Dataset)))
combined_plot_data$Dataset <- factor(combined_plot_data$Dataset, levels = dataset_order_labels)

# 4. Generate color palette with fallback to viridis
palette_colors <- tryCatch({
  cols <- wesanderson::wes_palette(name = palette_name, type = "discrete")
  n_cols <- length(unique(combined_plot_data$Dataset))
  if (length(cols) < n_cols) {
    rep(cols, length.out = n_cols)
  } else {
    cols[1:n_cols]
  }
}, error = function(e) {
  viridis::viridis(length(unique(combined_plot_data$Dataset)))
})

# 5. Build ggplot object
p <- ggplot2::ggplot(combined_plot_data, ggplot2::aes(x = Specificity, y = Sensitivity, color = Dataset)) +
  ggplot2::geom_line(size = 1.25, alpha = alpha) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  ggplot2::scale_color_manual(values = palette_colors) +
  ggplot2::labs(
    title = "ROC Curves Validation Comparison",
    subtitle = subtitle,
    x = "1 - Specificity",
    y = "Sensitivity",
    color = "Validation Cohort"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = 0.01)
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = 0.01)
  ) +
  ggprism::theme_prism(base_size = base_size) +
  ggplot2::theme(
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0.8)),
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8),
    plot.title = ggplot2::element_text(hjust = 0.5)
  )

print(p)

# 6. Save outputs if requested
if (save_plots) {
  if (is.null(save_dir)) stop("save_dir must be provided when save_plots = TRUE.")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  ggplot2::ggsave(filename = file.path(save_dir, "best_model_roc_plot.pdf"),
                  plot = p, width = plot_width, height = plot_height, device = "pdf")
  cat("Plot saved to:", file.path(save_dir, "best_model_roc_plot.pdf"), "\n")
}

if (save_data) {
  if (is.null(save_dir)) stop("save_dir must be provided when save_data = TRUE.")
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  data_path <- file.path(save_dir, csv_filename)
  write.csv(combined_plot_data, data_path, row.names = FALSE)
  cat("ROC curve data saved to: ", data_path, "\n")
}

return(plot_data_list)
}


#' Select Best Performing Model (Fixed)
#'
#' NOTE: this assumes `object@all.results` is a data.frame with a `Model`
#' column and the requested `metric` column, which is only true for BINARY
#' task_type. Do not call this on a `Train_Model` object produced by
#' `ModelTrainAnalysis()` in multiclass/regression mode -- that mode already
#' selects and stores the best model internally in `@best.model.result`; use
#' the downstream functions in `model_05b_nonbinary_downstream.R` instead.
#'
#' @param object A Train_Model object.
#' @param metric Primary metric (default "auc").
#' @param custom_selection Optional model name to force selection.
#' @return Updated Train_Model object.
#' @export
#' @examples
#' \dontrun{
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' set.seed(123)
#' idx <- sample(1:nrow(mtcars), 20)
#' model@filtered.set <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#' trained <- ModelTrainAnalysis(model, methods = c("glm", "rf"),
#'                               control = list(method = "cv", number = 3),
#'                               task_type = "binary", metric_selection = "ROC",
#'                               save_plots = FALSE)
#' trained <- SelectBestModel(trained, metric = "auc")
#' print(trained@best.model.result$model_type)
#' }
SelectBestModel <- function(object,
                            metric = "auc",
                            custom_selection = NULL) {

  if (!inherits(object, "Train_Model")) {
    stop("Input must be an object of class 'Train_Model'")
  }

  all.results <- slot(object, "all.results")
  model_list  <- slot(object, "train.models")

  if (!is.data.frame(all.results)) {
    stop("object@all.results is not a data.frame (this happens for multiclass/",
         "regression results from ModelTrainAnalysis()). SelectBestModel() only ",
         "supports binary task_type; the best model for other task types is ",
         "already stored in object@best.model.result. See ",
         "model_05b_nonbinary_downstream.R for matching downstream tools.")
  }

  available_metrics <- intersect(c("auc", "Sensitivity", "Specificity",
                                   "accuracy_score", "f1_score", "Precision",
                                   "recall_score", "Positive_predictive_value",
                                   "Negative_predictive_value"),
                                 colnames(all.results))

  if (!metric %in% available_metrics) {
    stop("Metric '", metric, "' not available.\n",
         "Available: ", paste(available_metrics, collapse = ", "))
  }

  if (!is.null(custom_selection)) {
    if (!custom_selection %in% all.results$Model) {
      stop("Model '", custom_selection, "' not found in trained models")
    }
    cat("User selected model: ", custom_selection, "\n")
    best.model.result <- all.results[all.results$Model == custom_selection, ]
    best_model_type   <- custom_selection
  } else {

    primary_vals <- all.results[[metric]]
    best_val <- max(primary_vals, na.rm = TRUE)

    candidates <- all.results[which(primary_vals == best_val), ]

    if (nrow(candidates) > 1 && "f1_score" %in% colnames(candidates)) {
      f1_vals <- candidates$f1_score
      f1_vals[is.na(f1_vals)] <- -Inf
      best_f1 <- max(f1_vals)
      candidates <- candidates[which(f1_vals == best_f1), ]
    }

    best.model.result <- candidates[1, , drop = FALSE]
    best_model_type   <- best.model.result$Model

    if (nrow(candidates) > 1) {
      cat("Multiple models with same metrics, using first one: ", best_model_type, "\n")
    } else {
      cat("Best model (", metric, "): ", best_model_type, "\n")
    }
  }

  best_model <- model_list[[best_model_type]]
  if (is.null(best_model)) {
    stop("Model '", best_model_type, "' not found in trained models")
  }

  object@best.model.result <- list(
    model              = best_model,
    model_type         = best_model_type,
    train_performance  = best.model.result
  )

  return(object)
}

#' Fit and Apply Feature Preprocessing (Reusable Core)
#'
#' Learns a \code{caret::preProcess} model on the TRAINING data only (never
#' on testing/validation data, to avoid leakage) and applies it to both
#' training and testing sets. This is the low-level, reusable building block
#' behind \code{PreprocessModelData()} and \code{ModelTrainAnalysis()}: given
#' any training/testing pair it returns the transformed data AND the fitted
#' model, so the exact same transformation can be re-applied later to new
#' data (e.g. an external validation cohort) via
#' \code{stats::predict(preProc_model, newdata)}.
#'
#' This function deliberately handles ONLY feature-level preprocessing
#' (centering, scaling, transformation, imputation, etc. -- anything
#' \code{caret::preProcess} supports). It does NOT perform class-imbalance
#' resampling (up/down/SMOTE/ROSE). Imbalance correction is intentionally
#' kept inside \code{train_and_evaluate_models()}'s CV loop (via
#' \code{trainControl(sampling = ...)}) rather than here, because resampling
#' rows for class balance -- unlike centering/scaling -- changes WHICH rows
#' exist. Doing that once on the whole training set before cross-validation
#' leaks synthetic/duplicated rows across folds and inflates CV performance
#' estimates; doing it inside each fold (as caret's `sampling` argument does)
#' keeps every held-out fold untouched by the resampling. See Details in
#' \code{train_and_evaluate_models()}.
#'
#' @param training Training data frame (features + outcome column).
#' @param testing Testing data frame with the same columns as `training`.
#' @param group_col Name of the outcome column to exclude from preprocessing.
#' @param method Character vector of preProcess methods, e.g.
#'   \code{c("center", "scale", "YeoJohnson", "bagImpute")}. Set to `NULL` or
#'   `character(0)` to skip preprocessing entirely (data is returned as-is).
#' @param columns Optional character vector of feature columns to preprocess.
#'   Defaults to all numeric feature columns (non-numeric/factor columns are
#'   left untouched, since methods like center/scale/YeoJohnson only apply to
#'   numeric data).
#' @return A list with:
#'   \item{training}{Preprocessed training data frame.}
#'   \item{testing}{Preprocessed testing data frame.}
#'   \item{preProc}{The fitted \code{caret::preProcess} object (or `NULL` if
#'     `method` was empty), reusable via \code{predict(preProc, newdata)}.}
#'   \item{columns}{The feature columns that were preprocessed.}
#' @importFrom caret preProcess
#' @importFrom stats predict
#' @export
FitApplyPreprocess <- function(training, testing, group_col,
                               method = c("center", "scale", "YeoJohnson", "bagImpute"),
                               columns = NULL) {
  if (!is.data.frame(training) || !is.data.frame(testing)) {
    stop("`training` and `testing` must both be data frames.")
  }
  if (is.null(group_col) || !group_col %in% colnames(training)) {
    stop("group_col ('", group_col, "') not found in `training`.")
  }
  if (!group_col %in% colnames(testing)) {
    stop("group_col ('", group_col, "') not found in `testing`.")
  }
  
  feature_cols <- setdiff(colnames(training), group_col)
  if (is.null(columns)) {
    columns <- feature_cols[vapply(training[feature_cols], is.numeric, logical(1))]
  } else {
    missing_cols <- setdiff(columns, feature_cols)
    if (length(missing_cols) > 0) {
      stop("`columns` not found among training feature columns: ",
           paste(missing_cols, collapse = ", "))
    }
  }
  
  if (is.null(method) || length(method) == 0 || length(columns) == 0) {
    if (length(columns) == 0 && !is.null(method) && length(method) > 0) {
      message("No numeric feature columns found to preprocess; returning data unchanged.")
    }
    return(list(training = training, testing = testing, preProc = NULL, columns = columns))
  }
  
  cat("Learning preprocessing (", paste(method, collapse = ", "), ") from training data (",
      length(columns), "feature(s))...\n")
  preProc <- caret::preProcess(training[, columns, drop = FALSE], method = method)
  
  training[, columns] <- stats::predict(preProc, training[, columns, drop = FALSE])
  testing[, columns]  <- stats::predict(preProc, testing[, columns, drop = FALSE])
  cat("Preprocessing applied to", length(columns), "feature column(s).\n")
  
  list(training = training, testing = testing, preProc = preProc, columns = columns)
}

#' Preprocess a Train_Model Object's Split Data (Reusable Pipeline Step)
#'
#' Convenience wrapper around \code{FitApplyPreprocess()} that operates
#' directly on a \code{Train_Model} object's \code{split.data} slot (or a
#' plain list with \code{training}/\code{testing} elements), and writes the
#' result to \code{split.scale.data} (and, mirroring the rest of this
#' pipeline, copies it into \code{filtered.set} so it is immediately usable
#' by \code{ModelTrainAnalysis()} without a separate feature-selection step).
#' The fitted preprocessing model is stored in
#' \code{process.info$preprocess_model} so it can be reused later --
#' e.g. reapplied to an external validation set at deployment time via
#' \code{predict(object@process.info$preprocess_model, newdata)} -- instead
#' of being re-fit (and potentially fit on the wrong data) every time.
#'
#' @param model_obj A `Train_Model` object with a populated `split.data` slot,
#'   or a list with `training`/`testing` data frames.
#' @param method Character vector of preProcess methods (see
#'   \code{caret::preProcess}). Default \code{c("center","scale","YeoJohnson","bagImpute")}.
#' @param group_col Outcome column name. If `NULL` and `model_obj` is a
#'   `Train_Model`, taken from `model_obj@group_col`.
#' @param columns Optional character vector restricting which feature columns
#'   are preprocessed (default: all numeric feature columns).
#' @return The updated `model_obj` (same input class), with
#'   `split.scale.data` and `filtered.set` populated and the fitted
#'   preprocessing model retained for reuse.
#' @export
#' @examples
#' \dontrun{
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' # Fake split (just for demo)
#' model@split.data <- list(training = mtcars[1:20, ], testing = mtcars[21:32, ])
#' model <- PreprocessModelData(model, method = c("center", "scale"))
#' print(dim(model@split.scale.data$training))
#' }
PreprocessModelData <- function(model_obj,
                                method = c("center", "scale", "YeoJohnson", "bagImpute"),
                                group_col = NULL,
                                columns = NULL) {
  is_S4 <- inherits(model_obj, "Train_Model")
  
  if (is_S4) {
    if (length(model_obj@split.data) == 0 ||
        !all(c("training", "testing") %in% names(model_obj@split.data))) {
      stop("model_obj@split.data must contain 'training' and 'testing'. ",
           "Run the data-splitting step first.")
    }
    training <- model_obj@split.data$training
    testing  <- model_obj@split.data$testing
    if (is.null(group_col)) group_col <- model_obj@group_col
  } else if (is.list(model_obj) && all(c("training", "testing") %in% names(model_obj))) {
    training <- model_obj$training
    testing  <- model_obj$testing
  } else {
    stop("model_obj must be a Train_Model object (with split.data) or a ",
         "list with 'training'/'testing' elements.")
  }
  
  res <- FitApplyPreprocess(training = training, testing = testing,
                            group_col = group_col, method = method, columns = columns)
  
  scaled_list <- list(training = res$training, testing = res$testing)
  
  if (is_S4) {
    model_obj@split.scale.data <- scaled_list
    model_obj@filtered.set     <- model_obj@split.scale.data
    model_obj@process.info$preprocess_model  <- res$preProc
    model_obj@process.info$preprocess_method <- method
    model_obj@process.info$preprocess_columns <- res$columns
    return(model_obj)
  } else {
    return(list(
      split.scale.data = scaled_list,
      filtered.set = scaled_list,
      preprocess_model = res$preProc,
      preprocess_method = method,
      preprocess_columns = res$columns
    ))
  }
}
