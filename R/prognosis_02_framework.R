#' Extract Survival Task from PrognosiX or TaskSurv Object
#'
#' This helper function extracts or returns a \code{TaskSurv} object from either
#' a \code{PrognosiX} S4 object or a \code{TaskSurv} object. It serves as a
#' unified interface for downstream functions that require a survival task.
#'
#' @param object An object of class \code{PrognosiX} or \code{TaskSurv}.
#'   If \code{PrognosiX}, the function extracts the survival data, time column,
#'   and event column to create a \code{TaskSurv} object.
#'
#' @details
#' When \code{object} is already a \code{TaskSurv}, this function also performs
#' a defensive sanity check on \code{object$row_roles$use}. Row indices produced
#' by functions such as \code{caret::createDataPartition(..., list = FALSE)} are
#' returned as a single-column matrix rather than a plain integer vector. If such
#' a matrix is passed directly to \code{TaskSurv$filter()}, mlr3's internal
#' resampling machinery silently misnames its row-id column (e.g.
#' \code{row_id.Resample1} instead of \code{row_id}), which later surfaces as an
#' opaque \code{"column not found: [row_id]"} error deep inside
#' \code{mlr3::Resampling$instantiate()}. To prevent this, \code{surv_extract_task}
#' coerces any matrix-typed \code{row_roles$use} back to a plain integer vector
#' before returning the task, so all downstream benchmarking, tuning, and
#' evaluation functions receive a well-formed task.
#'
#' @return A \code{\link[mlr3proba]{TaskSurv}} object with sanitized row roles.
#'
#' @seealso \code{\link{surv_create_surv_task}} for creating tasks from data frames
#' @export
#'
#' @examples
#' veteran <- survival::veteran
#' stat <- CreateStatObject(raw.data = veteran, clean.data = veteran,
#'                          group_col = "status", na.action = "allow")
#' prog <- Stat_to_PrognosiX(stat, "time", "status", na_action = "omit",
#'                           min_events = 10, verbose = FALSE)
#'
#' # From a PrognosiX object:
#' task <- surv_extract_task(prog)
#' task$nrow
#'
#' # Directly from a data frame:
#' task2 <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#' task2$nrow
surv_extract_task <- function(object) {
  .check_prognosis_packages()
  if (inherits(object, 'PrognosiX')) {
    data <- object@survival.data
    time_col <- object@time_col
    event_col <- object@status_col
    return(surv_create_surv_task(data, time_col, event_col))
  } else if (inherits(object, 'TaskSurv')) {
    if (is.matrix(object$row_roles$use)) {
      object$row_roles$use <- as.integer(as.vector(object$row_roles$use))
    }
    return(object)
  }
  stop('Input must be a PrognosiX object or TaskSurv')
}

# Internal helper: prepare binned predicted-vs-observed calibration data.
#
# @details Bin count is adaptively capped so each bin retains a statistically
#   meaningful sample size (default target: >= 10 observations/bin). Without
#   this guard, small task subsets (e.g. a held-out validation split) combined
#   with a fixed n_bins can produce bins with only 2-3 observations; the
#   per-bin KM estimate at time_point then frequently collapses to exactly 0
#   or 1, producing a jagged, misleading calibration curve driven by binning
#   noise rather than genuine miscalibration.
.prepare_cal_data <- function(learner, object, time_point, n_bins) {
  if (inherits(object, "TaskSurv")) {
    task <- object
  } else if (inherits(object, "PrognosiX")) {
    task <- surv_extract_task(object)
  } else {
    stop("object must be a TaskSurv or PrognosiX object")
  }
  
  n_bins <- max(n_bins, 5L)
  # Adaptive bin count: cap n_bins so each bin has a statistically meaningful
  # number of observations. With too few points per bin, the per-bin KM
  # estimate at time_point becomes unstable (often collapsing to exactly 0
  # or 1), producing a jagged, misleading calibration curve rather than a
  # true miscalibration signal -- this is especially visible when comparing
  # a small validation set against a larger training set.
  min_per_bin <- 10L
  adaptive_n_bins <- min(n_bins, max(3L, floor(task$nrow / min_per_bin)))
  if (adaptive_n_bins < n_bins) {
    message(sprintf(
      "[Calibration] Reduced n_bins from %d to %d (task has %d rows; targeting >= %d obs/bin for stable per-bin KM estimates).",
      n_bins, adaptive_n_bins, task$nrow, min_per_bin
    ))
  }
  n_bins <- adaptive_n_bins
  
  learner$predict_type <- "distr"
  pred <- learner$predict(task)
  
  surv_prob <- .extract_surv_prob(pred$distr, time_point, task)
  if (is.null(surv_prob)) return(NULL)
  surv_prob <- as.numeric(surv_prob)
  
  if (length(surv_prob) != task$nrow) {
    warning("Length of survival probabilities does not match task rows. Calibration skipped.")
    return(NULL)
  }
  
  data_df <- as.data.frame(task$data())
  time    <- data_df[[task$target_names[1L]]]
  status  <- data_df[[task$target_names[2L]]]
  
  df <- data.frame(pred = surv_prob, time = time, status = status)
  df <- df[order(df$pred), ]
  
  breaks <- unique(quantile(df$pred, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE))
  if (length(breaks) < 3L) {
    warning("Not enough unique predicted probabilities to form bins. Calibration skipped.")
    return(NULL)
  }
  df$bin <- cut(df$pred, breaks = breaks, include.lowest = TRUE)
  
  obs_surv <- sapply(split(df, df$bin), function(bin_data) {
    if (nrow(bin_data) < 2L) return(NA_real_)
    km <- survival::survfit(survival::Surv(time, status) ~ 1, data = bin_data)
    sp <- summary(km, times = time_point, extend = TRUE)$surv
    if (length(sp) == 0L) NA_real_ else sp[[1L]]
  })
  
  bin_centers <- tapply(df$pred, df$bin, mean, na.rm = TRUE)
  cal_df <- na.omit(data.frame(predicted = as.numeric(bin_centers),
                               observed  = as.numeric(obs_surv)))
  
  if (nrow(cal_df) < 2L) {
    warning("Not enough valid bins for calibration (need at least 2).")
    return(NULL)
  }
  return(cal_df)
}
# ==============================================================================

#' Available Survival Learners
#'
#' A character vector of all survival learner IDs currently available in the
#' \code{mlr3} environment. This is generated at package load time by filtering
#' \code{mlr3::mlr_learners$keys()} for learners with the \code{"surv."} prefix.
#'
#' @format A character vector of learner IDs, e.g., \code{"surv.coxph"},
#'   \code{"surv.ranger"}, etc. Returns an empty vector if \code{mlr3} is not installed.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("mlr3", quietly = TRUE)) {
#'   head(surv_keys, 5)}
#' }
surv_keys <- if (requireNamespace("mlr3", quietly = TRUE)) {
  mlr3::mlr_learners$keys()[grep("^surv\\.", mlr3::mlr_learners$keys())]
} else {
  character(0)
}

#' Flexible Search Space Manager for Survival Analysis
#'
#' Returns a predefined hyperparameter search space for a given survival learner.
#' The search space includes sensible ranges for tuning, with dynamic scaling of
#' parameters like \code{mtry} based on the number of features in the task.
#'
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.ranger"}).
#' @param object Optional. A \code{TaskSurv} or \code{PrognosiX} object used to
#'   determine the number of features for dynamic parameter scaling. If \code{NULL},
#'   defaults to 100 features.
#'
#' @return A \code{\link[paradox]{ParamSet}} object defining the hyperparameter
#'   search space. If no predefined space exists for the learner, attempts to
#'   fetch from \code{mlr3tuningspaces}; if that fails, returns an empty \code{ParamSet}.
#'
#' @details
#' The function maintains predefined search spaces for over 30 survival learners,
#' organized into categories:
#' \itemize{
#'   \item \strong{Random Forests & Ensemble Trees}: \code{surv.ranger}, \code{surv.rfsrc}, \code{surv.aorsf}, etc.
#'   \item \strong{Gradient Boosting}: \code{surv.xgboost.cox}, \code{surv.gbm}, \code{surv.mboost}, etc.
#'   \item \strong{Regularized Regression}: \code{surv.glmnet}, \code{surv.cv_glmnet}, \code{surv.penalized}, etc.
#'   \item \strong{Decision Trees & SVM}: \code{surv.rpart}, \code{surv.ctree}, \code{surv.svm}, etc.
#'   \item \strong{Splines & Flexible Models}: \code{surv.flexreg}, \code{surv.flexspline}, etc.
#'   \item \strong{Neural Networks}: \code{surv.survdnn}
#'   \item \strong{Non-parametric}: \code{surv.kaplan}, \code{surv.nelson} (empty ParamSet)
#' }
#'
#' @seealso \code{\link{surv_get_tuning_config}} for tuning strategy
#' @export
#' @importFrom paradox ps p_int p_dbl p_fct p_lgl
#' @importFrom mlr3tuningspaces lts
#' @examples
#' \dontrun{
#' data(pro_obj_test)
#' library(mlr3proba)
#' 
#' # Get search space for random forest
#' ps <- surv_get_search_space("surv.ranger")
#' print(ps)
#' 
#' # With task for dynamic scaling
#' data("veteran", package = "survival")
#' task <- TaskSurv$new("veteran", backend = veteran, time = "time", event = "status")
#' ps_dynamic <- surv_get_search_space("surv.rfsrc", object = task)
#' }
surv_get_search_space <- function(learner_id, object = NULL) {
  .check_prognosis_packages()
  task <- if (!is.null(object)) surv_extract_task(object) else NULL
  
  # Get the number of features for dynamic scaling of parameters like mtry
  p <- if (!is.null(task)) length(task$feature_names) else 100
  
  search_spaces <- list(
    
    # === 1. Random Forests & Ensemble Trees ===
    "surv.ranger"      = ps(num.trees = p_int(100, 1000), mtry.ratio = p_dbl(0.1, 0.8), min.node.size = p_int(1, 20), splitrule = p_fct(c("logrank", "extratrees"))),
    "surv.rfsrc"       = ps(ntree = p_int(100, 1000), mtry = p_int(1, max(2, floor(sqrt(p)*2))), nodesize = p_int(1, 20), splitrule = p_fct(c("logrank", "random"))),
    "surv.aorsf"       = ps(n_tree = p_int(100, 500), leaf_min_events = p_int(1, 10), split_min_events = p_int(5, 20)),
    "surv.blockforest" = ps(n_trees = p_int(100, 500), block.weights = p_fct(c("proportional", "equal"))),
    "surv.cforest"     = ps(ntree = p_int(100, 500), mtry = p_int(1, max(2, floor(sqrt(p)))), mincriterion = p_dbl(0.5, 0.99)),
    "surv.bart"        = ps(num_trees = p_int(20, 200), k = p_dbl(1, 3), power = p_dbl(1, 3)),
    
    # === 2. Gradient Boosting Machines (GBM) ===
    "surv.xgboost.cox" = ps(nrounds = p_int(50, 500), eta = p_dbl(1e-3, 0.3, logscale = TRUE), max_depth = p_int(2, 8), subsample = p_dbl(0.5, 1)),
    "surv.xgboost.aft" = ps(nrounds = p_int(50, 500), eta = p_dbl(1e-3, 0.3, logscale = TRUE), aft_loss_distribution = p_fct(c("normal", "logistic")), max_depth = p_int(2, 8)),
    "surv.gbm"         = ps(n.trees = p_int(100, 1000), interaction.depth = p_int(1, 5), shrinkage = p_dbl(1e-3, 0.1, logscale = TRUE), n.minobsinnode = p_int(2, 15)),
    "surv.mboost"      = ps(mstop = p_int(50, 500), nu = p_dbl(0.01, 0.2), baselearner = p_fct(c("bbs", "bols", "btree"))),
    "surv.blackboost"  = ps(mstop = p_int(50, 500), maxdepth = p_int(2, 8), nu = p_dbl(0.01, 0.2)),
    "surv.gamboost"    = ps(mstop = p_int(50, 500), nu = p_dbl(0.01, 0.2)),
    "surv.glmboost"    = ps(mstop = p_int(50, 500), nu = p_dbl(0.01, 0.2)),
    "surv.coxboost"    = ps(stepno = p_int(10, 200), penalty = p_dbl(1, 100, logscale = TRUE)),
    "surv.cv_coxboost" = ps(maxstepno = p_int(50, 200), penalty = p_dbl(1, 100, logscale = TRUE)),
    
    # === 3. Regularized & Penalized Regression ===
    "surv.glmnet"      = ps(alpha = p_dbl(0, 1), lambda = p_dbl(1e-4, 1, logscale = TRUE)),
    "surv.cv_glmnet"   = ps(alpha = p_dbl(0, 1)), 
    "surv.penalized"   = ps(lambda1 = p_dbl(0, 20), lambda2 = p_dbl(0, 20)),
    "surv.priority_lasso" = ps(block1.penalization = p_dbl(0, 1), lambda.type = p_fct(c("lambda.min", "lambda.1se"))),
    "surv.cv_ncvsurv"  = ps(penalty = p_fct(c("MCP", "SCAD", "lasso")), alpha = p_dbl(0.1, 1)),
    "surv.coxph"       = ps(ties = p_fct(c("efron", "breslow"))),
    
    # === 4. Decision Trees & Support Vector Machines ===
    "surv.rpart"       = ps(cp = p_dbl(1e-4, 0.1, logscale = TRUE), maxdepth = p_int(1, 30)),
    "surv.ctree"       = ps(mincriterion = p_dbl(0.5, 0.99), minsplit = p_int(2, 30), minbucket = p_int(1, 20)),
    "surv.svm"         = ps(type = p_fct(c("regression", "vanbelle1", "vanbelle2")), kernel = p_fct(c("lin_kernel", "rbf_kernel", "poly_kernel")), mu = p_dbl(0, 1)),
    
    # === 5. Splines & Flexible Parametric Models ===
    "surv.flexreg"     = ps(dist = p_fct(c("weibull", "gengamma", "genf", "gompertz"))),
    "surv.flexspline"  = ps(k = p_int(1, 10), scale = p_fct(c("hazard", "odds", "normal"))),
    "surv.gam.cox"     = ps(select = p_lgl()), 
    
    # === 6. Neural Networks & Baseline Estimators ===
    "surv.survdnn"     = ps(epochs = p_int(10, 100), lr = p_dbl(1e-4, 1e-2, logscale = TRUE), batch_size = p_int(16, 128)),
    "surv.kaplan"      = ps(), # Non-parametric (No Tuning)
    "surv.nelson"      = ps()  # Non-parametric (No Tuning)
  )
  
  # --- Selection Logic ---
  if (learner_id %in% names(search_spaces)) {
    return(search_spaces[[learner_id]])
  } else {
    # Attempt to fetch from mlr3tuningspaces (Expert Default) if not in the list
    t_space <- tryCatch({lts(learner_id)$values}, error = function(e) NULL)
    if (!is.null(t_space)) return(t_space)
    
    message(sprintf("[-] Info: No predefined space for '%s', returning empty ParamSet.", learner_id))
    return(ps())
  }
}

#' Get Recommended Tuning Configuration
#'
#' Provides a recommended tuning strategy (tuner and terminator) based on the
#' complexity of the specified survival learner. Complex models (e.g., random
#' forests, XGBoost) default to random search, while simpler models use grid search.
#'
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.coxph"}).
#' @param tuning_budget An integer specifying the number of evaluations allowed
#'   during tuning. Default is \code{50}.
#' @importFrom mlr3tuning tnr trm
#' @return A list 
#' @details
#' The function distinguishes between complex learners (requiring more
#' sophisticated tuning strategies) and simple learners:
#' \itemize{
#'   \item \strong{Complex learners}: \code{surv.ranger}, \code{surv.xgboost.*},
#'     \code{surv.gbm}, \code{surv.cforest} -> random search with budget.
#'   \item \strong{Simple learners}: For learners with parameters, grid search
#'     with adaptive resolution based on the number of parameters.
#'   \item \strong{No-tuning learners}: \code{surv.kaplan}, \code{surv.nelson}
#'     -> returns \code{NULL} for both components.
#' }
#'
#' @examples
#' \dontrun{
#' # Get tuning config for random forest
#' config <- surv_get_tuning_config("surv.ranger", tuning_budget = 100)
#' print(config$tuner$id)
#' 
#' # For Cox PH (simple learner)
#' config_cox <- surv_get_tuning_config("surv.coxph")
#' print(config_cox$tuner$id)
#' }
#'
#' @seealso \code{\link{surv_get_search_space}} for available search spaces
#' @export
surv_get_tuning_config <- function(learner_id, tuning_budget = 50) {
  .check_prognosis_packages()
  # Select appropriate tuning strategy based on algorithm characteristics
  complex_learners <- c("surv.ranger", "surv.xgboost", "surv.gbm", "surv.cforest")
  
  if (learner_id %in% complex_learners) {
    # Complex models use random search or Bayesian optimization
    tuner <- tnr("random_search")
    terminator <- trm("evals", n_evals = tuning_budget)
  } else {
    # Simple models use grid search
    search_space <- surv_get_search_space(learner_id)
    n_params <- length(search_space$params)
    
    if (n_params == 0) {
      tuner <- NULL
      terminator <- NULL
    } else {
      resolution <- ceiling(tuning_budget^(1/n_params))
      tuner <- tnr("grid_search", resolution = resolution)
      terminator <- trm("evals", n_evals = tuning_budget)
    }
  }
  
  list(tuner = tuner, terminator = terminator)
}

# ==============================================================================
# 2. Data Processing Module
# ==============================================================================

#' Create a Survival Analysis Task
#'
#' Creates a \code{TaskSurv} object from a data frame for use with \code{mlr3proba}
#' survival analysis workflows. The function automatically coerces the data to
#' \code{data.table} for optimized performance.
#'
#' @param data A data frame containing the dataset with survival time and event
#'   status columns.
#' @param time_col A character string specifying the name of the column containing
#'   survival times.
#' @param event_col A character string specifying the name of the column containing
#'   event status (typically 1 for event, 0 for censored).
#' @param id A character string specifying the task identifier. Default is
#'   \code{"survival_task"}.
#'
#' @return A \code{\link[mlr3proba]{TaskSurv}} object ready for use with
#'   \code{mlr3} survival learners.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' data("veteran", package = "survival")
#' 
#' task <- surv_create_surv_task(
#'   data = veteran,
#'   time_col = "time",
#'   event_col = "status",
#'   id = "veteran_task"
#' )
#' print(task)
#' }
#'
#' @seealso \code{\link{surv_extract_task}} for extracting tasks from other objects
#' @export
surv_create_surv_task <- function(data, time_col, event_col, id = "survival_task") {
  .check_prognosis_packages()
  # Coerce to data.table for optimized performance in mlr3
  data <- as.data.table(data)
  task <- TaskSurv$new(
    id = id,
    backend = data,
    time = time_col,
    event = event_col
  )
  # Ensure backend is preserved when the task is saved/loaded
  return(task)
}
# ==============================================================================
# 3. Model Training & Tuning Module
# ==============================================================================

#' Instantiate and Configure a Survival Learner
#'
#' Creates a survival learner instance with appropriate \code{predict_type} and
#' automatically handles categorical features by adding an encoding pipeline
#' if necessary.
#'
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.coxph"}).
#' @param task A \code{TaskSurv} object used to determine feature types and
#'   learner capabilities.
#' @return A configured \code{\link[mlr3]{Learner}} object, possibly wrapped in
#'   a \code{PipeOp} pipeline if encoding is required.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Instantiates the learner using \code{lrn(learner_id)}.
#'   \item Sets \code{predict_type} to \code{"distr"} if available, otherwise
#'     \code{"crank"}.
#'   \item Checks if the task contains factor/character features and if the
#'     learner supports them. If not, adds a \code{po("encode")} pipeline.
#' }
#' @import mlr3pipelines
#' @export
#' @examples
#' veteran <- survival::veteran
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#'
#' lrn_cox <- surv_get_learner("surv.coxph", task)
#' lrn_cox$id
surv_get_learner <- function(learner_id, task) {
  lrn_obj <- lrn(learner_id)
  
  if ("distr" %in% lrn_obj$predict_types) {
    lrn_obj$predict_type <- "distr"
  } else if ("crank" %in% lrn_obj$predict_types) {
    lrn_obj$predict_type <- "crank"
  }
  
  task_ftypes <- task$feature_types$type
  unsupported_factors <- ("factor" %in% task_ftypes || "character" %in% task_ftypes) && 
                         !("factor" %in% lrn_obj$feature_types)
  
  if (unsupported_factors) {
    lrn_obj <- po("encode", method = "treatment") %>>% lrn_obj
    lrn_obj <- mlr3::as_learner(lrn_obj)
    
    if ("distr" %in% lrn_obj$predict_types) {
      lrn_obj$predict_type <- "distr"
    } else if ("crank" %in% lrn_obj$predict_types) {
      lrn_obj$predict_type <- "crank"
    }
  }
  return(lrn_obj)
}

#' Train and Tune a Survival Learner with Hyperparameter Optimization
#'
#' Trains a survival learner on a given task, optionally tuning its
#' hyperparameters using a predefined or custom search space. Handles
#' \code{GraphLearner} objects by correctly prefixing parameter names when an
#' automatic encoding pipeline is added (e.g., for factor features).
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param learner_id Character string identifying the learner (e.g., "surv.coxph").
#' @param search_space A \code{\link[paradox]{ParamSet}} defining the hyperparameter
#'   search space. If \code{NULL}, a default is generated via
#'   \code{surv_get_search_space}.
#' @param resampling A \code{\link[mlr3]{Resampling}} object. Defaults to 5-fold CV.
#' @param measure A \code{\link[mlr3]{Measure}} object. Defaults to \code{surv.cindex}.
#' @param tuning_budget Integer number of evaluations allowed during tuning.
#' @param tuner A \code{\link[mlr3tuning]{Tuner}} object. Defaults to random search.
#' @param seed Integer seed for reproducibility.
#'
#' @return A list with four components:
#' \describe{
#'   \item{learner}{The trained learner with optimal parameters.}
#'   \item{best_params}{List of best hyperparameter values (names are prefixed
#'         if a \code{GraphLearner} was used).}
#'   \item{tuning_result}{Data frame of all evaluated parameter sets and performance.}
#'   \item{cv_performance}{Numeric cross-validated performance score.}
#' }
#'
#' @importFrom mlr3 lrn rsmp msr
#' @importFrom mlr3tuning tnr trm TuningInstanceSingleCrit
#' @importFrom paradox ps
#' @export
#'
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#'
#' tune <- surv_train_and_tune(task, "surv.ranger", tuning_budget = 5)
#' lrn <- tune$learner
#' round(tune$cv_performance, 3)
#'
#' lrn$predict_type <- "distr"
#' perf <- surv_evaluate_model(lrn, task,
#'                             measures = list(mlr3::msr("surv.cindex")))
#' perf$surv.cindex
#' }
surv_train_and_tune <- function(object,
                                learner_id,
                                search_space = NULL,
                                resampling = NULL,
                                measure = NULL,
                                tuning_budget = 50,
                                tuner = NULL,
                                seed = 123) {
  # ---- Package checks ----
  .check_prognosis_packages()
  
  # ---- Extract task ----
  task <- surv_extract_task(object)
  set.seed(seed)
  
  # ---- Instantiate learner (may return GraphLearner) ----
  learner <- tryCatch({
    surv_get_learner(learner_id, task)
  }, error = function(e) {
    stop(sprintf(
      "Failed to instantiate learner '%s': %s\nConsider running mlr3extralearners::install_learner('%s')",
      learner_id, e$message, learner_id
    ))
  })
  
  # ---- Get or validate search space ----
  if (is.null(search_space)) {
    search_space <- surv_get_search_space(learner_id, object = task)
  }
  
  # ---- No tuning needed if search space is empty ----
  if (length(search_space$params) == 0) {
    message(sprintf("[-] Learner '%s' requires no tuning. Training directly...", learner_id))
    learner$train(task)
    return(list(
      learner = learner,
      best_params = list(),
      tuning_result = NULL,
      cv_performance = NA_real_
    ))
  }
  
  # ============================================================
  # CRITICAL FIX: GraphLearner parameter prefixing
  # ============================================================
  # When surv_get_learner() wraps the learner in a GraphLearner
  # (e.g., to encode factors), the parameter names of the inner
  # learner are prefixed with the GraphLearner's set_id.
  # We must rename the search_space parameters accordingly.
  if (inherits(learner, "GraphLearner")) {
    prefix <- learner$param_set$set_id
    if (!is.null(prefix) && nzchar(prefix)) {
      old_names <- names(search_space$params)
      new_names <- paste0(prefix, ".", old_names)
      names(search_space$params) <- new_names
      message("[Tune] GraphLearner detected. Param names prefixed: ",
              paste(new_names, collapse = ", "))
    }
  }
  
  # ---- Default resampling / measure / tuner ----
  if (is.null(resampling)) resampling <- mlr3::rsmp("cv", folds = 5)
  if (is.null(measure)) measure <- mlr3::msr("surv.cindex")
  if (is.null(tuner)) tuner <- mlr3tuning::tnr("random_search")
  terminator <- mlr3tuning::trm("evals", n_evals = tuning_budget)
  
  # ---- Create tuning instance ----
  instance <- mlr3tuning::TuningInstanceSingleCrit$new(
    task = task,
    learner = learner,
    resampling = resampling,
    measure = measure,
    search_space = search_space,
    terminator = terminator
  )
  
  # ---- Run tuning ----
  message(sprintf("[*] Starting tuning for learner '%s' (Budget: %d evals)...",
                  learner_id, tuning_budget))
  tuner$optimize(instance)
  
  # ---- Train final model with best parameters ----
  learner$param_set$values <- instance$result_learner_param_vals
  learner$train(task)
  
  # ---- Return results ----
  list(
    learner = learner,
    best_params = instance$result_learner_param_vals,
    tuning_result = instance$result,
    cv_performance = instance$result_y[[1]]  # single-criterion score
  )
}

# ==============================================================================
# 4. Model Evaluation Module
# ==============================================================================

#' Evaluate Model Performance on a Survival Task
#'
#' Evaluates a trained survival learner on a given task using specified measures.
#' The function automatically selects appropriate measures based on the learner's
#' \code{predict_type}.
#'
#' @param learner A trained \code{\link[mlr3]{Learner}} object.
#' @param object A \code{TaskSurv} or \code{PrognosiX} object. The function
#'   extracts the task using \code{surv_extract_task()}.
#' @param measures A list of \code{\link[mlr3]{Measure}} objects. If \code{NULL},
#'   automatically selects:
#'   \itemize{
#'     \item \code{surv.cindex} (always included)
#'     \item \code{surv.graf} if learner supports \code{"distr"} predictions
#'   }
#'
#' @return A data frame with one row containing the performance metrics for each
#'   requested measure.
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' learner <- lrn("surv.coxph")$train(task)
#' 
#' # Auto-select measures
#' perf <- surv_evaluate_model(learner, task)
#' print(perf)
#' 
#' # Custom measures
#' perf_custom <- surv_evaluate_model(
#'   learner, task,
#'   measures = list(msr("surv.cindex"), msr("surv.logloss"))
#' )
#' }
#'
#' @seealso \code{\link{surv_train_and_tune}}, \code{\link{surv_benchmark_learners}}
#' @export
surv_evaluate_model <- function(learner, object, measures = NULL) {
  task <- surv_extract_task(object)
  
  # Smart measure selection based on learner capabilities
  if (is.null(measures)) {
    if (learner$predict_type == "distr") {
      measures <- list(msr("surv.cindex"), msr("surv.graf")) # C-index & Brier Score
    } else {
      measures <- list(mlr3::msr("surv.cindex")) # Fallback to C-index only
    }
  }
  
  # Generate predictions on the task (Note: usually done on test data)
  predictions <- learner$predict(task)
  
  # Calculate scores safely
  scores <- sapply(measures, function(m) {
    tryCatch({
      predictions$score(m)
    }, error = function(e) {
      NA_real_
    })
  })
  
  # Format results
  names(scores) <- sapply(measures, function(m) m$id)
  as.data.frame(t(scores))
}

# ==============================================================================
# 5. Batch Benchmark Module
# ==============================================================================

#' Batch Train and Benchmark Multiple Learners with CV Performance
#'
#' Trains and evaluates multiple survival learners using cross-validation.
#' Optionally performs hyperparameter tuning for each learner before benchmarking.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param learner_ids A character vector of learner IDs (e.g., \code{c("surv.coxph", "surv.ranger")}).
#' @param tune A logical value. Should hyperparameter tuning be performed?
#'   Default is \code{TRUE}.
#' @param resampling A \code{\link[mlr3]{Resampling}} strategy. If \code{NULL},
#'   defaults to 5-fold cross-validation.
#' @param measures A list of evaluation measures. If \code{NULL}, uses
#'   \code{surv.cindex} for CV and training evaluation.
#' @param tuning_budget An integer specifying the number of tuning evaluations
#'   when \code{tune = TRUE}. Default is \code{50}.
#'
#' @return A list where each element corresponds to a learner and contains:
#'   \describe{
#'     \item{learner}{The trained learner object.}
#'     \item{best_params}{A list of best hyperparameter values (if tuning was performed).}
#'     \item{cv_performance}{A numeric cross-validated C-index.}
#'     \item{performance}{A data frame of training set metrics.}
#'   }
#' @importFrom mlr3 rsmp msr resample 
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' veteran$celltype <- as.factor(veteran$celltype)
#' task <- surv_create_surv_task(veteran, "time", "status", "veteran_task")
#' 
#' # Without tuning (fast)
#' results <- surv_benchmark_learners(
#'   object = task,
#'   learner_ids = c("surv.coxph", "surv.ranger"),
#'   tune = FALSE
#' )
#' 
#' # With tuning (slower)
#' \dontrun{
#' results_tuned <- surv_benchmark_learners(
#'   object = task,
#'   learner_ids = c("surv.coxph", "surv.ranger"),
#'   tune = TRUE,
#'   tuning_budget = 20
#' )
#' }
#' }
#'
#' @seealso \code{\link{surv_summarize_benchmark}} for summarizing results
#' @export
surv_benchmark_learners <- function(object,
                                    learner_ids,
                                    tune = TRUE,
                                    resampling = NULL,
                                    measures = NULL,
                                    tuning_budget = 50) {
  task <- surv_extract_task(object)
  
  # Default resampling: 5-fold CV
  if (is.null(resampling)) {
    resampling <- mlr3::rsmp("cv", folds = 5)
  }
  
  # Default measure for tuning and CV
  if (is.null(measures)) {
    measures <- list(msr("surv.cindex"))
  }
  cv_measure <- measures[[1]]   # Use first measure for CV
  
  results <- list()
  
  for (learner_id in learner_ids) {
    message(sprintf("\n========== Processing Learner: %s ==========", learner_id))
    
    learner_result <- tryCatch({
      if (tune) {
        # Tuning workflow (includes CV performance from tuning)
        tune_res <- surv_train_and_tune(
          object = task,
          learner_id = learner_id,
          resampling = resampling,
          measure = cv_measure,
          tuning_budget = tuning_budget
        )
        learner <- tune_res$learner
        best_params <- tune_res$best_params
        cv_perf <- tune_res$cv_performance
      } else {
        # No tuning: use a FRESH clone for CV (the trained learner must not be
        # reused for CV, as it has already seen all data -- that would invalidate
        # the CV estimate).
        learner_for_cv <- surv_get_learner(learner_id, task)
        rr <- resample(task, learner_for_cv, resampling, store_models = FALSE)
        cv_perf <- rr$aggregate(cv_measure)
        # Now train the final model on the full dataset for downstream use
        learner <- surv_get_learner(learner_id, task)
        learner$train(task)
        best_params <- list()
      }
      
      # Training set (apparent) performance
      train_perf <- surv_evaluate_model(learner, task, measures)
      
      list(
        learner = learner,
        best_params = best_params,
        cv_performance = cv_perf,
        performance = train_perf
      )
    }, error = function(e) {
      message(sprintf("[X] [%s] Failed: %s", learner_id, e$message))
      NULL
    })
    
    if (!is.null(learner_result)) {
      results[[learner_id]] <- learner_result
      message(sprintf("[OK] [%s] Completed successfully (CV C-index = %.4f)", 
                      learner_id, learner_result$cv_performance))
    }
  }
  
  return(results)
}

#' Summarize Benchmark Results into a Leaderboard
#'
#' Converts the output from \code{surv_benchmark_learners} into a sorted
#' leaderboard data frame for easy comparison of model performance.
#'
#' @param benchmark_results The list output from \code{\link{surv_benchmark_learners}}.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{learner}{Learner ID.}
#'     \item{cv_score}{Cross-validated C-index score.}
#'     \item{...}{Additional performance metrics from training set evaluation.}
#'   }
#'   The data frame is sorted by \code{cv_score} in descending order.
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' 
#' # Run benchmark (without tuning for speed)
#' results <- surv_benchmark_learners(
#'   object = task,
#'   learner_ids = c("surv.coxph", "surv.ranger"),
#'   tune = FALSE
#' )
#' 
#' # Summarize
#' summary_df <- surv_summarize_benchmark(results)
#' print(summary_df)
#' }
#'
#' @seealso \code{\link{surv_benchmark_learners}}
#' @export
surv_summarize_benchmark <- function(benchmark_results) {
  perf_list <- lapply(names(benchmark_results), function(learner_id) {
    res <- benchmark_results[[learner_id]]
    if (!is.null(res)) {
      # Extract cv_score (if missing, use training C-index as fallback)
      cv_score <- if (!is.null(res$cv_performance) && !is.na(res$cv_performance)) {
        res$cv_performance
      } else if (!is.null(res$performance$surv.cindex)) {
        warning(sprintf("No CV score for %s, using training C-index as fallback.", learner_id))
        res$performance$surv.cindex
      } else {
        NA_real_
      }
      df <- data.frame(
        learner = learner_id,
        cv_score = cv_score,
        stringsAsFactors = FALSE
      )
      cbind(df, res$performance)
    }
  })
  perf_df <- data.table::rbindlist(perf_list, fill = TRUE)
  if (nrow(perf_df) > 0 && "cv_score" %in% colnames(perf_df)) {
    perf_df <- perf_df[order(-perf_df$cv_score), ]
  }
  return(as.data.frame(perf_df))
}

# ==============================================================================
# 6. Utility Functions Module
# ==============================================================================

#' List Available Survival Learners
#'
#' Returns a character vector of all survival analysis learner IDs currently
#' available in the \code{mlr3} environment.
#'
#' @return A character vector of learner IDs with the \code{"surv."} prefix.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("mlr3", quietly = TRUE)) {
#'   available <- surv_list_available_learners()
#'   print(head(available, 10))
#' surv_list_available_learners()
#' list_surv_feature_methods()
#' print_surv_feature_methods()
#' }
#' }
#'
#' @seealso \code{\link{surv_keys}} for the global list of survival learners
#' @export
surv_list_available_learners <- function() {
  surv_learners <- mlr_learners$keys()[grep("^surv\\.", mlr_learners$keys())]
  return(surv_learners)
}


# ==============================================================================
# 7. Advanced Visualization & Interpretability Module
# ==============================================================================
#' Plot Risk Stratification Kaplan-Meier Curve
#'
#' Generates Kaplan-Meier curves stratified by risk groups based on a model's
#' predicted risk scores or a clinical variable. Supports multiple cutoff
#' determination methods.
#'
#' @param learner A trained mlr3 learner (must support \code{"crank"} predictions
#'   if \code{group_col} is \code{NULL}).
#' @param object A \code{PrognosiX} or \code{TaskSurv} object.
#' @param cutoff_method Method for determining cutoffs:
#'   \itemize{
#'     \item \code{"median"}: split at median (binary).
#'     \item \code{"tertile"}: tertiles (3 groups).
#'     \item \code{"quartile"}: quartiles (4 groups).
#'     \item \code{"p_optimize"}: bootstrap resampling to find cutoff that
#'       minimises log‑rank p‑value (repeatedly, then take median).
#'       \strong{Warning:} This is exploratory; the resulting p‑value is
#'       over‐optimistic and must be validated in an independent set.
#'     \item \code{"maxstat"}: use \code{maxstat::maxstat.test()} to find
#'       the cutoff that maximises the log‑rank statistic (Hothorn & Lausen, 2003).
#'       Supports \code{n_boot} resampling to obtain a stable cutoff.
#'       \strong{Warning:} Same overfitting caveat as \code{p_optimize}.
#'     \item \code{"custom"}: user‐supplied cutoffs.
#'   }
#' @param custom_cutoffs Numeric vector for custom method.
#' @param n_boot Number of bootstrap samples for \code{"p_optimize"} and
#'   \code{"maxstat"}. Default 100.
#' @param fraction Subsample fraction for \code{"p_optimize"} (default 0.7).
#' @param minprop Minimum proportion of observations per group for
#'   \code{"maxstat"} (default 0.3).
#' @param conf_int Logical; draw confidence bands.
#' @param risk_table Logical; show risk table below plot.
#' @param palette_name Wes Anderson palette name.
#' @param show_cutoff Logical; show cutoff values in subtitle.
#' @param title Custom plot title.
#' @param group_col Optional character; name of a clinical variable in
#'   \code{info.data} to use instead of risk scores.
#'
#' @return A \code{ggsurvplot} object (with cutoffs stored as attribute).
#' @export
surv_plot_risk_km <- function(learner, object,
                              cutoff_method = c("median", "tertile", "quartile",
                                                "p_optimize", "maxstat", "custom"),
                              custom_cutoffs = NULL,
                              n_boot = 100,
                              fraction = 0.7,
                              minprop = 0.3,
                              conf_int = FALSE,
                              risk_table = FALSE,
                              palette_name = "AsteroidCity1",
                              show_cutoff = TRUE,
                              title = NULL,
                              group_col = NULL) {
  
  cutoff_method <- match.arg(cutoff_method)
  
  # ---- 1. Extract survival data and optional grouping variable ----
  if (inherits(object, "PrognosiX")) {
    surv_df <- as.data.frame(object@survival.data)
    info_df <- as.data.frame(object@info.data)
    if (!is.null(group_col) && group_col %in% colnames(info_df)) {
      common <- intersect(rownames(surv_df), rownames(info_df))
      if (length(common) == 0) stop("No matching rows between survival.data and info.data.")
      surv_df <- surv_df[common, , drop = FALSE]
      info_df <- info_df[common, , drop = FALSE]
      surv_df[[group_col]] <- info_df[[group_col]]
    }
    task <- surv_extract_task(object)
  } else if (inherits(object, "TaskSurv")) {
    surv_df <- as.data.frame(object$data())
    task <- object
  } else {
    stop("object must be a PrognosiX or TaskSurv.")
  }
  
  # ---- 2. Determine strata ----
  if (!is.null(group_col) && group_col %in% colnames(surv_df)) {
    strata_var <- surv_df[[group_col]]
    lp <- NULL
  } else {
    if (!("crank" %in% learner$predict_types)) learner$predict_type <- "crank"
    lp <- learner$predict(task)$crank
    if (length(unique(lp)) < 2) stop("All risk scores are identical.")
    strata_var <- NULL
  }
  
  # ---- 3. Determine cutoffs ----
  if (!is.null(strata_var)) {
    if (is.factor(strata_var)) {
      cutoffs_used <- levels(strata_var)
    } else {
      cutoffs_used <- sort(unique(strata_var))
    }
    surv_df$risk_group <- factor(strata_var)
  } else {
    # Use risk scores to define groups
    if (cutoff_method == "custom") {
      if (is.null(custom_cutoffs)) stop("custom_cutoffs required.")
      cutoffs_used <- sort(custom_cutoffs)
    } else if (cutoff_method %in% c("median", "tertile", "quartile")) {
      # Simple quantile-based methods
      if (cutoff_method == "median") {
        cutoffs_used <- median(lp, na.rm = TRUE)
      } else if (cutoff_method == "tertile") {
        q <- quantile(lp, probs = c(1/3, 2/3), na.rm = TRUE)
        if (length(unique(q)) < 2) { warning("Tertile cutoffs not unique, using median."); cutoffs_used <- median(lp) } else cutoffs_used <- as.numeric(q)
      } else if (cutoff_method == "quartile") {
        q <- quantile(lp, probs = c(1/4, 2/4, 3/4), na.rm = TRUE)
        if (length(unique(q)) < 3) { warning("Quartile cutoffs not unique, using median."); cutoffs_used <- median(lp) } else cutoffs_used <- as.numeric(q)
      }
    } else if (cutoff_method == "p_optimize") {
      # ---- p_optimize: bootstrap minimising log-rank p ----
      message("p_optimize: Bootstrap resampling to find cutoff that minimises log-rank p-value.")
      message("   This is an exploratory method; the resulting p-value is over-optimistic.")
      message("   Do not report the p-value from this split as confirmatory.")
      if (n_boot < 10) n_boot <- 10
      n_samples <- length(lp)
      sample_size <- max(10, floor(n_samples * fraction))
      boot_cutoffs <- numeric(n_boot)
      for (i in seq_len(n_boot)) {
        idx <- sample(seq_len(n_samples), size = sample_size, replace = TRUE)
        boot_lp <- lp[idx]
        boot_time <- surv_df$time[idx]
        boot_status <- surv_df$status[idx]
        unique_cuts <- sort(unique(boot_lp))
        if (length(unique_cuts) < 3) {
          boot_cutoffs[i] <- median(boot_lp, na.rm = TRUE)
          next
        }
        # Exclude extreme values to ensure both groups have events
        inner_cuts <- unique_cuts[-c(1, length(unique_cuts))]
        if (length(inner_cuts) < 1) {
          boot_cutoffs[i] <- median(boot_lp, na.rm = TRUE)
          next
        }
        p_vals <- sapply(inner_cuts, function(cut) {
          grp <- ifelse(boot_lp > cut, "High", "Low")
          if (length(unique(grp)) < 2) return(NA)
          tab <- table(grp, boot_status)
          if (any(tab < 2)) return(NA)
          fit <- tryCatch(survival::survdiff(survival::Surv(boot_time, boot_status) ~ grp),
                          error = function(e) NULL)
          if (is.null(fit)) return(NA)
          pchisq(fit$chisq, length(fit$n) - 1, lower.tail = FALSE)
        })
        valid <- !is.na(p_vals)
        if (sum(valid) == 0) {
          boot_cutoffs[i] <- median(boot_lp, na.rm = TRUE)
        } else {
          best_idx <- which.min(p_vals[valid])
          boot_cutoffs[i] <- inner_cuts[valid][best_idx]
        }
      }
      cutoffs_used <- median(boot_cutoffs, na.rm = TRUE)
      message(sprintf("   Final cutoff (median of %d bootstraps): %.4f", n_boot, cutoffs_used))
      
    } else if (cutoff_method == "maxstat") {
      # ---- maxstat: use maxstat::maxstat.test ----
      if (!requireNamespace("maxstat", quietly = TRUE)) {
        stop("Package 'maxstat' is required for cutoff_method = 'maxstat'.")
      }
      message("maxstat: Using maximally selected rank statistics (Hothorn & Lausen, 2003).")
      message("   This method finds the cutoff that maximises the log-rank statistic.")
      message("   The resulting cutoff is data-driven and must be validated externally.")
      if (n_boot < 10) n_boot <- 10
      # We'll run maxstat on bootstrap samples and take median cutoff
      n_samples <- length(lp)
      boot_cutoffs <- numeric(n_boot)
      for (i in seq_len(n_boot)) {
        idx <- sample(seq_len(n_samples), size = n_samples, replace = TRUE)
        boot_lp <- lp[idx]
        boot_time <- surv_df$time[idx]
        boot_status <- surv_df$status[idx]
        # Fit maxstat
        ms <- tryCatch(
          maxstat::maxstat.test(
            survival::Surv(boot_time, boot_status) ~ boot_lp,
            data = data.frame(boot_lp, boot_time, boot_status),
            smethod = "LogRank",
            minprop = minprop
          ),
          error = function(e) NULL
        )
        if (is.null(ms) || is.na(ms$estimate)) {
          boot_cutoffs[i] <- median(boot_lp, na.rm = TRUE)
        } else {
          boot_cutoffs[i] <- ms$estimate
        }
      }
      cutoffs_used <- median(boot_cutoffs, na.rm = TRUE)
      message(sprintf("   Final cutoff (median of %d bootstraps): %.4f", n_boot, cutoffs_used))
    }
    
    # Assign risk groups
    surv_df$risk_group <- .assign_risk_groups(lp, cutoffs_used)
  }
  
  # ---- 4. Fit KM curves ----
  fit <- survival::survfit(survival::Surv(time, status) ~ risk_group, data = surv_df)
  
  # ---- 5. Plot ----
  n_groups <- length(unique(surv_df$risk_group))
  cols <- .get_palette(palette_name, n_groups)
  
  if (is.null(title)) title <- if (is.null(group_col)) "Risk Stratification" else paste("KM by", group_col)
  
  p <- survminer::ggsurvplot(
    fit, data = surv_df,
    pval = (n_groups > 1), pval.method = (n_groups > 1),
    conf.int = conf_int, risk.table = risk_table,
    title = title,
    subtitle = if (show_cutoff && !is.null(cutoffs_used)) paste("Cutoffs:", paste(round(cutoffs_used, 3), collapse = ", ")) else "",
    palette = cols,
    ggtheme = ggprism::theme_prism() + ggplot2::theme(legend.title = ggplot2::element_blank()),
    legend = "right"
  )
  
  attr(p, "cutoffs_used") <- cutoffs_used
  return(p)
}

# Internal helper for risk group assignment (used by above)
.assign_risk_groups <- function(lp, cutoffs) {
  cutoffs <- sort(cutoffs)
  n_cuts <- length(cutoffs)
  if (n_cuts == 1) {
    factor(ifelse(lp > cutoffs, "High Risk", "Low Risk"),
           levels = c("Low Risk", "High Risk"))
  } else if (n_cuts == 2) {
    cut(lp, breaks = c(-Inf, cutoffs[1], cutoffs[2], Inf),
        labels = c("Low Risk", "Medium Risk", "High Risk"))
  } else {
    labels <- c("Low Risk", paste("Risk Group", 2:(n_cuts)), "High Risk")
    if (length(labels) != n_cuts + 1) 
      labels <- paste("Group", 1:(n_cuts+1))
    cut(lp, breaks = c(-Inf, cutoffs, Inf), labels = labels)
  }
}


#' Extract Cutoffs from a Risk Stratification Plot
#'
#' Retrieves the cutoffs used to create risk groups in a Kaplan-Meier plot
#' generated by \code{surv_plot_risk_km}.
#'
#' @param km_plot A \code{ggsurvplot} object returned by \code{\link{surv_plot_risk_km}}.
#'
#' @return A numeric vector of cutoffs used in the plot, or \code{NULL} if not available.
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' learner <- lrn("surv.coxph")$train(task)
#' 
#' p <- surv_plot_risk_km(learner, task, cutoff_method = "median")
#' cutoffs <- get_cf(p)
#' print(cutoffs)
#' }
#'
#' @seealso \code{\link{surv_plot_risk_km}}
#' @export
get_cf <- function(km_plot) {
  attr(km_plot, "cutoffs_used")
}


#' Generate Clinical Nomogram for Survival Model (with PH test)
#'
#' Creates a nomogram for predicting survival probabilities at specified time
#' points. Automatically selects time points from observed event quantiles if
#' the user-provided ones fail. Also performs a proportional hazards assumption
#' test using \code{cox.zph} and prints the result.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param selected_features Character vector of feature names to include.
#'   If \code{NULL}, uses the top 5 features from the task.
#' @param time_points Numeric vector of time points; if \code{NULL} or if
#'   the model fails at these points, falls back to quantiles of event times.
#' @param time_unit Optional label for the x‑axis (e.g., "days", "months").
#'
#' @return A \code{nomogram} object (invisibly) and prints the nomogram plot,
#'   or \code{NULL} if the model cannot be fitted.
#' @export
surv_generate_nomogram <- function(object,
                                   selected_features = NULL,
                                   time_points = NULL,
                                   time_unit = NULL) {
  
  if (!requireNamespace("rms", quietly = TRUE)) {
    stop("Package 'rms' is required for nomogram generation.")
  }
  
  task <- surv_extract_task(object)
  data <- as.data.frame(task$data())
  
  # ---- 1. Select features ----
  if (is.null(selected_features)) {
    selected_features <- head(task$feature_names, 5)
  }
  selected_features <- intersect(selected_features, colnames(data))
  if (length(selected_features) == 0) {
    stop("No valid features selected.")
  }
  
  # ---- 2. Prepare time points ----
  event_times <- data[[task$target_names[1]]][data[[task$target_names[2]]] == 1]
  if (length(event_times) < 2) {
    stop("Fewer than 2 event times; cannot estimate survival curve.")
  }
  
  if (!is.null(time_points)) {
    time_points <- time_points[time_points >= min(event_times) & time_points <= max(event_times)]
    if (length(time_points) < 2) {
      warning("User-provided time points outside event range or insufficient. Falling back to quantiles.")
      time_points <- NULL
    }
  }
  
  if (is.null(time_points)) {
    probs <- c(0.25, 0.5, 0.75)
    time_points <- round(quantile(event_times, probs = probs, na.rm = TRUE), 0)
    time_points <- unique(time_points)
    if (length(time_points) < 2) {
      time_points <- unique(c(min(event_times), median(event_times), max(event_times)))
    }
    message(sprintf("Using automatically selected time points: %s",
                    paste(time_points, collapse = ", ")))
  }
  
  # ---- 3. Fit Cox model with rms ----
  dd <- rms::datadist(data)
  old_dd <- options()$datadist
  on.exit(options(datadist = old_dd), add = TRUE)
  options(datadist = dd)
  
  formula_str <- paste("survival::Surv(", task$target_names[1], ", ",
                       task$target_names[2], ") ~ ",
                       paste(selected_features, collapse = " + "))
  
  fit <- tryCatch({
    rms::cph(as.formula(formula_str),
             data = data, x = TRUE, y = TRUE, surv = TRUE)
  }, error = function(e) {
    warning("cph model fitting failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(fit)) {
    message("Nomogram generation aborted due to model fitting error.")
    return(NULL)
  }
  
  # ---- 4. Proportional Hazards Assumption Test ----
  if (requireNamespace("survival", quietly = TRUE)) {
    ph_test <- tryCatch(
      survival::cox.zph(fit),
      error = function(e) {
        warning("cox.zph failed: ", e$message)
        return(NULL)
      }
    )
    if (!is.null(ph_test)) {
      cat("\n========== Proportional Hazards Assumption Test ==========\n")
      print(ph_test)
      cat("  If any p-value < 0.05, the PH assumption may be violated.\n")
      cat("  Consider including time-dependent covariates or stratification.\n")
      cat("===========================================================\n\n")
    }
  }
  
  # ---- 5. Create survival function and nomogram ----
  surv_obj <- rms::Survival(fit)
  
  surv_funcs <- lapply(time_points, function(t) {
    function(x) surv_obj(t, x)
  })
  
  if (is.null(time_unit)) {
    funlabel <- paste0(time_points, " Survival")
  } else {
    unit_label <- paste0(toupper(substring(time_unit, 1, 1)),
                         substring(time_unit, 2))
    funlabel <- paste0(time_points, "-", unit_label, " Survival")
  }
  
  nom <- tryCatch({
    rms::nomogram(fit,
                  fun = surv_funcs,
                  funlabel = funlabel,
                  lp = FALSE)
  }, error = function(e) {
    if (grepl("approx", e$message, fixed = TRUE)) {
      warning("Interpolation error in nomogram. This often means the survival ",
              "curve has only one distinct value at the chosen time points. ",
              "Try different time points or reduce the number of features.")
    } else {
      warning("Nomogram creation failed: ", e$message)
    }
    return(NULL)
  })
  
  if (is.null(nom)) {
    return(NULL)
  }
  
  plot(nom)
  return(invisible(nom))
}

#' SurvSHAP(t) Explanations for Survival Models -- Production Version
#'
#' Computes time-dependent SurvSHAP(t) explanations for survival predictions.
#' Supports both local (individual patient) and global (population-level)
#' interpretations using Kernel SHAP.
#'
#' @param learner A trained \code{mlr3} \code{LearnerSurv} object.
#' @param task A \code{TaskSurv} object.
#' @param type A character string specifying the explanation type. Must be one of
#'   \code{"local"} (explain specific patients) or \code{"global"} (population-level
#'   importance averaged over many patients).
#' @param n_explain An integer specifying the number of observations to explain.
#'   For \code{type = "local"}: number of patients to explain. For
#'   \code{type = "global"}: number of observations to aggregate over.
#'   Default is \code{NULL} (auto-set to 20 for local, 50 for global).
#' @param n_background An integer specifying the background data size for Kernel SHAP.
#'   Default is \code{50L}.
#' @param n_timepoints An integer specifying how many evaluation time points to use.
#'   If \code{NULL}, uses all available time points. Default is \code{NULL}.
#' @param n_top_features An integer specifying the top features to show in plots.
#'   Default is \code{6L}.
#' @param bar_color A character string specifying the hex color for bar plots.
#'   Default is \code{"#2980b9"}.
#' @param seed An integer seed for reproducibility. Default is \code{123L}.
#' @param verbose A logical value. Print progress messages. Default is \code{TRUE}.
#'
#' @return A list with six components:
#'   \describe{
#'     \item{shap_long}{A tidy data frame with columns \code{feature}, \code{time},
#'       \code{shap_value}, and \code{observation}.}
#'     \item{explainer}{A \code{survex} explainer object.}
#'     \item{eval_times}{The evaluation time points actually used.}
#'     \item{plots}{A list with \code{bar_plot} (feature importance) and
#'       \code{line_plot} (SHAP dynamics over time).}
#'     \item{original_features}{A data frame with the original feature values.}
#'   }
#'
#' @details
#' The function uses the \code{survex} package for efficient SHAP computation.
#' Key features:
#' \itemize{
#'   \item \strong{Time-dependent}: SHAP values are computed at multiple time points.
#'   \item \strong{Sampling acceleration}: Sub-samples time points and observations for speed.
#'   \item \strong{Native mlr3 support}: Handles factor/character columns without errors.
#' }
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' learner <- lrn("surv.coxph")$train(task)
#' 
#' # Global explanation (population-level)
#' shap_result <- surv_explain_shap(
#'   learner = learner,
#'   task = task,
#'   type = "global",
#'   n_explain = 20,
#'   n_background = 20,
#'   n_timepoints = 5
#' )
#' 
#' # View bar plot
#' print(shap_result$plots$bar_plot)
#' 
#' # Local explanation (single patient)
#' shap_local <- surv_explain_shap(
#'   learner = learner,
#'   task = task,
#'   type = "local",
#'   n_explain = 1,
#'   n_background = 20
#' )
#' }
#'
#' @seealso \code{\link{surv_plot_shap_beeswarm}} for visualizing SHAP values
#' @export
surv_explain_shap <- function(
    learner,
    task,
    type               = c("local", "global"),
    n_explain          = NULL,
    n_background       = 50L,
    n_timepoints       = NULL,
    n_top_features     = 6L,
    bar_color          = "#2980b9",
    seed               = 123L,
    verbose            = TRUE) {
  
  type <- match.arg(type)
  set.seed(seed)
  
  # ---- Input validation ----
  if (!inherits(learner, "LearnerSurv"))
    stop("learner must be mlr3::LearnerSurv. Got: ", class(learner)[1])
  if (!inherits(task, "TaskSurv"))
    stop("task must be mlr3::TaskSurv. Got: ", class(task)[1])
  
  if (n_background < 10L)
    warning("n_background < 10: estimates may be very noisy. Recommended >= 50.")
  
  # ---- Set defaults ----
  if (is.null(n_explain)) {
    n_explain <- if (type == "local") 20L else 50L
  }
  n_explain    <- min(as.integer(n_explain), nrow(task$data()))
  n_background <- min(as.integer(n_background), nrow(task$data()))
  
  if (verbose) {
    cat(sprintf("[SurvSHAP] mode=%s, n_explain=%d, n_background=%d\n",
                type, n_explain, n_background))
  }
  
  # ---- Build explainer ----
  data     <- as.data.frame(task$data())
  features <- data[, task$feature_names, drop = FALSE]
  
  target   <- survival::Surv(data[[task$target_names[1L]]],
                             data[[task$target_names[2L]]])
  
  predict_surv_fn <- function(model, newdata, times) {
    old_type <- model$predict_type
    model$predict_type <- "distr"
    on.exit(model$predict_type <- old_type, add = TRUE)
    
    time_col <- task$target_names[1L]
    status_col <- task$target_names[2L]
    
    tmp_df <- newdata
    tmp_df[[time_col]] <- 1
    tmp_df[[status_col]] <- 0
    
    tmp_task <- mlr3proba::TaskSurv$new("tmp", as.data.frame(tmp_df),
                                        time = time_col, event = status_col)
    pred <- model$predict(tmp_task)
    t(as.matrix(pred$distr$survival(times)))
  }
  
  explainer <- survex::explain_survival(
    learner, data = features, y = target,
    predict_survival_function = predict_surv_fn,
    label = learner$id, verbose = FALSE
  )
  
  # ---- Time-point subsampling for speed ----
  eval_times_full <- explainer$times
  if (!is.null(n_timepoints) && length(eval_times_full) > n_timepoints) {
    idx_subsample <- seq(1L, length(eval_times_full),
                         length.out = as.integer(n_timepoints))
    eval_times_use <- eval_times_full[idx_subsample]
    if (verbose) {
      cat(sprintf("[SurvSHAP] Subsampling times: %d -> %d points (%.1fx speedup)\n",
                  length(eval_times_full), length(eval_times_use),
                  length(eval_times_full) / length(eval_times_use)))
    }
    explainer$times <- eval_times_use
  }
  
  # ---- Sample observations to explain ----
  obs_idx <- sample(nrow(features), size = n_explain, replace = FALSE)
  new_obs <- features[obs_idx, , drop = FALSE]
  
  # ---- Compute SHAP (either local or global mode) ----
  if (type == "global") {
    if (verbose) cat("[SurvSHAP] Computing GLOBAL SurvSHAP(t)...\n")
    shap_obj <- tryCatch(
      survex::model_survshap(
        explainer,
        new_observation    = new_obs,
        N                  = n_background,
        calculation_method = "kernelshap",
        aggregation_method = "mean_absolute"
      ),
      error = function(e) {
        stop("model_survshap failed. Diagnostic:\n  ", conditionMessage(e))
      }
    )
  } else {
    if (verbose) cat("[SurvSHAP] Computing LOCAL SurvSHAP(t)...\n")
    shap_obj <- tryCatch(
      survex::predict_parts(
        explainer,
        new_observation = new_obs,
        type            = "survshap",
        output_type     = "survival",
        N               = n_background
      ),
      error = function(e) {
        stop("predict_parts failed. Diagnostic:\n  ", conditionMessage(e))
      }
    )
  }
  
  # ---- Convert to tidy long format ----
  shap_long <- .survshap_to_long(shap_obj, obs_idx, features)
  
  if (is.null(shap_long) || nrow(shap_long) == 0L)
    stop("No SHAP rows returned. Check task/learner compatibility.")
  
  # ---- Generate plots ----
  plots <- .generate_shap_plots(shap_long, n_top_features, bar_color,
                                learner$id, type)
  
  if (verbose) {
    cat(sprintf("[SurvSHAP] Complete. %d observations x %d features x %d times.\n",
                length(unique(shap_long$observation)),
                length(unique(shap_long$feature)),
                length(unique(shap_long$time[!is.na(shap_long$time)]))))
  }
  
  list(
    shap_long         = shap_long,
    explainer         = explainer,
    eval_times        = if (!is.null(n_timepoints)) eval_times_use else eval_times_full,
    plots             = plots,
    original_features = features
  )
}


# =========================================================================
# Beeswarm / Violin Plot 
# =========================================================================

#' Survival SHAP Beeswarm/Violin Summary Plot
#'
#' This function creates a beeswarm or violin summary plot for SHAP (SHapley 
#' Additive exPlanations) values from survival models. It aggregates SHAP values 
#' across observations and optionally at a specific time point, displaying the 
#' most important features based on mean absolute SHAP values.
#'
#' @param shap_result A list object returned by a survival SHAP calculation function,
#'   typically containing `shap_long` (long-format SHAP values) and optionally
#'   `original_features` (original feature values for coloring).
#' @param time_point A numeric value specifying the time point at which to aggregate
#'   SHAP values. If `NULL` or all SHAP values have missing time, SHAP values are
#'   averaged across all time points. Default is `NULL`.
#' @param top_n An integer specifying the number of top features to display based
#'   on mean absolute SHAP values. Default is `8L`.
#' @param method A character string specifying the plot type. Must be one of
#'   `"beeswarm"` (default) or `"violin"`.
#' @param color_low A character string specifying the color for low feature values
#'   in the color gradient. Default is `"#2c7bb6"` (blue).
#' @param color_high A character string specifying the color for high feature values
#'   in the color gradient. Default is `"#d7191c"` (red).
#' @param title A character string for the plot title. If `NULL`, a default title
#'   is generated with time point and sample size information. Default is `NULL`.
#'
#' @return A \code{ggplot} object representing the SHAP summary plot. Points are
#'   colored by feature values if `original_features` is available in `shap_result`.
#'
#' @details
#' The function first aggregates SHAP values from the long-format input. If a
#' `time_point` is specified, SHAP values are filtered to the closest available
#' evaluation time before aggregation; otherwise, values are averaged across all
#' time points. Features are ranked by mean absolute SHAP value, and the top
#' \code{top_n} features are displayed.
#'
#' Two visualization methods are supported:
#' \itemize{
#'   \item \code{"beeswarm"}: Points are distributed along the y-axis to avoid
#'     overlap, with optional coloring by feature values.
#'   \item \code{"violin"}: Violin plots show the distribution of SHAP values,
#'     with quasi-random points overlaid for individual observations.
#' }
#'
#' When `original_features` is provided in `shap_result`, points are colored by
#' the corresponding feature values. Numeric features use a continuous color
#' gradient, while categorical features use discrete colors. The function uses
#' the \pkg{ggbeeswarm} package for beeswarm and quasirandom geometries.
#'
#' @note
#' The `shap_result` object must contain a `shap_long` data frame with columns
#' `observation`, `feature`, `shap_value`, and optionally `time`. If
#' `original_features` is provided, it should be a data frame with observations
#' as row names and features as columns.
#'
#' This function requires the following packages: \pkg{ggplot2}, \pkg{ggbeeswarm},
#' \pkg{ggprism}, \pkg{dplyr}, and \pkg{tidyr}.
#'
#' @importFrom dplyr group_by summarise arrange desc slice rename left_join
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_vline theme element_text
#'   labs scale_color_gradient geom_violin
#' @importFrom ggbeeswarm geom_beeswarm geom_quasirandom
#' @importFrom ggprism theme_prism
#' @export
#'
#' @examples
#' \dontrun{
#' # Load required libraries
#' library(dplyr)
#' library(ggplot2)
#'
#' # Example SHAP result structure (simulated)
#' set.seed(123)
#' shap_long <- data.frame(
#'   observation = rep(1:50, each = 5),
#'   feature = rep(c("age", "sex", "bmi", "stage", "treatment"), 50),
#'   shap_value = rnorm(250, 0, 1),
#'   time = rep(c(12, 24, 36), length.out = 250)
#' )
#'
#' original_features <- data.frame(
#'   age = rnorm(50, 65, 10),
#'   sex = factor(sample(c("M", "F"), 50, replace = TRUE)),
#'   bmi = rnorm(50, 28, 5),
#'   stage = factor(sample(1:4, 50, replace = TRUE)),
#'   treatment = factor(sample(c("A", "B"), 50, replace = TRUE))
#' )
#' rownames(original_features) <- 1:50
#'
#' shap_result <- list(
#'   shap_long = shap_long,
#'   original_features = original_features
#' )
#'
#' # Basic beeswarm plot (averaged across time)
#' surv_plot_shap_beeswarm(shap_result, top_n = 5)
#'
#' # Plot at specific time point with custom colors
#' surv_plot_shap_beeswarm(
#'   shap_result,
#'   time_point = 24,
#'   top_n = 6,
#'   color_low = "darkblue",
#'   color_high = "darkred"
#' )
#'
#' # Violin plot method
#' surv_plot_shap_beeswarm(
#'   shap_result,
#'   method = "violin",
#'   top_n = 4,
#'   title = "SHAP Summary - Violin Plot"
#' )
#'
#' # Without feature coloring
#' shap_result_no_color <- list(shap_long = shap_long)
#' surv_plot_shap_beeswarm(shap_result_no_color, top_n = 5)
#' }
surv_plot_shap_beeswarm <- function(shap_result,
                                    time_point = NULL,
                                    top_n      = 8L,
                                    method     = c("beeswarm", "violin"),
                                    color_low  = "#2c7bb6",
                                    color_high = "#d7191c",
                                    title      = NULL) {
  
  method    <- match.arg(method)
  shap_long <- shap_result$shap_long
  
  if (is.null(shap_long) || nrow(shap_long) == 0L)
    stop("shap_long is empty.")
  
  n_obs <- length(unique(shap_long$observation))
  
  if (is.null(time_point) || all(is.na(shap_long$time))) {
    shap_agg   <- shap_long %>%
      dplyr::group_by(observation, feature) %>%
      dplyr::summarise(shap = mean(shap_value, na.rm = TRUE), .groups = "drop")
    time_label <- ""
  } else {
    eval_times <- sort(unique(shap_long$time[!is.na(shap_long$time)]))
    closest_t  <- eval_times[which.min(abs(eval_times - time_point))]
    shap_agg   <- shap_long %>%
      dplyr::filter(abs(time - closest_t) < 1e-9) %>%
      dplyr::rename(shap = shap_value)
    time_label <- sprintf(" (t = %.1f)", closest_t)
  }
  
  feat_imp <- shap_agg %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(mean_abs = mean(abs(shap), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice(seq_len(min(as.integer(top_n), dplyr::n())))
  
  top_features <- feat_imp$feature
  plot_df      <- shap_agg %>% dplyr::filter(feature %in% top_features)
  plot_df$feature <- factor(plot_df$feature, levels = rev(top_features))
  
  color_by_value <- FALSE
  is_numeric_scale <- FALSE
  
  if (!is.null(shap_result$original_features)) {
    var_vals <- shap_result$original_features
    
    var_vals[] <- lapply(var_vals, as.character)
    
    var_long <- tidyr::pivot_longer(
      cbind(observation = rownames(var_vals), var_vals),
      cols = -observation, names_to = "feature", values_to = "feature_value"
    )
    plot_df <- dplyr::left_join(plot_df, var_long, by = c("observation", "feature"))
    color_by_value <- !all(is.na(plot_df$feature_value))

    plot_df$feature_value_num <- suppressWarnings(as.numeric(as.character(plot_df$feature_value)))
    is_numeric_scale <- !all(is.na(plot_df$feature_value_num))
  }
  
  if (!requireNamespace("ggbeeswarm", quietly = TRUE))
    stop("Package 'ggbeeswarm' required.")
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = shap, y = feature))
  
  if (method == "beeswarm") {
    if (color_by_value) {
      if (is_numeric_scale) {
        p <- p + ggbeeswarm::geom_beeswarm(ggplot2::aes(color = feature_value_num), size = 2, cex = 2)
      } else {
        p <- p + ggbeeswarm::geom_beeswarm(ggplot2::aes(color = factor(feature_value)), size = 2, cex = 2)
      }
    } else {
      p <- p + ggbeeswarm::geom_beeswarm(ggplot2::aes(color = feature), size = 2, cex = 2)
    }
  } else {
    fill_aes <- if (color_by_value) {
      if (is_numeric_scale) ggplot2::aes(color = feature_value_num) else ggplot2::aes(color = factor(feature_value))
    } else {
      ggplot2::aes(color = feature)
    }
    p <- p + ggplot2::geom_violin(fill = "gray90", alpha = 0.5) + ggbeeswarm::geom_quasirandom(fill_aes, size = 2, width = 0.3)
  }
  
  p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggprism::theme_prism(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  
  if (color_by_value) {
    if (is_numeric_scale) {
      p <- p + ggplot2::scale_color_gradient(low = color_low, high = color_high, name = "Feature Value")
    } else {
      p <- p + ggplot2::labs(color = "Category")
    }
  }
  
  if (is.null(title))
    title <- sprintf("SurvSHAP Summary%s (%d patients)", time_label, n_obs)
  p + ggplot2::labs(title = title, x = "SHAP value", y = NULL)
}
# ==============================================================================
# 9. Clinical Reporting & Subgroup Analysis
# ==============================================================================

#' Subgroup Forest Plot for Survival Models (with info.data Support)
#'
#' Creates a publication‑quality forest plot showing hazard ratios for a risk
#' score within subgroups defined by categorical variables. The function now
#' retrieves subgroup variables from both \code{survival.data} and \code{info.data}
#' of the \code{PrognosiX} object, so clinical metadata are fully accessible.
#'
#' @param learner A trained \code{mlr3} \code{LearnerSurv} object.
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param subgroup_vars Character vector of categorical column names.
#' @param prog Optional \code{PrognosiX} object; if provided, used for data extraction.
#' @param var_labels Optional named vector to relabel section headers.
#' @param level_order Optional named list of level order per variable.
#' @param ref_line Numeric. Reference vertical line (default 1).
#' @param xlim Numeric length‑2. If \code{NULL}, automatically computed.
#' @param ticks_at Numeric vector. If \code{NULL}, automatically computed.
#' @param ticks_digits Integer. Decimal digits for x-axis tick labels.
#' @param x_trans Passed to \code{forestploter::forest()}; \code{"log"} (default).
#' @param digits Integer for HR/CI decimals.
#' @param p_digits Integer for p‑value decimals.
#' @param box_col Color for the point‑estimate box.
#' @param ci_col Color for CI whiskers.
#' @param theme_args Named list passed to \code{forestploter::forest_theme()}.
#' @param arrow_labels Character length‑2 or NULL.
#' @param footnote Optional footnote string.
#' @param save_plot Logical; save plot to PDF.
#' @param save_dir Character; output directory.
#' @param file_name Character; output filename.
#' @param width Numeric; plot width (inches).
#' @param height Numeric; plot height (inches).
#'
#' @return Invisibly, a list with components \code{plot}, \code{data}, and \code{fits}.
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#' lrn <- surv_get_learner("surv.coxph", task)
#' lrn$train(task)
#'
#' forest <- surv_plot_subgroup_forest(lrn, task, subgroup_vars = c("celltype", "trt"))
#' forest
#' }
surv_plot_subgroup_forest <- function(learner, object, subgroup_vars, prog = NULL,
                                      var_labels = NULL, level_order = NULL,
                                      ref_line = 1,
                                      xlim = NULL, ticks_at = NULL,
                                      ticks_digits = 2,
                                      x_trans = "log", digits = 2, p_digits = 3,
                                      box_col = "#377eb8", ci_col = "black",
                                      theme_args = list(), arrow_labels = NULL,
                                      footnote = NULL, save_plot = FALSE,
                                      save_dir = NULL, file_name = NULL,
                                      width = 8, height = 6) {
  
  # ---- package checks ----
  if (!requireNamespace("forestploter", quietly = TRUE)) {
    stop("Package 'forestploter' is required for this plot.\n",
         "Install it with: install.packages('forestploter')")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for Cox models.")
  }
  
  # ---- Phase 1: Resolve data sources ----
  if (inherits(object, "TaskSurv")) {
    task <- object
    row_filter <- task$row_ids
    cohort_name <- task$id
    if (is.null(prog)) {
      prog <- get0("prog", envir = parent.frame(), inherits = TRUE)
      if (is.null(prog)) prog <- get0("prog", envir = .GlobalEnv)
    }
    missing_in_task <- setdiff(subgroup_vars, task$backend$cols)
    if (length(missing_in_task) > 0L && is.null(prog)) {
      stop(sprintf(
        "Subgroup variables [%s] not found in Task backend. Provide 'prog'.",
        paste(missing_in_task, collapse = ", ")
      ))
    }
  } else if (inherits(object, "PrognosiX")) {
    prog <- object
    cohort_name <- "Full Cohort"
    task <- tryCatch(surv_extract_task(prog), error = function(e) NULL)
    if (!is.null(task)) {
      row_filter <- task$row_ids
    } else {
      row_filter <- 1:nrow(prog@survival.data)
    }
  } else {
    stop("'object' must be a TaskSurv or PrognosiX object.")
  }
  
  if (!is.null(task)) {
    time_col   <- task$target_names[1L]
    status_col <- task$target_names[2L]
  } else {
    time_col   <- prog@time_col %||% "time"
    status_col <- prog@status_col %||% "status"
  }
  
  # ---- Phase 2: Build a unified data frame (features + clinical vars) ----
  if (!is.null(prog) && inherits(prog, "PrognosiX")) {
    # Combine survival.data (features + time/status) with info.data (clinical)
    surv_df <- as.data.frame(prog@survival.data)
    info_df <- as.data.frame(prog@info.data)
    
    # Ensure row names align
    common_ids <- intersect(rownames(surv_df), rownames(info_df))
    if (length(common_ids) == 0) {
      stop("No common sample IDs between survival.data and info.data.")
    }
    surv_df <- surv_df[common_ids, , drop = FALSE]
    info_df <- info_df[common_ids, , drop = FALSE]
    
    # Start with surv_df, then add missing subgroup variables from info_df
    data <- surv_df
    for (v in subgroup_vars) {
      if (v %in% colnames(info_df) && !(v %in% colnames(data))) {
        data[[v]] <- info_df[[v]]
      }
    }
    # Check for still-missing variables
    missing_vars <- setdiff(subgroup_vars, colnames(data))
    if (length(missing_vars) > 0L) {
      warning("Variables not found in prog (neither survival.data nor info.data): ",
              paste(missing_vars, collapse = ", "))
    }
    
    # Keep only the variables we need (time, status, subgroup vars)
    needed_cols <- unique(c(time_col, status_col, subgroup_vars))
    data <- data[, intersect(needed_cols, colnames(data)), drop = FALSE]
    
  } else {
    # Fallback: use task data only (if prog not available)
    cols_to_fetch <- unique(c(time_col, status_col, subgroup_vars))
    valid_cols <- intersect(cols_to_fetch, task$backend$cols)
    data <- as.data.frame(task$backend$data(rows = row_filter, cols = valid_cols))
  }
  
  # Compute risk scores from the learner
  if (!is.null(task)) {
    predictions <- learner$predict(task)
    data$risk_score <- predictions$crank
  } else {
    stop("Unable to compute risk scores.")
  }
  
  # ---- Phase 3: Run Cox regressions per subgroup ----
  valid_vars <- intersect(subgroup_vars, colnames(data))
  if (length(valid_vars) == 0L) {
    stop("No valid subgroup variables found.")
  }
  
  results <- list()
  fit_list <- list()
  
  for (var in valid_vars) {
    if (all(is.na(data[[var]]))) next
    
    levels_var <- unique(na.omit(data[[var]]))
    if (!is.null(level_order) && !is.null(level_order[[var]])) {
      levels_var <- level_order[[var]][level_order[[var]] %in% levels_var]
    }
    
    for (lev in levels_var) {
      subset_data <- data[data[[var]] == lev & !is.na(data[[var]]), ]
      n_total  <- nrow(subset_data)
      n_events <- sum(subset_data[[status_col]], na.rm = TRUE)
      
      if (n_events < 5L) {
        message(sprintf("Skipping %s - %s (events=%d < 5).", var, lev, n_events))
        next
      }
      
      fit <- tryCatch({
        survival::coxph(
          survival::Surv(subset_data[[time_col]], subset_data[[status_col]]) ~ risk_score,
          data = subset_data
        )
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        hr <- exp(stats::coef(fit))
        ci <- exp(confint(fit))
        if (is.na(hr) || is.infinite(hr) || any(is.na(ci))) next
        p_val <- summary(fit)$coefficients[1, "Pr(>|z|)"]
        
        results[[paste0(var, "_", lev)]] <- data.frame(
          Variable = var,
          Subgroup = as.character(lev),
          N        = n_total,
          HR       = as.numeric(hr),
          Lower    = as.numeric(ci[1L]),
          Upper    = as.numeric(ci[2L]),
          P        = p_val,
          stringsAsFactors = FALSE
        )
        fit_list[[paste0(var, "_", lev)]] <- fit
      }
    }
  }
  
  if (length(results) == 0L) {
    warning("No subgroups met the criteria. Returning NULL.")
    return(NULL)
  }
  
  res_df <- do.call(rbind, results)
  rownames(res_df) <- NULL
  
  # ---- Phase 4: Build forestploter display table ----
  rows <- list()
  
  for (var in unique(res_df$Variable)) {
    sub_df <- res_df[res_df$Variable == var, ]
    if (!is.null(level_order) && !is.null(level_order[[var]])) {
      sub_df <- sub_df[order(match(sub_df$Subgroup, level_order[[var]])), ]
    } else {
      sub_df <- sub_df[order(sub_df$Subgroup), ]
    }
    
    header_label <- if (!is.null(var_labels) && var %in% names(var_labels)) var_labels[var] else var
    rows[[length(rows) + 1]] <- data.frame(
      Label = header_label,
      N     = NA_integer_,
      HR    = NA_real_,
      Lower = NA_real_,
      Upper = NA_real_,
      P     = NA_real_,
      RowType = "header",
      stringsAsFactors = FALSE
    )
    
    for (i in seq_len(nrow(sub_df))) {
      rows[[length(rows) + 1]] <- data.frame(
        Label = paste0("    ", sub_df$Subgroup[i]),
        N     = sub_df$N[i],
        HR    = sub_df$HR[i],
        Lower = sub_df$Lower[i],
        Upper = sub_df$Upper[i],
        P     = sub_df$P[i],
        RowType = "subgroup",
        stringsAsFactors = FALSE
      )
    }
  }
  
  forest_df <- do.call(rbind, rows)
  if (is.null(forest_df) || nrow(forest_df) == 0) {
    stop("No rows to plot.")
  }
  
  # ---- Automated xlim and ticks ----
  if (is.null(xlim) || is.null(ticks_at)) {
    hr_vals <- forest_df$HR[!is.na(forest_df$HR) & is.finite(forest_df$HR)]
    lower_vals <- forest_df$Lower[!is.na(forest_df$Lower) & is.finite(forest_df$Lower)]
    upper_vals <- forest_df$Upper[!is.na(forest_df$Upper) & is.finite(forest_df$Upper)]
    all_vals <- c(hr_vals, lower_vals, upper_vals)
    all_vals <- all_vals[all_vals > 0]
    
    if (length(all_vals) > 0) {
      min_val <- min(all_vals, na.rm = TRUE)
      max_val <- max(all_vals, na.rm = TRUE)
      
      log_min <- log(min_val)
      log_max <- log(max_val)
      span <- log_max - log_min
      expand <- max(0.08, span * 0.6)
      new_min <- exp(log_min - expand)
      new_max <- exp(log_max + expand)
      
      new_min <- min(new_min, 1)
      new_max <- max(new_max, 1)
      
      if (new_max / new_min < 1.15) {
        mid <- sqrt(new_min * new_max)
        new_min <- mid / 1.15
        new_max <- mid * 1.15
      }
      
      if (is.null(xlim)) {
        xlim <- c(new_min, new_max)
      }
      
      if (is.null(ticks_at)) {
        log_range <- log(c(new_min, new_max))
        log_ticks <- pretty(log_range, n = 4)
        ticks <- exp(log_ticks)
        ticks <- ticks[ticks >= new_min * 0.99 & ticks <= new_max * 1.01]
        
        if (!any(abs(ticks - 1) < 0.0001)) {
          ticks <- sort(c(ticks, 1))
        }
        if (length(ticks) > 4) {
          ticks <- sort(ticks)
          keep_idx <- c(1, which.min(abs(ticks - 1)), length(ticks))
          if (length(ticks) >= 3) {
            mid_idx <- round(median(seq_along(ticks)))
            keep_idx <- sort(unique(c(keep_idx, mid_idx)))
          }
          ticks <- ticks[keep_idx]
        }
        ticks_at <- unique(round(ticks, 8))
        ticks_at <- ticks_at[ticks_at > 0]
      }
    } else {
      if (is.null(xlim)) xlim <- c(0.5, 2)
      if (is.null(ticks_at)) ticks_at <- c(0.5, 1, 2)
    }
  }
  
  # Safety checks
  if (is.null(xlim) || length(xlim) != 2 || any(!is.finite(xlim))) {
    xlim <- c(0.5, 2)
  }
  xlim <- sort(xlim)
  if (xlim[1] <= 0) xlim[1] <- 0.01
  
  if (!is.null(ticks_at)) {
    ticks_at <- ticks_at[ticks_at >= xlim[1] & ticks_at <= xlim[2]]
    ticks_at <- ticks_at[ticks_at > 0]
    if (length(ticks_at) > 5) {
      ticks_at <- sort(ticks_at)
      keep <- c(1, round(length(ticks_at)/2), length(ticks_at))
      one_idx <- which.min(abs(ticks_at - 1))
      keep <- sort(unique(c(keep, one_idx)))
      ticks_at <- ticks_at[keep]
    }
    if (length(ticks_at) < 2) {
      log_ticks <- pretty(log(c(xlim[1], xlim[2])), n = 3)
      ticks_at <- exp(log_ticks)
      ticks_at <- ticks_at[ticks_at >= xlim[1] & ticks_at <= xlim[2]]
    }
    ticks_at <- unique(sort(ticks_at))
  } else {
    log_ticks <- pretty(log(c(xlim[1], xlim[2])), n = 3)
    ticks_at <- exp(log_ticks)
    ticks_at <- ticks_at[ticks_at >= xlim[1] & ticks_at <= xlim[2]]
    ticks_at <- unique(sort(ticks_at))
  }
  
  # ---- Format columns ----
  fmt <- paste0("%.", digits, "f")
  forest_df$disp_hr <- ifelse(
    forest_df$RowType == "header", "",
    sprintf(paste0(fmt, " (", fmt, ", ", fmt, ")"),
            forest_df$HR, forest_df$Lower, forest_df$Upper)
  )
  forest_df$disp_p <- ifelse(
    forest_df$RowType == "header" | is.na(forest_df$P), "",
    formatC(forest_df$P, digits = p_digits, format = "f")
  )
  forest_df$disp_n <- ifelse(is.na(forest_df$N), "", forest_df$N)
  forest_df[[" "]] <- paste(rep(" ", 20), collapse = " ")
  
  display_df <- forest_df[, c("Label", "disp_n", " ", "disp_hr", "disp_p")]
  names(display_df) <- c(" ", "N", " ", "HR (95% CI)", "p value")
  
  is_summary <- forest_df$RowType == "header"
  
  # ---- Theme ----
  base_theme <- list(
    base_size    = 10,
    ci_pch       = 15,
    ci_col       = ci_col,
    ci_fill      = box_col,
    ci_lty       = 1,
    ci_lwd       = 1.5,
    ci_Theight   = 0.2,
    refline_lwd  = 1,
    refline_lty  = "dashed",
    refline_col  = "grey40",
    footnote_col = "grey50",
    xaxis_gp     = grid::gpar(fontsize = 9, fontfamily = "", cex = 0.9)
  )
  
  if (packageVersion("forestploter") >= "0.3.0") {
    base_theme$refline_gp <- grid::gpar(lwd = base_theme$refline_lwd,
                                        lty = base_theme$refline_lty,
                                        col = base_theme$refline_col)
    base_theme$footnote_gp <- grid::gpar(col = base_theme$footnote_col)
    base_theme$refline_lwd <- NULL
    base_theme$refline_lty <- NULL
    base_theme$refline_col <- NULL
    base_theme$footnote_col <- NULL
  }
  
  if (!is.null(arrow_labels)) {
    base_theme$arrow_lab  <- arrow_labels
    base_theme$arrow_type <- "open"
    base_theme$arrow_col  <- "black"
  }
  theme_args_final <- utils::modifyList(base_theme, theme_args)
  tm <- do.call(forestploter::forest_theme, theme_args_final)
  
  if (is.null(footnote)) {
    footnote <- ''
  }
  
  # ---- Draw ----
  p <- forestploter::forest(
    display_df,
    est          = forest_df$HR,
    lower        = forest_df$Lower,
    upper        = forest_df$Upper,
    ci_column    = 3,
    is_summary   = is_summary,
    ref_line     = ref_line,
    xlim         = xlim,
    ticks_at     = ticks_at,
    ticks_digits = ticks_digits,
    x_trans      = x_trans,
    footnote     = footnote,
    theme        = tm
  )
  print(p)
  
  # ---- Save ----
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- getwd()
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    if (is.null(file_name)) file_name <- "subgroup_forest.pdf"
    grDevices::pdf(file.path(save_dir, file_name), width = width, height = height)
    plot(p)
    grDevices::dev.off()
  }
  
  invisible(list(plot = p, data = forest_df, fits = fit_list,
                 xlim = xlim, ticks_at = ticks_at, ticks_digits = ticks_digits))
}
# ==============================================================================
# 9. Stability & Sensitivity Analysis Module
# ==============================================================================

#' Robust Feature Stability Analysis for Sparse Survival Models
#'
#' Evaluates the stability of feature selection under repeated subsampling using
#' \code{glmnet::cv.glmnet} with Cox regression. Provides the Jaccard stability
#' index and selection frequencies for each feature.
#'
#' @param object A \code{PrognosiX} object, data frame, or \code{TaskSurv}.
#' @param time_col A character string specifying the name of the time column
#'   (required if \code{object} is a data frame). Default is \code{"time"}.
#' @param status_col A character string specifying the name of the status column
#'   (required if \code{object} is a data frame). Default is \code{"status"}.
#' @param n_repeat An integer specifying the number of subsampling iterations.
#'   Default is \code{30}.
#' @param train_ratio A numeric value specifying the proportion of data to sample
#'   each iteration. Default is \code{0.8}.
#' @param alpha Elastic net mixing parameter: 1 = LASSO, 0 = Ridge, 0.5 = elastic net.
#'   Default is \code{1}.
#' @param palette_name A character string specifying the name of the Wes Anderson
#'   palette. Default is \code{"AsteroidCity1"}.
#' @param seed An integer seed for reproducibility. Default is \code{2025}.
#' @param verbose A logical value. Print progress messages. Default is \code{TRUE}.
#'
#' @return A list with four components:
#'   \describe{
#'     \item{stability_index}{Jaccard stability index (mean pairwise Jaccard similarity).}
#'     \item{frequencies}{A data frame of features and their selection frequencies.}
#'     \item{plot}{A \code{ggplot} object of the top 15 features.}
#'     \item{success_rate}{The proportion of successful iterations.}
#'   }
#'
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#'
#' stab <- surv_analyze_feature_stability(task, "time", "status", n_repeat = 5, alpha = 1)
#' stab$stability_index
#'
#' abl <- surv_analyze_feature_ablation(task, "surv.coxph",
#'                                      features_to_test = c("karno", "celltype", "trt"))
#' abl$results
#'
#' sens <- surv_analyze_model_sensitivity(task, "surv.coxph", analysis_type = "sample_size")
#' sens$results
#' }
#'
#' @seealso \code{\link{surv_feature_selection_multi}} for multi-method selection
#' @export
surv_analyze_feature_stability <- function(object,
                                           time_col = "time",
                                           status_col = "status",
                                           n_repeat = 30,
                                           train_ratio = 0.8,
                                           alpha = 1,
                                           palette_name = "AsteroidCity1",
                                           seed = 2025,
                                           verbose = TRUE) {
  
  # Extract data
  if (inherits(object, "PrognosiX")) {
    data <- object@survival.data
    time_col <- object@time_col
    status_col <- object@status_col
  } else if (is.data.frame(object)) {
    data <- object
  } else if (inherits(object, "TaskSurv")) {
    data <- as.data.frame(object$data())
    time_col <- object$target_names[1]
    status_col <- object$target_names[2]
  } else {
    stop("object must be a PrognosiX object, data frame, or TaskSurv")
  }
  
  # Required packages
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package 'glmnet' is required for stability analysis.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")
  
  set.seed(seed)
  total_samples <- nrow(data)
  feature_cols <- setdiff(colnames(data), c(time_col, status_col))
  
  # One-hot encode factor features (robust, avoids mlr3 encoding issues)
  encode_data <- function(df) {
    df_encoded <- df[, c(time_col, status_col), drop = FALSE]
    for (col in feature_cols) {
      x <- df[[col]]
      if (is.factor(x) || is.character(x)) {
        x <- as.factor(x)
        dummies <- model.matrix(~ x - 1)
        colnames(dummies) <- paste0(col, "_", levels(x))
        df_encoded <- cbind(df_encoded, dummies)
      } else {
        df_encoded[[col]] <- x
      }
    }
    df_encoded <- df_encoded[, !duplicated(colnames(df_encoded))]
    colnames(df_encoded) <- make.names(colnames(df_encoded))
    return(df_encoded)
  }
  
  data_encoded <- encode_data(data)
  encoded_features <- setdiff(colnames(data_encoded), c(time_col, status_col))
  
  if (length(encoded_features) == 0) {
    stop("No features available after encoding. Check input data.")
  }
  
  if (verbose) {
    message(sprintf("\n[*] Starting Feature Stability Analysis (%d iterations)...", n_repeat))
    message(sprintf("    Subsample ratio: %.0f%% | Elastic net alpha = %.1f", train_ratio * 100, alpha))
    message(sprintf("    Original features: %d | Encoded features: %d", 
                    length(feature_cols), length(encoded_features)))
  }
  
  selected_sets <- list()
  success_count <- 0
  
  # Progress bar if interactive
  if (verbose && interactive() && requireNamespace("utils", quietly = TRUE)) {
    pb <- utils::txtProgressBar(min = 0, max = n_repeat, style = 3)
  } else {
    pb <- NULL
  }
  
  for (i in seq_len(n_repeat)) {
    idx <- sample(total_samples, size = floor(total_samples * train_ratio))
    sub_data <- data_encoded[idx, , drop = FALSE]
    
    y <- survival::Surv(sub_data[[time_col]], sub_data[[status_col]])
    x <- as.matrix(sub_data[, encoded_features, drop = FALSE])
    
    fit <- tryCatch({
      glmnet::cv.glmnet(x = x, y = y, family = "cox", alpha = alpha, nfolds = 5)
    }, error = function(e) {
      if (verbose && i <= 3) warning(sprintf("Iter %d failed: %s", i, e$message))
      NULL
    })
    
    if (is.null(fit)) {
      selected_sets[[i]] <- character(0)
      next
    }
    
    coefs <- as.matrix(stats::coef(fit, s = "lambda.min"))
    selected <- rownames(coefs)[coefs[, 1] != 0]
    selected <- setdiff(selected, "(Intercept)")
    
    if (length(selected) == 0) selected <- "none"
    selected_sets[[i]] <- selected
    success_count <- success_count + 1
    
    if (!is.null(pb)) utils::setTxtProgressBar(pb, i)
  }
  
  if (!is.null(pb)) close(pb)
  
  if (verbose) {
    message(sprintf("\n    Successful iterations: %d / %d", success_count, n_repeat))
  }
  
  if (success_count == 0) {
    warning("No successful iterations. Returning empty result.")
    return(list(stability_index = NA, frequencies = data.frame(), plot = NULL, success_rate = 0))
  }
  
  valid_sets <- selected_sets[sapply(selected_sets, length) > 0 & 
                                sapply(selected_sets, function(x) !all(x == "none"))]
  
  if (length(valid_sets) == 0) {
    warning("No features selected in any iteration. Try increasing sample size or changing alpha.")
    return(list(stability_index = NA, frequencies = data.frame(), plot = NULL, success_rate = 0))
  }
  
  # Selection frequencies
  all_selected <- unlist(valid_sets)
  freq_tab <- sort(table(all_selected), decreasing = TRUE)
  freq_df <- data.frame(
    Feature = names(freq_tab),
    Frequency = as.numeric(freq_tab) / length(valid_sets)
  )
  
  # Jaccard stability index
  stab_index <- NA
  if (length(valid_sets) > 1) {
    jaccard_vals <- c()
    for (i in seq_len(length(valid_sets) - 1)) {
      for (j in (i + 1):length(valid_sets)) {
        inter <- length(intersect(valid_sets[[i]], valid_sets[[j]]))
        union <- length(union(valid_sets[[i]], valid_sets[[j]]))
        if (union > 0) jaccard_vals <- c(jaccard_vals, inter / union)
      }
    }
    stab_index <- mean(jaccard_vals, na.rm = TRUE)
  }
  
  # Plot top 15 features with optional wesanderson color
  plot_df <- head(freq_df, 15)
  p <- NULL
  if (nrow(plot_df) > 0) {
    # Determine bar color
    if (requireNamespace("wesanderson", quietly = TRUE)) {
      bar_color <- wesanderson::wes_palette(palette_name, n = 1, type = "continuous")
    } else {
      bar_color <- "#2980b9"  # default blue
      if (verbose) message("    Note: Install 'wesanderson' for more colorful plots.")
    }
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = stats::reorder(Feature, Frequency), y = Frequency)) +
      ggplot2::geom_bar(stat = "identity", fill = bar_color, width = 0.6) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::labs(
        x = "Feature",
        y = "Selection Frequency",
        title = "Feature Selection Stability (LASSO/Elastic Net)",
        subtitle = sprintf(
          "Subsampling ratio: %.0f%% | Iterations: %d | Jaccard Index: %.3f",
          train_ratio * 100, length(valid_sets), stab_index
        )
      ) +
      ggprism::theme_prism(base_size = 12) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    print(p)
  }
  
  if (verbose) {
    message(sprintf("  [OK] Stability index (Jaccard): %.3f", stab_index))
  }
  
  list(
    stability_index = stab_index,
    frequencies = freq_df,
    plot = p,
    success_rate = success_count / n_repeat
  )
}

#' Analyze Model Performance Sensitivity
#'
#' Evaluates how model performance (C-index) changes under varying conditions
#' such as sample size or censoring rate. This helps assess the robustness of
#' the model to data perturbations.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.coxph"}).
#' @param analysis_type A character string specifying the type of sensitivity
#'   analysis. Must be one of \code{"sample_size"} or \code{"censoring"}.
#' @param param_values A numeric vector of parameter values to test. For
#'   \code{"sample_size"}, these are proportions (e.g., \code{c(0.3, 0.5, 0.7, 0.9, 1.0)}).
#'   For \code{"censoring"}, these are additional censoring proportions
#'   (e.g., \code{c(0.1, 0.2, 0.3, 0.5)}). If \code{NULL}, uses reasonable defaults.
#' @param palette_name A character string (kept for compatibility, no longer used).
#'
#' @return A list with two components:
#'   \describe{
#'     \item{results}{A data frame with columns: \code{Parameter}, \code{C_Index}, and \code{SE}.}
#'     \item{plot}{A \code{ggplot} object showing the sensitivity trajectory.}
#'   }
#'
#' @export
surv_analyze_model_sensitivity <- function(object, learner_id, analysis_type = c("sample_size", "censoring"),
                                           param_values = NULL, palette_name = "AsteroidCity1") {

  `%>%` <- magrittr::`%>%`  
  
  task <- surv_extract_task(object)  
  
  analysis_type <- match.arg(analysis_type)
  original_data <- data.table::as.data.table(task$data())
  
  base::message(sprintf("\n[*] Starting Sensitivity Analysis: %s...", analysis_type))
  
  results <- base::data.frame(Parameter = base::numeric(), C_Index = base::numeric(), SE = base::numeric())
  
  # Default parameters if not provided
  if (base::is.null(param_values)) {
    if (analysis_type == "sample_size") param_values <- c(0.3, 0.5, 0.7, 0.9, 1.0)
    if (analysis_type == "censoring") param_values <- c(0.1, 0.2, 0.3, 0.5)
  }
  
  for (val in param_values) {
    
    cv_scores <- base::numeric(5L)
    
    for (fold in seq_len(5L)) {
      temp_data <- data.table::copy(original_data)
      
      if (analysis_type == "sample_size") {
        keep_idx  <- base::sample(base::nrow(temp_data), size = base::max(20L, base::floor(base::nrow(temp_data) * val)))
        temp_data <- temp_data[keep_idx, ]
      } else if (analysis_type == "censoring") {
        time_col   <- task$target_names[1L]
        status_col <- task$target_names[2L]
        event_idx  <- which(temp_data[[status_col]] == 1)
        n_to_censor <- base::floor(base::length(event_idx) * val)
        if (n_to_censor > 0L) {
          censor_idx <- base::sample(event_idx, size = n_to_censor, replace = FALSE)
          new_times  <- stats::runif(n_to_censor, min = 0,
                                     max = temp_data[[time_col]][censor_idx])
          new_times  <- base::pmax(new_times, .Machine$double.eps)
          temp_data[[time_col]][censor_idx]   <- new_times
          temp_data[[status_col]][censor_idx] <- 0L
        }
      }

      temp_task <- surv_create_surv_task(temp_data, task$target_names[1], task$target_names[2], id = "temp_sens")
      
      res <- base::tryCatch({
        rr <- mlr3::resample(temp_task, mlr3::lrn(learner_id), mlr3::rsmp("cv", folds = 3), store_models = FALSE)
        rr$aggregate(mlr3::msr("surv.cindex"))
      }, error = function(e) NA)
      
      cv_scores[fold] <- res
    }
    
    results <- base::rbind(results, base::data.frame(
      Parameter = val,
      C_Index = base::mean(cv_scores, na.rm = TRUE),
      SE = stats::sd(cv_scores, na.rm = TRUE) / base::sqrt(base::sum(!base::is.na(cv_scores)))
    ))
  }
  
  x_label <- ifelse(analysis_type == "sample_size",
                    "Proportion of Total Sample Size",
                    "Additional Censoring Proportion (fraction of events re-censored)")
  
  p <- ggplot2::ggplot(results, ggplot2::aes(x = Parameter, y = C_Index)) +
    ggplot2::geom_line(color = "black", size = 1) +
    ggplot2::geom_point(color = "black", size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = C_Index - SE, ymax = C_Index + SE),
                           width = 0.01, color = "grey50") +
    ggplot2::labs(
      x = x_label, 
      y = "Cross-Validated C-Index (+/- SE)", 
      title = "Model Sensitivity Analysis",
      subtitle = sprintf("Model: %s | Perturbation: %s", learner_id, analysis_type)
    ) +
    ggprism::theme_prism(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  base::print(p)
  base::message("  [OK] Sensitivity analysis complete.")
  base::list(results = results, plot = p)
}
#' Feature Ablation Sensitivity Analysis
#'
#' Evaluates the impact of removing individual features on model performance.
#' For each feature, the model is retrained without that feature and the change
#' in C-index is recorded.
#'
#' @param object A \code{TaskSurv} or \code{PrognosiX} object.
#' @param learner_id A character string specifying the learner ID (e.g., \code{"surv.coxph"}).
#' @param features_to_test A character vector of feature names to test. If \code{NULL},
#'   tests all features in the task. Default is \code{NULL}.
#'
#' @return A list with three components:
#'   \describe{
#'     \item{results}{A data frame with columns: \code{Feature_Removed}, \code{New_CIndex},
#'       and \code{Performance_Drop}.}
#'     \item{plot}{A \code{ggplot} object showing the performance drop for each feature.}
#'     \item{baseline}{The baseline C-index with all features.}
#'   }
#' @importFrom mlr3 rsmp msr resample 
#' @details
#' The function:
#' \enumerate{
#'   \item Calculates baseline performance using all features.
#'   \item For each feature, creates a task without that feature and evaluates performance.
#'   \item Reports the drop in C-index for each feature.
#'   \item Creates a bar plot sorted by impact magnitude.
#' }
#'
#' @examples
#' \dontrun{
#' library(mlr3proba)
#' library(survival)
#' 
#' data("veteran", package = "survival")
#' task <- surv_create_surv_task(veteran, "time", "status")
#' 
#' # Ablation analysis on all features
#' ablation <- surv_analyze_feature_ablation(
#'   object = task,
#'   learner_id = "surv.coxph"
#' )
#' print(ablation$plot)
#' 
#' # Test specific features
#' ablation_subset <- surv_analyze_feature_ablation(
#'   object = task,
#'   learner_id = "surv.coxph",
#'   features_to_test = c("age", "karno", "celltype")
#' )
#' }
#'
#' @export
surv_analyze_feature_ablation <- function(object, learner_id, features_to_test = NULL) {
  task <- surv_extract_task(object)
  
  if (is.null(features_to_test)) {
    features_to_test <- task$feature_names
  }
  
  message(sprintf("\n[*] Starting Feature Ablation Analysis for %d features...", length(features_to_test)))
  
  # 1. Calculate Baseline Performance (using all features)
  learner <- surv_get_learner(learner_id, task)
  
  # Use a small 3-fold CV to get a stable baseline
  resampling <- mlr3::rsmp("cv", folds = 3)
  baseline_rr <- resample(task, learner, resampling, store_models = FALSE)
  baseline_cindex <- baseline_rr$aggregate(msr("surv.cindex"))
  
  results <- list()
  
  # 2. Iterate and Remove Each Feature
  for (feat in features_to_test) {
    # Create a task without the specific feature
    remaining_feats <- setdiff(task$feature_names, feat)
    temp_task <- task$clone()$select(remaining_feats)
    
    # Evaluate model performance without this feature
    temp_rr <- tryCatch({
      resample(temp_task, learner, resampling, store_models = FALSE)
    }, error = function(e) NULL)
    
    if (!is.null(temp_rr)) {
      new_cindex <- temp_rr$aggregate(msr("surv.cindex"))
      drop <- baseline_cindex - new_cindex
      
      results[[feat]] <- data.frame(
        Feature_Removed = feat,
        New_CIndex = new_cindex,
        Performance_Drop = drop
      )
    }
  }
  
  res_df <- do.call(rbind, results)
  res_df <- res_df[order(-res_df$Performance_Drop), ] # Sort by impact
  
  # 3. Generate Visualization
  p <- ggplot2::ggplot(res_df, ggplot2::aes(x = stats::reorder(Feature_Removed, Performance_Drop), y = Performance_Drop)) +
    ggplot2::geom_bar(stat = "identity", fill = ifelse(res_df$Performance_Drop > 0, "#E67E22", "#95A5A6")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = "Feature Removed",
      y = "Drop in C-Index (Baseline - New)",
      title = "Feature Ablation Sensitivity Analysis",
      subtitle = sprintf("Baseline C-Index (Full Model): %.4f", baseline_cindex)
    ) +
    ggprism::theme_prism()
  
  print(p)
  message("  [OK] Feature ablation analysis complete.")
  
  return(list(results = res_df, plot = p, baseline = baseline_cindex))
}

# ==============================================================================
# 11. Clinical Validation: Calibration & Time-dependent AUC
#' Plot Calibration Curve (Apparent)
#'
#' Plots predicted vs. observed (KM-based) survival probability at a fixed
#' evaluation time, using adaptively-sized quantile bins.
#'
#' @param learner Trained mlr3 learner (must support "distr" predict_type)
#' @param object A \code{TaskSurv} or \code{PrognosiX} object
#' @param time_point Numeric. Evaluation time point for calibration.
#' @param n_bins Integer. Requested number of bins (default 10; automatically
#'   reduced for small datasets, see \code{.prepare_cal_data}).
#' @param apparent Logical. Whether this is apparent (training-set,
#'   optimistic) calibration (default TRUE).
#' @param print_metrics Logical. Print calibration metrics to console (default TRUE).
#' @param show_ici Logical. Show ICI in the plot subtitle (default TRUE).
#' @return A ggplot object (invisibly carrying \code{calibration_metrics} and
#'   \code{calibration_data} as attributes), or \code{NULL} if calibration
#'   could not be computed.
#' @export
#'
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' veteran$celltype <- as.character(veteran$celltype)
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#' lrn <- surv_get_learner("surv.ranger", task)
#' lrn$train(task)
#'
#' cal <- surv_plot_calibration(lrn, task, time_point = 90)
#'
#' set.seed(1)
#' idx <- sample(seq_len(task$nrow), size = floor(0.7 * task$nrow))
#' train_task <- task$clone()$filter(idx)
#' val_task   <- task$clone()$filter(setdiff(seq_len(task$nrow), idx))
#' lrn2 <- surv_get_learner("surv.ranger", train_task)
#' lrn2$train(train_task)
#' cal_cmp <- surv_plot_comparison_calibration(lrn2, train_task, val_task, time_point = 90)
#' }
surv_plot_calibration <- function(learner, object, time_point, n_bins = 10,
                                  apparent = TRUE,
                                  print_metrics = TRUE, show_ici = TRUE) {
  if (!("distr" %in% learner$predict_types)) {
    warning("Learner does not support 'distr' predict_type. Calibration skipped.")
    return(NULL)
  }
  old_predict_type <- learner$predict_type
  on.exit({ learner$predict_type <- old_predict_type }, add = TRUE)
  
  cal_df <- .prepare_cal_data(learner = learner, object = object,
                              time_point = time_point, n_bins = n_bins)
  if (is.null(cal_df)) return(NULL)
  
  if (apparent) {
    message(paste(
      "[Calibration] Evaluating on the TRAINING SET.",
      "This is APPARENT (optimistic) calibration.",
      "For valid calibration use an independent test set or CV out-of-fold predictions."
    ))
  }
  
  metrics <- .compute_cal_metrics(cal_df)
  
  if (print_metrics) {
    cat("\n========== Calibration Metrics at t =", time_point, "==========\n")
    if (apparent) cat("** APPARENT (training-set) calibration -- interpret with caution **\n")
    cat(sprintf("Calibration slope    (ideal 1): %.4f\n", metrics$slope))
    cat(sprintf("Calibration intercept (ideal 0): %.4f\n", metrics$intercept))
    cat(sprintf("R-squared:                       %.4f\n", metrics$r_squared))
    cat(sprintf("Mean Absolute Error (MAE):       %.4f\n", metrics$mae))
    cat(sprintf("Integrated Calibration Index:    %s\n",
                ifelse(is.na(metrics$ici), "NA", sprintf("%.4f", metrics$ici))))
    cat(sprintf("E50 (median abs error):          %.4f\n", metrics$e50))
    cat(sprintf("E90 (90th pct abs error):        %.4f\n", metrics$e90))
    cat("====================================================\n")
  }
  
  slope_label <- round(metrics$slope, 3)
  ici_label   <- ifelse(is.na(metrics$ici), "NA", round(metrics$ici, 3))
  plot_title  <- ifelse(apparent,
                        paste("Apparent Calibration at t =", time_point, "(training set)"),
                        paste("Calibration Curve at t =", time_point))
  sub_title   <- ifelse(show_ici,
                        paste(learner$id, "| Slope =", slope_label, "| ICI =", ici_label),
                        learner$id)
  
  p <- ggplot2::ggplot(cal_df, ggplot2::aes(x = predicted, y = observed)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "red", linewidth = 0.8) +
    ggplot2::geom_line(color = "#2980b9", linewidth = 0.9) +
    ggplot2::geom_point(size = 3, color = "#2980b9") +
    ggplot2::labs(x = "Predicted Survival Probability",
                  y = "Observed Survival Probability (KM)",
                  title = plot_title, subtitle = sub_title) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) # Set display range without clipping geom_line paths; no forced 1:1 aspect
  
  if (requireNamespace("ggprism", quietly = TRUE)) {
    p <- p + ggprism::theme_prism()
  } else {
    p <- p + ggplot2::theme_minimal()
  }
  
  attr(p, "calibration_metrics") <- metrics
  attr(p, "calibration_data")    <- cal_df
  return(p)
}


#' Plot Comparison Calibration (Apparent Training vs Validation)
#'
#' Overlays calibration curves for the training set (apparent/optimistic) and
#' a held-out validation set on the same plot.
#'
#' @param learner Trained mlr3 learner (must support "distr" predict_type)
#' @param train_task Training TaskSurv object
#' @param val_task Validation TaskSurv object
#' @param time_point Numeric. Evaluation time point for calibration.
#' @param n_bins Integer. Requested number of bins (default 10; each of
#'   \code{train_task} and \code{val_task} is independently adjusted for its
#'   own sample size, see \code{.prepare_cal_data} -- a message is printed
#'   when the two datasets end up with different effective bin counts).
#' @param print_metrics Logical. Print metrics for both datasets (default TRUE).
#' @return A ggplot object showing calibration curves for training and validation.
#' @export
#'
#' @examples
#' \dontrun{
#' comp_plot <- surv_plot_comparison_calibration(learner, train_task, val_task, time_point = 365)
#' }
surv_plot_comparison_calibration <- function(learner, train_task, val_task,
                                             time_point, n_bins = 10,
                                             print_metrics = TRUE) {
  if (!("distr" %in% learner$predict_types)) {
    warning("Learner does not support 'distr' predict_type. Calibration skipped.")
    return(NULL)
  }
  old_predict_type <- learner$predict_type
  on.exit({ learner$predict_type <- old_predict_type }, add = TRUE)
  
  cal_train <- .prepare_cal_data(learner, train_task, time_point, n_bins)
  cal_val   <- .prepare_cal_data(learner, val_task, time_point, n_bins)
  
  if (is.null(cal_train) || is.null(cal_val)) {
    warning("Insufficient bins or missing columns for calibration comparison.")
    return(NULL)
  }
  
  m_train <- .compute_cal_metrics(cal_train)
  m_val   <- .compute_cal_metrics(cal_val)
  
  if (print_metrics) {
    cat("\n========== Calibration Metrics at t =", time_point, "==========\n")
    cat("** TRAINING (apparent -- optimistic for flexible learners): **\n")
    cat(sprintf("  Slope=%.4f | Intercept=%.4f | R2=%.4f | MAE=%.4f | ICI=%.4f\n",
                m_train$slope, m_train$intercept, m_train$r_squared, m_train$mae, m_train$ici))
    cat("** VALIDATION (unbiased estimate): **\n")
    cat(sprintf("  Slope=%.4f | Intercept=%.4f | R2=%.4f | MAE=%.4f | ICI=%.4f\n",
                m_val$slope, m_val$intercept, m_val$r_squared, m_val$mae, m_val$ici))
    cat("  Note: Training calibration is always expected to be more optimistic.\n")
    if (nrow(cal_train) != nrow(cal_val)) {
      cat(sprintf("  Note: Effective bin counts differ (train=%d, val=%d) due to per-set adaptive binning.\n",
                  nrow(cal_train), nrow(cal_val)))
    }
    cat("=================================================================\n")
  }
  
  cal_train$Dataset <- "Training (apparent)"
  cal_val$Dataset   <- "Validation"
  cal_all <- rbind(cal_train, cal_val)
  
  ici_tr  <- ifelse(is.na(m_train$ici), "NA", round(m_train$ici, 3))
  ici_val <- ifelse(is.na(m_val$ici),   "NA", round(m_val$ici,   3))
  
  p <- ggplot2::ggplot(cal_all, ggplot2::aes(x = predicted, y = observed, color = Dataset)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "gray40", linewidth = 0.8) +
    ggplot2::geom_line(linewidth = 1.0) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = c("Training (apparent)" = "#2980b9", "Validation" = "#e74c3c")) +
    ggplot2::labs(x = "Predicted Survival Probability",
                  y = "Observed Survival Probability (KM)",
                  title = paste("Calibration Comparison at t =", time_point),
                  subtitle = paste0(learner$id, " | Train ICI = ", ici_tr, " (apparent) | Val ICI = ", ici_val)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) # Set display range without clipping geom_line paths; no forced 1:1 aspect
  
  if (requireNamespace("ggprism", quietly = TRUE)) {
    p <- p + ggprism::theme_prism()
  } else {
    p <- p + ggplot2::theme_minimal()
  }
  
  return(p)
}

#' Plot Continuous Time-dependent AUC (Dynamic AUC)
#' @param learner Trained mlr3 learner
#' @param object A \code{TaskSurv} or \code{PrognosiX} object. The function
#'   extracts the task using \code{surv_extract_task()}.
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' veteran$celltype <- as.character(veteran$celltype)
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#' lrn <- surv_get_learner("surv.ranger", task)
#' lrn$train(task)
#'
#' auc_plot <- surv_plot_time_dependent_auc(lrn, task)
#'
#' # Comparison needs two tasks (e.g. train vs validation split):
#' set.seed(1)
#' idx <- sample(seq_len(task$nrow), size = floor(0.7 * task$nrow))
#' train_task <- task$clone()$filter(idx)
#' val_task   <- task$clone()$filter(setdiff(seq_len(task$nrow), idx))
#' lrn2 <- surv_get_learner("surv.ranger", train_task)
#' lrn2$train(train_task)
#' auc_cmp <- surv_plot_comparison_auc(lrn2, train_task, val_task)
#' }
surv_plot_time_dependent_auc <- function(learner, object) {
  task <- surv_extract_task(object)
  if (!requireNamespace("risksetROC", quietly = TRUE)) stop("Please install 'risksetROC'")
  
  data <- as.data.frame(task$data())
  pred <- learner$predict(task)
  
  # Extract Time, Status, and Marker (linear predictor or crank)
  # risksetROC needs 'marker' where higher value means higher risk
  stime <- data[[task$target_names[1]]]
  status <- data[[task$target_names[2]]]
  marker <- pred$crank 
  
  # Get unique event times for plotting (as in your snippet)
  utimes <- sort(unique(stime[status == 1]))
  
  message("[*] Calculating Time-dependent AUC across all event times...")
  
  # Use risksetROC to calculate AUC at each time point
  out <- risksetROC::risksetAUC(
    Stime = stime,
    status = status,
    marker = marker,
    tmax = max(stime) * 0.95, # Avoid instability at the very end
    plot = FALSE
  )
  
  # Prepare data for ggplot
  auc_df <- data.frame(
    times = utimes,
    tAUC = out$AUC[match(utimes, out$utimes)]
  )
  # Remove NAs
  auc_df <- na.omit(auc_df)
  
  # Create the plot similar to your provided ggplot snippet
  p <- ggplot2::ggplot(auc_df, ggplot2::aes(x = times, y = tAUC)) +
    ggplot2::geom_step(direction = "vh", color = "#2C3E50", size = 1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0.5, ymax = tAUC), fill = "#3498DB", alpha = 0.1) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ggplot2::labs(
      x = "Evaluation Time Points",
      y = "AUC",
      title = "Time-Dependent AUC (Dynamic Accuracy)",
      subtitle = sprintf("Model: %s | Integrated AUC: %.3f", learner$id, mean(auc_df$tAUC))
    ) +
    ggplot2::ylim(0.4, 1.0) +
    ggprism::theme_prism()
  
  print(p)
  return(auc_df)
}

#' Univariate Cox Filtering with Support for Clinical Metadata
#'
#' Performs univariate Cox regression on each feature (or clinical variable) and
#' selects those with p‑value below a threshold. When a \code{PrognosiX} object
#' is given, both \code{survival.data} and \code{info.data} are merged to allow
#' filtering on clinical variables.
#'
#' @param object A \code{PrognosiX} or \code{TaskSurv} object, or a data frame.
#' @param p_threshold Numeric; p‑value cutoff (default 0.05).
#' @param features Character vector of feature names to test. If \code{NULL},
#'   all numeric columns (plus any specified clinical variables) are used.
#' @param include_clinical Logical; if \code{TRUE} and \code{object} is
#'   \code{PrognosiX}, also test all columns in \code{info.data} except
#'   \code{time} and \code{status}. Default \code{FALSE}.
#'
#' @return A list with components:
#'   \item{task}{A \code{TaskSurv} object containing only significant features.}
#'   \item{table}{Data frame of univariate results for all tested variables.}
#'   \item{plot}{A ggplot object showing -log10(p) for each variable.}
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' task <- surv_create_surv_task(veteran, "time", "status")
#' uni <- surv_filter_features_clinical(task, p_threshold = 0.1)
#' }
surv_filter_features_clinical <- function(object,
                                          p_threshold = 0.05,
                                          features = NULL,
                                          include_clinical = FALSE) {
  
  message("[*] Starting Feature Selection: Univariate Cox Filtering...")
  
  # ---- 1. Extract data ----
  if (inherits(object, "PrognosiX")) {
    surv_df <- as.data.frame(object@survival.data)
    info_df <- as.data.frame(object@info.data)
    # Merge by row names
    common <- intersect(rownames(surv_df), rownames(info_df))
    if (length(common) == 0) stop("No matching rows between survival.data and info.data.")
    surv_df <- surv_df[common, , drop = FALSE]
    info_df <- info_df[common, , drop = FALSE]
    # Combine, but avoid duplicating time/status if already in info
    # Remove time/status columns from info if they are already in surv_df
    time_col <- object@time_col %||% "time"
    status_col <- object@status_col %||% "status"
    info_cols <- setdiff(colnames(info_df), c(time_col, status_col))
    data <- cbind(surv_df, info_df[, info_cols, drop = FALSE])
    target_time <- time_col
    target_status <- status_col
  } else if (inherits(object, "TaskSurv")) {
    data <- as.data.frame(object$data())
    target_time <- object$target_names[1]
    target_status <- object$target_names[2]
  } else if (is.data.frame(object)) {
    data <- object
    target_time <- "time"
    target_status <- "status"
  } else {
    stop("object must be a PrognosiX, TaskSurv, or data.frame.")
  }
  
  # Determine candidate features
  all_cols <- colnames(data)
  if (is.null(features)) {
    # Use all numeric columns except time/status
    numeric_cols <- names(which(sapply(data, is.numeric)))
    features <- setdiff(numeric_cols, c(target_time, target_status))
    if (include_clinical && inherits(object, "PrognosiX")) {
      # Also include all columns from info.data (non-numeric allowed)
      clinical_vars <- setdiff(colnames(info_df), c(target_time, target_status, features))
      # Convert character/factor to numeric? We'll test as-is; coxph can handle factors.
      features <- c(features, clinical_vars)
    }
  }
  
  # Remove any that are not in data
  features <- intersect(features, colnames(data))
  if (length(features) == 0) stop("No valid features to test.")
  
  # ---- 2. Run univariate Cox ----
  unv_results <- list()
  for (feat in features) {
    formula_str <- sprintf("survival::Surv(%s, %s) ~ %s", target_time, target_status, feat)
    fit <- tryCatch({
      survival::coxph(as.formula(formula_str), data = data)
    }, error = function(e) NULL)
    if (!is.null(fit)) {
      s <- summary(fit)
      p_val <- s$coefficients[1, "Pr(>|z|)"]
      unv_results[[feat]] <- data.frame(
        Feature = feat,
        HR = s$conf.int[1, "exp(coef)"],
        P_Value = p_val,
        Lower_CI = s$conf.int[1, "lower .95"],
        Upper_CI = s$conf.int[1, "upper .95"],
        stringsAsFactors = FALSE
      )
    }
  }
  unv_df <- do.call(rbind, unv_results)
  
  # ---- 3. Select significant features ----
  significant_feats <- unv_df$Feature[unv_df$P_Value < p_threshold]
  message(sprintf("  [[OK]] Univariate filter: %d -> %d features", length(features), length(significant_feats)))
  
  # ---- 4. Create filtered task if object is a TaskSurv or PrognosiX ----
  if (inherits(object, "PrognosiX") || inherits(object, "TaskSurv")) {
    # We need to produce a new TaskSurv with only significant features
    # For PrognosiX, we can create a temporary data frame
    if (inherits(object, "PrognosiX")) {
      # Use the merged data (surv_df only) but keep only significant features
      surv_data <- object@survival.data
      # We need to select columns from surv_data (which has all features and time/status)
      # But time/status are not in features, we keep them
      keep_cols <- c(target_time, target_status, significant_feats)
      keep_cols <- intersect(keep_cols, colnames(surv_data))
      new_surv <- surv_data[, keep_cols, drop = FALSE]
      # Create a new PrognosiX? Better to return the task directly.
      task <- surv_create_surv_task(new_surv, target_time, target_status)
    } else {
      task <- object$clone()$select(significant_feats)
    }
  } else {
    task <- NULL
  }
  
  # ---- 5. Plot ----
  p <- ggplot2::ggplot(unv_df, ggplot2::aes(x = stats::reorder(Feature, -P_Value), y = -log10(P_Value))) +
    ggplot2::geom_segment(ggplot2::aes(xend = Feature, yend = 0), color = "grey") +
    ggplot2::geom_point(ggplot2::aes(color = P_Value < p_threshold), size = 3) +
    ggplot2::geom_hline(yintercept = -log10(p_threshold), linetype = "dashed", color = "red") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Univariate Feature Filtering", y = "-log10(P-value)", x = "Features") +
    ggprism::theme_prism()
  
  print(p)
  
  return(list(task = task, table = unv_df, plot = p))
}

#' Run Algorithm Benchmark for Survival Analysis
#'
#' Compares multiple survival learners using cross-validation and visualizes
#' the C-index performance with a boxplot.
#'
#' @param object A \code{PrognosiX} object or a \code{TaskSurv}.
#' @param learners_list A list of \code{mlr3} learners. If \code{NULL}, a default
#'   set of learners is used (CoxPH, Lasso, RandomForest, XGBoost).
#' @param resampling An \code{mlr3} resampling strategy. Defaults to 5-fold CV.
#' @param palette Character or character vector. Either a ColorBrewer palette name
#'   (e.g. \code{"Set2"}, \code{"Dark2"}, \code{"Paired"}) or a vector of colors.
#'   If \code{NULL} (default), a custom clean palette is used.
#'
#' @return A list containing:
#'   \item{bmr}{The \code{BenchmarkResult} object.}
#'   \item{table}{Aggregated performance table.}
#'   \item{plot}{A \code{ggplot} object of the performance comparison.}
#'
#' @export
surv_run_algorithm_benchmark <- function(object,
                                         learners_list = NULL,
                                         resampling = NULL,
                                         palette = NULL) {
  task <- surv_extract_task(object)
  requireNamespace("mlr3viz", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("ggprism", quietly = TRUE)
  
  if (is.null(learners_list)) {
    learners_list <- list(
      lrn("surv.coxph",     id = "CoxPH"),
      lrn("surv.cv_glmnet", id = "Lasso"),
      lrn("surv.ranger",    id = "RandomForest"),
      lrn("surv.xgboost",   id = "XGBoost")
    )
  }
  
  # Ensure all learners use the same predict_type for fair comparison
  for (l in learners_list) {
    if ("distr" %in% l$predict_types) l$predict_type <- "distr"
  }
  
  if (is.null(resampling)) resampling <- mlr3::rsmp("cv", folds = 5)
  
  message("[*] Running Benchmark: Comparing Algorithms via 5-fold CV...")
  
  design <- benchmark_grid(task, learners_list, resampling)
  bmr    <- benchmark(design)
  
  # Measure performance
  measures <- list(msr("surv.cindex"))
  perf_tab <- bmr$aggregate(measures)
  
  # ---------- Color handling ----------
  n_learners <- length(learners_list)
  
  if (is.null(palette)) {
    # Default clean palette
    cols <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860")
  } else if (length(palette) == 1L &&
             requireNamespace("RColorBrewer", quietly = TRUE) &&
             palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    # ColorBrewer palette name
    max_col <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
    n_col   <- min(max(3L, n_learners), max_col)
    cols    <- RColorBrewer::brewer.pal(n_col, palette)
  } else {
    # User-supplied color vector
    cols <- palette
  }
  
  # Recycle if needed
  if (length(cols) < n_learners) {
    cols <- rep(cols, length.out = n_learners)
  }
  
  # Visualization
  p <- autoplot(bmr, measure = msr("surv.cindex")) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      title    = "Algorithm Performance Comparison",
      subtitle = "5-Fold Cross-Validation C-Index",
      x        = NULL,
      y        = "C-Index"
    ) +
    ggprism::theme_prism() +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1, vjust = 1),
      legend.position = "none",
      plot.margin     = ggplot2::margin(t = 10, r = 10, b = 25, l = 10)
    )
  
  print(p)
  
  return(list(bmr = bmr, table = perf_tab, plot = p))
}

# ==============================================================================
# 12. Pipeline Helper Functions
# ==============================================================================

#' Check Data Quality for Survival Analysis
#' @param data Data frame to check
#' @param time_col Time column name
#' @param event_col Event/status column name
#' @return Invisible NULL, stops on error
#' @keywords internal
check_data_quality <- function(data, time_col, event_col) {
  # Check if data is a data frame
  if (!is.data.frame(data)) {
    stop("train_data must be a data frame.")
  }
  
  # Check required columns exist
  if (!time_col %in% colnames(data)) {
    stop(sprintf("Time column '%s' not found in data.", time_col))
  }
  if (!event_col %in% colnames(data)) {
    stop(sprintf("Event column '%s' not found in data.", event_col))
  }
  
  # Check for missing values in target columns
  time_missing <- sum(is.na(data[[time_col]]))
  event_missing <- sum(is.na(data[[event_col]]))
  
  if (time_missing > 0) {
    warning(sprintf("Time column has %d missing values.", time_missing))
  }
  if (event_missing > 0) {
    warning(sprintf("Event column has %d missing values.", event_missing))
  }
  
  # Check time is numeric
  if (!is.numeric(data[[time_col]])) {
    stop("Time column must be numeric.")
  }
  
  # Check event is binary (0/1)
  event_vals <- unique(na.omit(data[[event_col]]))
  if (!all(event_vals %in% c(0, 1))) {
    warning("Event column should be binary (0/1). Found values: ", paste(event_vals, collapse = ", "))
  }
  
  # Check for non-finite times
  non_finite <- sum(!is.finite(data[[time_col]]))
  if (non_finite > 0) {
    stop(sprintf("Time column has %d non-finite values (Inf, -Inf, or NaN).", non_finite))
  }
  
  # Check for negative times
  neg_times <- sum(data[[time_col]] < 0, na.rm = TRUE)
  if (neg_times > 0) {
    warning(sprintf("Time column has %d negative values.", neg_times))
  }
  
  # Check censoring rate
  event_rate <- mean(data[[event_col]] == 1, na.rm = TRUE)
  if (event_rate < 0.05) {
    warning(sprintf("Very low event rate (%.1f%%). Model may be unstable.", event_rate * 100))
  }
  if (event_rate > 0.95) {
    warning(sprintf("Very high event rate (%.1f%%). Check data coding.", event_rate * 100))
  }
  
  message(sprintf("[[OK]] Data quality check passed. N=%d, Events=%.1f%%", 
                  nrow(data), event_rate * 100))
  invisible(NULL)
}

#' Create Step Directory for Pipeline
#' @param base_dir Base output directory
#' @param step_num Step number
#' @param step_name Step name
#' @return Full path to created directory
#' @keywords internal
create_step_dir <- function(base_dir, step_num, step_name) {
  dir_name <- sprintf("Step%02d_%s", step_num, step_name)
  full_path <- file.path(base_dir, dir_name)
  if (!dir.exists(full_path)) {
    dir.create(full_path, recursive = TRUE, showWarnings = FALSE)
  }
  message(sprintf("  -> Step %d: %s", step_num, step_name))
  return(full_path)
}


# ==============================================================================
# 13. Multi-Dataset Validation & Comparison Module
# ==============================================================================

#' Evaluate Model on a New Dataset (Validation/Test)
#' @param learner Trained mlr3 learner
#' @param test_data New data frame (validation/test set)
#' @param task_ref Reference task to ensure column consistency
#' @return Prediction object
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#'
#' set.seed(1)
#' idx <- sample(seq_len(task$nrow), size = floor(0.7 * task$nrow))
#' train_task <- task$clone()$filter(idx)
#' val_data   <- as.data.frame(task$data())[setdiff(seq_len(task$nrow), idx), ]
#'
#' lrn <- surv_get_learner("surv.coxph", train_task)
#' lrn$train(train_task)
#'
#' res <- surv_predict_on_validation(lrn, val_data, train_task)
#' res$prediction$score(mlr3::msr("surv.cindex"))
#' }
surv_predict_on_validation <- function(learner, test_data, task_ref) {
  # Ensure the test data has the same structure/encoding as training data
  test_task <- surv_create_surv_task(
    data = test_data, 
    time_col = task_ref$target_names[1], 
    event_col = task_ref$target_names[2],
    id = "validation_set"
  )
  
  # Perform prediction
  pred <- learner$predict(test_task)
  return(list(task = test_task, prediction = pred))
}

#' Plot Comparison Time-dependent AUC (Train vs Validation)
#' @param learner Trained learner
#' @param train_task Training task
#' @param val_task Validation task
#' @return ggplot object
#' @export
surv_plot_comparison_auc <- function(learner, train_task, val_task) {
  if (!requireNamespace("risksetROC", quietly = TRUE)) stop("Please install 'risksetROC'")
  
  # Calculate for Training
  pred_train <- learner$predict(train_task)
  auc_train <- risksetROC::risksetAUC(
    Stime = train_task$data()[[train_task$target_names[1]]],
    status = train_task$data()[[train_task$target_names[2]]],
    marker = pred_train$crank,
    tmax = max(train_task$data()[[train_task$target_names[1]]]) * 0.9,
    plot = FALSE
  )
  
  # Calculate for Validation
  pred_val <- learner$predict(val_task)
  auc_val <- risksetROC::risksetAUC(
    Stime = val_task$data()[[val_task$target_names[1]]],
    status = val_task$data()[[val_task$target_names[2]]],
    marker = pred_val$crank,
    tmax = max(val_task$data()[[val_task$target_names[1]]]) * 0.9,
    plot = FALSE
  )
  
  # Combine Data
  df_train <- data.frame(Time = auc_train$utimes, AUC = auc_train$AUC, Group = "Training Set")
  df_val <- data.frame(Time = auc_val$utimes, AUC = auc_val$AUC, Group = "Validation Set")
  df_all <- rbind(df_train, df_val)
  
  p <- ggplot2::ggplot(df_all, ggplot2::aes(x = Time, y = AUC, color = Group)) +
    ggplot2::geom_step(size = 1) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
    ggplot2::scale_color_manual(values = c("Training Set" = "#2980b9", "Validation Set" = "#e74c3c")) +
    ggplot2::labs(title = "Time-Dependent AUC Comparison", subtitle = learner$id) +
    ggprism::theme_prism()
  
  print(p)
  return(p)
}

#' Multi-strategy feature selection for survival analysis
#'
#' Applies multiple feature selection methods to a survival task and returns a
#' consensus set of features. Supports 12+ algorithms including univariate Cox,
#' penalized Cox (LASSO, Ridge, Elastic Net), random forest importance,
#' XGBoost importance, VIMP, stepwise, stability selection,  and Boruta.
#'
#' @param object A \code{PrognosiX} object or a \code{TaskSurv} (mlr3 survival task).
#' @param methods Character vector of method names to apply. Available methods:
#'   \describe{
#'     \item{\code{"uni_cox"}}{Univariate Cox regression; keeps features with p < \code{p_threshold}.}
#'     \item{\code{"lasso"}}{LASSO penalized Cox (alpha = 1, lambda.min).}
#'     \item{\code{"ridge"}}{Ridge penalized Cox (alpha = 0, lambda.min).}
#'     \item{\code{"enet"}}{Elastic net (alpha = 0.5, lambda.min).}
#'     \item{\code{"rf_imp"}}{Random forest (ranger) permutation importance; keeps top \code{top_ratio} features.}
#'     \item{\code{"rfsrc_imp"}}{Random survival forest (randomForestSRC) importance; keeps top \code{top_ratio}.}
#'     \item{\code{"xgb_imp"}}{XGBoost gain importance; keeps top \code{top_ratio}.}
#'     \item{\code{"vimp"}}{VIMP variable importance from \code{randomForestSRC::vimp} (recommended).}
#'     \item{\code{"boruta"}}{Boruta wrapper algorithm (requires \code{Boruta} package; disabled by default).}
#'     \item{\code{"stepwise"}}{Stepwise Cox regression (both directions, AIC). Only for low-dimensional data (p < 30).}
#'     \item{\code{"stab_sel"}}{Stability selection using \code{c060::stabpath} with Lasso.}
#'   }
#' @param p_threshold Numeric. P-value threshold for univariate Cox (default 0.05).
#' @param top_ratio Numeric. For importance-based methods (RF, XGBoost, VIMP), keep this proportion of top features (default 0.5).
#' @param combine Character. How to combine results from different methods:
#'   \itemize{
#'     \item \code{"union"} - take union of all selected feature sets.
#'     \item \code{"intersection"} - take intersection (common features).
#'     \item \code{"freq"} - keep features selected by at least \code{freq_cutoff} methods.
#'   }
#' @param freq_cutoff Integer. Minimum number of methods that must select a feature when \code{combine = "freq"} (default 2).
#' @param verbose Logical. Print progress messages (default TRUE).
#' @param use_boruta Logical. Enable Boruta (can be slow and may fail on survival data). Default FALSE.
#'
#' @return A list with three components:
#'   \describe{
#'     \item{\code{selected}}{Character vector of finally selected feature names.}
#'     \item{\code{method_table}}{Data frame with rows = all features, columns = each method, indicating selection status (TRUE/FALSE).}
#'     \item{\code{method_results}}{List of raw outputs from each method (e.g., fitted models, importance vectors) for further inspection.}
#'   }
#'
#' @importFrom stats as.formula
#' @importFrom utils head
#' @importFrom survival Surv coxph
#' @importFrom MASS stepAIC
#' @importFrom mlr3 lrn
#' @importFrom mlr3proba TaskSurv
#' @importFrom randomForestSRC vimp
#' @importFrom c060 stabpath
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#'
#' feat_sel <- surv_feature_selection_multi(
#'   object = task, methods = c("uni_cox", "lasso"),
#'   p_threshold = 0.1, combine = "union", verbose = FALSE
#' )
#' feat_sel$selected
#' }
surv_feature_selection_multi <- function(object,
                                         methods = c("uni_cox", "lasso", "rf_imp", "vimp"),
                                         p_threshold = 0.05,
                                         top_ratio = 0.5,
                                         combine = c("union", "intersection", "freq"),
                                         freq_cutoff = 2,
                                         verbose = TRUE,
                                         use_boruta = FALSE) {
  
  combine <- match.arg(combine)
  
  # freq_cutoff only affects combine = "freq"; warn if the user set it
  # explicitly under a different combine mode, since it would otherwise
  # be silently ignored.
  if (!missing(freq_cutoff) && combine != "freq") {
    warning(sprintf(
      "`freq_cutoff` was supplied but has no effect when combine = '%s' (it only applies when combine = 'freq'). The value will be ignored.",
      combine
    ), call. = FALSE)
  }
  
  if (!requireNamespace("mlr3", quietly = TRUE))
    stop("Package 'mlr3' required.")
  
  # Extract task and data
  if (inherits(object, "PrognosiX")) {
    task <- surv_extract_task(object)
  } else if (inherits(object, "TaskSurv")) {
    task <- object
  } else {
    stop("object must be a PrognosiX or TaskSurv object.")
  }
  
  data <- as.data.frame(task$data())
  time_var <- task$target_names[1]
  status_var <- task$target_names[2]
  all_features <- task$feature_names
  
  if (verbose) cat("\n[Multi-Feature Selection] Methods:", paste(methods, collapse=", "), "\n")
  
  # Helper: keep top proportion of features based on importance vector
  keep_top <- function(imp_vec, ratio) {
    if (length(imp_vec) == 0) return(character(0))
    n_keep <- max(1, floor(length(imp_vec) * ratio))
    names(sort(imp_vec, decreasing = TRUE)[1:n_keep])
  }
  
  selection_list <- list()
  raw_results <- list()
  
  # ----------------------------------------------------------------------------
  # 1. Univariate Cox
  # ----------------------------------------------------------------------------
  if ("uni_cox" %in% methods) {
    if (verbose) cat("  - Running univariate Cox (p <", p_threshold, ")...\n")
    uni <- tryCatch({
      surv_filter_features_clinical(object, p_threshold = p_threshold)
    }, error = function(e) NULL)
    if (!is.null(uni)) {
      selection_list$uni_cox <- uni$task$feature_names
      raw_results$uni_cox <- uni$table
    } else {
      selection_list$uni_cox <- character(0)
    }
  }
  
  # ----------------------------------------------------------------------------
  # 2-4. glmnet family: LASSO, Ridge, Elastic Net
  # ----------------------------------------------------------------------------
  if ("lasso" %in% methods || "ridge" %in% methods || "enet" %in% methods) {
    if (verbose) cat("  - Running glmnet family (LASSO, Ridge, Elastic Net)...\n")
    
    # Helper to extract non-zero coefficients from glmnet model
    extract_glmnet_features <- function(lrn_obj, alpha_val) {
      selected <- character(0)
      if (!is.null(lrn_obj)) {
        coef_mat <- tryCatch(as.matrix(stats::coef(lrn_obj$model, s = "lambda.min")), error = function(e) NULL)
        if (!is.null(coef_mat)) {
          selected <- rownames(coef_mat)[abs(coef_mat[,1]) > 1e-6]
          selected <- setdiff(selected, "(Intercept)")
        }
      }
      return(selected)
    }
    
    if ("lasso" %in% methods) {
      lasso_lrn <- tryCatch(lrn("surv.cv_glmnet", alpha = 1, s = "lambda.min")$train(task), error = function(e) NULL)
      selection_list$lasso <- extract_glmnet_features(lasso_lrn, 1)
      raw_results$lasso <- lasso_lrn
    }
    
    if ("ridge" %in% methods) {
      ridge_lrn <- tryCatch(lrn("surv.cv_glmnet", alpha = 0, s = "lambda.min")$train(task), error = function(e) NULL)
      selection_list$ridge <- extract_glmnet_features(ridge_lrn, 0)
      raw_results$ridge <- ridge_lrn
    }
    
    if ("enet" %in% methods) {
      enet_lrn <- tryCatch(lrn("surv.cv_glmnet", alpha = 0.5, s = "lambda.min")$train(task), error = function(e) NULL)
      selection_list$enet <- extract_glmnet_features(enet_lrn, 0.5)
      raw_results$enet <- enet_lrn
    }
  }
  
  # ----------------------------------------------------------------------------
  # 5. Random Forest (ranger) permutation importance
  # ----------------------------------------------------------------------------
  if ("rf_imp" %in% methods) {
    if (verbose) cat("  - Running Random Forest (ranger) importance...\n")
    rf_lrn <- tryCatch({
      lrn("surv.ranger", importance = "permutation")$train(task)
    }, error = function(e) NULL)
    selected <- character(0)
    if (!is.null(rf_lrn) && !is.null(rf_lrn$importance())) {
      imp <- rf_lrn$importance()
      selected <- keep_top(imp, top_ratio)
    }
    selection_list$rf_imp <- selected
    raw_results$rf_imp <- if(!is.null(rf_lrn)) rf_lrn$importance() else NULL
  }
  
  # ----------------------------------------------------------------------------
  # 6. Random Survival Forest (rfsrc) VIMP importance
  # ----------------------------------------------------------------------------
  if ("rfsrc_imp" %in% methods) {
    if (verbose) cat("  - Running Random Survival Forest (rfsrc) importance...\n")
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
      warning("Package 'randomForestSRC' not installed. Skipping rfsrc_imp.")
    } else {
      rfsrc_lrn <- tryCatch({
        lrn("surv.rfsrc", importance = "permute")$train(task)
      }, error = function(e) NULL)
      selected <- character(0)
      if (!is.null(rfsrc_lrn) && !is.null(rfsrc_lrn$importance())) {
        imp <- rfsrc_lrn$importance()
        selected <- keep_top(imp, top_ratio)
      }
      selection_list$rfsrc_imp <- selected
      raw_results$rfsrc_imp <- if(!is.null(rfsrc_lrn)) rfsrc_lrn$importance() else NULL
    }
  }
  
  # ----------------------------------------------------------------------------
  # 7. XGBoost gain importance
  # ----------------------------------------------------------------------------
  if ("xgb_imp" %in% methods) {
    if (verbose) cat("  - Running XGBoost importance (gain)...\n")
    xgb_lrn <- tryCatch({
      lrn("surv.xgboost.cox", nrounds = 50)$train(task)
    }, error = function(e) NULL)
    selected <- character(0)
    if (!is.null(xgb_lrn) && !is.null(xgb_lrn$importance())) {
      imp_df <- xgb_lrn$importance()
      if (nrow(imp_df) > 0) {
        imp_vec <- setNames(imp_df$Gain, imp_df$Feature)
        selected <- keep_top(imp_vec, top_ratio)
      }
    }
    selection_list$xgb_imp <- selected
    raw_results$xgb_imp <- if(!is.null(xgb_lrn)) xgb_lrn$importance() else NULL
  }
  
  # ----------------------------------------------------------------------------
  # 8. VIMP (Variable Importance) from randomForestSRC -- Standalone algorithm
  # ----------------------------------------------------------------------------
  if ("vimp" %in% methods) {
    if (verbose) cat("  - Running VIMP variable importance...\n")
    if (!requireNamespace("randomForestSRC", quietly = TRUE)) {
      warning("Package 'randomForestSRC' not installed. Skipping vimp.")
    } else {
      vimp_fit <- tryCatch({
        randomForestSRC::vimp(
          formula = as.formula(paste("Surv(", time_var, ",", status_var, ") ~ .")),
          data = data[, c(all_features, time_var, status_var)],
          importance = "permute"
        )
      }, error = function(e) NULL)
      selected <- character(0)
      if (!is.null(vimp_fit)) {
        imp <- vimp_fit$importance
        selected <- keep_top(imp, top_ratio)
      }
      selection_list$vimp <- selected
      raw_results$vimp <- if(!is.null(vimp_fit)) vimp_fit else NULL
    }
  }
  
  # ----------------------------------------------------------------------------
  # 9. Boruta (optional, may fail on survival data)
  # ----------------------------------------------------------------------------
  if (use_boruta && "boruta" %in% methods) {
    if (verbose) cat("  - Running Boruta (optional, may be slow)...\n")
    if (!requireNamespace("Boruta", quietly = TRUE)) {
      warning("Package 'Boruta' not installed. Skipping boruta.")
    } else {
      selected <- character(0)
      tryCatch({
        cox_lrn <- lrn("surv.coxph")$train(task)
        risk <- cox_lrn$predict(task)$crank
        boruta_data <- data[, all_features, drop = FALSE]
        boruta_data$.risk <- risk
        boruta_fit <- Boruta::Boruta(.risk ~ ., data = boruta_data, doTrace = 0)
        selected <- names(boruta_fit$finalDecision[boruta_fit$finalDecision == "Confirmed"])
      }, error = function(e) {
        warning("Boruta failed: ", e$message)
      })
      selection_list$boruta <- selected
      raw_results$boruta <- if(exists("boruta_fit")) boruta_fit else NULL
    }
  }
  
  # ----------------------------------------------------------------------------
  # 10. Stepwise Cox (low-dimensional only)
  # ----------------------------------------------------------------------------
  if ("stepwise" %in% methods) {
    if (verbose) cat("  - Running stepwise Cox (direction = both, AIC)...\n")
    selected <- character(0)
    tryCatch({
      # Prevent high-dimensional explosion: take top 20 features from uni_cox filter first
      pre_feats <- if (length(selection_list$uni_cox) > 0) {
        head(selection_list$uni_cox, min(20, length(selection_list$uni_cox)))
      } else {
        head(all_features, min(20, length(all_features)))
      }
      if (length(pre_feats) > 0) {
        full_form <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~", 
                                      paste(pre_feats, collapse = " + ")))
        cox_full <- survival::coxph(full_form, data = data)
        step_mod <- MASS::stepAIC(cox_full, direction = "both", trace = 0)
        selected <- names(stats::coef(step_mod))
      }
    }, error = function(e) {
      warning("Stepwise failed: ", e$message)
    })
    selection_list$stepwise <- selected
    raw_results$stepwise <- if(exists("step_mod")) step_mod else NULL
  }
  
  # ----------------------------------------------------------------------------
  # 11. Stability Selection (via c060::stabpath)
  # ----------------------------------------------------------------------------
  if ("stab_sel" %in% methods) {
    if (verbose) cat("  - Running stability selection (c060::stabpath)...\n")
    selected <- character(0)
    if (!requireNamespace("c060", quietly = TRUE)) {
      warning("Package 'c060' not installed. Skipping stab_sel.")
    } else {
      tryCatch({
        x_mat <- as.matrix(data[, all_features, drop = FALSE])
        y_mat <- survival::Surv(data[[time_var]], data[[status_var]])
        stab_path <- c060::stabpath(y = y_mat, x = x_mat, steps = 50)
        selection_prob <- stab_path$stabpath
        if (!is.null(selection_prob) && ncol(selection_prob) > 0) {
          prob_avg <- colMeans(selection_prob, na.rm = TRUE)
          selected <- names(prob_avg)[prob_avg > 0.6]
        }
      }, error = function(e) {
        warning("Stability selection failed: ", e$message)
      })
      selection_list$stab_sel <- selected
      raw_results$stab_sel <- if(exists("stab_path")) stab_path else NULL
    }
  }
  
  
  # ----------------------------------------------------------------------------
  # 12. Combine results
  # ----------------------------------------------------------------------------
  method_names <- names(selection_list)
  if (length(method_names) == 0) {
    warning("No method succeeded. Returning all features.")
    selected_final <- all_features
  } else {
    method_table <- data.frame(Feature = all_features)
    for (m in method_names) {
      method_table[[m]] <- all_features %in% selection_list[[m]]
    }
    
    if (combine == "union") {
      selected_final <- unique(unlist(selection_list))
    } else if (combine == "intersection") {
      selected_final <- Reduce(intersect, selection_list)
    } else if (combine == "freq") {
      counts <- rowSums(method_table[, -1])
      selected_final <- all_features[counts >= freq_cutoff]
    }
    
    if (verbose) {
      cat(sprintf("\n[Combined] %s selection: %d features out of %d\n",
                  toupper(combine), length(selected_final), length(all_features)))
    }
  }
  
  return(list(
    selected = selected_final,
    method_table = method_table,
    method_results = raw_results
  ))
}

#' Decision Curve Analysis for One or More Survival Models
#'
#' Computes and plots decision curves at a specified time point for one or more
#' survival models using standard Kaplan-Meier corrections via the 'dcurves' package.
#' Supports highly flexible aesthetic configurations for colors and linetypes.
#'
#' @param learners A **named list** of trained `mlr3` survival learners.
#'   Each element must be a learner that supports `"distr"` predictions.
#'   Example: `list("Ranger" = learner_ranger)`.
#' @param object A `TaskSurv` or `PrognosiX` object containing the validation data.
#' @param eval_time Numeric. The time point at which to evaluate event probabilities.
#' @param thresholds Numeric vector. Risk thresholds (probabilities of event) at
#'   which net benefit is calculated. Defaults to `seq(0.01, 0.99, length.out = 50)`.
#' @param colors Character vector or single string. Can be:
#'   - A built-in palette keyword: `"default"`, `"clinical"`, `"vibrant"`, or `"jama"`.
#'   - A fully/partially named vector, e.g. `c(Ranger = "#2c7fb8")`. Missing reference
#'     strategies (`TreatAll`, `TreatNone`) will be filled automatically with high-contrast distinct colors.
#' @param linetypes Named character vector of line types for the strategies.
#'   Can be partially named (e.g., `c(Ranger = "solid")`); missing reference styles
#'   will default to distinct `"dashed"` and `"dotted"` profiles automatically.
#' @param include_reference Logical. Should the "Treat All" and "Treat None" curves be added? Default is `TRUE`.
#' @param ylim Numeric vector of length 2. Y-axis limits for net benefit. Defaults to `c(-0.05, NA)`.
#' @param title Character. Custom plot title.
#' @param subtitle Character. Custom plot subtitle.
#' @param print_stats Logical. Should summary statistics be printed to the console? Default `TRUE`.
#' @param clin_range Numeric vector of length 2. Default is `c(0.05, 0.5)`.
#'
#' @return A list with three components: `plot`, `table`, and `summary`.
#'
#' @importFrom survival Surv
#' @importFrom dcurves dca
#' @importFrom tibble as_tibble
#' @importFrom ggplot2 ggplot aes geom_line labs scale_color_manual scale_linetype_manual coord_cartesian theme_minimal
#' @importFrom tidyr pivot_wider
#' @export
#' @examples
#' \dontrun{
#' veteran <- survival::veteran
#' task <- surv_create_surv_task(veteran, time_col = "time", event_col = "status")
#' lrn <- surv_get_learner("surv.coxph", task)
#' lrn$train(task)
#' lrn$predict_type <- "distr"
#'
#' dca_res <- plot_dca_survival(
#'   learners = list("Cox" = lrn), object = task,
#'   eval_time = 90, clin_range = c(0.05, 0.5), print_stats = FALSE
#' )
#' }
plot_dca_survival <- function(learners,
                                object,
                                eval_time,
                                thresholds = seq(0.01, 0.99, length.out = 50),
                                colors = NULL,
                                linetypes = NULL,
                                include_reference = TRUE,
                                ylim = c(-0.05, NA),
                                title = NULL,
                                subtitle = NULL,
                                print_stats = TRUE,
                                clin_range = c(0.05, 0.5)) {
    
    # ---- 1. Input validation and Task extraction  ------------------------------
    if (!is.list(learners) || is.null(names(learners))) {
      stop("'learners' must be a named list (e.g., list(Model1 = lrn1, Model2 = lrn2))")
    }
    if (length(learners) == 0) stop("At least one learner is required.")
    if (!is.numeric(eval_time) || length(eval_time) != 1) stop("'eval_time' must be a single number.")
    if (!is.numeric(clin_range) || length(clin_range) != 2 || clin_range[1] >= clin_range[2]) {
      stop("'clin_range' must be a numeric vector of length 2 with min < max.")
    }
    
    if (inherits(object, "TaskSurv")) {
      task <- object
    } else if (inherits(object, "PrognosiX")) {
      task <- surv_extract_task(object)
    } else {
      stop("object must be a TaskSurv or PrognosiX object.")
    }
    
    data_df <- as.data.frame(task$data())
    time_var <- task$target_names[1]
    status_var <- task$target_names[2]
    
    # ---- 2. Safely extract predicted event probability from mlr3 learners [1 - S(t)] ---------
    dca_data <- data_df[, c(time_var, status_var), drop = FALSE]
    
    for (nm in names(learners)) {
      lrn <- learners[[nm]]
      if (!"distr" %in% lrn$predict_types) {
        stop(sprintf("Learner '%s' does not support 'distr' predictions.", lrn$id))
      }
      lrn$predict_type <- "distr"
      pred <- lrn$predict(task)
      
      surv_prob <- as.numeric(pred$distr$survival(eval_time))
      dca_data[[nm]] <- 1 - surv_prob
    }
    
    # ---- 3. Standard survival DCA calculation using dcurves -----------------------
    formula_str <- sprintf("survival::Surv(%s, %s) ~ %s", 
                           time_var, status_var, 
                           paste(names(learners), collapse = " + "))
    
    dca_obj <- dcurves::dca(
      formula = as.formula(formula_str),
      data = dca_data,
      time = eval_time,
      thresholds = thresholds
    )
    
    nb_table <- as.data.frame(tibble::as_tibble(dca_obj))
    
    colnames(nb_table)[colnames(nb_table) == "variable"]    <- "Strategy"
    colnames(nb_table)[colnames(nb_table) == "threshold"]   <- "Threshold"
    colnames(nb_table)[colnames(nb_table) == "net_benefit"] <- "NetBenefit"
    
    if (!include_reference) {
      nb_table <- nb_table[!nb_table$Strategy %in% c("all", "none"), ]
    } else {
      nb_table$Strategy <- ifelse(nb_table$Strategy == "all", "TreatAll",
                                  ifelse(nb_table$Strategy == "none", "TreatNone", nb_table$Strategy))
    }
    
    all_strategies <- unique(nb_table$Strategy)
    model_names <- names(learners)
    
    # ---- 4. Automated aesthetic settings (high-contrast colors and line types) ----------------------------
    default_ref_colors <- c("TreatAll" = "#E69F00", "TreatNone" = "#000000") 
    default_ref_lts    <- c("TreatAll" = "dashed", "TreatNone" = "dotted")
    
    if (is.null(colors)) {
      hues <- seq(15, 375, length.out = length(model_names) + 1)
      mod_cols <- grDevices::hcl(hues[1:length(model_names)], l = 55, c = 90)
      names(mod_cols) <- model_names
      colors <- c(mod_cols, default_ref_colors)
    } else if (is.character(colors) && length(colors) == 1) {
      pal_choice <- tolower(colors)
      if (pal_choice == "clinical") {
        mod_cols <- c("#00A087FF", "#3C5488FF", "#4DBBD5FF")[1:length(model_names)]
        names(mod_cols) <- model_names
        colors <- c(mod_cols, "TreatAll" = "#E64B35FF", "TreatNone" = "#111111")
      } else if (pal_choice == "vibrant") {
        mod_cols <- c("#0073C2FF", "#CD534CFF", "#7AA6C2FF")[1:length(model_names)]
        names(mod_cols) <- model_names
        colors <- c(mod_cols, "TreatAll" = "#EFC000FF", "TreatNone" = "#000000")
      } else if (pal_choice == "jama") {
        mod_cols <- c("#374E55FF", "#DF8F44FF", "#00A1D5FF")[1:length(model_names)]
        names(mod_cols) <- model_names
        colors <- c(mod_cols, "TreatAll" = "#B24745FF", "TreatNone" = "#79AF97FF")
      } else {
        hues <- seq(15, 375, length.out = length(model_names) + 1)
        mod_cols <- grDevices::hcl(hues[1:length(model_names)], l = 55, c = 90)
        names(mod_cols) <- model_names
        colors <- c(mod_cols, default_ref_colors)
      }
    } else {
      missing_refs <- setdiff(c("TreatAll", "TreatNone"), names(colors))
      if (length(missing_refs) > 0) {
        colors <- c(colors, default_ref_colors[missing_refs])
      }
    }
    
    if (is.null(linetypes)) {
      lt <- rep("solid", length(all_strategies))
      names(lt) <- all_strategies
      if ("TreatAll" %in% all_strategies) lt["TreatAll"] <- "dashed"
      if ("TreatNone" %in% all_strategies) lt["TreatNone"] <- "dotted"
      linetypes <- lt
    } else {
      missing_lts <- setdiff(c("TreatAll", "TreatNone"), names(linetypes))
      if (length(missing_lts) > 0) {
        linetypes <- c(linetypes, default_ref_lts[missing_lts])
      }
      unmapped_models <- setdiff(model_names, names(linetypes))
      if (length(unmapped_models) > 0) {
        extra_lts <- rep("solid", length(unmapped_models))
        names(extra_lts) <- unmapped_models
        linetypes <- c(linetypes, extra_lts)
      }
    }
    
    # ---- 5. Build ggplot (using coord_cartesian to avoid data clipping warnings) -----------------
    if (is.null(title)) title <- "Decision Curve Analysis"
    if (is.null(subtitle)) subtitle <- paste("Time =", eval_time, "| Validation set")
    
    p <- ggplot2::ggplot(nb_table, ggplot2::aes(x = Threshold, y = NetBenefit, 
                                                color = Strategy, linetype = Strategy)) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::scale_color_manual(values = colors) +
      ggplot2::scale_linetype_manual(values = linetypes) +
      ggplot2::labs(title = title, subtitle = subtitle,
                    x = "Risk Threshold (Probability of Event)",
                    y = "Net Benefit",
                    color = "Strategy", linetype = "Strategy") +
      ggplot2::coord_cartesian(ylim = ylim) # Core fix: avoid warnings and preserve full curve continuity
    
    if (requireNamespace("ggprism", quietly = TRUE)) {
      p <- p + ggprism::theme_prism()
    } else {
      p <- p + ggplot2::theme_minimal()
    }
    
    # ---- 6. Metric evaluation within clinically relevant range (smartly adapts to user-defined risk threshold) ---------------------
    max_avail_thr <- max(nb_table$Threshold, na.rm = TRUE)
    if (clin_range[2] > max_avail_thr) {
      clin_range[2] <- max_avail_thr # Auto-align evaluation upper bound to user input max threshold, prevent out-of-bounds
    }
    
    summary_stats <- data.frame()
    
    for (nm in model_names) {
      model_nb <- nb_table[nb_table$Strategy == nm, ]
      idx <- which(model_nb$Threshold >= clin_range[1] & model_nb$Threshold <= clin_range[2])
      if (length(idx) == 0) idx <- seq_len(nrow(model_nb))
      
      nb_clin <- model_nb$NetBenefit[idx]
      thr_clin <- model_nb$Threshold[idx]
      
      max_nb   <- max(nb_clin, na.rm = TRUE)
      best_thr <- thr_clin[which.max(nb_clin)]
      avg_nb   <- mean(nb_clin, na.rm = TRUE)
      
      treat_all_clin <- nb_table[nb_table$Strategy == "TreatAll" & nb_table$Threshold %in% thr_clin, "NetBenefit"]
      avg_gain <- if(length(treat_all_clin) == length(nb_clin)) mean(nb_clin - treat_all_clin, na.rm = TRUE) else NA
      
      summary_stats <- rbind(summary_stats, data.frame(
        Model = nm,
        ClinRange_Min = clin_range[1],
        ClinRange_Max = clin_range[2],
        Max_NetBenefit = round(max_nb, 4),
        Threshold_at_Max = round(best_thr, 4),
        Avg_NetBenefit = round(avg_nb, 4),
        Avg_NetBenefit_Gain = round(avg_gain, 4)  # Removed invalid and non-standard AUC_NetBenefit metric
      ))
    }
    
    if (print_stats) {
      cat("\n========== DCA Summary (Clinical Range: [", clin_range[1], ", ", clin_range[2], "]) ==========\n", sep = "")
      cat("  Calculated using standard Kaplan-Meier survival adjustments via 'dcurves'.\n")
      print(summary_stats, row.names = FALSE)
      cat("========================================================================\n")
    }
    
    # ---- 7. Matrix form output (Core fix: remove unique attribute columns, ensure wide table perfect alignment) -------
    nb_table_core <- nb_table[, c("Threshold", "Strategy", "NetBenefit")]
    dca_wide <- tidyr::pivot_wider(nb_table_core, names_from = Strategy, values_from = NetBenefit)
    
    invisible(list(
      plot = p,
      table = as.data.frame(dca_wide),
      summary = summary_stats
    ))
  }
                           
#' List all available feature selection methods in surv_feature_selection_multi
#'
#' @param verbose If TRUE, prints the table to console. If FALSE, returns the data frame.
#' @return A data frame with columns: Method, Description, RequiredPackages, Recommendation.
#' @export
#'
#' @examples
#' \dontrun{
#' list_surv_feature_methods()
#' }
list_surv_feature_methods <- function(verbose = TRUE) {
  methods_df <- data.frame(
    Method = c(
      "uni_cox", "lasso", "ridge", "enet", "rf_imp", "rfsrc_imp",
      "xgb_imp", "vimp", "boruta", "stepwise", "stab_sel"
    ),
    Description = c(
      "Univariate Cox regression (p < threshold)",
      "LASSO penalized Cox (lambda.min)",
      "Ridge penalized Cox",
      "Elastic net (alpha = 0.5, cross-validated)",
      "Random forest (ranger) permutation importance, keep top ratio",
      "Random survival forest (randomForestSRC) importance, keep top ratio",
      "XGBoost gain importance, keep top ratio",
      "VIMP variable importance from randomForestSRC (robust, recommended)",
      "Boruta wrapper algorithm (default OFF, may fail on survival data)",
      "Stepwise Cox regression (both directions, AIC) - low-dimensional only",
      "Stability selection with Lasso via c060::stabpath"
    ),
    RequiredPackages = c(
      "survival (built-in)",
      "glmnet (via mlr3learners)",
      "glmnet",
      "glmnet",
      "ranger (mlr3learners)",
      "randomForestSRC",
      "xgboost (mlr3extralearners)",
      "randomForestSRC",
      "Boruta, randomForest (optional)",
      "MASS",
      "c060"
    ),
    Recommendation = c(
      "* * * * * (must-have)",
      "* * * * * (top choice)",
      "* * * * * (high collinearity)",
      "* * * * * (often best glmnet)",
      "* * * * * (nonlinear effects)",
      "* * * * * (survival-specialized)",
      "* * * * * (handles missing data)",
      "* * * * * (stable, official RF method)",
      "* * * * * (use with caution, OFF by default)",
      "* * * * * (only for low dimension, p < 30)",
      "* * * * * (robust for high-dim)"
    ),
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat("\n============================================================\n")
    cat("Survival Feature Selection Methods (updated)\n")
    cat("============================================================\n\n")
    print(methods_df, row.names = FALSE)
    cat("\n[Note] Use combine = 'union', 'intersection', or 'freq'.\n")
    cat("Recommended production set: c('uni_cox', 'lasso', 'vimp', 'rf_imp')\n")
  }
  
  invisible(methods_df)
}

# Short alias
print_surv_feature_methods <- function() {
  list_surv_feature_methods(verbose = TRUE)
}

# =========================================================================
# Internal Helpers
# =========================================================================
#' Convert SurvSHAP Result to Long Format
#'
#' Internal helper function to convert the output from \code{survex} SHAP
#' calculations into a tidy long-format data frame.
#'
#' @param shap_obj The SHAP object returned by \code{survex}.
#' @param obs_idx The indices of observations explained.
#' @param features The original feature data frame.
#'
#' @return A data frame with columns: \code{time}, \code{feature}, \code{shap_value},
#'   and \code{observation}.
#'
#' @keywords internal
#' @noRd
.survshap_to_long <- function(shap_obj, obs_idx, features) {
  result <- shap_obj$result
  times  <- shap_obj$eval_times
  var_vals <- features[obs_idx, , drop = FALSE]
  obs_names <- rownames(var_vals)
  if (is.null(obs_names)) obs_names <- as.character(obs_idx)
  
  if (is.data.frame(result)) {
    df_long <- .df_to_long(result, times, obs_names[1L])
    return(df_long)
  }
  
  if (is.list(result)) {
    n <- length(result)
    long_list <- lapply(seq_len(n), function(i) {
      .df_to_long(as.data.frame(result[[i]]), times, obs_names[i])
    })
    return(do.call(rbind, long_list))
  }
  
  stop("Cannot parse SurvSHAP result. Unexpected structure: ", class(result)[1])
}

.df_to_long <- function(df, times, obs_label) {
  df <- df[, !names(df) %in% c("times", "time", "_times_", "id", "B"),
           drop = FALSE]
  
  n_rows <- nrow(df)
  t_vec <- if (!is.null(times) && length(times) == n_rows) {
    times
  } else {
    rn <- rownames(df)
    t_parsed <- suppressWarnings(as.numeric(sub("^t=", "", rn)))
    if (all(!is.na(t_parsed))) t_parsed else seq_len(n_rows)
  }
  
  df_long <- tidyr::pivot_longer(
    cbind(data.frame(.time = t_vec, stringsAsFactors = FALSE), df),
    cols = -".time", names_to = "feature", values_to = "shap_value"
  )
  names(df_long)[names(df_long) == ".time"] <- "time"
  df_long$observation <- obs_label
  as.data.frame(df_long)
}

.generate_shap_plots <- function(shap_long, n_top, bar_col, model_id, type) {
  feat_imp <- shap_long %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(mean_abs = mean(abs(shap_value), na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice(seq_len(min(n_top, dplyr::n())))
  
  feat_imp$feature <- factor(feat_imp$feature, levels = rev(feat_imp$feature))
  
  n_obs_used <- length(unique(shap_long$observation))
  title_text <- sprintf("SurvSHAP Importance (%s mode, n=%d): %s",
                        toupper(type), n_obs_used, model_id)
  
  p_bar <- ggplot2::ggplot(feat_imp, ggplot2::aes(x = feature, y = mean_abs)) +
    ggplot2::geom_bar(stat = "identity", fill = bar_col, width = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title_text, x = "Feature", y = "Mean |SHAP|") +
    ggprism::theme_prism()
  
  if (!all(is.na(shap_long$time))) {
    top_feats  <- feat_imp$feature
    line_data  <- shap_long %>%
      dplyr::filter(feature %in% top_feats) %>%
      dplyr::group_by(feature, time) %>%
      dplyr::summarise(mean_shap = mean(shap_value, na.rm = TRUE),
                       .groups = "drop")
    
    p_line <- ggplot2::ggplot(line_data,
                              ggplot2::aes(x = time, y = mean_shap, color = feature)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::labs(title = "SHAP Dynamics over Time",
                    x = "Time", y = "Average SHAP Value") +
      ggprism::theme_prism()
  } else {
    p_line <- NULL
  }
  
  list(bar_plot = p_bar, line_plot = p_line)
}

.prognosis_optional_packages <- c(
  "mlr3", "mlr3proba", "mlr3tuning", "mlr3learners",
  "mlr3extralearners", "survival", "tidyverse", "paradox", "data.table"
)

.check_prognosis_packages <- function() {
  missing <- .prognosis_optional_packages[
    !vapply(.prognosis_optional_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing) > 0) {
    stop(
      "Missing packages required for PrognosiX framework: ",
      paste(missing, collapse = ", "),
      ". Install them before using prognosis-related functions."
    )
  }
  invisible(TRUE)
}


# Unified interface to calculate calibration metrics
.compute_cal_metrics <- function(cal_df) {
  metrics <- list()
  lm_fit            <- lm(observed ~ predicted, data = cal_df)
  metrics$slope     <- unname(stats::coef(lm_fit)[2L])
  metrics$intercept <- unname(stats::coef(lm_fit)[1L])
  metrics$r_squared <- summary(lm_fit)$r.squared
  
  errors            <- abs(cal_df$observed - cal_df$predicted)
  metrics$mae       <- mean(errors, na.rm = TRUE)
  metrics$ici       <- .compute_ici(cal_df)
  metrics$e50       <- unname(quantile(errors, 0.5, na.rm = TRUE))
  metrics$e90       <- unname(quantile(errors, 0.9, na.rm = TRUE))
  return(metrics)
}

# Extractor for predicted survival probabilities from mlr3 distribution objects
.extract_surv_prob <- function(distr, time_point, task) {
  if (inherits(distr, "Matdist") || is.environment(distr)) {
    if (!is.null(distr$survival)) {
      sp <- distr$survival(time_point)
      return(if (is.matrix(sp)) as.numeric(sp) else sp)
    } else if (!is.null(distr$cdf)) {
      sp <- 1 - distr$cdf(time_point)
      return(if (is.matrix(sp)) as.numeric(sp) else sp)
    }
    warning("No survival or cdf method in distr object.")
    return(NULL)
  } else if (is.matrix(distr) || is.array(distr)) {
    times <- attr(distr, "times")
    if (is.null(times)) {
      data_tmp <- as.data.frame(task$data())
      times <- sort(unique(data_tmp[[task$target_names[1L]]][data_tmp[[task$target_names[2L]]] == 1]))
    }
    t_idx <- which.min(abs(times - time_point))
    sp <- if (length(dim(distr)) == 2L) distr[, t_idx] else distr[, t_idx, 1L]
    return(as.numeric(sp))
  }
  warning("Unknown distr type")
  NULL
}

# Integrated Calibration Index (ICI) computation engine
.compute_ici <- function(cal_df) {
  if (nrow(cal_df) < 2L) return(NA_real_)
  tryCatch({
    if (nrow(cal_df) >= 4L) {
      lo    <- loess(observed ~ predicted, data = cal_df,
                     span = 1.0, degree = 1L, surface = "direct")
      x_seq <- seq(min(cal_df$predicted), max(cal_df$predicted), length.out = 200L)
      y_hat <- predict(lo, newdata = data.frame(predicted = x_seq))
      valid <- !is.na(y_hat)
      if (sum(valid) < 2L) return(NA_real_)
      x_v <- x_seq[valid]; d_v <- abs(y_hat[valid] - x_v)
    } else {
      x_v <- cal_df$predicted; d_v <- abs(cal_df$observed - x_v)
    }
    rng <- diff(range(x_v))
    if (rng < 1e-6) return(NA_real_)
    sum(diff(x_v) * (d_v[-length(d_v)] + d_v[-1L]) / 2L) / rng
  }, error = function(e) NA_real_)
}

