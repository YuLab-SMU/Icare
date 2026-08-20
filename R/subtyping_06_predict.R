# =============================================================================
# Module 3 - Validation Set Subtype Prediction
# =============================================================================
# Functions:
#   Sub_extract_norm_params  Extract per-variable normalisation params from
#                            a trained Subtyping object (supports all methods)
#   Sub_apply_norm_params    Apply training normalisation params to new data
#   Sub_predict_subtypes     Predict subtypes with a trained model (kmeans/lpa/nmf)
#   Sub_create_val_object    One-step: create + normalise + predict for a val set
# =============================================================================


# -----------------------------------------------------------------------------
#' Extract Normalisation Parameters from a Trained Subtyping Object
#'
#' Reads \code{object@clean.data} and the normalisation method stored in
#' \code{object@scale.data}'s \code{"normalization_info"} attribute (set by
#' \code{normalize_data} in module 1) to build a per-variable parameter list
#' that can later be applied to a validation set via \code{Sub_apply_norm_params}.
#'
#' The returned list has one entry per numeric variable:
#' \preformatted{
#'   list(
#'     var1 = list(method = "min_max", min = ..., max = ...),
#'     var2 = list(method = "z_score", mean = ..., sd = ...),
#'     var3 = list(method = "log",     min_raw = ...),
#'     ...
#'   )
#' }
#'
#' @param object  A \code{Subtyping} object whose \code{scale.data} slot has
#'   been filled by \code{Sub_normalize_process}.
#' @param normalize_method  Fallback method string used when the
#'   \code{"normalization_info"} attribute is absent (default \code{"min_max"}).
#' @param verbose  Print a summary (default \code{TRUE}).
#' @return Named list of per-variable normalisation parameters.
#' @examples
#' \dontrun{
#'   # Assuming 'obj' is a trained Subtyping object with scale.data
#'   norm_params <- Sub_extract_norm_params(obj, verbose = TRUE)
#'   print(norm_params[[1]])
#' }
#' @export
Sub_extract_norm_params <- function(object,
                                    normalize_method = "min_max",
                                    verbose          = TRUE) {

  if (!inherits(object, "Subtyping"))
    stop("'object' must be a 'Subtyping' object.")

  raw_data <- object@clean.data
  if (is.null(raw_data) || nrow(raw_data) == 0)
    stop("'object@clean.data' is empty.")

  # Retrieve per-column method map (produced by normalize_data in module 1)
  norm_info <- attr(object@scale.data, "normalization_info")

  num_cols <- names(which(vapply(raw_data, is.numeric, logical(1))))

  params <- lapply(num_cols, function(col) {
    x      <- raw_data[[col]]
    method <- if (!is.null(norm_info) && !is.null(norm_info[[col]]))
                norm_info[[col]]
              else
                normalize_method
    if (method %in% c("box_cox", "yeo_johnson")) {
      stop("Normalization method '", method, "' requires parameter estimation on new data ",
           "and cannot be reliably applied to validation sets. Please re-run normalization ",
           "with a simpler method (e.g., 'min_max', 'z_score', 'log').")
    }
    p <- list(method = method)

    switch(method,
      min_max = ,
      min_max_scale = {
        p$min <- min(x, na.rm = TRUE)
        p$max <- max(x, na.rm = TRUE)
      },
      z_score = ,
      z_score_standardize = {
        p$mean <- mean(x, na.rm = TRUE)
        p$sd   <- sd(x,   na.rm = TRUE)
      },
      center = {
        p$mean <- mean(x, na.rm = TRUE)
      },
      scale = {
        p$sd <- sd(x, na.rm = TRUE)
      },
      max_abs = {
        p$max_abs <- max(abs(x), na.rm = TRUE)
      },
      log = {
        p$min_raw <- min(x, na.rm = TRUE)
      },
      box_cox = ,
      yeo_johnson = {
        # These transforms depend on a fitted object; we store the raw stats
        # as a fallback and warn.
        p$mean    <- mean(x, na.rm = TRUE)
        p$sd      <- sd(x,   na.rm = TRUE)
        p$min_raw <- min(x,   na.rm = TRUE)
        warning(sprintf(
          "Column '%s': method '%s' is not invertible without the fitted caret object. ",
          col, method),
          "Validation data will be re-transformed independently using the same method.",
          call. = FALSE)
      },
      {
        warning(sprintf("Unknown normalisation method '%s' for column '%s'. Skipping.", method, col),
                call. = FALSE)
      }
    )
    p
  })
  names(params) <- num_cols

  if (verbose) {
    method_tbl <- table(vapply(params, `[[`, character(1), "method"))
    cat("Sub_extract_norm_params: extracted parameters for", length(params), "variables.\n")
    cat("  Methods used:", paste(names(method_tbl), method_tbl, sep = "x", collapse = ", "), "\n")
  }

  params
}


# -----------------------------------------------------------------------------
#' Apply Training Normalisation Parameters to a Validation Set
#'
#' Accepts the extended parameter list produced by \code{Sub_extract_norm_params}
#' and applies the appropriate transformation to each variable in
#' \code{object@clean.data}, writing the result to \code{object@scale.data}.
#'
#' Supported methods (matched from the \code{method} field of each parameter
#' entry): \code{min_max}, \code{z_score}, \code{center}, \code{scale},
#' \code{max_abs}, \code{log}, \code{box_cox}, \code{yeo_johnson}.
#'
#' For backwards compatibility the function also accepts the old format
#' \code{list(var = list(min=, max=))} (without a \code{method} field), which
#' is treated as \code{min_max}.
#'
#' @param object      A \code{Subtyping} object whose \code{clean.data} slot
#'                    contains the raw (un-scaled) validation features.
#' @param norm_params Named list produced by \code{Sub_extract_norm_params}, or
#'   the legacy \code{list(var = list(min=, max=))} format.
#' @param verbose     Print progress messages (default \code{TRUE}).
#' @return The updated \code{Subtyping} object with \code{scale.data} filled.
#' @examples
#' \dontrun{
#'   # Assuming 'val_obj' is a validation Subtyping object and 'norm_params' from training
#'   val_obj <- Sub_apply_norm_params(val_obj, norm_params = norm_params)
#' }
#' @export
Sub_apply_norm_params <- function(object,
                                  norm_params,
                                  verbose = TRUE) {

  if (!inherits(object, "Subtyping"))
    stop("'object' must be a 'Subtyping' object.")
  if (!is.list(norm_params) || is.null(names(norm_params)))
    stop("'norm_params' must be a named list returned by Sub_extract_norm_params.")

  scaled     <- object@clean.data
  n_applied  <- 0L
  n_skipped  <- 0L

  for (col in names(norm_params)) {
    if (!col %in% colnames(scaled) || !is.numeric(scaled[[col]])) {
      n_skipped <- n_skipped + 1L
      next
    }

    p      <- norm_params[[col]]
    x      <- scaled[[col]]

    # Backwards-compatibility: old format had no 'method' field
    method <- if (!is.null(p$method)) p$method else "min_max"

    scaled[[col]] <- switch(method,

      min_max = ,
      min_max_scale = {
        rng <- p$max - p$min
        if (is.na(rng) || rng == 0) rep(0, length(x))
        else pmin(pmax((x - p$min) / rng, 0), 1)
      },

      z_score = ,
      z_score_standardize = {
        if (is.na(p$sd) || p$sd == 0) rep(0, length(x))
        else (x - p$mean) / p$sd
      },

      center = {
        x - p$mean
      },

      scale = {
        if (is.na(p$sd) || p$sd == 0) rep(0, length(x))
        else x / p$sd
      },

      max_abs = {
        if (is.na(p$max_abs) || p$max_abs == 0) rep(0, length(x))
        else x / p$max_abs
      },

      log = {
        # Shift using training minimum to guarantee positivity
        shift <- if (p$min_raw <= 0) abs(p$min_raw) + 1 else 0
        log(x + shift)
      },

      box_cox = ,
      yeo_johnson = {
        # Cannot replay a fitted caret object; re-transform independently
        if (method == "box_cox") {
          boxcox_transform(x)
        } else {
          yeojohnson_transform(x)
        }
      },

      {
        # Unknown method - leave unchanged
        warning(sprintf("Unknown method '%s' for column '%s'; left unchanged.", method, col),
                call. = FALSE)
        x
      }
    )

    n_applied <- n_applied + 1L
  }

  if (verbose) {
    cat(sprintf(
      "Normalisation applied: %d variables transformed, %d skipped (absent or non-numeric).\n",
      n_applied, n_skipped
    ))
    cat("scale.data dimensions:", dim(scaled), "\n")
  }

  object@scale.data <- scaled
  return(object)
}


# -----------------------------------------------------------------------------
#' Predict Subtypes for a Validation \code{Subtyping} Object
#'
#' Dispatches to the appropriate prediction method depending on \code{method}:
#' \describe{
#'   \item{\code{"kmeans"}}{Assigns each validation sample to the nearest
#'     training centroid (Euclidean distance in the scaled feature space).}
#'   \item{\code{"lpa"}}{Uses \code{mclust::predict.Mclust} to classify
#'     samples under the fitted Gaussian-mixture model.}
#'   \item{\code{"nmf"}}{Projects samples via NNLS onto the training W matrix
#'     (features x components); assignment is the component with the highest
#'     coefficient.}
#' }
#'
#' In all cases the predicted label is written to
#' \code{object@info.data$cluster_<method>}.
#'
#' @param object       Validation \code{Subtyping} object (must have
#'                     \code{scale.data} filled, e.g. via
#'                     \code{Sub_apply_norm_params}).
#' @param train_object Trained \code{Subtyping} object that carries the fitted
#'                     model inside \code{cluster.results}.
#' @param method       One of \code{"kmeans"}, \code{"lpa"}, \code{"nmf"}.
#' @param prefix       Subtype label prefix for NMF (default \code{"S"}).
#' @param verbose      Print progress (default \code{TRUE}).
#' @return The updated validation \code{Subtyping} object.
#' @importFrom stats dist
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#'   # Assuming 'val_obj' is a validation Subtyping object and 'train_obj' is trained
#'   val_obj <- Sub_predict_subtypes(val_obj, train_obj, method = "kmeans")
#'   val_obj <- Sub_predict_subtypes(val_obj, train_obj, method = "lpa")
#'   val_obj <- Sub_predict_subtypes(val_obj, train_obj, method = "nmf")
#' set.seed(1)
#' require('NMF')
#' demo_df <- data.frame(
#'   group = rep(c(0, 1), each = 30),
#'   feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'   feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'   feat3 = rnorm(60, 10, 2),
#'   feat4 = c(rnorm(30, 1, 0.3), rnorm(30, 3, 0.3))
#'  )
#' rownames(demo_df) <- paste0("S", 1:60)  
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#' split <- SplitSubtypingObject(sub_obj, p = 0.7, stratify_by = "group")
#' sub_train <- Sub_normalize_process(split$train, normalize_method = "min_max")
#' norm_params <- Sub_extract_norm_params(sub_train, verbose = FALSE)
#' sub_test <- Sub_apply_norm_params(split$test, norm_params = norm_params, verbose = FALSE)
#'
#' sub_train <- Sub_nmf_estimate(sub_train, rank_range = 2:3, nrun = 2, seed = 1)
#' sub_train <- Sub_nmf_best_rank(sub_train, nrun = 2)
#' sub_train <- Sub_nmf_assign_subtypes(sub_train)
#' sub_train <- Sub_nmf_train_model(
#'   sub_train, best_k = sub_train@Optimal.cluster, nrun = 2,
#'   model_name = "demo_nmf_model",save_dir   ="./"
#' )
#'
#' # Predict subtypes on the held-out validation split.
#' sub_test <- Sub_predict_subtypes(sub_test, train_object = sub_train,
#'                                  method = "nmf", verbose = FALSE)
#' table(sub_test@info.data$cluster_nmf)
#' }
Sub_predict_subtypes <- function(object,
                                 train_object,
                                 method  = c("nmf", "kmeans", "lpa"),
                                 prefix  = "S",
                                 verbose = TRUE) {

  method <- match.arg(method)

  if (!inherits(object,       "Subtyping")) stop("'object' must be a 'Subtyping' object.")
  if (!inherits(train_object, "Subtyping")) stop("'train_object' must be a 'Subtyping' object.")

  val_scaled <- object@scale.data
  if (nrow(val_scaled) == 0)
    stop("'object@scale.data' is empty. Run Sub_apply_norm_params first.")

  # Ensure info.data rows match clean.data rows
  if (nrow(object@info.data) == 0)
    object@info.data <- data.frame(row.names = rownames(object@clean.data))
  object@info.data <- object@info.data[rownames(object@clean.data), , drop = FALSE]

  col_name <- paste0("cluster_", method)

  # ---------------------------- K-means ------------------------------------
  if (method == "kmeans") {

    km_res <- train_object@cluster.results[["kmeans.result"]]
    if (is.null(km_res))
      stop("No kmeans.result found in train_object. Run Sub_kmeans_with_optimal_k first.")

    centers <- km_res$model$centers   # K x features matrix
    common  <- intersect(colnames(val_scaled), colnames(centers))
    if (length(common) == 0) stop("No overlapping features between val data and kmeans centroids.")
    if (length(common) < ncol(centers) && verbose)
      warning(ncol(centers) - length(common), " centroid features absent in val data - using common features only.")

    val_mat <- as.matrix(val_scaled[, common, drop = FALSE])
    ctr_mat <- centers[, common, drop = FALSE]

    # Assign each sample to the nearest centroid
    dist_mat <- apply(val_mat, 1, function(x) {
      sqrt(rowSums(sweep(ctr_mat, 2, x)^2))
    })                                       # K x n_val
    nearest  <- apply(dist_mat, 2, which.min)
    labels   <- factor(nearest)

    if (verbose) {
      cat("K-means prediction (nearest centroid):\n")
      print(table(labels))
    }
    object@info.data[[col_name]] <- labels[match(rownames(object@info.data),
                                                  names(nearest))]
  }

  # ---------------------------- LPA (mclust) -------------------------------
  else if (method == "lpa") {

    lpa_res <- train_object@cluster.results[["lpa.result"]]
    if (is.null(lpa_res))
      stop("No lpa.result found in train_object. Run Sub_lpa_with_optimal_k first.")

    mclust_obj <- lpa_res$lpa_object
    if (!inherits(mclust_obj, "Mclust"))
      stop("lpa.result$lpa_object is not a Mclust object.")

    # align features
    train_vars <- names(mclust_obj$parameters$mean)  # works for univariate too
    if (is.null(train_vars)) {
      # multivariate: mean is a matrix (features x components)
      train_vars <- rownames(mclust_obj$parameters$mean)
    }
    if (is.null(train_vars)) train_vars <- colnames(mclust_obj$data)

    common <- intersect(colnames(val_scaled), train_vars)
    if (length(common) == 0) stop("No overlapping features between val data and LPA model.")

    val_sub <- val_scaled[, common, drop = FALSE]

    pred    <- mclust::predict.Mclust(mclust_obj, newdata = val_sub)
    labels  <- factor(pred$classification)

    if (verbose) {
      cat("LPA prediction (Mclust MAP classification):\n")
      print(table(labels))
    }
    object@info.data[[col_name]] <- labels[match(rownames(object@info.data),
                                                  rownames(val_sub))]
  }

  # ---------------------------- NMF (NNLS projection) ----------------------
  else if (method == "nmf") {

    fit <- train_object@cluster.results[["nmf.result"]][["best_estimate"]]
    if (is.null(fit))
      stop("No nmf best_estimate found in train_object. Run Sub_nmf_best_rank first.")

    W <- NMF::basis(fit)             # features x K  (NMF convention)
    feature_names <- rownames(W)

    common <- intersect(feature_names, colnames(val_scaled))
    if (length(common) == 0)
      stop("No overlapping features between val data and NMF W matrix.")
    if (length(common) < nrow(W) && verbose)
      warning(nrow(W) - length(common),
              " NMF features absent in val data - filled with 0.")

    # Build input matrix: features x samples (fill missing features with 0)
    val_mat <- matrix(0,
                      nrow = nrow(W),
                      ncol = nrow(val_scaled),
                      dimnames = list(feature_names, rownames(val_scaled)))
    val_mat[common, ] <- t(as.matrix(val_scaled[, common, drop = FALSE]))
    val_mat[val_mat < 0] <- 0   # NMF requires non-negative input

    # NNLS per sample: solve  W * h ~ v,  h >= 0
    val_H <- apply(val_mat, 2, function(v) nnls::nnls(W, v)$x)  # K x n_val
    group_idx   <- apply(val_H, 2, which.max)
    labels      <- factor(paste0(prefix, group_idx))
    names(labels) <- colnames(val_mat)

    if (verbose) {
      cat("NMF prediction (NNLS H-matrix projection):\n")
      print(table(labels))
    }
    object@info.data[[col_name]] <- labels[match(rownames(object@info.data),
                                                  names(labels))]
  }

  if (verbose)
    cat(sprintf("Subtype labels written to info.data$%s (%d samples).\n",
                col_name, sum(!is.na(object@info.data[[col_name]]))))

  return(object)
}


# -----------------------------------------------------------------------------
#' Create, Normalise, and Predict Subtypes for a Validation Set
#'
#' Convenience wrapper that chains:
#' \enumerate{
#'   \item \code{CreateSubtypingObject} -- wraps raw validation data.
#'   \item \code{Sub_apply_norm_params} -- applies training min-max parameters.
#'   \item \code{Sub_predict_subtypes}  -- predicts subtypes for each requested
#'     method.
#' }
#'
#' @param clean.data  Data frame of raw validation features
#'                    (samples x variables).
#' @param info.data   Data frame of clinical/meta variables for the validation
#'                    set (optional; same row order as \code{clean.data}).
#' @param norm_params Named list produced by \code{Sub_extract_norm_params},
#'   or the legacy \code{list(var = list(min=, max=))} format.
#'   Typically stored at \code{stat_obj@process.info$norm_params}.
#' @param train_object Trained \code{Subtyping} object.
#' @param methods     Character vector of methods to predict; any subset of
#'                    \code{c("nmf", "kmeans", "lpa")} (default: all three,
#'                    but only methods with a fitted model are actually run).
#' @param prefix      Subtype label prefix for NMF (default \code{"S"}).
#' @param verbose     Print progress (default \code{TRUE}).
#' @return A fully populated validation \code{Subtyping} object.
#' @examples
#' \dontrun{
#'   # Assuming 'train_obj' is a trained Subtyping object and 'norm_params' extracted
#'   val_obj <- Sub_create_val_object(
#'     clean.data = val_clean,
#'     info.data = val_info,
#'     norm_params = norm_params,
#'     train_object = train_obj,
#'     methods = c("kmeans", "lpa", "nmf")
#'   )
#' }
#' @export
Sub_create_val_object <- function(clean.data,
                                  info.data    = data.frame(),
                                  norm_params,
                                  train_object,
                                  methods  = c("nmf", "kmeans", "lpa"),
                                  prefix   = "S",
                                  verbose  = TRUE) {

  if (!is.data.frame(clean.data)) stop("'clean.data' must be a data frame.")
  if (!inherits(train_object, "Subtyping"))
    stop("'train_object' must be a 'Subtyping' object.")

  # ---- 1. Create Subtyping object ----
  if (verbose) cat("-- Step 1: Creating validation Subtyping object...\n")
  val_obj <- CreateSubtypingObject(
    clean.data = clean.data,
    info.data  = if (nrow(info.data) > 0) info.data else data.frame()
  )

  # ---- 2. Apply training normalisation parameters ----
  if (verbose) cat("-- Step 2: Applying training min-max parameters...\n")
  val_obj <- Sub_apply_norm_params(val_obj, norm_params = norm_params,
                                   verbose = verbose)

  # ---- 3. Predict subtypes for each requested method ----
  available <- list(
    kmeans = !is.null(train_object@cluster.results[["kmeans.result"]]),
    lpa    = !is.null(train_object@cluster.results[["lpa.result"]]),
    nmf    = !is.null(train_object@cluster.results[["nmf.result"]][["best_estimate"]])
  )

  for (m in methods) {
    if (!m %in% c("nmf", "kmeans", "lpa")) {
      warning("Unknown method '", m, "' - skipped."); next
    }
    if (!available[[m]]) {
      if (verbose)
        cat(sprintf("-- Step 3 [%s]: No fitted model found in train_object - skipped.\n", m))
      next
    }
    if (verbose) cat(sprintf("-- Step 3 [%s]: Predicting subtypes...\n", m))
    val_obj <- Sub_predict_subtypes(val_obj, train_object,
                                    method  = m,
                                    prefix  = prefix,
                                    verbose = verbose)
  }

  if (verbose) {
    cat("\nValidation Subtyping object ready.\n")
    cat("  scale.data :", dim(val_obj@scale.data),  "\n")
    cat("  info.data  :", dim(val_obj@info.data),   "\n")
    cluster_cols <- grep("^cluster_", colnames(val_obj@info.data), value = TRUE)
    if (length(cluster_cols)) {
      cat("  Subtype columns:", paste(cluster_cols, collapse = ", "), "\n")
      for (cc in cluster_cols) {
        cat(sprintf("    %s: ", cc)); print(table(val_obj@info.data[[cc]]))
      }
    }
  }

  return(val_obj)
}
