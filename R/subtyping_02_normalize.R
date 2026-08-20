#' Normalize Process for Subtyping
#'
#' @param object Subtyping object or data frame.
#' @param normalize_method Normalization method.
#' @param group_col Group column.
#' @param max_unique_values Max unique values.
#' @return If input is a Subtyping object, returns the modified object with
#'   normalized data stored in the 'scale.data' slot. If input is a data frame,
#'   returns the normalized data frame.
#' @export
#' @examples
#' set.seed(1)
#' demo_df <- data.frame(
#'   id    = paste0("S", 1:60),
#'   group = rep(c(0, 1), each = 30),
#'   feat1 = c(rnorm(30, 5, 1), rnorm(30, 8, 1)),
#'   feat2 = c(rnorm(30, 2, 0.5), rnorm(30, 4, 0.5)),
#'   feat3 = rnorm(60, 10, 2)
#' )
#' stat_obj <- CreateStatObject(raw.data = demo_df, clean.data = demo_df,
#'                              group_col = "group", na.action = "allow")
#' sub_obj <- ConvertObject(stat_obj, to = "Subtyping")
#'
#' set.seed(1)
#' split <- SplitSubtypingObject(sub_obj, p = 0.7, stratify_by = "group")
#'
#' # Fit normalization on the training split only.
#' sub_train <- Sub_normalize_process(split$train, normalize_method = "min_max")
#' dim(sub_train@scale.data)
#'
#' # Extract training-derived parameters and apply them to the validation split
#' # (avoids leaking validation-set values into their own scaling).
#' norm_params <- Sub_extract_norm_params(sub_train, verbose = FALSE)
#' sub_test <- Sub_apply_norm_params(split$test, norm_params = norm_params, verbose = FALSE)
#' dim(sub_test@scale.data)
Sub_normalize_process <- function(object,
                                   normalize_method = "min_max_scale",
                                   group_col = "group",
                                   max_unique_values = 5) {
  cat("Input object class:", class(object), "\n")

  if (inherits(object, 'Subtyping')) {
    data <- slot(object, "clean.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'sub' or a data frame")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }

  cat("Starting normalization process...\n")

  nm_result <- normalize_data(data,
                              method = normalize_method,
                              group_col = group_col,  
                              max_unique_values = max_unique_values)
  # Note: normalize_data returns data frame, not list with 'scaled_data' in current implementation in module 1.
  # Checking module 1 normalize_data implementation...
  # It returns normalized_data directly.
  # So nm_result IS the data.
  nmdat <- nm_result
  
  if (inherits(object, 'Subtyping')) {
    if (!is.null(slotNames(object)) && "scale.data" %in% slotNames(object)) {
      object@scale.data <- nmdat
      cat("Normalized data stored in 'scale.data' slot.\n")
    } else {
      stop("The 'sub' object does not have a 'scale.data' slot.")
    }
    return(object)
  }

  cat("Normalization complete, returning normalized data frame.\n")
  return(nmdat)
}
