
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
                              normalize_method = normalize_method,
                              max_unique_values = max_unique_values)
  nmdat<-nm_result$scaled_data
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

