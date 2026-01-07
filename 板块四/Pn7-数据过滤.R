Prognos_filter_features <- function(
    object,
    time_col = "time",
    status_col = "status",
    use_filtered_features = TRUE 
) {
  
  if (!inherits(object, 'PrognosiX')) {
    stop("Input must be an object of class 'PrognosiX'.")
  }
  
  if (!is.list(object@split.data) ||
      !all(c("training", "testing") %in% names(object@split.data))) {
    stop("The 'split.data' slot in the 'PrognosiX' object must contain 'train' and 'test' datasets.")
  }
  
  train_data <- slot(object, "split.data")[["training"]]
  test_data <- slot(object, "split.data")[["testing"]]
  
  if (is.null(train_data) || !is.data.frame(train_data) || nrow(train_data) == 0) {
    stop("No valid training data found in the 'PrognosiX' object.")
  }
  if (is.null(test_data) || !is.data.frame(test_data) || nrow(test_data) == 0) {
    stop("No valid test data found in the 'PrognosiX' object.")
  }
  
  if (!use_filtered_features) {
    cat("Using original training and testing datasets without filtering features.\n")
    
    object@filtered.set <- list(
      training = train_data,
      testing = test_data
    )
    cat("Original data stored in the 'filtered.set' slot.\n")
    return(object)
  }else{
    if (!"lasso_results" %in% names(object@feature.result) ||
      is.null(object@feature.result[["lasso_results"]][["important_vars"]])) {
    stop("No important features found in the LASSO results.")
  }
  
  best_features <- object@feature.result[["lasso_results"]][["important_vars"]]
  filtered_features <- c(best_features, time_col, status_col) 
  
  filtered_train <- train_data[, filtered_features, drop = FALSE]
  filtered_test <- test_data[, filtered_features, drop = FALSE]
  
  cat("Data filtered to retain best features.\n")
  cat("Total features retained:", length(filtered_features), "\n")
  
  object@filtered.set <- list(
    training = filtered_train,
    testing = filtered_test
  )
  cat("Filtered data stored in the 'filtered.set' slot.\n")
  
  return(object)
  }
}
