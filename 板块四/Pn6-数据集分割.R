SplitData <- function(
    data,
    unwanted_var =NULL,
    train_ratio = 0.7,
    test_ratio = 0.3,
    seed = 123
) {
  cat("Starting data split process...\n")

  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
    stop("Invalid input: 'data' must be a non-empty data frame.")
  }

  total_ratio <- train_ratio + test_ratio
  if (abs(total_ratio - 1) > .Machine$double.eps ^ 0.5) {
    stop("The sum of train_ratio and test_ratio must be equal to 1.")
  }

  set.seed(seed)
  data <- data[, !names(data) %in% unwanted_var]
  shuffled_data <- data[sample(nrow(data)), ]
  original_row_names <- rownames(shuffled_data)
  n <- nrow(shuffled_data)

  n_train <- round(n * train_ratio)
  n_test <- n - n_train

  if (n_train < 1 || n_test < 1) {
    stop("Invalid ratio values: insufficient data for training or testing.")
  }

  train_data <- shuffled_data[1:n_train, ]
  test_data <- shuffled_data[(n_train + 1):n, ]

  rownames(train_data) <- original_row_names[1:n_train]
  rownames(test_data) <- original_row_names[(n_train + 1):n]

  cat("Data split completed successfully.\n")
  return(list(
    training = train_data,
    testing = test_data
  ))
}

SplitDataPrognosiX <- function(
    object = NULL,
    unwanted_var =NULL,
    train_ratio = 0.7,
    test_ratio = 0.3,
    seed = 123,
    use_subgroup_data=FALSE
) {

  if (inherits(object, 'PrognosiX')) {
    if (use_subgroup_data) {
      data <- slot(object, "sub.data")
      cat("Using subgroup analysis data...\n")
    } else {
      data <- slot(object, "survival.data")
      cat("Using original survival data...\n")
    }
  } else if (is.data.frame(object)) {
    cat("Using provided data frame...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }

  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input.")
  }

  cat("Proceeding with data splitting...\n")
  datalist <- SplitData(
    data = data,
    train_ratio = train_ratio,
    test_ratio = test_ratio,
    seed = seed
  )

  if (inherits(object, 'PrognosiX')) {
    object@split.data <- datalist
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'split.data' slot updated.\n")
    return(object)
  }

  cat("Data split process for 'PrognosiX' completed successfully.\n")
  return(datalist)
}
