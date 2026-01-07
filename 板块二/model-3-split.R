#' Split Data into Training and Testing Sets
#'
#' Splits a dataset into training and testing subsets, with options for:
#' - Random sampling (default)
#' - Stratified sampling by group variable
#' - Custom ratio allocation
#' - Reproducible random splits via seed setting
#'
#' @param data A data frame containing the data to be split. Must be non-empty.
#' @param train_ratio Numeric value between 0 and 1 specifying the proportion of 
#'   data for training (default: 0.7). Must sum with test_ratio to equal 1.
#' @param test_ratio Numeric value between 0 and 1 specifying the proportion of 
#'   data for testing (default: 0.3). Must sum with train_ratio to equal 1.
#' @param seed Integer specifying random seed for reproducibility (default: 123).
#' @param group_col Optional character string specifying column name for stratified 
#'   sampling. If provided, sampling will preserve class distribution (default: NULL).
#'
#' @return A list containing two elements:
#'   - `training`: Data frame containing training subset
#'   - `testing`: Data frame containing testing subset
#'   Both subsets maintain original row names and structure.
#' @import caret
#' @import stats 
#' @export
#'
#' @examples
#' # Basic random split
#' data(mtcars)
#' split1 <- split_data(mtcars)
#' 
#' # Stratified split by cylinder count
#' split2 <- split_data(mtcars, group_col = "cyl")
#' 
#' # Custom 80/20 split with seed
#' split3 <- split_data(mtcars, train_ratio = 0.8, test_ratio = 0.2, seed = 42)
#' 
#' # Verify stratified sampling
#' table(mtcars$cyl) / nrow(mtcars)
#' table(split2$training$cyl) / nrow(split2$training)
#' @import caret
split_data <- function(
    data,
    train_ratio = 0.7,
    test_ratio = 0.3,
    seed = 123,
    group_col = NULL
) {
  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
    stop("Invalid input: 'data' must be a non-empty data frame.")
  }

  total_ratio <- train_ratio + test_ratio
  if (abs(total_ratio - 1) > .Machine$double.eps ^ 0.5) {
    stop("The sum of train_ratio and test_ratio must be equal to 1.")
  }

  set.seed(seed)

  if (!is.null(group_col)) {
    if (!group_col %in% colnames(data)) {
      stop(paste("Column", group_col, "not found in the data frame."))
    }
    train_index <- createDataPartition(data[[group_col]], p = train_ratio, list = FALSE)
  } else {
    train_index <- sample(nrow(data), size = round(nrow(data) * train_ratio))
  }

  train_data <- data[train_index, ]
  test_data <- data[-train_index, ]

  n_train <- nrow(train_data)
  n_test <- nrow(test_data)
  if (n_train < 1 || n_test < 1) {
    stop("Invalid ratio values: insufficient data for training or testing.")
  }

  cat("Data split completed.\n")
  cat("Training data rows:", n_train, "\n")
  cat("Testing data rows:", n_test, "\n")

  return(list(training = train_data, testing = test_data))
}


#' Split Data into Training and Test Sets
#'
#' This function splits data into training and test sets, with optional CSV export capabilities.
#' It handles both data frames and Train_Model objects, and supports stratified sampling.
#'
#' @param object An input object of class 'Train_Model' or a data frame
#' @param train_ratio Proportion of data for training set (default: 0.7)
#' @param test_ratio Proportion of data for test set (default: 0.3)
#' @param seed Random seed for reproducibility (default: 123)
#' @param group_col Name of the grouping column for stratified sampling (optional)
#' @param save_data Logical indicating whether to save splits as CSV files (default: FALSE)
#' @param train_filename Name for training set CSV file (default: "train_data.csv")
#' @param test_filename Name for test set CSV file (default: "test_data.csv")
#' @param save_dir Directory path for saving CSV files (default: here::here("ModelData","Data"))
#'
#' @return Depending on input:
#'   - For 'Train_Model' objects: Returns updated object with split data
#'   - For data frames: Returns list containing train and test data frames
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with data frame
#' object_model <- SplitDatModel(object_model, train_ratio = 0.8, test_ratio = 0.2)
#' 
#' # With CSV export
#' object_model <- SplitDatModel(object_model, save_data = TRUE)

#' }
SplitDatModel <- function(
    object = NULL,
    train_ratio = 0.7,
    test_ratio = 0.3,
    seed = 123,
    group_col = NULL,
    save_data = TRUE,
    train_filename = "train_data.csv",
    test_filename = "test_data.csv",
    save_dir = here::here("ModelData","Data")
) {
  if (inherits(object, 'Train_Model')) {
    data <- slot(object, "clean.df")
    if (is.null(group_col)) {
      group_col <- object@group_col
    }
    return_object_modelect <- TRUE
  } else if (is.data.frame(object)) {
    data <- object
    return_object_modelect <- FALSE
  } else {
    stop("Input must be an object of class 'Train_Model' or a data frame.")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input.")
  }
  
  total_ratio <- train_ratio + test_ratio
  if (abs(total_ratio - 1) > .Machine$double.eps ^ 0.5) {
    stop("The sum of train_ratio and test_ratio must be equal to 1.")
  }
  cat("Train ratio:", train_ratio, "Test ratio:", test_ratio, "\n")
  
  datalist <- split_data(
    data = data,
    train_ratio = train_ratio,
    test_ratio = test_ratio,
    seed = seed,
    group_col = group_col
  )
  
  if (save_data) {
    tryCatch({
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("Created directory:", save_dir, "\n")
      }
      
      train_path <- file.path(save_dir, train_filename)
      write.csv(datalist$train, file = train_path, row.names = FALSE)
      cat("Training set saved to:", train_path, "\n")
      
      test_path <- file.path(save_dir, test_filename)
      write.csv(datalist$test, file = test_path, row.names = FALSE)
      cat("Test set saved to:", test_path, "\n")
    }, error = function(e) {
      warning("Failed to save CSV files: ", e$message)
    })
  }
  
  if (return_object_modelect) {
    object@split.data <- datalist
    cat("Updating 'Train_Model' object...\n")
    cat("The 'Train_Model' object has been updated with the following slots:\n")
    cat("- 'split.data' slot updated.\n")
    return(object)
  }
  
  return(datalist)
}
