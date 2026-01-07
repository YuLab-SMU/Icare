#' Balance Dataset Using Various Sampling Methods
#'
#' This function balances imbalanced datasets using oversampling, undersampling,
#' combined sampling, or SMOTE (Synthetic Minority Oversampling Technique).
#'
#' @param data A data frame containing the dataset to be balanced
#' @param group_col Character string specifying the name of the grouping variable (default: "group")
#' @param method Balancing method to use: "over" (oversampling), "under" (undersampling),
#'               "both" (combined sampling), or "smote" (SMOTE algorithm) (default: "both")
#' @param N Desired sample size after balancing (NULL for automatic calculation)
#' @param K Number of nearest neighbors to consider for SMOTE (default: 5)
#' @param dup_size Oversampling rate for SMOTE (0 = balanced, 1 = double minority class, etc.)
#' @param seed Random seed for reproducibility (default: 1234)
#'
#' @return A data frame with balanced class distribution
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data <- data.frame(
#'   group = factor(c(rep(0, 100), rep(1, 30))),
#'   var1 = rnorm(130),
#'   var2 = rnorm(130)
#' )
#'
#' # Oversample minority class
#' balanced_over <- balance_data_process(data, method = "over")
#'
#' # Undersample majority class
#' balanced_under <- balance_data_process(data, method = "under")
#'
#' # Use SMOTE
#' balanced_smote <- balance_data_process(data, method = "smote")
#' }
#'
#'
#' @importFrom ROSE ovun.sample
#' @importFrom smotefamily SMOTE
#' @export
balance_data_process <- function(data,
                                 group_col = "group",
                                 method = "both",
                                 N = NULL,
                                 K = 5,
                                 dup_size = 0,
                                 seed = 1234) {


  data[[group_col]] <- as.factor(data[[group_col]])
  valid_methods <- c("over", "under", "both", "smote")
  if (!method %in% valid_methods) {
    stop("Invalid method. Choose one of: ", paste(valid_methods, collapse = ", "))
  }

  if (is.null(N)) {
    class_counts <- table(data[[group_col]])
    if (method == "over") {
      N <- 2 * max(class_counts)
    } else if (method == "under") {
      N <- 2 * min(class_counts)
    } else if (method == "both") {
      N <- sum(class_counts)
    } else {
      N <- NULL
    }
  }

  set.seed(seed)
  if (method %in% c("over", "under", "both")) {
    formula <- as.formula(paste(group_col, "~ ."))
    balanced_data <- ROSE::ovun.sample(
      formula,
      data = data,
      method = method,
      N = N,
      seed = seed
    )$data
  } else {
    data<-convert_to_numeric(data)
    data<-convert_two_class_to_binary(data)

    X <- data[, !colnames(data) %in% group_col, drop = FALSE]
    y <- as.numeric(data[[group_col]])

    if (!all(sapply(X, is.numeric))) {
      stop("SMOTE requires all features to be numeric. Convert factors to dummy variables first.")
    }

    smote_result <- smotefamily::SMOTE(
      X = X,
      target = y,
      K = K,
      dup_size = dup_size
    )

    balanced_data <- data.frame(smote_result$data)
    colnames(balanced_data)[colnames(balanced_data) == "class"] <- group_col
  }


  cat("\nBalance Result:\n")
  cat("Original counts:", paste(names(table(data[[group_col]])),
                                table(data[[group_col]]), sep = "=", collapse = ", "), "\n")
  cat("Balanced counts:", paste(names(table(balanced_data[[group_col]])),
                                table(balanced_data[[group_col]]), sep = "=", collapse = ", "), "\n")
  cat("Method used:", method, "\n")


  return(balanced_data)
}



#' Visualize the balance of the dataset before and after balancing.
#'
#' This function generates a bar plot to visualize the distribution of class labels in the original and balanced datasets. The plot compares the class distribution between the two datasets and saves it as a PDF if required.
#'
#' @import ggplot2
#' @import wesanderson
#' @importFrom ggprism theme_prism
#' @importFrom scales percent
#' @import here
#' @param original_data A data frame containing the original (imbalanced) dataset.
#' @param balanced_data A data frame containing the balanced dataset.
#' @param group_col A string specifying the column name that represents the class labels in both datasets. Default is "group".
#' @param palette_name A string specifying the name of the color palette for the plot. Default is "Royal1".
#' @param base_size An integer specifying the base font size for the plot. Default is 14.
#' @param save_plots A logical value indicating whether to save the plot as a PDF. Default is `TRUE`.
#' @param save_dir A string specifying the directory where the plot will be saved. Default is "ModelData/balacing_info" using `here()`.
#' @param plot_width A numeric value specifying the width of the plot when saved. Default is 5.
#' @param plot_height A numeric value specifying the height of the plot when saved. Default is 5.
#'
#' @returns A ggplot object representing the class distribution before and after balancing.
#' @export
#'
#' @examples
#' plot <- visualize_balance(original_data = original_data, balanced_data = balanced_data, group_col = "group")
#' plot <- visualize_balance(original_data = original_data, balanced_data = balanced_data, save_plots = FALSE)

visualize_balance <- function(
    original_data,
    balanced_data,
    group_col = "group",
    palette_name = "Royal1",
    base_size = 14,
    save_plots = TRUE,
    save_dir = here("ModelData", "balacing_info"),
    plot_width = 5,
    plot_height = 5) {


  original_data$dataset <- "Original"
  balanced_data$dataset <- "Balanced"
  combined_data <- rbind(original_data, balanced_data)

  plot <- ggplot(combined_data, aes(x = .data[[group_col]], fill = dataset)) +
    geom_bar(position = "dodge") +
    geom_text(stat = 'count', aes(label = ..count..), position = position_dodge(width = 0.9), vjust = -0.5, size = 5) +
    geom_text(stat = 'count', aes(label = scales::percent(..count../sum(..count..), accuracy = 0.1)),
              position = position_dodge(width = 0.9), vjust = 1.5, size = 4) +
    labs(title = "Before and After Balancing",
         x = group_col,
         y = "Count",
         fill = "Dataset") +
    scale_fill_manual(values = wes_palette(palette_name)) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12)) +
    theme_prism(base_size = base_size)

  if (save_plots) {
    ggsave(filename = file.path(save_dir, "class_distribution_balance.pdf"),
           plot = plot,
           width = plot_width, height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "class_distribution_balance.pdf"), "\n")
  }

  return(plot)
}


#' Balance and visualize data distribution before and after balancing.
#'
#' This function detects class imbalance in the dataset and balances it if necessary using specified methods. It also visualizes the class distribution before and after balancing, and optionally saves the plot as a PDF.
#'
#' @import ggplot2
#' @import wesanderson
#' @importFrom ggprism theme_prism
#' @importFrom scales percent
#' @import here
#' @import stats
#' @param data A data frame containing the data to be balanced.
#' @param group_col A string specifying the column name that represents the class labels.
#' @param method A string specifying the balancing method. Can be "over", "under", or "both". Default is "both".
#' @param N A numeric value specifying the number of samples to generate after balancing. If `NULL`, it is automatically calculated based on the chosen method.
#' @param seed An integer to set the random seed for reproducibility. Default is 123.
#' @param palette_name A string specifying the color palette name for the plot. Default is "Royal1".
#' @param imbalance_threshold A numeric value representing the threshold for class imbalance. If the class ratio is below this threshold, balancing is triggered. Default is 0.15.
#' @param sample_size_threshold A numeric value representing the threshold for sample size. If the number of samples is below this threshold, balancing is triggered. Default is 1500.
#' @param force_balance A logical value indicating whether to force balancing regardless of the class distribution. Default is `FALSE`.
#' @param save_plots A logical value indicating whether to save the plot as a PDF. Default is `TRUE`.
#' @param save_dir A string specifying the directory where the plot will be saved. Default is "ModelData/balacing_info" using `here()`.
#' @param plot_width A numeric value specifying the width of the plot when saved. Default is 5.
#' @param plot_height A numeric value specifying the height of the plot when saved. Default is 5.
#' @param base_size An integer specifying the base font size for the plot. Default is 14.
#'
#' @returns A list containing the following elements:
#' - `original_data`: The original (imbalanced) data.
#' - `balanced_data`: The balanced dataset.
#' - `method`: The method used for balancing.
#' - `balanced_plot`: The plot visualizing the class distribution before and after balancing.
#'
#' @export
#'
#' @examples
#' result <- balance_and_visualize_data(data = dataset, group_col = "group")
#' result <- balance_and_visualize_data(data = dataset, group_col = "group", force_balance = TRUE)

balance_and_visualize_data <- function(
    data,
    group_col,
    method = "both",
    N = NULL,
    K = 5,
    dup_size = 0,
    seed = 123,
    palette_name = "Royal1",
    imbalance_threshold = 0.15,
    sample_size_threshold = 1500,
    force_balance = FALSE,
    save_plots = TRUE,
    save_dir = here("ModelData", "balacing_info"),
    plot_width = 5,
    plot_height = 5,
    base_size = 14) {

  class_counts <- table(data[[group_col]])
  class_ratio <- min(class_counts) / sum(class_counts)

  cat("Class counts before balancing:\n")
  print(class_counts)

  cat("\nforce_balance:", force_balance, "\n")
  cat("class_ratio:", class_ratio, "\n")
  cat("class_ratio < imbalance_threshold:", class_ratio < imbalance_threshold, "\n")
  cat("nrow(data):", nrow(data), "\n")
  cat("nrow(data) < sample_size_threshold:", nrow(data) < sample_size_threshold, "\n\n")
  method=method
  if (force_balance || (class_ratio < imbalance_threshold && nrow(data) < sample_size_threshold)) {
    cat("Data imbalance detected or force_balance is TRUE, or sample size is below the threshold. Balancing data.\n")

    balanced_data <- balance_data_process(
      data = data,
      group_col = group_col,
      method = method,
      N = N,
      K = K,
      dup_size = dup_size,
      seed = seed
    )

    class_counts_after <- table(balanced_data[[group_col]])
    class_ratio_after <- min(class_counts_after) / sum(class_counts_after)

    cat("Class ratio after balancing:", round(class_ratio_after, 4), "\n")
    cat("Class counts after balancing:\n")
    print(class_counts_after)

    plot <- visualize_balance(original_data = data,
                              balanced_data = balanced_data,
                              group_col = group_col,
                              palette_name = palette_name,
                              base_size = base_size,
                              save_plots = save_plots,
                              save_dir = save_dir,
                              plot_width = plot_width,
                              plot_height = plot_height)
    print(plot)

    balance_result <- list(original_data = data,
                           balanced_data = balanced_data,
                           method = method,
                           balanced_plot = plot)
    return(balance_result)
  } else {
    cat("Data is balanced within threshold or force_balance is FALSE. No balancing performed.\n")
    cat("Class ratio (below threshold):", round(class_ratio, 4), "\n\n")

    balance_result <- list(balanced_data = data,
                           method = NULL,
                           balanced_plot = NULL)

    return(balance_result)
  }
}


#' Balance Training Dataset Using Various Sampling Methods
#'
#' This function performs class balancing on the training dataset only (while keeping test/validation
#' sets unchanged) using various sampling techniques including oversampling, undersampling,
#' combined sampling, or SMOTE (Synthetic Minority Oversampling Technique). It handles both
#' Train_Model objects and raw data list inputs, with optional CSV export of the balanced training data.
#'
#' @param object Either a Train_Model object or a list containing training data and optionally
#'              testing/validation data
#' @param group_col Character string specifying the name of the grouping variable (default: "group")
#' @param method Balancing method to use: "over" (oversampling), "under" (undersampling),
#'               "both" (combined sampling), or "smote" (SMOTE algorithm) (default: "both")
#' @param N Desired sample size after balancing (NULL for automatic calculation based on method)
#' @param palette_name Color palette name for visualization plots (default: "Royal1")
#' @param seed Random seed for reproducibility (default: 123)
#' @param imbalance_threshold Minimum class ratio threshold to trigger balancing
#'                           (default: 0.15, range 0-1)
#' @param sample_size_threshold Maximum sample size threshold to apply balancing
#'                             (default: 1500 observations)
#' @param force_balance Logical, whether to force balancing regardless of thresholds
#'                     (default: FALSE)
#' @param save_plots Logical, whether to save balance visualization plots (default: TRUE)
#' @param save_dir Directory path to save plots and CSV files (default: here::here("Data","ModelData"))
#' @param plot_width Plot width in inches (default: 5)
#' @param plot_height Plot height in inches (default: 5)
#' @param base_size Base font size for plots (default: 14)
#' @param K Number of nearest neighbors to consider for SMOTE (default: 5)
#' @param dup_size Oversampling rate for SMOTE (0 = balanced, 1 = double minority class, etc.)
#' @param save_data Logical, whether to save balanced training data as CSV file (default: TRUE)
#' @param train_filename Name for the output CSV file containing balanced training data
#'                      (default: "train_data.csv")
#'
#' @return Depending on input type:
#'   - For Train_Model objects: Returns updated object with balanced training data in split.data slot
#'     and balancing information in process.info slot
#'   - For lists: Returns a list containing balanced training data and unchanged other datasets,
#'     along with balancing parameters
#'
#' @examples
#' \dontrun{
#' # Example with Train_Model object
#' object_model <- BalanceData(
#'   object = object_model,
#'   method = "smote",
#'   imbalance_threshold = 0.2
#' )
#'
#' }
#' @export
BalanceData <- function(
    object,
    group_col = "group",
    method = "both",
    N = NULL,
    palette_name = "Royal1",
    seed = 123,
    imbalance_threshold = 0.15,
    sample_size_threshold = 1500,
    force_balance = FALSE,
    save_plots = TRUE,
    save_dir = here::here("ModelData","Data"),  # Combined into one parameter
    plot_width = 5,
    plot_height = 5,
    base_size = 14,
    K = 5,
    dup_size = 0,
    save_data = TRUE,
    train_filename = "train_data.csv") {

  set.seed(seed)

  cat("\n=== Starting Training Data Balancing Process ===\n")
  cat("Balancing method:", method, "\n")
  if (!is.null(N)) cat("Target sample size per group:", N, "\n")

  if (inherits(object, 'Train_Model')) {
    group_col <- object@group_col %||% group_col
    training <- slot(object, "split.data")[["training"]]
    testing <- slot(object, "split.data")[["testing"]]
    validation <- if ("validation" %in% names(slot(object, "split.data")))
      slot(object, "split.data")[["validation"]] else NULL
    external_validation <- if ("external_validation" %in% names(slot(object, "split.data")))
      slot(object, "split.data")[["external_validation"]] else NULL
    return_object_modelect <- TRUE
  } else if (is.list(object)) {
    training <- object$training
    testing <- object$testing
    validation <- if ("validation" %in% names(object)) object$validation else NULL
    external_validation <- if ("external_validation" %in% names(object)) object$external_validation else NULL
    return_object_modelect <- FALSE
  } else {
    stop("Input must be an object of class 'Train_Model' or a list")
  }

  if (is.null(training) || nrow(training) == 0) {
    stop("No valid training data found in the input")
  }

  if (length(group_col) == 0 || is.null(group_col) || !group_col %in% colnames(training)) {
    cat("-> Group column is not valid, setting to NULL.\n")
    group_col <- NULL
  }

  needs_balancing <- function(data, group_col, force_balance, imbalance_threshold, sample_size_threshold) {
    if (is.null(group_col)) return(FALSE)
    class_counts <- table(data[[group_col]])
    class_ratio <- min(class_counts) / sum(class_counts)
    return(force_balance || (class_ratio < imbalance_threshold && nrow(data) < sample_size_threshold))
  }

  cat("\n Balancing TRAINING data...\n")
  if (needs_balancing(training, group_col, force_balance, imbalance_threshold, sample_size_threshold)) {
    plot_save_dir <- file.path(save_dir, "balancing_info", "train_balance")
    if (!dir.exists(plot_save_dir)) dir.create(plot_save_dir, recursive = TRUE)
  } else {
    plot_save_dir <- NULL
  }

  results_train <- balance_and_visualize_data(
    data = training,
    group_col = group_col,
    method = method,
    N = N,
    seed = seed,
    palette_name = palette_name,
    imbalance_threshold = imbalance_threshold,
    sample_size_threshold = sample_size_threshold,
    force_balance = force_balance,
    save_plots = save_plots,
    save_dir = plot_save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    K = K,
    dup_size = dup_size
  )

  balanced_train_data <- as.data.frame(results_train$balanced_data)
  cat("-> Training data balanced. Original rows:", nrow(training),
      "| Balanced rows:", nrow(balanced_train_data), "\n")

  if (save_data) {
    tryCatch({
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("Created directory:", save_dir, "\n")
      }

      full_path <- file.path(save_dir, train_filename)
      write.csv(balanced_train_data, file = full_path, row.names = FALSE)
      cat("-> Saved balanced training data to:", full_path, "\n")
    }, error = function(e) {
      warning("Failed to save CSV file: ", e$message)
    })
  }

  if (return_object_modelect) {
    object@split.data[["training"]] <- balanced_train_data

    balance_result <- list(
      balance_result_train = results_train,
      balancing_method = method,
      balancing_params = list(
        N = N,
        K = K,
        dup_size = dup_size,
        method = method,
        imbalance_threshold = imbalance_threshold,
        sample_size_threshold = sample_size_threshold,
        force_balance = force_balance,
        seed = seed
      )
    )

    object@process.info[["balance_result"]] <- balance_result

    cat("\nThe 'Train_Model' object has been updated with:\n")
    cat("- Balanced training data in split.data slot\n")
    cat("- Balancing results in process.info slot\n")

    return(object)
  } else {
    result <- list(
      training = balanced_train_data,
      testing = testing,
      validation = validation,
      external_validation = external_validation,
      balance_info = list(
        method = method,
        params = list(
          N = N,
          K = K,
          dup_size = dup_size,
          imbalance_threshold = imbalance_threshold,
          sample_size_threshold = sample_size_threshold,
          force_balance = force_balance,
          seed = seed
        )
      )
    )

    cat("\nReturning list with balanced training data (other datasets unchanged)\n")
    return(result)
  }
}
