#' Calculate Maximal Information Coefficients (MIC) Between Features and Target
#'
#' This function computes the Maximal Information Coefficient (MIC) between each column in the feature matrix \code{x}
#' and the target vector \code{y}, using the \code{mine} function from the \code{minerva} package.
#' @importFrom minerva mine
#' @param x A data frame or matrix containing the features (predictor variables).
#' @param y A numeric vector representing the target variable (response).
#'
#' @returns A data frame with two columns:
#' \describe{
#'   \item{\code{Feature}}{The names of the features (columns of \code{x}).}
#'   \item{\code{Maximal Information Coefficient}}{The corresponding MIC values between each feature and \code{y}.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' x <- data.frame(
#'   var1 = rnorm(100),
#'   var2 = runif(100),
#'   var3 = rpois(100, lambda = 3)
#' )
#' y <- rnorm(100)
#' mic_info <- mic_information(x, y)
#' print(mic_info)
#' }
mic_information <- function(x, y) {
  colnames(x) <- make.names(colnames(x))

  result <- minerva::mine(x, y)
  mic_values <- result$MIC
  feature_names <- colnames(x)
  mic_data <- data.frame(Feature = feature_names, "Maximal Information Coefficient" = mic_values)
  colnames(mic_data) <- c("Feature", "Maximal Information Coefficient")

  return(mic_data)
}

#' Compute Information Value (IV) for Each Feature
#'
#' This function calculates the Information Value (IV) of each feature with respect to the binary target variable.
#' IV is a commonly used metric in feature selection for binary classification tasks, especially in credit scoring.
#'
#' @param features A data frame or matrix containing the predictor variables.
#' @param labels A binary vector representing the target variable (0 or 1).
#'
#' @importFrom Information create_infotables
#' @returns A data frame containing:
#' \describe{
#'   \item{\code{Feature}}{The names of the predictor variables.}
#'   \item{\code{IV}}{The computed Information Value for each variable.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(Information)
#' features <- data.frame(
#'   age = c(25, 45, 35, 50, 23),
#'   income = c(50000, 80000, 62000, 75000, 40000)
#' )
#' labels <- c(0, 1, 0, 1, 0)
#' iv_result <- Information_Value(features, labels)
#' print(iv_result)
#' }
Information_Value <- function(features, labels) {
  dt_s <- cbind(features, label = labels)
  colnames(dt_s) <- make.names(colnames(dt_s))

  info_table <- Information::create_infotables(dt_s, y = "label")
  iv_values <- info_table[["Summary"]][["IV"]]

  iv_data <- data.frame(Feature = info_table[["Summary"]][["Variable"]], IV = iv_values)
  iv_data <- iv_data[order(-iv_data$IV), ]

  return(iv_data)
}


#' Compute Mutual Information Between Features and Label
#'
#' This function calculates the mutual information (MI) between each feature in the input dataset
#' and a binary classification label. Mutual information is a measure of the dependency between
#' two variables and is commonly used in feature selection for classification tasks.
#' @importFrom infotheo mutinformation
#' @param dataSet A data frame or matrix containing the input features. Each column represents a feature.
#' @param label A binary vector or factor representing the class labels (e.g., 0 and 1).
#'
#' @returns A data frame containing:
#' \describe{
#'   \item{\code{Feature}}{The names of the input features.}
#'   \item{\code{MI}}{The mutual information value for each feature with respect to the label.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(infotheo)
#' dataSet <- data.frame(
#'   var1 = rnorm(100),
#'   var2 = runif(100),
#'   var3 = rnorm(100, mean = 5)
#' )
#' label <- sample(0:1, 100, replace = TRUE)
#' result <- mutual_information(dataSet, label)
#' print(result)
#' }
mutual_information <- function(dataSet, label) {
  colnames(dataSet) <- make.names(colnames(dataSet))

  if(is.numeric(label) && all(label %in% c(0,1))) {
    label <- factor(label, levels = c(0,1))
  }

  num <- ncol(dataSet)
  mi_values <- numeric(num)

  for (i in 1:num) {
    discrete_feature <- cut(dataSet[, i], breaks = 10)
    mi_values[i] <- infotheo::mutinformation(discrete_feature, label)
  }

  col_names <- colnames(dataSet)
  mi_data <- data.frame(Feature = col_names, MI = mi_values)
  mi_data <- mi_data[order(-mi_data$MI), ]

  return(mi_data)
}

#' Feature Selection Based on Information Value and AUC
#'
#' This function performs feature selection by ranking features using Information Value (IV)
#' and evaluating model performance using AUC (Area Under the ROC Curve). It incrementally
#' adds top IV-ranked features to a linear model and selects the subset that yields the highest AUC.
#' @importFrom pROC roc auc
#' @importFrom Information create_infotables
#' @param X A data frame or matrix containing the predictor variables (features).
#' @param Y A numeric or factor vector representing the binary outcome variable (0/1).
#' @param max_features An integer specifying the maximum number of features to evaluate.
#'
#' @returns A list containing:
#' \describe{
#'   \item{\code{SelectedFeatures}}{The names of the selected features that yielded the highest AUC.}
#'   \item{\code{AUC}}{A numeric vector of AUC values for each subset of features evaluated.}
#'   \item{\code{Method}}{A string indicating the feature selection method used: "Information Value".}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(pROC)
#' library(Information)
#' X <- data.frame(a = rnorm(100), b = runif(100), c = rnorm(100, 2))
#' Y <- sample(0:1, 100, replace = TRUE)
#' result <- feature_selection_iv(X, Y, max_features = 3)
#' print(result)
#' }
feature_selection_iv <- function(X, Y, max_features) {
  colnames(X) <- make.names(colnames(X))

  iv_data <- Information_Value(X, Y)

  if (any(is.na(iv_data$Feature))) {
    warning("NA values found in iv_data$Feature. Removing these rows.")
    iv_data <- iv_data[complete.cases(iv_data), ]
  }

  iv_data$Feature <- make.names(iv_data$Feature)

  if (nrow(iv_data) < max_features) {
    warning("Number of features in iv_data is less than max_features. Adjusting max_features.")
    max_features <- nrow(iv_data)
  }

  if (max_features == 0) {
    stop("No valid features found in iv_data.")
  }

  auc_values <- numeric(max_features)
  selected_features <- character(max_features)

  for (i in 1:max_features) {
    current_features <- iv_data$Feature[1:i]
    if (any(is.na(current_features))) {
      stop("Encountered NA in feature list.")
    }

    formula_str <- paste("Y ~", paste(current_features, collapse = "+"))
    formula_model <- as.formula(formula_str)

    missing_features <- setdiff(current_features, colnames(X))
    if (length(missing_features) > 0) {
      stop("The following features are not found in X: ", paste(missing_features, collapse=", "))
    }

    stepwise_model <- lm(formula_model, data = data.frame(Y = Y, X))
    predictions <- predict(stepwise_model, newdata = data.frame(Y = Y, X), type = "response")
    roc_auc <- pROC::auc(pROC::roc(Y, predictions))

    auc_values[i] <- roc_auc
    selected_features[i] <- paste(names(coef(stepwise_model))[-1], collapse = ", ")
  }

  return(list(SelectedFeatures = selected_features[which.max(auc_values)],
              AUC = auc_values,
              Method = "Information Value"))
}

#' Feature Selection Based on Mutual Information and AUC
#'
#' This function performs feature selection by ranking features using Mutual Information (MI)
#' and evaluating model performance using AUC (Area Under the ROC Curve). It incrementally
#' adds the top MI-ranked features to a linear model and selects the feature subset with the highest AUC.
#' @importFrom pROC roc auc
#' @importFrom infotheo mutinformation
#' @param X A data frame or matrix containing the predictor variables (features).
#' @param Y A numeric or factor vector representing the binary outcome variable (0/1).
#' @param max_features An integer specifying the maximum number of features to evaluate.
#'
#' @returns A list containing:
#' \describe{
#'   \item{\code{SelectedFeatures}}{The names of the selected features that yielded the highest AUC.}
#'   \item{\code{AUC}}{A numeric vector of AUC values corresponding to the evaluated feature subsets.}
#'   \item{\code{Method}}{A string indicating the feature selection method used: "Mutual Information".}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(pROC)
#' library(infotheo)
#' X <- data.frame(a = rnorm(100), b = runif(100), c = rnorm(100, 2))
#' Y <- sample(0:1, 100, replace = TRUE)
#' result <- feature_selection_mi(X, Y, max_features = 3)
#' print(result)
#' }
feature_selection_mi <- function(X, Y, max_features) {
  colnames(X) <- make.names(colnames(X))



  mi_data <- mutual_information(X, Y)

  if (any(is.na(mi_data$Feature))) {
    warning("NA values found in mi_data$Feature. Removing these rows.")
    mi_data <- mi_data[complete.cases(mi_data), ]
  }

  mi_data$Feature <- make.names(mi_data$Feature)

  if (nrow(mi_data) < max_features) {
    warning("Number of features in mi_data is less than max_features. Adjusting max_features.")
    max_features <- nrow(mi_data)
  }

  if (max_features == 0) {
    stop("No valid features found in mi_data.")
  }

  auc_values <- numeric(max_features)
  selected_features <- character(max_features)

  for (i in 1:max_features) {
    current_features <- mi_data$Feature[1:i]
    if (any(is.na(current_features))) {
      stop("Encountered NA in feature list from mi_data.")
    }

    formula_str <- paste("Y ~", paste(current_features, collapse = "+"))
    formula_model <- as.formula(formula_str)

    missing_features <- setdiff(current_features, colnames(X))
    if (length(missing_features) > 0) {
      stop("The following features are not found in X: ", paste(missing_features, collapse=", "))
    }

    stepwise_model <- lm(formula_model, data = data.frame(Y = Y, X))
    predictions <- predict(stepwise_model, newdata = data.frame(Y = Y, X), type = "response")
    roc_auc <- pROC::auc(pROC::roc(Y, predictions))

    auc_values[i] <- roc_auc
    selected_features[i] <- paste(names(coef(stepwise_model))[-1], collapse = ", ")
  }

  return(list(SelectedFeatures = selected_features[which.max(auc_values)],
              AUC = auc_values,
              Method = "Mutual Information"))
}


#' Feature Selection Based on Maximal Information Coefficient and AUC
#'
#' This function selects predictive features by ranking them using the Maximal Information Coefficient (MIC),
#' and evaluates model performance with incrementally added features using AUC (Area Under the ROC Curve).
#' It returns the optimal feature subset that yields the highest AUC when used in a linear model.
#'
#' @importFrom pROC roc auc
#' @importFrom minerva mine
#' @param X A data frame or matrix of predictor variables.
#' @param Y A numeric vector or factor representing the binary response variable (0/1).
#' @param max_features An integer specifying the maximum number of top-ranked features to evaluate.
#'
#' @returns A list containing:
#' \describe{
#'   \item{\code{SelectedFeatures}}{The feature names that yield the highest AUC when used in a linear model.}
#'   \item{\code{AUC}}{A numeric vector of AUC values corresponding to each subset of top features.}
#'   \item{\code{Method}}{The method used for feature ranking, i.e., "Maximal Information Coefficient".}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(minerva)
#' library(pROC)
#' X <- data.frame(x1 = rnorm(100), x2 = rnorm(100), x3 = runif(100))
#' Y <- sample(0:1, 100, replace = TRUE)
#' result <- feature_selection_mic_stepwise(X, Y, max_features = 3)
#' print(result)
#' }
feature_selection_mic_stepwise <- function(X, Y, max_features) {
  colnames(X) <- make.names(colnames(X))

  mic_data <- mic_information(X, Y)
  mic_data <- mic_data[order(-mic_data$`Maximal Information Coefficient`), ]

  if (any(is.na(mic_data$Feature))) {
    warning("NA values found in mic_data$Feature. Removing these rows.")
    mic_data <- mic_data[complete.cases(mic_data), ]
  }

  mic_data$Feature <- make.names(mic_data$Feature)

  if (nrow(mic_data) < max_features) {
    warning("Number of features in mic_data is less than max_features. Adjusting max_features.")
    max_features <- nrow(mic_data)
  }

  if (max_features == 0) {
    stop("No valid features found in mic_data.")
  }

  auc_values <- numeric(max_features)
  selected_features <- character(max_features)

  for (i in 1:max_features) {
    current_features <- mic_data$Feature[1:i]
    if (any(is.na(current_features))) {
      stop("Encountered NA in feature list from mic_data.")
    }

    formula_str <- paste("Y ~", paste(current_features, collapse = "+"))
    formula_model <- as.formula(formula_str)

    missing_features <- setdiff(current_features, colnames(X))
    if (length(missing_features) > 0) {
      stop("The following features are not found in X: ", paste(missing_features, collapse=", "))
    }

    stepwise_model <- lm(formula_model, data = data.frame(Y = Y, X))
    predictions <- predict(stepwise_model, newdata = data.frame(Y = Y, X), type = "response")
    roc_auc <- pROC::auc(pROC::roc(Y, predictions))

    auc_values[i] <- roc_auc
    selected_features[i] <- paste(names(coef(stepwise_model))[-1], collapse = ", ")
  }

  return(list(SelectedFeatures = selected_features[which.max(auc_values)],
              AUC = auc_values,
              Method = "Maximal Information Coefficient"))
}

#' Plot AUC vs Number of Features
#'
#' Creates a visualization of model performance (AUC) against the number of features
#' used, highlighting the optimal number of features identified during feature selection.
#' The plot compares multiple feature selection methods (IV, MIC, MI) and can save
#' the output to PDF.
#'
#' @param combined_plot_data A data.frame containing plotting data with columns:
#'   \itemize{
#'     \item Method: Feature selection method (IV/MIC/MI)
#'     \item Num_Features: Number of features used
#'     \item AUC: Corresponding AUC value
#'     \item dAUC: AUC change from previous feature count (optional)
#'   }
#' @param best_features_num Integer specifying the optimal number of features identified
#' @param palette_name Character name of color palette (default: "AsteroidCity1").
#'   Must be a valid palette from `wesanderson` package.
#' @param base_size Base font size for plot elements (default: 14)
#' @param save_plots Logical indicating whether to save plot to file (default: TRUE)
#' @param save_dir Directory path to save plot (default: here::here("ModelData", "sel_feature"))
#' @param plot_width Plot width in inches (default: 5)
#' @param plot_height Plot height in inches (default: 5)
#'
#' @return A ggplot object containing the performance comparison plot. When `save_plots=TRUE`,
#'   also saves a PDF file to the specified directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After running SelFeatureSet
#' plot_data <- object@feature.selection$plot_auc_features$data
#' best_n <- object@feature.selection$best_features$best_features_num
#' 
#' plot_num_auc(
#'   combined_plot_data = plot_data,
#'   best_features_num = best_n,
#'   palette_name = "Royal1",
#'   save_dir = "results/feature_selection"
#' )
#' }
plot_num_auc <- function(combined_plot_data,
                         best_features_num,
                         palette_name = "AsteroidCity1",
                         base_size = 14,
                         save_plots = TRUE,
                         save_dir = here::here("ModelData", "sel_feature"),
                         plot_width = 5,
                         plot_height = 5) {
  
  
  
  combined_plot <- ggplot(combined_plot_data, aes(x = Num_Features, y = AUC, color = Method)) +
    geom_line(size = 1.2, alpha = 0.8) +
    geom_point(size = 2.5, shape = 19) +
    geom_vline(xintercept = best_features_num, 
               linetype = "dashed", 
               color = "black", 
               show.legend = FALSE)  +
    scale_color_manual(values = wes_palette(palette_name)) +
    labs(
      title = "AUC vs Number of Features",
      subtitle = paste("Optimal features:", best_features_num),
      x = "Number of Features",
      y = "AUC",
      color = "Method"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.8)),
      legend.title = element_text(face = "bold", size = base_size * 0.7),
      legend.text = element_text(size = base_size * 0.6),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )
  
  print(combined_plot)
  if (save_plots) {
    if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    ggplot2::ggsave(
      filename = file.path(save_dir, "Feature_Selection_Performance.pdf"),
      plot = combined_plot,
      width = plot_width,
      height = plot_height,
      device = "pdf"
    )
    cat("Plot saved to:", file.path(save_dir, "Feature_Selection_Performance.pdf"), "\n")
  }
  
  return(combined_plot)
}


#' Feature Selection and Performance Evaluation
#'
#' Performs feature selection using three different methods (Information Value, Maximal 
#' Information Coefficient, and Mutual Information) and evaluates their performance based 
#' on AUC metrics. Identifies the optimal number of features and generates visualization.
#'
#' @param object A `Train_Model` object containing training data or a list with 
#'   `features` and `label` components. If NULL, `features` and `label` must be provided.
#' @param max_features Maximum number of features to consider (default: NULL, uses all features).
#' @param xintercept X-axis intercept for plot (default: NULL, uses max_features).
#' @param AUC_change_threshold Minimum AUC improvement threshold for feature selection 
#'   (default: 0.01).
#' @param group_col Name of the response variable column (default: "group").
#' @param palette_name Name of color palette for plotting (default: "AsteroidCity1").
#' @param data_type Type of data to use: "clean" for raw or "scale" for standardized 
#'   (default: "clean").
#' @param save_plots Logical indicating whether to save plots (default: TRUE).
#' @param save_dir Directory to save results (default: here::here("ModelData", "sel_feature")).
#' @param plot_width Plot width in inches (default: 5).
#' @param plot_height Plot height in inches (default: 5).
#' @param base_size Base font size for plots (default: 14).
#' @param seed Random seed for reproducibility (default: 123).
#'
#' @return Depending on input type:
#'   - For `Train_Model` object: Returns updated object with feature selection results stored in:
#'     * `@feature.selection`: Detailed selection results
#'     * `@feature.result`: Best feature subset
#'   - For other inputs: Returns list with:
#'     * `plot`: ggplot object of performance comparison
#'     * `feature_data`: Combined performance data for all methods
#'
#' @details
#' The function performs the following steps:
#' 1. Validates input data and parameters
#' 2. Runs three feature selection methods (IV, MIC, MI)
#' 3. Calculates AUC performance for each method across feature counts
#' 4. Identifies optimal feature set based on AUC improvement threshold
#' 5. Generates performance comparison plot
#' 6. Updates input object or returns results
#'
#' @section Output Files:
#' - `feature_selection_results.csv`: CSV file with performance metrics
#' - `Combined_AUC_vs_Features.pdf`: Plot of AUC vs feature count (if save_plots=TRUE)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # With Train_Model object
#' # Custom parameters
#' object_model<-SelFeatureSet(
#'   object = object_model,
#'   max_features = 50,
#'   AUC_change_threshold = 0.005,
#'   palette_name = "Royal1"
#' )
#' }
SelFeatureSet <- function(
    object = NULL,
    max_features = NULL,
    xintercept = NULL,
    AUC_change_threshold = 0.01,
    group_col = "group",
    palette_name = "AsteroidCity1",
    data_type = "clean",
    save_plots = TRUE,
    save_dir = here::here("ModelData", "sel_feature"),
    plot_width = 5,
    plot_height = 5,
    base_size = 14,
    seed = 123,
    csv_filename="feature_selection_results.csv",
    save_data =TRUE) {
  
  set.seed(seed)
  
  if (missing(object) && (is.null(features) || is.null(label))) {
    stop("Either 'object' must be provided or both 'features' and 'label' must be provided.")
  }
  
  cat("Validating input...\n")
  if (!is.null(object)) {
    if (inherits(object, 'Train_Model')) {
      group_col <- slot(object, "group_col")
      
      cat("'object' is of type 'Train_Model'. Extracting training data...\n")
      if (data_type == "scale") {
        train <- slot(object, "split.scale.data")$training
      } else if (data_type == "clean") {
        train <- slot(object, "split.data")[["training"]]
      } else {
        stop("Invalid 'data_type'. Use 'clean' or 'scale'.")
      }
      
      if (!inherits(train, "data.frame")) {
        stop("'train' slot in 'object' must be a data frame.")
      }
      if (ncol(train) < 2) {
        stop("'train' data frame must have at least one feature column and one response column.")
      }
      
      X_train <- train[, -ncol(train)]
      y_train <- train[[group_col]]
      if (is.null(y_train)) {
        stop(sprintf("Response column '%s' not found in 'train' data frame.", group_col))
      }
      
      cat("Extracted training features and response variable.\n")
      
      X_train <- as.data.frame(lapply(X_train, function(col) as.numeric(as.character(col))))
      y_train <- as.numeric(as.character(y_train))
      
    } else if (is.list(object)) {
      # Handle list input
      cat("'object' is a list.\n")
      if (is.null(features) || is.null(label)) {
        if (!("features" %in% names(object)) || !("label" %in% names(object))) {
          stop("If 'object' is a list, it must contain 'features' and 'label'.")
        }
        X_train <- object$features
        y_train <- object$label
      } else {
        stop("Invalid type for 'object'. Must be 'Train_Model' or a list.")
      }
    } else {
      stop("Invalid type for 'object'. Must be 'Train_Model' or a list.")
    }
  }
  
  if (is.null(max_features)) {
    max_features <- ncol(X_train)
  }
  
  if (is.null(xintercept)) {
    xintercept <- max_features
  }
  
  X_train <- as.data.frame(X_train)
  cat("Starting feature selection...\n")
  
  selected_features_iv <- feature_selection_iv(
    X = X_train,
    Y = y_train,
    max_features = max_features
  )
  
  selected_features_mic <- feature_selection_mic_stepwise(
    X = X_train,
    Y = y_train,
    max_features = max_features
  )
  
  selected_features_mi <- feature_selection_mi(
    X = X_train,
    Y = y_train,
    max_features = max_features
  )
  
  cat("Feature selection completed for IV, MIC, and MI methods.\n")
  
  selected_features_iv$SelectedFeatures <- unlist(strsplit(selected_features_iv$SelectedFeatures, ", "))
  selected_features_mic$SelectedFeatures <- unlist(strsplit(selected_features_mic$SelectedFeatures, ", "))
  selected_features_mi$SelectedFeatures <- unlist(strsplit(selected_features_mi$SelectedFeatures, ", "))
  
  selected_features_iv_data <- data.frame(
    Method = "Information Value (IV)",  # Information Value method
    Num_Features = 1:max_features, 
    AUC = selected_features_iv$AUC
  )
  
  selected_features_mic_data <- data.frame(
    Method = "Maximal Information Coefficient (MIC)",  
    Num_Features = 1:max_features, 
    AUC = selected_features_mic$AUC
  )
  
  selected_features_mi_data <- data.frame(
    Method = "Mutual Information (MI)",  
    Num_Features = 1:max_features, 
    AUC = selected_features_mi$AUC
  )
  
  combined_plot_data <- rbind(
    selected_features_iv_data, 
    selected_features_mic_data, 
    selected_features_mi_data
  )
  
  combined_plot_data <- combined_plot_data %>%
    dplyr::group_by(Method) %>%
    dplyr::mutate(
      max_AUC = max(AUC, na.rm = TRUE),
      dAUC = (max_AUC - AUC) / max_AUC 
    ) %>%
    dplyr::ungroup()
  
  cat("AUC changes calculated.\n")
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    cat("Created directory:", save_dir, "\n")
  }

  if (save_data) {
      full_path <- file.path(save_dir, csv_filename)
      write.csv(combined_plot_data, file = full_path, row.names = FALSE)
      cat("Saved feature selection results to:", full_path, "\n")
  }
  change_points <- combined_plot_data %>%
    dplyr::filter(dAUC <= AUC_change_threshold) %>%
    dplyr::group_by(Method) %>%
    dplyr::arrange(Num_Features) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  if (nrow(change_points) == 0) {
    cat("No change points found. Selecting features with maximum AUC.\n")
    change_points <- combined_plot_data %>%
      dplyr::group_by(Method) %>%
      dplyr::filter(AUC == max(AUC)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }
  
  best_features <- change_points %>%
    dplyr::arrange(desc(AUC)) %>%
    dplyr::slice(1)
  
  best_features_num <- best_features$Num_Features
  best_method <- best_features$Method
  
  cat(sprintf("Best features selected. Method: %s, Number of Features: %d\n", 
              best_method, best_features_num))
  
  best_features_subset <- switch(
    best_method,
    "Information Value (IV)" = {
      if (length(selected_features_iv$SelectedFeatures) < best_features_num) {
        stop("Not enough features in 'IV' selection.")
      }
      head(selected_features_iv$SelectedFeatures, best_features_num)
    },
    "Maximal Information Coefficient (MIC)" = {
      if (length(selected_features_mic$SelectedFeatures) < best_features_num) {
        stop("Not enough features in 'MIC' selection.")
      }
      head(selected_features_mic$SelectedFeatures, best_features_num)
    },
    "Mutual Information (MI)" = {
      if (length(selected_features_mi$SelectedFeatures) < best_features_num) {
        warning("Not enough features in 'MI' selection. Reducing the number of features.")
        best_features_num <- length(selected_features_mi$SelectedFeatures)
      }
      head(selected_features_mi$SelectedFeatures, best_features_num)
    }
  )
  
  combined_plot <- plot_num_auc(
    combined_plot_data = combined_plot_data,
    best_features_num = best_features_num,
    palette_name = palette_name,
    base_size = base_size,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height
  )
  
  cat("Plotting completed.\n")
  
  if (inherits(object, 'Train_Model')) {
    if (!is.null(object@feature.selection)) {
      warning("Overwriting existing feature selection results.")
    }
    
    object@feature.selection <- list(
      selected_features_iv = selected_features_iv,
      selected_features_mic = selected_features_mic,
      selected_features_mi = selected_features_mi,
      best_features = list(
        best_features_num = best_features_num,
        best_method = best_method
      ),
      plot_auc_features = combined_plot
    )
    
    object@feature.result <- list(best_features_subset = best_features_subset)
    
    cat("Updating 'Train_Model' object...\n")
    cat("The 'Train_Model' object has been updated with the following slots:\n")
    cat("- 'feature.selection' slot updated.\n")
    cat("- 'feature.result' slot updated.\n")
    
    return(object)
  }
  
  return(list(
    plot = combined_plot,
    feature_data = combined_plot_data
  ))
}
