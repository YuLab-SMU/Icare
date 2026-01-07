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
#'
#' @examples
#' # Example usage of Extract_filtered.set
#' # Assuming 'model_data' is an existing Train_Model object with a 'filtered.set' slot
#' data_sets <- Extract_filtered.set(object = model_data)
#' training_data <- data_sets$training
#' testing_data <- data_sets$testing
Extract_filtered.set <- function(object) {
  train <- tryCatch(object@filtered.set$training,
                    error = function(e) NULL)
  test <- tryCatch(object@filtered.set$testing,
                   error = function(e) NULL)
  return(list(training = train, testing = test))
}


#' Train and Evaluate Multiple Machine Learning Models
#'
#' This function trains and evaluates multiple machine learning models using caret,
#' with automatic selection of cross-validation strategy based on sample size.
#' Supports both Leave-One-Out Cross-Validation (LOOCV) for small datasets and
#' K-fold cross-validation for larger datasets.
#'
#' @param data A data frame containing both features and the target variable
#' @param methods Character vector of caret model methods to train
#' @param control List of control parameters for trainControl (method, number, repeats)
#' @param tune_grids List of tuning grids for each model method
#' @param classProbs Logical indicating whether class probabilities should be computed
#' @param allowParallel Logical indicating whether to allow parallel processing
#' @param group_col Name of the target variable column (default: "group")
#' @param loocv_threshold Sample size threshold below which LOOCV is used (default: 100)
#'
#' @return A list of trained caret models, one for each method that was successfully trained
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data <- data.frame(
#'   group = factor(c(rep("Class1", 50), rep("Class2", 50))),
#'   feature1 = rnorm(100),
#'   feature2 = rnorm(100)
#' )
#'
#' control <- list(method = "repeatedcv", number = 10, repeats = 3)
#' methods <- c("rf", "glm", "svmRadial")
#' tune_grids <- list(
#'   rf = data.frame(mtry = c(2, 5)),
#'   svmRadial = data.frame(sigma = 0.1, C = 1)
#' )
#'
#' results <- train_and_evaluate_models(
#'   data = data,
#'   methods = methods,
#'   control = control,
#'   tune_grids = tune_grids
#' )
#' }
#' @importFrom doParallel registerDoParallel
#' @import parallel
#' @importFrom foreach registerDoSEQ
train_and_evaluate_models <- function(data,
                                      methods,
                                      control,
                                      tune_grids,
                                      classProbs = TRUE,
                                      allowParallel = TRUE,
                                      group_col = "group",
                                      loocv_threshold = 100) {
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)

  data[[group_col]] <- factor(data[[group_col]])
  levels(data[[group_col]]) <- make.names(levels(data[[group_col]]))

  n_samples <- nrow(data)
  if (n_samples <= loocv_threshold) {
    fitControl <- trainControl(
      method = "LOOCV",
      allowParallel = allowParallel,
      classProbs = classProbs,
      summaryFunction = twoClassSummary,
      savePredictions = "final"
    )
    cat(paste("Using LOOCV (n =", n_samples, "samples)"))
  } else {
    fitControl <- trainControl(
      method = control$method,
      number = control$number,
      repeats = control$repeats,
      allowParallel = allowParallel,
      classProbs = classProbs,
      summaryFunction = twoClassSummary,
      savePredictions = "final"
    )
    cat(paste("Using", control$method, "with", control$number, "folds and",
                  control$repeats, "repeats (n =", n_samples, "samples)"))
  }

  results <- list()

  for (method in methods) {
    if (method %in% caret::modelLookup()$model) {
      tune_grid <- tune_grids[[method]]

      model <- tryCatch({
        if (method == "glm") {
          train(
            as.formula(paste(group_col, "~ .")),
            data = data,
            method = method,
            family = "binomial",
            trControl = fitControl,
            metric = "ROC"
          )
        } else {
          train(
            as.formula(paste(group_col, "~ .")),
            data = data,
            method = method,
            trControl = fitControl,
            tuneGrid = tune_grid,
            metric = "ROC"
          )
        }
      }, error = function(e) {
        warning(paste("Error training model", method, ":", e$message))
        return(NULL)
      })

      if (!is.null(model)) {
        results[[method]] <- model
      }
    } else {
      warning(paste("Model", method, "is not in caret's built-in library"))
    }
  }

  stopCluster(cl)
  registerDoSEQ()

  return(results)
}




#' Evaluate Model Performance Metrics
#'
#' This function evaluates the performance of one or more classification models using
#' various metrics including sensitivity, specificity, accuracy, AUC, etc. It supports
#' both single model evaluation and multiple model comparison.
#'
#' @param data A data frame containing the test data with both features and target variable
#' @param model_result Either a single trained model object or a list of trained models
#' @param group_col Name of the target variable column (default: "group")
#' @param use_youden Logical indicating whether to use Youden's J statistic to determine optimal cutoff (default: FALSE)
#' @param custom_cutoff Optional custom probability cutoff (overrides use_youden if provided)
#'
#' @return A data frame containing performance metrics for each model. For multiple models,
#'         returns a combined data frame with one row per model.
#'
#' @examples
#' \dontrun{
#' # Example with single model
#' model <- train(...)  # trained model
#' test_data <- data.frame(...)  # test data
#' performance <- evaluate_model_performance(
#'   data = test_data,
#'   model_result = model,
#'   use_youden = TRUE
#' )
#'
#' # Example with multiple models
#' model_list <- list(
#'   model1 = train(...),
#'   model2 = train(...)
#' )
#' performance <- evaluate_model_performance(
#'   data = test_data,
#'   model_result = model_list,
#'   custom_cutoff = 0.4
#' )
#' }
#' @export
#' @import stats
#' @import caret
evaluate_model_performance<- function(data, model_result,
                                       group_col = "group",
                                       custom_cutoff = NULL) {
  results <- list()
  set.seed(1234)



  if (!("method" %in% names(model_result))) {
    for (model_name in names(model_result)) {
      model <- model_result[[model_name]]

      pred_prob <- predict(model, newdata = data, type = "prob")[, 2]

   if (!is.null(custom_cutoff)) {
        cutoff <- custom_cutoff
      } else {
        cutoff <- 0.5
      }

      pred <- ifelse(pred_prob > cutoff, 1, 0)
      cf <- table(factor(data[[group_col]], levels = c(0,1)),
                  factor(pred, levels = c(0,1)))

      if (dim(cf)[1] == 2 && dim(cf)[2] == 2) {
        TN <- cf[1, 1]
        FP <- cf[1, 2]
        FN <- cf[2, 1]
        TP <- cf[2, 2]

        Sensitivity <- TP / (TP + FN)
        Specificity <- TN / (TN + FP)
        Positive_predictive_value <- TP / (TP + FP)
        Negative_predictive_value <- TN / (TN + FN)
        accuracy_score <- sum(diag(cf)) / sum(cf)
        Precision <- Positive_predictive_value
        recall_score <- Sensitivity
        f1_score <- 2 * (Precision * recall_score) / (Precision + recall_score)

        roc_data <- pROC::roc(data[[group_col]], pred_prob)
        auc_value <- pROC::auc(roc_data)

        result <- data.frame(
          Model = model_name,
          Cutoff_Used = cutoff,
          Sensitivity = Sensitivity,
          Specificity = Specificity,
          Positive_predictive_value = Positive_predictive_value,
          Negative_predictive_value = Negative_predictive_value,
          accuracy_score = accuracy_score,
          Precision = Precision,
          recall_score = recall_score,
          f1_score = f1_score,
          auc = auc_value
        )

        results[[model_name]] <- result
      } else {
        cat(paste("Model", model_name, "produced invalid confusion matrix. Skipping."))
      }
    }
    return(do.call(rbind, results))

  } else {
    model <- model_result
    pred_prob <- predict(model, newdata = data, type = "prob")[, 2]

   if (!is.null(custom_cutoff)) {
      cutoff <- custom_cutoff
    } else {
      cutoff <- 0.5
    }

    pred <- ifelse(pred_prob > cutoff, 1, 0)
    cf <- table(factor(data[[group_col]], levels = c(0,1)),
                factor(pred, levels = c(0,1)))

    if (dim(cf)[1] == 2 && dim(cf)[2] == 2) {
      TN <- cf[1, 1]
      FP <- cf[1, 2]
      FN <- cf[2, 1]
      TP <- cf[2, 2]

      Sensitivity <- TP / (TP + FN)
      Specificity <- TN / (TN + FP)
      Positive_predictive_value <- TP / (TP + FP)
      Negative_predictive_value <- TN / (TN + FN)
      accuracy_score <- sum(diag(cf)) / sum(cf)
      Precision <- Positive_predictive_value
      recall_score <- Sensitivity
      f1_score <- 2 * (Precision * recall_score) / (Precision + recall_score)

      roc_data <- pROC::roc(data[[group_col]], pred_prob)
      auc_value <- pROC::auc(roc_data)

      return(data.frame(
        Model = model$method,
        Cutoff_Used = cutoff,
        Sensitivity = Sensitivity,
        Specificity = Specificity,
        Positive_predictive_value = Positive_predictive_value,
        Negative_predictive_value = Negative_predictive_value,
        accuracy_score = accuracy_score,
        Precision = Precision,
        recall_score = recall_score,
        f1_score = f1_score,
        auc = auc_value
      ))
    } else {
      cat("Invalid confusion matrix produced by the model.")
      return(NULL)
    }
  }
}


#' Plot ROC Curves for Multiple Models
#'
#' Generates ROC curves for multiple models and saves the plotting data in CSV format.
#'
#' @param model_list List of models to evaluate
#' @param validation_data Data frame containing validation data
#' @param group_col Name of the response variable column (default: "group")
#' @param palette_name Name of color palette (default: "AsteroidCity1")
#' @param base_size Base font size (default: 14)
#' @param save_plots Logical to save plots (default: TRUE)
#' @param save_data  Logical to save data (default: TRUE)
#' @param save_dir Directory to save outputs (default: here::here("ModelData", "best.model.result"))
#' @param plot_width Plot width in inches (default: 5)
#' @param plot_height Plot height in inches (default: 5)
#' @param alpha Transparency level (default: 1)
#'
#' @return List of ROC objects and saves CSV data if requested
#' @export
#'
#' @examples
#' \dontrun{
#' roc_results <- plot_roc_curve(model_list = my_models,
#'                              validation_data = test_data)
#' }
plot_roc_curve <- function(model_list,
                           validation_data,
                           group_col = "group",
                           palette_name = "AsteroidCity1",
                           base_size = 14,
                           save_plots = TRUE,
                           save_dir = here::here("ModelData", "best.model.result"),
                           plot_width = 5,
                           plot_height = 5,
                           alpha = 1,
                           save_data =TRUE) {

  roc_list <- list()
  plot_list <- list()
  auc_results <- numeric()

  for (model_name in names(model_list)) {
    model <- model_list[[model_name]]
    predictions <- predict(model, validation_data, type = "prob")[, 2]
    roc_obj <- roc(validation_data[[group_col]], predictions, levels = c("0", "1"), direction = "<")
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
    cat("Failed to use specified palette, falling back to viridis default colors.")
    viridis::viridis(n_colors)
  })

  # Create plot
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
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.8)),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90")
    )

  print(p)

  if (save_data) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

    ggsave(filename = file.path(save_dir, "roc_curves.pdf"),
           plot = p,
           width = plot_width,
           height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "roc_curves.pdf"), "\n")

  }
  
  if(save_data){
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
#' This function checks if the specified column in the provided dataset has at least two distinct levels. It is used to ensure that the validation data contains more than one class, which is a requirement for classification tasks.
#'
#' @param data A data frame or tibble containing the dataset.
#' @param group_col The name of the column to check, which should represent class labels.
#'
#' @returns TRUE if the column contains at least two levels, otherwise stops with an error message.
#' @export
#'
#' @examples
#' # Check if the 'group' column in the dataset 'my_data' contains at least two levels
#' check_factor_level(data = my_data, group_col = "group")
#'
#' @details
#' This function ensures that the column representing class labels (`group_col`) has at least two levels, which is essential for classification problems. If the column has fewer than two levels, the function stops and throws an error indicating that the validation data is not suitable for classification.
check_factor_level <- function(data, group_col) {
  levels_present <- levels(as.factor(data[[group_col]]))

  if (length(levels_present) < 2) {
    stop("Validation data must contain at least two class levels.")
  }

  return(TRUE)
}


#' Comprehensive Model Training and Analysis Pipeline
#'
#' This function performs a complete machine learning workflow including:
#' 1. Data preparation and validation
#' 2. Model training with cross-validation
#' 3. Performance evaluation on training data
#' 4. ROC curve generation
#' 5. Results organization and visualization
#'
#' @param object Either a Train_Model object or a list containing training and test data
#' @param methods Character vector of caret model methods to train (default: common classifiers)
#' @param control List of control parameters for trainControl (default: 10-fold repeated CV)
#' @param tune_grids List of tuning grids for each model method
#' @param loocv_threshold Sample size threshold below which LOOCV is used (default: 100)
#' @param classProbs Logical indicating whether class probabilities should be computed
#' @param allowParallel Logical indicating whether to allow parallel processing
#' @param group_col Name of the target variable column (default: "group")
#' @param palette_name Name of color palette for ROC curves (default: "AsteroidCity1")
#' @param base_size Base font size for plots (default: 14)
#' @param save_plots Logical indicating whether to save plots (default: TRUE)
#' @param save_dir Directory path to save plots (default: here("ModelData", "best_model_result"))
#' @param plot_width Width of saved plots in inches (default: 5)
#' @param plot_height Height of saved plots in inches (default: 5)
#' @param seed Random seed for reproducibility (default: 123)
#' @param alpha Transparency level for ROC curves (default: 0.8)
#' @param use_youden Logical indicating whether to use Youden's J statistic for cutoff
#' @param custom_cutoff Optional custom probability cutoff (overrides use_youden if provided)
#'
#' @return If input is Train_Model object, returns updated object with results in slots.
#'         If input is list, returns list with performance metrics and ROC data.
#'
#' @examples
#' \dontrun{
#'
#' # Example with list input
#' object_model <- ModelTrainAnalysis(
#'   object = object_model,
#'   methods = c("gbm", "xgbTree"),
#'   palette_name = "Royal1"
#' )
#' }
#' @export
#' @import caret
#' @import ggplot2
#' @importFrom wesanderson wes_palette
#' @importFrom viridis viridis
#' @import here
#' @import stats
#' @importFrom doParallel registerDoParallel
#' @import parallel
#' @importFrom foreach registerDoSEQ
ModelTrainAnalysis <- function(object,
                               methods = c("glm", "rpart", "naive_bayes", "bayesglm", "rf",
                                           "xgbTree", "svmRadial", "svmLinear", "gbm", "earth", "glmnet"),
                               control = list(method = "repeatedcv", number = 10, repeats = 5),
                               tune_grids = list(
                                 glm = NULL,
                                 rpart = expand.grid(cp = seq(0.0001, 0.01, length.out = 10)),
                                 naive_bayes = NULL,
                                 bayesglm = NULL,
                                 rf = expand.grid(mtry = 1:5),
                                 xgbTree = expand.grid(
                                   nrounds = 100,
                                   max_depth = c(2, 4, 6),
                                   eta = c(0.01, 0.1),
                                   gamma = 0,
                                   colsample_bytree = 1,
                                   min_child_weight = 1,
                                   subsample = 1
                                 ),
                                 svmRadial = expand.grid(sigma = 0.01, C = 2^(-1:2)),
                                 svmLinear = expand.grid(C = c(0.01, 0.1, 1)),
                                 gbm = expand.grid(
                                   n.trees = c(50, 100),
                                   interaction.depth = c(2, 3),
                                   shrinkage = c(0.001, 0.01),
                                   n.minobsinnode = c(10, 20)
                                 ),
                                 earth = expand.grid(degree = 1:2, nprune = 2:10),
                                 glmnet = expand.grid(
                                   alpha = c(0.1, 0.5, 0.9),
                                   lambda = 10^seq(-4, -1, 1)
                                 )
                               ),
                               loocv_threshold = 100,
                               classProbs = TRUE,
                               allowParallel = FALSE,
                               group_col = "group",
                               palette_name = "AsteroidCity1",
                               base_size = 14,
                               save_plots = TRUE,
                               save_dir = here::here("ModelData", "best_model_result"),
                               plot_width = 5,
                               plot_height = 5,
                               seed = 123,
                               alpha = 0.8,
                               custom_cutoff = NULL) {

  set.seed(seed)

  if (inherits(object, "Train_Model")) {
    cat("Input is of class 'Train_Model'. Extracting datasets...\n")
    data_sets <- Extract_filtered.set(object)
    train_data <- data_sets$training
    test_data <- data_sets$testing
    group_col <- object@group_col
  } else if (is.list(object) && all(c("train", "test") %in% names(object))) {
    cat("Input is a list with 'train' and 'test' elements.\n")
    train_data <- object$train
    test_data <- object$test
  } else {
    stop("Input must be an object of class 'Train_Model' or a list with 'train' and 'test' elements")
  }

  if (nrow(train_data) < loocv_threshold) {
    cat(paste0("Sample size <", loocv_threshold, ". Using LOOCV strategy.\n"))
    control <- list(method = "LOOCV")
  }

  cat("Data extracted. Checking factor levels in the training data...\n")
  check_factor_level(data = train_data, group_col = group_col)

  cat("Training and evaluating models...\n")
  model_list <- train_and_evaluate_models(
    data = train_data,
    methods = methods,
    control = control,
    tune_grids = tune_grids,
    classProbs = classProbs,
    allowParallel = allowParallel,
    group_col = group_col,
    loocv_threshold=loocv_threshold
  )

  cat("Evaluating models on the training dataset...\n")
  train_performance <- evaluate_model_performance(
    data = train_data,
    model_result = model_list,
    group_col = group_col,
    custom_cutoff = custom_cutoff
  )

  cat("Sorting results by accuracy score...\n")
  train_performance <- train_performance[order(train_performance$accuracy_score, decreasing = TRUE), ]

  cat("Generating ROC curves on training data...\n")
  roc_list <- plot_roc_curve(
    model_list = model_list,
    validation_data = train_data,
    group_col = group_col,
    palette_name = palette_name,
    base_size = base_size,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    alpha = alpha
  )

  if (inherits(object, "Train_Model")) {
    object@all.results <- train_performance
    object@train.models <- model_list
    cat("Updating 'Train_Model' object...\n")
    return(object)
  } else {
    cat("Returning results as a list...\n")
    return(list(train_performance = train_performance, roc_list = roc_list))
  }
}

#' Plot ROC Curves for Best Model Across Multiple Datasets
#'
#' Generates and optionally saves ROC curves comparing model performance across:
#' - Training set
#' - Testing set
#' - Validation set (optional)
#' - External validation set (optional)
#' Includes AUC values with confidence intervals in legend labels.
#'
#' @param best_model A trained classification model object with predict method
#' @param train_data Data frame containing training data with response variable
#' @param test_data Data frame containing testing data with response variable
#' @param validation_data Data frame containing validation data (optional)
#' @param external_validation Data frame containing external validation data (optional)
#' @param group_col Character string specifying the response variable column name (default: "group")
#' @param palette_name Character string specifying Wes Anderson palette name (default: "AsteroidCity1")
#' @param base_size Numeric base font size for plot (default: 14)
#' @param save_plots Logical indicating whether to save plots (default: TRUE)
#' @param save_dir Character string specifying directory to save outputs (default: here("ModelData", "best_model_result"))
#' @param plot_width Numeric plot width in inches (default: 5)
#' @param plot_height Numeric plot height in inches (default: 5)
#' @param alpha Numeric transparency value for ROC curves (default: 1)
#' @param subtitle Character string for plot subtitle (default: "Training and Testing Datasets")
#' @param save_data Logical indicating whether to save ROC curve data (default: TRUE)
#' @param csv_filename Character string for saved data filename (default: "best_model_roc_data.csv")
#'
#' @return A list of data frames containing ROC curve coordinates and AUC information for each dataset:
#' \itemize{
#'   \item training - Data frame with training set ROC coordinates
#'   \item testing - Data frame with testing set ROC coordinates
#'   \item validation - Data frame with validation set ROC coordinates (if provided)
#'   \item external_validation - Data frame with external validation ROC coordinates (if provided)
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After training a model (e.g., randomForest)
#' model <- randomForest(group ~ ., data = train_data)
#'
#' # Basic usage with training and test sets
#' roc_data <- plot_best_model_roc(
#'   best_model = model,
#'   train_data = train_data,
#'   test_data = test_data
#' )
#'
#' # Full usage with all datasets
#' roc_data <- plot_best_model_roc(
#'   best_model = model,
#'   train_data = train_data,
#'   test_data = test_data,
#'   validation_data = valid_data,
#'   external_validation = ext_data,
#'   palette_name = "GrandBudapest1",
#'   save_dir = "results/roc_plots"
#' )
#' }
#'
plot_best_model_roc <- function(best_model,
                                train_data ,
                                test_data ,
                                validation_data ,
                                external_validation,
                                group_col = "group",
                                palette_name = "AsteroidCity1",
                                base_size = 14,
                                save_plots = TRUE,
                                save_dir = here("ModelData", "best_model_result"),
                                plot_width = 5,
                                plot_height = 5,
                                alpha = 1,
                                subtitle = "Training and Testing Datasets",
                                save_data = TRUE,
                                csv_filename = "best_model_roc_data.csv") {

  plot_data_list <- list()

  if (!is.null(train_data)) {
    training_predictions <- predict(best_model, newdata = train_data, type = "prob")[, 2]
    roc_training <- roc(train_data[[group_col]], training_predictions, levels = c("0", "1"), direction = "<")

    auc_training <- auc(roc_training)
    auc_ci_training <- ci.auc(roc_training)
    training_plot_data <- data.frame(
      Specificity = 1 - roc_training$specificities,
      Sensitivity = roc_training$sensitivities,
      Dataset = paste0("Training Set (AUC = ", sprintf("%.3f", auc_training),
                       " ± ", sprintf("%.3f", (auc_ci_training[3]-auc_ci_training[1])/2), ")")
    )

    plot_data_list$training <- training_plot_data
  }

  if (!is.null(test_data)) {
    testing_predictions <- predict(best_model, newdata = test_data, type = "prob")[, 2]
    roc_testing <- roc(test_data[[group_col]], testing_predictions, levels = c("0", "1"), direction = "<")

    auc_testing <- auc(roc_testing)
    auc_ci_testing <- ci.auc(roc_testing)
    testing_plot_data <- data.frame(
      Specificity = 1 - roc_testing$specificities,
      Sensitivity = roc_testing$sensitivities,
      Dataset = paste0("Testing Set (AUC = ", sprintf("%.3f", auc_testing),
                       " ± ", sprintf("%.3f", (auc_ci_testing[3]-auc_ci_testing[1])/2), ")")
    )

    plot_data_list$testing <- testing_plot_data
  }



  if (!is.null(validation_data)) {
    validation_data <- match_factor_levels(validation_data, train_data)
    validation_data[[group_col]] <- factor(validation_data[[group_col]], levels = c("0", "1"))
    validation_predictions <- predict(best_model, newdata = validation_data, type = "prob")[, 2]
    roc_validation <- roc(validation_data[[group_col]], validation_predictions, levels = c("0", "1"), direction = "<")

    if (auc(roc_validation) < 0.5) {
      cat("Warning: Model predictions are inverted. Reversing prediction probabilities.\n")
      validation_predictions <- 1 - validation_predictions
      roc_validation <- roc(validation_data[[group_col]], validation_predictions, levels = c("0", "1"), direction = ">")
    }

    auc_validation <- auc(roc_validation)
    auc_ci_validation <- ci.auc(roc_validation)
    validation_plot_data <- data.frame(
      Specificity = 1 - roc_validation$specificities,
      Sensitivity = roc_validation$sensitivities,
      Dataset = paste0("Validation Set (AUC = ", sprintf("%.3f", auc_validation),
                       " ± ", sprintf("%.3f", (auc_ci_validation[3] - auc_ci_validation[1]) / 2), ")")
    )
    plot_data_list$validation <- validation_plot_data
  }

  if (!is.null(external_validation)) {
    external_validation <- na.omit(external_validation)

    external_validation <- match_factor_levels(external_validation, train_data)

    external_validation[[group_col]] <- factor(external_validation[[group_col]], levels = c("0", "1"))
    external_validation_predictions <- predict(best_model, newdata = external_validation, type = "prob")[, 2]
    roc_external_validation <- roc(external_validation[[group_col]], external_validation_predictions, levels = c("0", "1"), direction = "<")

    if (auc(roc_external_validation) < 0.5) {
      cat("Warning: Model predictions are inverted. Reversing prediction probabilities.\n")
      external_validation_predictions <- 1 - external_validation_predictions
      roc_external_validation <- roc(external_validation[[group_col]], external_validation_predictions, levels = c("0", "1"), direction = ">")
    }

    auc_external_validation <- auc(roc_external_validation)
    auc_ci_external_validation <- ci.auc(roc_external_validation)
    external_validation_plot_data <- data.frame(
      Specificity = 1 - roc_external_validation$specificities,
      Sensitivity = roc_external_validation$sensitivities,
      Dataset = paste0("External Validation Set (AUC = ", sprintf("%.3f", auc_external_validation),
                       " ± ", sprintf("%.3f", (auc_ci_external_validation[3] - auc_ci_external_validation[1]) / 2), ")")
    )
    plot_data_list$external_validation <- external_validation_plot_data
  }


  combined_plot_data <- do.call(rbind, plot_data_list)
  combined_plot_data$Dataset <- factor(combined_plot_data$Dataset,
                                       levels = unlist(lapply(plot_data_list, function(x) unique(x$Dataset))))


  palette_colors <- wes_palette(name = palette_name, n = length(unique(combined_plot_data$Dataset)), type = "discrete")

  p <- ggplot(combined_plot_data, aes(x = Specificity, y = Sensitivity, color = Dataset)) +
    geom_line(size = 1.25, alpha = 0.8) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = palette_colors) +
    labs(
      title = "ROC Curves Validation Comparison",
      x = "1 - Specificity",
      y = "Sensitivity",
      color = "Validation Cohort"
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1),
      expand = expansion(mult = 0.01)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1),
      expand = expansion(mult = 0.01)
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = scales::alpha("white", 0.8)),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  print(p)

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    ggsave(filename = file.path(save_dir, "best_model_roc_plot.pdf"),
           plot = p,
           width = plot_width,
           height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "best_model_roc_plot.pdf"), "\n")
  }

  if (save_data) {
    data_path <- file.path(save_dir, csv_filename)
    write.csv(combined_plot_data, data_path, row.names = FALSE)
    cat("ROC curve data saved to: ", data_path, "\n")
  }

  return(plot_data_list)
}


#' Select Best Performing Model from Training Results
#'
#' This function selects the best performing model from a Train_Model object based on
#' AUC (Area Under the Curve) metric by default. It also supports custom model selection.
#'
#' @param object A Train_Model object containing training results
#' @param metric Evaluation metric to use (default: "auc")
#' @param custom_selection Name of specific model to select (bypasses metric selection)
#'
#' @return Returns the input Train_Model object with best model information updated in
#'         the best.model.result slot
#'
#' @examples
#' \dontrun{
#' # Automatic selection based on best AUC
#' object_model <- SelectBestModel(object_model)
#'
#' # Manually select specific model
#' object_model <- SelectBestModel(object_model, custom_selection = "Random Forest")
#' }
#' @export
SelectBestModel <- function(object,
                            metric = "auc",
                            custom_selection = NULL) {

  if (!inherits(object, "Train_Model")) {
    stop("Input must be an object of class 'Train_Model'")
  }

  all.results <- slot(object, "all.results")
  model_list <- slot(object, "train.models")

  available_metrics <- c("Sensitivity", "Specificity", "accuracy_score",
                         "Precision", "recall_score", "f1_score", "auc",
                         "Positive_predictive_value",
                         "Negative_predictive_value")
  available_metrics <- available_metrics[available_metrics %in% colnames(all.results)]

  if (!metric %in% colnames(all.results)) {
    stop("Specified metric '", metric, "' not available.\n",
         "Available metrics: ", paste(available_metrics, collapse = ", "))
  }

  if (!is.null(custom_selection)) {
    if (!custom_selection %in% all.results$Model) {
      stop("Model '", custom_selection, "' not found in trained models")
    }

    cat("User selected model: ", custom_selection)
    best_model_type <- custom_selection
    best.model.result <- all.results[all.results$Model == custom_selection, ]
  } else {
    best_index <- which(all.results[[metric]] == max(all.results[[metric]]))

    if (length(best_index) > 1) {
      cat("Multiple models with same ", metric, ", using f1_score as tie-breaker")
      best_index <- which(all.results[[metric]] == max(all.results[[metric]]) &
                            all.results$f1_score == max(all.results$f1_score[best_index]))
    }

    best.model.result <- all.results[best_index[1], ]
    best_model_type <- best.model.result$Model

    cat("Best model (", metric, "): ", best_model_type)
  }

  best_model <- model_list[[best_model_type]]
  if (is.null(best_model)) {
    stop("Model '", best_model_type, "' not found in trained models")
  }

  object@best.model.result <- list(
    model = best_model,
    model_type = best_model_type,
    train_performance = best.model.result
  )

  return(object)
}


