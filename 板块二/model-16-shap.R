#' Calculate SHAP (SHapley Additive exPlanations) Values for Classification Model
#'
#' This function computes SHAP values for a LightGBM classification model, providing
#' interpretable feature importance measures. It handles both model training and
#' SHAP value calculation in one pipeline.
#'
#' @param train_data A data.frame containing the training data with features and target variable.
#' @param best_model A parsnip model specification object (for LightGBM). If NULL, 
#'                  will use default LightGBM classification settings.
#' @param group_col Character. Name of the column containing the target/grouping variable 
#'                 (default: "group").
#' @param bg_sample_size Integer. Size of background sample for SHAP calculation 
#'                      (default: 100). Smaller values speed up computation but may 
#'                      reduce accuracy.
#' @param seed Integer. Random seed for reproducibility (default: 123).
#' @param save_data Logical. Whether to save SHAP results to disk (default: TRUE).
#' @param save_dir Character. Directory path for saving results (default: 
#'                "ModelData/best_model_result/shap_results" using here::here()).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{shap_values} - A shapviz object containing SHAP values and feature data
#'   \item \code{metadata} - A list with:
#'     \itemize{
#'       \item \code{bg_sample_size} - Actual background sample size used
#'       \item \code{non_converging} - Indices of rows that failed to converge (if any)
#'     }
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Prepares the data and recipe for modeling
#' 2. Trains a LightGBM classification model
#' 3. Computes SHAP values using the kernelshap algorithm
#' 4. Optionally saves two matrices to CSV:
#'    - SHAP values matrix (feature contributions per observation)
#'    - Feature values matrix (original feature values)
#'
#' @note
#' - Requires the following packages: recipes, parsnip, workflows, lightgbm, kernelshap, shapviz
#' - For large datasets, consider reducing bg_sample_size for faster computation
#' - Non-converging samples will trigger a warning but not stop execution
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' library(recipes)
#' library(parsnip)
#' 
#' # Prepare sample data
#' data <- iris
#' data$group <- ifelse(data$Species == "setosa", 1, 0)
#' 
#' # Calculate SHAP values
#' shap_results <- calculate_shap(
#'   train_data = data,
#'   group_col = "group",
#'   bg_sample_size = 50
#' )
#' 
#' # Visualize results
#' shapviz::sv_importance(shap_results$shap_values)
#' }
#'
#' @importFrom recipes recipe
#' @importFrom parsnip boost_tree set_engine set_mode
#' @importFrom workflows workflow add_model add_recipe fit
#' @import stats
#' @importFrom kernelshap kernelshap
#' @importFrom shapviz shapviz
#' @importFrom here here
calculate_shap <- function(train_data,
                           best_model,
                           group_col = "group",
                           bg_sample_size = 100,
                           seed = 123,
                           save_data = TRUE,
                           save_dir = here("ModelData", "best_model_result","shap_results")) {
  
  set.seed(seed)
  train_data[[group_col]] <- as.factor(train_data[[group_col]])
  
  cat("Creating recipe for model...\n")
  rec <- recipe(as.formula(paste(group_col, "~ .")), data = train_data)
  
  cat("Setting up LightGBM model...\n")
  lgb_model <- boost_tree() %>%
    set_engine('lightgbm') %>%
    set_mode('classification')
  
  lgb_wflow <- workflow() %>%
    add_model(lgb_model) %>%
    add_recipe(rec)
  
  cat("Fitting the model...\n")
  lgb_fit <- workflows::fit(lgb_wflow, data = train_data)
  
  pred_fun <- function(model, new_data) {
    predict(model, new_data = new_data, type = "prob")$.pred_1
  }
  
  bg_sample <- train_data[sample(1:nrow(train_data), min(100, nrow(train_data))), ]
  
  cat("Calculating SHAP values...\n")
  ks <- kernelshap(pred_fun, 
                   object = lgb_fit, 
                   X = train_data[, -which(names(train_data) == group_col)], 
                   bg_X = bg_sample)
  
  non_converging <- attr(ks, "non_converging")
  if (!is.null(non_converging) && length(non_converging) > 0) {
    warning("Non-converging rows: ", paste(non_converging, collapse = ", "))
  }
  
  shp <- shapviz(ks)
  shap_values_matrix<-shp[["S"]]
  if (save_data) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat("Created directory:", save_dir, "\n")
    }
    shap_values_matrix <- as.data.frame(shp[["S"]])
    feature_value_matrix <- as.data.frame(shp[["X"]])
    
    csv_path1 <- file.path(save_dir, "shap_values_matrix.csv")
    write.csv(shap_values_matrix, csv_path1, row.names = TRUE)
    cat("SHAP values matrix saved to:", csv_path1, "\n")
    
    csv_path2 <- file.path(save_dir, "feature_value_matrix.csv")
    write.csv(feature_value_matrix, csv_path2, row.names = TRUE)
    cat("Feature values matrix saved to:", csv_path2, "\n")
  }
  
  return(list(
    shap_values = shp,
    metadata = list(
      bg_sample_size = bg_sample_size,
      non_converging = non_converging
    )
  ))
}


#' Calculate and Store SHAP Values for Best_Model Object
#'
#' This function computes SHAP (SHapley Additive exPlanations) values for a Best_Model object,
#' storing the results in the object's `shap.result` slot. It serves as a wrapper for the
#' `calculate_shap()` function, specifically designed to work with Best_Model class objects.
#'
#' @param object A Best_Model class object containing trained model and datasets.
#' @param best_model Optional. A parsnip model specification object. If NULL (default),
#'                  will use the model stored in the Best_Model object.
#' @param group_col Character. Name of the target/grouping variable column. If NULL,
#'                 will use the value stored in the Best_Model object (default: "group").
#' @param bg_sample_size Integer. Size of background sample for SHAP calculation
#'                      (default: 100). Smaller values speed up computation but may
#'                      reduce accuracy.
#' @param seed Integer. Random seed for reproducibility (default: 123).
#' @param save_data Logical. Whether to save SHAP results to disk (default: TRUE).
#' @param save_dir Character. Directory path for saving results (default:
#'                "ModelData/best_model_result/shap_results" using here::here()).
#'
#' @return The input Best_Model object with updated `shap.result` slot containing:
#' \itemize{
#'   \item `shap_values`: A shapviz object with SHAP values and feature data
#'   \item `metadata`: List with background sample size and non-converging rows info
#' }
#'
#' @details
#' The function performs the following steps:
#' 1. Extracts training data and model specifications from the Best_Model object
#' 2. Computes SHAP values using `calculate_shap()`
#' 3. Stores results in the object's `shap.result` slot
#' 4. Optionally saves SHAP matrices to CSV files
#'
#' @note
#' Requires the object to be of class Best_Model with the following slots:
#' - `best.model`: The trained model object
#' - `filtered.set`: List containing training and testing datasets
#' - `group_col`: Name of the target variable column
#'
#' @examples
#' \dontrun{
#' # Assuming 'bm' is a valid Best_Model object
#' bm_with_shap <- ModelShapValue(
#'   object = bm,
#'   bg_sample_size = 50,
#'   seed = 42
#' )
#' 
#' # Access results
#' shap_values <- bm_with_shap@shap.result$shap_values
#' }
#'
#' @seealso
#' \code{\link{calculate_shap}} for the underlying SHAP computation function
#'
#' @importFrom here here
#' @export
ModelShapValue <- function(object,
                           best_model,
                           group_col = "group",
                           bg_sample_size = 100,
                           seed = 123,
                           save_data = TRUE,
                           save_dir = here("ModelData", "best_model_result","shap_results")) {
  
  if (inherits(object, "Best_Model")) {
    cat("Input is of class 'Best_Model'. Extracting datasets...\n")
    best_model<-object@best.model[[1]]
    best_model_type<-object@best.model.type
    group_col  <-  object@group_col
    data_sets <- object@filtered.set
    train_data <- data_sets$training
    test_data <- data_sets$testing
  }
  
  cat("Generating SHAP value...\n")
  shap_results <- calculate_shap(train_data=train_data,
                                 best_model=best_model,
                                 group_col = group_col,
                                 bg_sample_size = bg_sample_size,
                                 seed = seed,
                                 save_data =save_data,
                                 save_dir =save_dir)
  
  object@shap.result <- shap_results
  cat("Updating 'Best_Model' object...\n")
  cat("The 'Best_Model' object has been updated with the following slots:\n")
  cat("- 'shap.result ' slot updated.\n")
  
  return(object)
  
}




#' Generate SHAP Dependence Plots for Model Interpretation
#'
#' Creates SHAP (SHapley Additive exPlanations) dependence plots to visualize how individual
#' features affect model predictions. Works with either Best_Model objects or direct
#' shapviz objects.
#'
#' @param object A Best_Model class object containing SHAP results. Alternatively,
#'               can be NULL if shp parameter is provided directly.
#' @param shp Optional. A shapviz object containing pre-computed SHAP values.
#'            Required if object is not a Best_Model.
#' @param features Character vector or numeric. Specific features to plot. Can be:
#'                - NULL (plots top n_top features)
#'                - Numeric (top n features)
#'                - Character vector (specific feature names)
#' @param n_top Integer. Number of top features to plot when features=NULL (default: 3).
#' @param palette_name Character. Name of color palette for plots (default: "AsteroidCity1").
#'                    Falls back to viridis if specified palette not found.
#' @param base_size Numeric. Base font size for plots (default: 14).
#' @param save_plots Logical. Whether to save plots to disk (default: TRUE).
#' @param plot_width Numeric. Width of saved plots in inches (default: 8).
#' @param plot_height Numeric. Height of saved plots in inches (default: 6).
#' @param save_dir Character. Directory to save plots and data (default:
#'                "ModelData/best_model_result/shap_results" using here::here()).
#'
#' @return Depending on input:
#' \itemize{
#'   \item For Best_Model input: Updated Best_Model object with plots stored in
#'         @shap.result[["shap_dependences"]]
#'   \item For shapviz input: List of ggplot objects or single plot if only one feature
#' }
#'
#' @details
#' For each selected feature, the function:
#' 1. Creates a dependence plot showing SHAP values vs feature values
#' 2. Adds a smoothed trend line for numeric features
#' 3. Optionally saves plots as PDF files
#' 4. Returns either updated Best_Model object or plot list
#'
#' @examples
#' \dontrun{
#' # Using Best_Model object
#' bm <- ModelShapDependence(best_model_object,
#'                          features = c("age", "income"),
#'                          save_plots = TRUE)
#'
#' # Using shapviz object directly
#' plots <- ModelShapDependence(object = NULL,
#'                            shp = shap_values,
#'                            n_top = 5)
#' }
#'
#' @seealso
#' \code{\link[shapviz]{sv_dependence}} for the underlying plotting function,
#' \code{\link{ModelShapValue}} for SHAP value calculation
#'
#' @import ggplot2 
#' @importFrom wesanderson wes_palette
#' @importFrom viridis viridis
#' @importFrom here here
#' @export
ModelShapDependence <- function(object,
                                shp = NULL,
                                features = NULL,
                                n_top = 3,
                                palette_name = "AsteroidCity1",
                                base_size = 14,
                                save_plots = TRUE,
                                plot_width = 8,
                                plot_height = 6,
                                save_dir = here("ModelData", "best_model_result","shap_results")) {
  
  cat("=== Starting SHAP Dependence Analysis ===\n")
  
  if (!inherits(object, "Best_Model") && is.null(shp)) {
    stop("Must provide either a Best_Model object or shapviz object")
  }
  
  if (inherits(object, "Best_Model")) {
    cat("Extracting SHAP values from Best_Model object...\n")
    shp <- object@shap.result[["shap_values"]]
    if (is.null(shp)) stop("No SHAP values found in Best_Model object")
  } else {
    cat("Using provided shapviz object...\n")
  }
  
  all_features <- colnames(shp$S)
  cat("Available features:", paste(all_features, collapse = ", "), "\n")
  
  if (is.null(features)) {
    cat("Selecting top", n_top, "features by importance...\n")
    features <- get_top_features(shp, n = n_top)
  } else if (is.numeric(features)) {
    cat("Selecting top", features, "features by importance...\n")
    features <- get_top_features(shp, n = features)
  } else if (is.character(features)) {
    cat("Using user-specified features...\n")
    invalid_features <- setdiff(features, all_features)
    if (length(invalid_features) > 0) {
      stop("Invalid features: ", paste(invalid_features, collapse = ", "),
           "\nAvailable features: ", paste(all_features, collapse = ", "))
    }
  }
  cat("Final features to analyze:", paste(features, collapse = ", "), "\n")
  
  colors <-wes_palette(palette_name, length(features), type = "continuous")
  
  cat("Generating dependence plots...\n")
  plots <- lapply(seq_along(features), function(i) {
    p <- sv_dependence(shp, v = features[i]) +
      theme_classic(base_size = base_size) +
      labs(title = paste("SHAP Dependence on", features[i]),
           x = features[i],
           y = "SHAP Value") +
      scale_color_gradientn(colors = colors,
                            name = "Feature Value",
                            guide = guide_colorbar(barheight = unit(2, "cm")))
    
    if (is.numeric(shp$X[[features[i]]])) {
      p <- p + geom_smooth(color = "red", se = FALSE, method = "loess")
    }
    return(p)
  })
  names(plots) <- features
  
  
  print(plots)
  
  if (save_plots) {
    cat("Saving plots to", save_dir, "...\n")
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat("Created directory:", save_dir, "\n")
    }
    lapply(names(plots), function(feat) {
      ggsave(
        filename = file.path(save_dir, paste0("shap_dependence_", feat, ".pdf")),
        plot = plots[[feat]], 
        width = plot_width, 
        height = plot_height
      )
      cat("Plot saved to:", file.path(save_dir), "\n")
      cat("Saved plot for feature:", feat, "\n")
    })
  }
  
  if (inherits(object, "Best_Model")) {
    cat("Updating Best_Model object with plots...\n")
    object@shap.result[["shap_dependences"]] <- plots
    return(object)
  } else {
    cat("Returning plot list...\n")
    if (length(plots) == 1) plots[[1]] else plots
  }
  
  cat("=== Analysis Completed ===\n")
}





#' Extract Top Important Features from SHAP Values
#'
#' Calculates and returns the top n most important features based on mean absolute SHAP values.
#' This is a helper function typically used to identify the most influential features
#' for visualization or further analysis.
#'
#' @param shp A shapviz object containing SHAP values (must have $S matrix component).
#' @param n Integer. Number of top features to return (default: 3).
#'
#' @return Character vector of feature names, ordered by decreasing importance.
#'         If n exceeds the number of available features, returns all features.
#'
#' @details
#' Feature importance is calculated as the mean absolute SHAP value across all observations:
#' \deqn{Importance_j = mean(|SHAP_j|)}
#' where \eqn{SHAP_j} are the SHAP values for feature j.
#'
#' @examples
#' \dontrun{
#' # Assuming 'shap_result' is a valid shapviz object
#' top_features <- get_top_features(shap_result, n = 5)
#' 
#' # Use in visualization
#' ModelShapDependence(object = NULL, shp = shap_result, features = top_features)
#' }
#'
#' @seealso
#' \code{\link{ModelShapDependence}} for a function that uses these top features,
#' \code{\link[shapviz]{shapviz}} for SHAP value object structure
#'
#' @export
get_top_features <- function(shp, n = 3) {
  imp <- colMeans(abs(shp$S))
  names(sort(imp, decreasing = TRUE)[1:min(n, length(imp))])
}





#' Get Top N Most Important Features from SHAP Values
#'
#' This function calculates and returns the names of the top N most important features
#' based on mean absolute SHAP values. Feature importance is determined by averaging
#' the absolute SHAP values across all observations for each feature.
#'
#' @param shp A shapviz object containing SHAP values in its `S` matrix component.
#'            The object should have a `$S` matrix where columns represent features
#'            and rows represent observations.
#' @param n Integer specifying the number of top features to return (default: 1).
#'          If n exceeds the number of available features, all features will be returned.
#'
#' @return A character vector containing the names of the top N most important features,
#'         ordered by decreasing importance. Returns all feature names if n > number of features.
#'
#' @details
#' The importance score for each feature is calculated as:
#' \deqn{Importance_j = \frac{1}{m}\sum_{i=1}^{m}|SHAP_{i,j}|}
#' where:
#' \itemize{
#'   \item \eqn{m} is the number of observations
#'   \item \eqn{SHAP_{i,j}} is the SHAP value for feature j and observation i
#' }
#'
#' @examples
#' \dontrun{
#' # Get top 3 most important features
#' top_features <- get_top_features(shap_values, n = 3)
#' 
#' # Use the result to create dependence plots
#' sv_dependence(shap_values, v = top_features[1])
#' 
#' # When n > number of features
#' get_top_features(shap_values, n = 10)  # returns all features
#' }
#'
#' @seealso
#' \code{\link[shapviz]{shapviz}} for creating SHAP objects,
#' \code{\link{ModelShapDependence}} which uses this function internally
#'
#' @export
get_top_features <- function(shp, n = 1) {
  importance <- colMeans(abs(shp$S))
  names(sort(importance, decreasing = TRUE)[1:min(n, length(importance))])
}



#' Generate SHAP Force Plots for Model Interpretation
#'
#' Creates SHAP force plots to visualize individual prediction explanations.
#' Works with either Best_Model objects or direct shapviz objects, allowing
#' visualization of how features contribute to specific predictions.
#'
#' @param object A Best_Model class object containing SHAP results. Alternatively,
#'               can be NULL if shp parameter is provided directly.
#' @param shp Optional. A shapviz object containing pre-computed SHAP values.
#'            Required if object is not a Best_Model.
#' @param samples Selection of samples to plot. Can be:
#'                - NULL (auto-selects n_samples extreme cases)
#'                - "typical" (selects high/median/low predictions)
#'                - Numeric vector (specific sample indices)
#' @param n_samples Integer. Number of samples to select when samples=NULL (default: 3).
#' @param save_plots Logical. Whether to save plots to disk (default: TRUE).
#' @param palette_name Character. Name of color palette for plots (default: "AsteroidCity1").
#'                    Falls back to viridis if specified palette not found.
#' @param save_dir Character. Directory path for saving results (default:
#'                "ModelData/best_model_result/shap_results" using here::here()).
#' @param plot_width Numeric. Width of saved plots in inches (default: 5).
#' @param plot_height Numeric. Height of saved plots in inches (default: 5).
#' @param base_size Numeric. Base font size for plots (default: 14).
#'
#' @return Depending on input:
#' \itemize{
#'   \item For Best_Model input: Updated Best_Model object with plots stored in
#'         @shap.result[["force_plots"]] and sample indices in @shap.result[["last_force_samples"]]
#'   \item For shapviz input: List of ggplot objects or single plot if only one sample
#' }
#'
#' @details
#' The function provides three sample selection modes:
#' 1. Automatic extreme samples (largest absolute SHAP contributions)
#' 2. Typical samples (high/median/low predictions)
#' 3. User-specified sample indices
#'
#' Each force plot shows how features push the prediction from the base value
#' to the final model output for a single observation.
#'
#' @examples
#' \dontrun{
#' # Using Best_Model object with automatic sample selection
#' bm <- ModelShapForce(best_model_object, n_samples = 5)
#'
#' # Using shapviz object with specific samples
#' plots <- ModelShapForce(object = NULL,
#'                       shp = shap_values,
#'                       samples = c(42, 99))
#'
#' # Using typical samples
#' typical_plots <- ModelShapForce(best_model_object, samples = "typical")
#' }
#'
#' @seealso
#' \code{\link[shapviz]{sv_force}} for the underlying plotting function,
#' \code{\link{ModelShapValue}} for SHAP value calculation,
#' \code{\link{get_extreme_samples}} and \code{\link{get_typical_samples}} for sample selection
#'
#' @import ggplot2 
#' @importFrom wesanderson wes_palette
#' @importFrom here here
#' @export
ModelShapForce <- function(object,
                           shp = NULL,
                           samples = NULL,
                           n_samples = 3,
                           save_plots =TRUE,
                           palette_name = "AsteroidCity1",
                           save_dir = here("ModelData", "best_model_result","shap_results"),
                           plot_width = 5,
                           plot_height = 5,
                           base_size = 14) {
  
  if (inherits(object, "Best_Model")) {
    cat("Input is of class 'Best_Model'. Extracting SHAP data...\n")
    shp <- object@shap.result[["shap_values"]]
  } else if (is.null(shp)) {
    stop("Must provide either a Best_Model object or shapviz object")
  }
  
  all_samples <- 1:nrow(shp$S)
  
  if (is.null(samples)) {
    samples <- get_extreme_samples(shp, n = n_samples)
    cat("Auto-selected samples with largest |SHAP|:", paste(samples, collapse = ", "), "\n")
  } else if (identical(samples, "typical")) {
    samples <- get_typical_samples(shp)
    cat("Selected typical samples (high/median/low prediction):", paste(samples, collapse = ", "), "\n")
  } else if (is.numeric(samples)) {
    invalid <- setdiff(samples, all_samples)
    if (length(invalid) > 0) {
      stop("Invalid sample indices: ", paste(invalid, collapse = ", "),
           "\nValid range: 1 to ", length(all_samples))
    }
  } else {
    stop("`samples` must be NULL, 'typical', or numeric vector")
  }
  colors <- wes_palette(palette_name, length(samples), type = "continuous")
  
  plots <- lapply(seq_along(samples), function(i) {
    sv_force(shp, row_id = samples[i]) +
      theme_minimal(base_size = base_size) +
      labs(title = paste("SHAP Force Plot (Sample", samples[i], ")"),
           x = "Feature Contribution",
           subtitle = paste("Prediction:", round(get_prediction(shp, samples[i]), 3))) +
      theme(plot.title = element_text(color = colors[i]))
  })
  names(plots) <- paste0("sample_", samples)
  print(plots)
  if (save_plots) {
    if (!dir.exists(save_dir)) {
      cat("Creating directory:", save_dir, "\n")
      dir.create(save_dir, recursive = TRUE)
    }
    
    cat("Saving plots to:", save_dir, "\n")
    for (feat in names(plots)) {
      filename <- file.path(save_dir, paste0("shap_dependence_", feat, ".pdf"))
      ggsave(filename, plot = plots[[feat]], width = plot_width, height = plot_height)
    }
  }
  if (inherits(object, "Best_Model")) {
    cat("Updating 'Best_Model' object...\n")
    object@shap.result[["force_plots"]] <- plots
    object@shap.result[["last_force_samples"]] <- samples
    return(object)
  } else {
    if (length(plots) == 1) plots[[1]] else plots
  }
}


#' Calculate Final Prediction from SHAP Values
#'
#' Computes the final model prediction for a given sample by summing its SHAP values
#' with the baseline (expected) value. This reconstructs the model's output from
#' the SHAP decomposition.
#'
#' @param shp A shapviz object containing:
#' \itemize{
#'   \item `$S`: Matrix of SHAP values (samples × features)
#'   \item `$baseline`: Baseline prediction value (intercept)
#' }
#' @param sample_id Integer. Index of the sample to calculate prediction for
#'                 (row index in shp$S matrix).
#'
#' @return Numeric. The reconstructed prediction value for the specified sample,
#'         calculated as:
#'         \deqn{prediction = baseline + \sum_{j=1}^{m} SHAP_{i,j}}
#'         where:
#'         \itemize{
#'           \item \eqn{baseline} is the average prediction
#'           \item \eqn{SHAP_{i,j}} are the SHAP values for sample i
#'           \item \eqn{m} is the number of features
#'         }
#'
#' @details
#' This function implements the SHAP additive property:
#' \deqn{f(x) = \phi_0 + \sum_{j=1}^{M}\phi_j}
#' where \eqn{\phi_0} is the baseline and \eqn{\phi_j} are the SHAP values.
#'
#' @examples
#' \dontrun{
#' # Get prediction for sample 42
#' pred <- get_prediction(shap_values, sample_id = 42)
#'
#' # Verify against original predictions
#' all.equal(pred, original_predictions[42])  # Should be TRUE
#' }
#'
#' @seealso
#' \code{\link[shapviz]{shapviz}} for SHAP object structure,
#' \code{\link{ModelShapForce}} which uses this function internally
#'
#' @export
get_prediction <- function(shp, sample_id) {
  shp$baseline + sum(shp$S[sample_id, ])
}

#' Identify Samples with Most Extreme SHAP Values
#'
#' Finds samples with the largest absolute SHAP contributions, indicating
#' observations where features had the strongest impact on predictions
#' (either positive or negative).
#'
#' @param shp A shapviz object containing SHAP values in its `S` matrix component.
#'            The object should have:
#'            \itemize{
#'              \item `$S`: Matrix of SHAP values (samples × features)
#'            }
#' @param n Integer. Number of extreme samples to return (default: 3).
#'          If n exceeds number of samples, returns all samples sorted by importance.
#'
#' @return Integer vector of sample indices/row numbers, ordered by decreasing
#'         total absolute SHAP contribution. The most "extreme" sample is first.
#'
#' @details
#' Extreme samples are identified by calculating the sum of absolute SHAP values
#' for each observation:
#' \deqn{TotalSHAP_i = \sum_{j=1}^{m} |SHAP_{i,j}|}
#' where:
#' \itemize{
#'   \item \eqn{m} is the number of features
#'   \item \eqn{SHAP_{i,j}} is the SHAP value for sample i and feature j
#' }
#'
#' These samples typically represent:
#' - Highly influential predictions
#' - Potential outliers
#' - Interesting edge cases for model interpretation
#'
#' @examples
#' \dontrun{
#' # Get top 5 most extreme samples
#' extreme_samples <- get_extreme_samples(shap_values, n = 5)
#'
#' # View their force plots
#' sv_force(shap_values, row_id = extreme_samples[1])
#' }
#'
#' @seealso
#' \code{\link{get_typical_samples}} for finding representative samples,
#' \code{\link{ModelShapForce}} which uses this function for visualization
#'
#' @export
get_extreme_samples <- function(shp, n = 3) {
  total_shap <- rowSums(abs(shp$S))
  order(total_shap, decreasing = TRUE)[1:min(n, length(total_shap))]
}

#' Identify Representative Samples Across Prediction Range
#'
#' Selects three characteristic samples representing different prediction ranges:
#' high, median, and low predicted values. Useful for understanding model behavior
#' across typical prediction scenarios.
#'
#' @param shp A shapviz object containing:
#'            \itemize{
#'              \item `$S`: Matrix of SHAP values (samples × features)
#'              \item `$baseline`: Baseline prediction value
#'            }
#'
#' @return Integer vector of length 3 containing sample indices for:
#' \itemize{
#'   \item Highest prediction (most positive outcome)
#'   \item Median prediction (most typical outcome)
#'   \item Lowest prediction (most negative outcome)
#' }
#' 
#' @details
#' The function calculates final predictions and selects samples by:
#' \enumerate{
#'   \item Calculating predictions: \eqn{pred = baseline + \sum SHAP}
#'   \item Finding:
#'     \itemize{
#'       \item which.max(pred) - sample with highest prediction
#'       \item which.min(abs(pred - median(pred))) - sample closest to median
#'       \item which.min(pred) - sample with lowest prediction
#'     }
#' }
#' 
#' These samples represent the full spectrum of model predictions and are
#' particularly useful for:
#' \itemize{
#'   \item Understanding feature contributions across different prediction levels
#'   \interleave{}{}
#'   \item Identifying how the model behaves for typical vs. extreme cases
#' }
#'
#' @examples
#' \dontrun{
#' # Get typical samples from SHAP results
#' typical_samples <- get_typical_samples(shap_values)
#' 
#' # View their force plots
#' sv_force(shap_values, row_id = typical_samples[2])  # Median prediction
#' }
#'
#' @seealso 
#' \code{\link{get_extreme_samples}} for samples with largest SHAP contributions,
#' \code{\link{ModelShapForce}} which uses this function when samples="typical"
#'
#' @export
get_typical_samples <- function(shp) {
  pred <- shp$baseline + rowSums(shp$S)
  c(
    which.max(pred),                     
    which.min(abs(pred - median(pred))), 
    which.min(pred)                      
  )
}







#' Visualize SHAP Feature Importance
#'
#' Generates SHAP feature importance plots (beeswarm or bar plot) to show the impact
#' of each feature on model predictions. Works with both Best_Model objects and
#' direct shapviz inputs.
#'
#' @param object A Best_Model class object containing SHAP results. Alternatively,
#'               can be NULL if shp parameter is provided directly.
#' @param shp Optional shapviz object. Required if object is not a Best_Model.
#' @param plot_type Character. Type of importance plot: "beeswarm" (default) or "bar".
#' @param n_features Integer or NULL. Number of top features to display. If NULL,
#'                  shows all features (default: NULL).
#' @param palette_name Character. Color palette name (default: "AsteroidCity1").
#'                    Falls back to viridis if specified palette not available.
#' @param base_size Numeric. Base font size for plots (default: 14).
#' @param save_plots Logical. Whether to save plots to disk (default: TRUE).
#' @param save_dir Character. Directory path for saving results (default:
#'                "ModelData/best_model_result/shap_results" using here::here()).
#' @param plot_width Numeric. Width of saved plots in inches (default: 5).
#' @param plot_height Numeric. Height of saved plots in inches (default: 5).
#'
#' @return Depending on input:
#' \itemize{
#'   \item For Best_Model input: Updated object with plot stored in
#'         @shap.result[["importance_<plot_type>"]]
#'   \item For shapviz input: ggplot object of the importance plot
#' }
#'
#' @details
#' The function provides two visualization types:
#' \describe{
#'   \item{beeswarm}{Shows distribution of SHAP values for each feature, with
#'                  points colored by feature value}
#'   \item{bar}{Displays mean absolute SHAP values as a bar plot}
#' }
#'
#' @examples
#' \dontrun{
#' # Using Best_Model object
#' bm <- ModelShapImportance(best_model_object,
#'                          plot_type = "beeswarm",
#'                          n_features = 10)
#'
#' # Using shapviz object directly
#' beeswarm_plot <- ModelShapImportance(object = NULL,
#'                                     shp = shap_values,
#'                                     plot_type = "beeswarm")
#' }
#'
#' @seealso
#' \code{\link[shapviz]{sv_importance}} for the underlying plotting function,
#' \code{\link{get_top_features}} for feature selection
#'
#' @import ggplot2 
#' @importFrom wesanderson wes_palette
#' @importFrom here here
#' @export
ModelShapImportance <- function(object,
                                shp = NULL,
                                plot_type = "beeswarm",
                                n_features = NULL,
                                palette_name = "AsteroidCity1",
                                base_size = 14,
                                save_plots = TRUE,
                                save_dir = here("ModelData", "best_model_result","shap_results"),
                                plot_width = 5,
                                plot_height = 5) {
  
  if (!inherits(object, "Best_Model") && is.null(shp)) {
    stop("Must provide either a Best_Model object or shapviz object")
  }
  cat("=== Starting SHAP Importance Analysis ===\n")
  if (!plot_type %in% c("beeswarm", "bar")) {
    stop("plot_type must be either 'beeswarm' or 'bar'")
  }
  
  if (inherits(object, "Best_Model")) {
    shp <- object@shap.result[["shap_values"]]
    if (is.null(shp)) stop("No SHAP values found in Best_Model object")
  }
  
  all_features <- colnames(shp$S)
  if (is.null(n_features)) {
    features <- all_features
  } else {
    features <- get_top_features(shp, n = n_features)
  }
  
  colors <- tryCatch(
    wes_palette(palette_name, length(features), type = "continuous"),
    error = function(e) viridis::viridis(length(features))
  )
  
  p <- sv_importance(shp, kind = plot_type) +
    theme_classic(base_size = base_size) +
    labs(title = paste("SHAP Feature Importance (", plot_type, ")"),
         x = "SHAP Value",
         y = "Feature") +
    scale_color_gradientn(colors = colors,
                          name = "Feature Value",
                          guide = guide_colorbar(barheight = unit(2, "cm")))
  
  if (!is.null(n_features)) {
    p <- p + scale_y_discrete(limits = rev(features))
  }
  
  if (save_plots) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    
    filename <- file.path(save_dir, 
                          paste0("shap_importance_", plot_type, ".pdf"))
    ggsave(filename, plot = p, 
           width = plot_width, height = plot_height)
    cat("Saving plots to:", save_dir, "\n")
    
  }
  
  if (inherits(object, "Best_Model")) {
    object@shap.result[[paste0("importance_", plot_type)]] <- p
    return(object)
  } else {
    return(p)
  }
}

#' Select Top Important Features by SHAP Values
#'
#' Identifies the most influential features based on mean absolute SHAP values.
#' Useful for focusing analysis on key predictors.
#'
#' @param shp A shapviz object containing SHAP values in its `S` matrix.
#' @param n Integer. Number of top features to return (default: 10).
#'
#' @return Character vector of feature names, ordered by decreasing importance.
#'         If n exceeds number of features, returns all features sorted.
#'
#' @details
#' Feature importance is calculated as:
#' \deqn{Importance_j = \frac{1}{m}\sum_{i=1}^{m}|SHAP_{i,j}|}
#' where \eqn{m} is the number of samples and \eqn{SHAP_{i,j}} is the SHAP value
#' for feature j and sample i.
#'
#' @examples
#' \dontrun{
#' # Get top 5 most important features
#' top_features <- get_top_features(shap_values, n = 5)
#'
#' # Use in importance plot
#' ModelShapImportance(shap_values, n_features = 5)
#' }
#'
#' @seealso
#' \code{\link{ModelShapImportance}} which uses this function internally
#'
#' @export
get_top_features <- function(shp, n = 10) {
  imp <- colMeans(abs(shp$S))
  names(sort(imp, decreasing = TRUE)[1:min(n, length(imp))])
}


