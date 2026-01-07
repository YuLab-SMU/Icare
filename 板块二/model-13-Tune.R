#' Visualize Hyperparameter Tuning Results
#'
#' Generates and optionally saves ggplot2 visualization of model hyperparameter tuning results.
#' Works with both raw tuning objects and `Train_Model` S4 objects. For `Train_Model` objects,
#' automatically updates the object with the plot and returns the modified object.
#'
#' @param object Either:
#'   - A model tuning object from caret (train) containing hyperparameter results
#'   - A `Train_Model` S4 object containing tuning results in @best.model.result slot
#' @param group_col Character specifying the grouping variable column name (default: "group")
#' @param palette_name Character specifying Wes Anderson palette name for coloring
#'        (default: "AsteroidCity1")
#' @param save_plots Logical indicating whether to save plots to file (default: TRUE)
#' @param save_dir Character specifying directory path to save outputs
#'        (default: here("ModelData", "tuning_results"))
#' @param plot_width Numeric plot width in inches (default: 6)
#' @param plot_height Numeric plot height in inches (default: 6)
#' @param base_size Numeric base font size for plot (default: 14)
#'
#' @return Depending on input:
#'   - For regular tuning objects: Invisibly returns the ggplot object
#'   - For `Train_Model` objects: Returns the modified S4 object with plot stored in
#'     @best.model.result[["hyperparameter_plot"]] slot
#'
#' @export
#'
#' @examples
#' \dontrun{
#' object_model<-ModelHyperparameterPlot(object_model)
#' }
#'
#' @importFrom wesanderson wes_palette
#'
#' @import here
ModelHyperparameterPlot<-function(object,
         group_col = "group",
         palette_name = "AsteroidCity1",
         save_plots = TRUE,
         save_dir = here("ModelData", "tuning_results"),
         plot_width = 6,
         plot_height = 6,
         base_size=14) {

  if (inherits(object, "Train_Model")) {
    cat("Input is of class 'Train_Model'. Extracting datasets...\n")
    group_col <-object@group_col
    best_model <- object@best.model.result[["model"]]
    model_type<-object@best.model.result[["model_type"]]
  }
  trellis.par.set(caretTheme())
  p<-ggplot(best_model)+
    theme_minimal(base_size = base_size)+
    scale_color_manual(values = wes_palette(palette_name)) +
    theme(
      axis.text.x = element_text(size = 12, color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95")
    )
  print(p)

  if (save_plots) {
    ggsave(filename = file.path(save_dir, "hyperparameter_tuning_curve.pdf"), plot = p,
           width = plot_width, height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "hyperparameter_tuning_curve.pdf"), "\n")
  }
  if (inherits(object, "Train_Model")) {
    object@best.model.result[["hyperparameter_plot"]] <- p
    cat("Updating 'Train_Model' object...\n")
    cat("The 'Train_Model' object has been updated with the following slots:\n")
    cat("- 'best.model.result' slot updated.\n")
    return(object)
  }
}

#' Generate and Evaluate Hyperparameter Tuning Suggestions
#'
#' This function generates intelligent hyperparameter tuning suggestions based on the best performing model
#' in a Train_Model object. It uses a combination of methods:
#' 1. For continuous parameters: Expands the search space around the current best values using logarithmic scaling
#' 2. For discrete parameters: Explores neighboring integer values
#' 3. For parameters with constraints: Ensures values stay within valid ranges
#'
#' The function implements a "local search" strategy by focusing the parameter search around the current best values,
#' while expanding the range by a specified factor. This approach is more efficient than grid search when you have
#' a good starting point (the current best model).
#'
#' @param object A Train_Model object containing the best model to tune
#' @param expand_factor Numeric factor to expand parameter ranges around current best values (default: 1.5)
#' @param suggestions Optional list of manual tuning suggestions (bypasses automatic generation)
#' @param control trainControl object specifying resampling method (default: 5-fold CV)
#' @param classProbs Logical indicating whether to compute class probabilities (default: TRUE)
#' @param allowParallel Logical indicating whether to allow parallel processing (default: TRUE)
#' @param group_col Name of the grouping variable column (default: "group")
#' @param tune_grids List of default tuning grids for reference ranges
#'
#' @return Returns the input Train_Model object with updated tuning results, or a list containing:
#'         - model: The tuned model
#'         - train_performance: Performance metrics on training data
#'         - suggestions: The parameter grid used for tuning
#'
#' @examples
#' \dontrun{
#' # Manual tuning suggestions
#' custom_tune <- list(
#'   rf = expand.grid(mtry = c(3, 5, 7)),
#'   gbm = expand.grid(n.trees = c(100, 150), interaction.depth = c(3, 5))
#' )
#' object_model <- ModelTuneSuggestion(object_model, suggestions = custom_tune)
#' }
#' @export
#' @import caret
#' @import stats
ModelTuneSuggestion <- function(object,
                                expand_factor = 1.5,
                                suggestions = NULL,
                                control = trainControl(method = "cv", number = 5),
                                classProbs = TRUE,
                                allowParallel = TRUE,
                                group_col = "group",
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
                                )
)
{

  if (!inherits(object, "Train_Model")) {
    stop("Input object must be of class 'Train_Model'")
  }

  if (is.null(object@best.model.result)) {
    stop("No best model result found in the object")
  }

  if (expand_factor <= 0) {
    stop("expand_factor must be a positive number")
  }

  model_type <- object@best.model.result[["model_type"]]
  best_model <- object@best.model.result[["model"]]
  best_params <- best_model$bestTune
  data_sets <- Extract_filtered.set(object)
  train_data <- data_sets$training
  test_data <- data_sets$testing

  if (is.null(suggestions)) {
    suggestions <- list()

    if (model_type == "gbm") {
      current_n_trees <- best_params$n.trees
      current_depth <- best_params$interaction.depth
      current_shrinkage <- best_params$shrinkage
      current_minobs <- best_params$n.minobsinnode

      suggestions$gbm <- expand.grid(
        n.trees = seq(max(50, current_n_trees/expand_factor),
                      current_n_trees*expand_factor,
                      length.out = 3) %>% round(),
        interaction.depth = seq(max(1, current_depth-1),
                                current_depth+1,
                                by = 1) %>% unique(),
        shrinkage = pmax(0.0001, pmin(0.1, exp(seq(log(max(0.0005, current_shrinkage/expand_factor)),
                                                   log(min(0.1, current_shrinkage*expand_factor)),
                                                   length.out = 3)))),
        n.minobsinnode = seq(max(5, current_minobs/expand_factor),
                             current_minobs*expand_factor,
                             length.out = 3) %>% round()
      )

    } else if (model_type == "rf") {
      current_mtry <- best_params$mtry
      max_features <- ncol(train_data) - 1

      suggestions$rf <- expand.grid(
        mtry = unique(pmax(1, pmin(max_features,
                                   seq(floor(current_mtry/expand_factor),
                                       ceiling(current_mtry*expand_factor),
                                       by = 1))))
      )

    } else if (model_type == "svmLinear") {
      current_C <- best_params$C

      suggestions$svmLinear <- expand.grid(
        C = pmax(0.001, pmin(10, exp(seq(log(max(0.005, current_C/expand_factor)),
                                         log(min(5, current_C*expand_factor)),
                                         length.out = 5))))
      )

    } else if (model_type == "svmRadial") {
      current_C <- best_params$C
      current_sigma <- best_params$sigma

      new_sigma <- if (current_sigma == max(tune_grids$svmRadial$sigma)) {
        seq(current_sigma, current_sigma*2, length.out = 3)
      } else if (current_sigma == min(tune_grids$svmRadial$sigma)) {
        seq(max(0.005, current_sigma/2), current_sigma*1.5, length.out = 3)
      } else {
        seq(current_sigma/1.5, current_sigma*1.5, length.out = 3)
      }

      suggestions$svmRadial <- expand.grid(
        sigma = pmax(0.001, pmin(1, new_sigma)),
        C = pmax(0.05, pmin(10, exp(seq(log(max(0.05, current_C/expand_factor)),
                                        log(min(10, current_C*expand_factor)),
                                        length.out = 3))))
      )

    } else if (model_type == "glmnet") {
      current_alpha <- best_params$alpha
      current_lambda <- best_params$lambda

      new_alpha <- if (current_alpha == max(tune_grids$glmnet$alpha)) {
        seq(current_alpha, 1, length.out = 3)
      } else if (current_alpha == min(tune_grids$glmnet$alpha)) {
        seq(0, current_alpha*1.5, length.out = 3)
      } else {
        seq(max(0, current_alpha-0.2), min(1, current_alpha+0.2), length.out = 3)
      }

      suggestions$glmnet <- expand.grid(
        alpha = pmax(0, pmin(1, new_alpha)),
        lambda = 10^seq(log10(max(1e-5, current_lambda/expand_factor)),
                        log10(min(1, current_lambda*expand_factor)),
                        length.out = 5)
      )

    } else if (model_type == "xgbTree") {
      current_nrounds <- best_params$nrounds
      current_max_depth <- best_params$max_depth
      current_eta <- best_params$eta

      suggestions$xgbTree <- expand.grid(
        nrounds = seq(max(50, current_nrounds/expand_factor),
                      current_nrounds*expand_factor,
                      length.out = 3) %>% round(),
        max_depth = seq(max(1, current_max_depth-1),
                        current_max_depth+1,
                        by = 1) %>% unique(),
        eta = pmax(0.001, pmin(0.3, exp(seq(log(max(0.001, current_eta/expand_factor)),
                                            log(min(0.3, current_eta*expand_factor)),
                                            length.out = 3)))),
        gamma = 0,
        colsample_bytree = 1,
        min_child_weight = 1,
        subsample = 1
      )

    } else if (model_type == "rpart") {
      current_cp <- best_params$cp

      suggestions$rpart <- expand.grid(
        cp = pmax(0.0001, pmin(0.1, exp(seq(log(max(0.0001, current_cp/expand_factor)),
                                            log(min(0.1, current_cp*expand_factor)),
                                            length.out = 10))))
      )

    } else if (model_type == "earth") {
      current_degree <- best_params$degree
      current_nprune <- best_params$nprune

      suggestions$earth <- expand.grid(
        degree = seq(max(1, current_degree-1), current_degree+1, by = 1) %>% unique(),
        nprune = seq(max(2, current_nprune/expand_factor),
                     current_nprune*expand_factor,
                     length.out = 5) %>% round()
      )

    } else {
      warning(paste("Automatic tuning suggestions not implemented for model type:", model_type,
                    "\nUsing default tuning grid or user-provided suggestions."))
    }

    if (!is.null(suggestions[[model_type]])) {
      cat("Parameter suggestions generated for:\n")
      cat(paste0("  - ", model_type, "\n"))
      print(suggestions[[model_type]])
    }
  }

  cat("\n--- Training Model with Suggested Parameters ---\n")

  if (object.size(train_data) > 1e8) {
    warning("Large training dataset detected (", format(object.size(train_data), units = "MB"),
            "). Consider subsampling for faster tuning.")
  }

  model_fit <- train_and_evaluate_models(
    data = train_data,
    methods = model_type,
    control = control,
    tune_grids = suggestions,
    classProbs = classProbs,
    allowParallel = allowParallel,
    group_col = group_col
  )

  closeAllConnections()
  best_params <- model_fit$bestTune

  if (model_type == "gbm") {
    current_n_trees <- ifelse(is.null(best_params$n.trees) || is.na(best_params$n.trees),
                              100,
                              as.numeric(best_params$n.trees)[1])

    current_depth <- ifelse(is.null(best_params$interaction.depth) || is.na(best_params$interaction.depth),
                            2,
                            as.numeric(best_params$interaction.depth)[1])

    current_shrinkage <- ifelse(is.null(best_params$shrinkage) || is.na(best_params$shrinkage),
                                0.01,
                                as.numeric(best_params$shrinkage)[1])

    current_minobs <- ifelse(is.null(best_params$n.minobsinnode) || is.na(best_params$n.minobsinnode),
                             10,
                             as.numeric(best_params$n.minobsinnode)[1])

    validate_parameter <- function(param, name, default) {
      if (is.na(param) || !is.numeric(param)) {
        warning(paste("Invalid", name, "value, using default:", default))
        return(default)
      }
      return(param)
    }

    current_n_trees <- validate_parameter(current_n_trees, "n.trees", 100)
    current_depth <- validate_parameter(current_depth, "interaction.depth", 2)
    current_shrinkage <- validate_parameter(current_shrinkage, "shrinkage", 0.01)
    current_minobs <- validate_parameter(current_minobs, "n.minobsinnode", 10)

    cat("\n=== GBM Parameter Debug Info ===\n")
    cat("n.trees:", current_n_trees, "\n")
    cat("interaction.depth:", current_depth, "\n")
    cat("shrinkage:", current_shrinkage, "\n")
    cat("n.minobsinnode:", current_minobs, "\n")

    suggestions$gbm <- tryCatch({
      expand.grid(
        n.trees = {
          from <- max(50, current_n_trees/expand_factor)
          to <- current_n_trees * expand_factor
          cat("Generating n.trees sequence: from", from, "to", to, "\n")
          seq(from, to, length.out = 3) %>% round()
        },
        interaction.depth = {
          depths <- seq(max(1, current_depth-1), current_depth+1, by = 1)
          cat("Generating interaction.depth sequence:", paste(depths, collapse = ", "), "\n")
          unique(depths)
        },
        shrinkage = {
          shrinkage_seq <- exp(seq(log(max(0.0005, current_shrinkage/expand_factor)),
                                   log(min(0.1, current_shrinkage*expand_factor)),
                                   length.out = 3))
          cat("Generating shrinkage sequence:", paste(round(shrinkage_seq, 5), collapse = ", "), "\n")
          pmax(0.0001, pmin(0.1, shrinkage_seq))
        },
        n.minobsinnode = {
          minobs_seq <- seq(max(5, current_minobs/expand_factor),
                            current_minobs*expand_factor,
                            length.out = 3)
          cat("Generating n.minobsinnode sequence:", paste(round(minobs_seq), collapse = ", "), "\n")
          round(minobs_seq)
        }
      )
    }, error = function(e) {
      warning("Error generating GBM parameter grid, using defaults: ", e$message)
      expand.grid(
        n.trees = c(50, 100, 150),
        interaction.depth = c(1, 2, 3),
        shrinkage = c(0.001, 0.01, 0.1),
        n.minobsinnode = c(5, 10, 20)
      )
    })

    if (any(is.na(suggestions$gbm))) {
      warning("Generated parameter grid contains NA values, checking parameter ranges")
      print(suggestions$gbm)
    }
  }


  cat("Evaluating models on the train dataset...\n")
  train_performance <- evaluate_model_performance(
    data = train_data,
    model_result = model_fit,
    group_col = group_col
  )
  print(train_performance)

  if (inherits(object, "Train_Model")) {
    object@best.model.result[["model"]] <- model_fit
    object@best.model.result[["train_performance"]] <- train_performance
    object@best.model.result[["suggestions"]] <- suggestions
    if (!is.null(model_fit$bestTune)) {
      object@best.model.result[["bestTune"]] <- model_fit$bestTune
      best_params <- model_fit$bestTune
    } else if (!is.null(model_fit[[1]]$bestTune)) {
      object@best.model.result[["bestTune"]] <- model_fit[[1]]$bestTune
      best_params <- model_fit[[1]]$bestTune
    } else {
      warning("No bestTune parameters found in the model fit")
      best_params <- NULL
    }
    cat("Updating 'Train_Model' object...\n")
    cat("The 'Train_Model' object has been updated with the following slots:\n")
    cat("- Training performance metrics in 'best.model.result$train_performance'\n")
    cat("- Tuning suggestions in 'best.model.result$suggestions'\n")
    if (!is.null(best_params)) {
      cat("- Optimal tuning parameters in 'best.model.result$bestTune':\n")
      print(best_params)
    }
    return(object)
  } else {
    return(list(
      model = model_fit,
      train_performance = train_performance,
      suggestions = suggestions
    ))
  }
}

#' Plot ROC Curve Comparison Between Original and Tuned Models
#'
#' This function generates ROC curves comparing the performance of an original model
#' versus a tuned model on test data, including AUC values and confidence intervals.
#' The plot can be customized and optionally saved to disk.
#'
#' @param original_model The original model object (must have predict method for type = "prob")
#' @param tuned_model The tuned model object (must have predict method for type = "prob")
#' @param test_data Data frame containing test data with response variable
#' @param group_col Character string specifying the column name of the binary response variable (default = "group")
#' @param palette_name Character string specifying the Wes Anderson palette name (default = "AsteroidCity1")
#' @param base_size Base font size for the plot (default = 14)
#' @param save_plots Logical indicating whether to save plots (default = TRUE)
#' @param save_dir Directory path where plots should be saved (default = here("ModelData", "tuning_results"))
#' @param plot_width Width of the plot in inches (default = 5)
#' @param plot_height Height of the plot in inches (default = 5)
#' @param alpha Transparency level for plot lines (default = 1)
#' @param plot_title Title for the plot (default = "ROC Curve Comparison: Before vs After Tuning")
#'
#' @return Invisibly returns a list containing:
#' \itemize{
#'   \item roc_objects - List of ROC curve objects for both models
#'   \item auc_values - Named vector of AUC values
#' }
#'
#' @import ggplot2
#' @importFrom wesanderson wes_palette
#' @import here
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming you have original_model, tuned_model, and test_data prepared
#' results <- plot_tuning_comparison(
#'   original_model = original_rf,
#'   tuned_model = tuned_rf,
#'   test_data = test_data
#' )
#' }
plot_tuning_comparison <- function(original_model,
                                   tuned_model,
                                   test_data,
                                   group_col = "group",
                                   palette_name = "AsteroidCity1",
                                   base_size = 14,
                                   save_plots = TRUE,
                                   save_dir = here("ModelData", "tuning_results"),
                                   plot_width = 5,
                                   plot_height = 5,
                                   alpha = 1,
                                   plot_title = "ROC Curve Comparison: Before vs After Tuning") {

  roc_list <- list()
  auc_values <- numeric()

  cat("Generating predictions and ROC curves...\n")

  original_pred <- predict(original_model, newdata = test_data, type = "prob")[,2]
  original_roc <- roc(test_data[[group_col]], original_pred, levels = c("0", "1"), direction = "<")
  original_auc <- auc(original_roc)
  original_ci <- ci.auc(original_roc)

  roc_list$original <- original_roc
  auc_values["Original"] <- original_auc

  tuned_pred <- predict(tuned_model, newdata = test_data, type = "prob")[,2]
  tuned_roc <- roc(test_data[[group_col]], tuned_pred, levels = c("0", "1"), direction = "<")
  tuned_auc <- auc(tuned_roc)
  tuned_ci <- ci.auc(tuned_roc)

  roc_list$tuned <- tuned_roc
  auc_values["Tuned"] <- tuned_auc

  plot_data <- rbind(
    data.frame(
      Specificity = 1 - original_roc$specificities,
      Sensitivity = original_roc$sensitivities,
      Model = paste0("Original (AUC = ", round(original_auc, 3),
                     ", CI = [", round(original_ci[1], 3), ", ", round(original_ci[3], 3), "])")
    ),
    data.frame(
      Specificity = 1 - tuned_roc$specificities,
      Sensitivity = tuned_roc$sensitivities,
      Model = paste0("Tuned (AUC = ", round(tuned_auc, 3),
                     ", CI = [", round(tuned_ci[1], 3), ", ", round(tuned_ci[3], 3), "])")
    )
  )

  cat("Creating comparison plot...\n")
  p <- ggplot(plot_data, aes(x = Specificity, y = Sensitivity, color = Model)) +
    geom_line(size = 1.25, alpha = alpha) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    scale_color_manual(values = wes_palette(palette_name)) +
    labs(title = plot_title,
         subtitle = "Including AUC and 95% Confidence Intervals",
         x = "1 - Specificity",
         y = "Sensitivity",
         color = "Model Version (AUC and CI)") +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = c(0.95, 0.05),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.8)),
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey40")
    )

  print(p)

  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    ggsave(filename = file.path(save_dir, "tuning_comparison.pdf"),
           plot = p,
           width = plot_width,
           height = plot_height,
           device = "pdf")
    cat("Plot saved to:", file.path(save_dir, "tuning_comparison.pdf"), "\n")
  }

  cat("\nAUC Values:\n")
  print(auc_values)

  invisible(list(
    roc_objects = roc_list,
    auc_values = auc_values
  ))
}

#' Model Tuning Comparison and Visualization
#'
#' Compares original and tuned model performance using ROC curves and updates Train_Model object.
#'
#' @param object A Train_Model object containing original and tuned models
#' @param group_col Name of the outcome/group column (default: "group")
#' @param palette_name Color palette name (default: "AsteroidCity1")
#' @param plot_title Plot title (default: "Model Performance: Before vs After Hyperparameter Tuning")
#' @param base_size Base font size (default: 14)
#' @param save_plots Whether to save plots (default: TRUE)
#' @param save_dir Directory to save plots (default: here("ModelData", "tuning_results"))
#' @param plot_width Plot width in inches (default: 5)
#' @param plot_height Plot height in inches (default: 5)
#' @param alpha Line transparency (default: 1)
#'
#' @return The updated Train_Model object with comparison results stored in best.model.result slot
#' @export
ModelTuneComparison <- function(object,
                                group_col = "group",
                                palette_name = "AsteroidCity1",
                                plot_title = "Model Performance: Before vs After Hyperparameter Tuning",
                                base_size = 14,
                                save_plots = TRUE,
                                save_dir = here("ModelData", "tuning_results"),
                                plot_width = 5,
                                plot_height = 5,
                                alpha = 1) {

  if (!inherits(object, "Train_Model")) {
    stop("Input object must be of class 'Train_Model'")
  }

  cat("Extracting model components...\n")
  model_type <- object@best.model.result[["model_type"]]
  tuned_model <- object@best.model.result[["model"]][[1]]
  original_model <- object@train.models[[model_type]]
  data_sets <- Extract_filtered.set(object)
  train_data <- data_sets$training
  test_data <- data_sets$testing

  cat("Generating ROC comparison plot...\n")
  roc_comparison <- plot_tuning_comparison(
    original_model = original_model,
    tuned_model = tuned_model,
    test_data = train_data,
    group_col = group_col,
    palette_name = palette_name,
    plot_title = plot_title,
    base_size = base_size,
    save_plots = save_plots,
    save_dir = save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    alpha = alpha
  )

  if (inherits(object, "Train_Model")) {
    object@best.model.result[["roc_comparison"]] <- roc_comparison
    cat("\nUpdating 'Train_Model' object...\n")
    cat("The following slots have been updated:\n")
    cat("- 'best.model.result$roc_comparison' (ROC comparison results)\n")

  }

  return(object)
}

