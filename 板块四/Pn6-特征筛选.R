library(glmnet)
library(survival)
lasso_feature_selection <- function(data,
                                    time_col,
                                    status_col,
                                    alpha = 1,
                                    seed = 123) {

  data[[status_col]] <- as.factor(data[[status_col]])
  x <- as.matrix(data[, !(names(data) %in% c(time_col, status_col))])
  y <- data[[status_col]]

  set.seed(seed)

  lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = alpha)
  best_lambda <- lasso_model$lambda.1se


  lasso_coef <- as.matrix(coef(lasso_model, s = best_lambda))
  important_vars <- rownames(lasso_coef)[lasso_coef != 0]
  important_vars <- important_vars[important_vars != "(Intercept)"]

  return(list(
    model = lasso_model,
    best_lambda = best_lambda,
    important_vars = important_vars,
    cv_fit = lasso_model,
    lasso_coef = lasso_coef
  ))
}



run_lasso_feature_selection <- function(object,
                                        time_col = "time",
                                        status_col = "status",
                                        group_col ="group",
                                        alpha = 1) {

  cat("Starting LASSO feature selection...\n")

  if (inherits(object, 'PrognosiX')) {
    split_data <- slot(object, "split.data")
    train_data <- split_data$train
    time_col <- slot(object, "time_col")
    status_col <- slot(object, "status_col")

  } else if (is.data.frame(object)) {
    train_data <- object
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }

  if (is.null(train_data) || nrow(train_data) == 0) {
    stop("Training data is empty. Please ensure it contains the necessary columns.")
  }
  if(!is.null(group_col) && group_col %in% colnames(train_data)){
    train_data[[group_col]] <- NULL
  }
  results <- lasso_feature_selection(data=train_data, 
                                     time_col=time_col, 
                                     status_col=status_col, 
                                     alpha=alpha)

  cat("LASSO feature selection completed.\n")
  cat("Identified important variables (total: ", length(results$important_vars), "):\n")
  print(results$important_vars)

  if (inherits(object, 'PrognosiX')) {
    object@feature.result <- list(lasso_results = results)
    object@feature.result[["important_vars"]]<-results$important_vars
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'feature.result' slot updated.\n")

    return(object)
  }

  return(results)
}


library(ggplot2)
library(here)

lasso_cv_visualization <- function(object,
                                   palette_name = "AsteroidCity1",
                                   base_size = 14,
                                   save_plots = TRUE,
                                   save_dir = here("PrognosiX", "Sel_feature"),
                                   plot_width = 5,
                                   plot_height = 5) {

  cat("Checking the input object...\n")
  if (inherits(object, 'PrognosiX')) {
    lasso_model <- slot(object, "feature.result")[["lasso_results"]][["model"]]
    cat("Input is a 'PrognosiX' object.\n")
  } else if (is.list(object)) {
    lasso_model <- object
    cat("Input is a list object.\n")
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }

  if (is.null(lasso_model$lambda) || is.null(lasso_model$cvm) || is.null(lasso_model$cvsd)) {
    stop("Invalid LASSO model: missing lambda, cvm, or cvsd.")
  }
  cat("LASSO model contains lambda, cvm, and cvsd values.\n")

  cv_data <- data.frame(Lambda = lasso_model$lambda,
                        MeanSquaredError = lasso_model$cvm,
                        StandardError = lasso_model$cvsd)

  colors <- as.vector(wes_palette(palette_name))
  primary_color <- colors[1]
  secondary_color <- colors[2]

  p <- ggplot(cv_data, aes(log(Lambda), MeanSquaredError)) +
    geom_line() +
    geom_ribbon(aes(ymin = MeanSquaredError - StandardError, ymax = MeanSquaredError + StandardError), alpha = 0.2) +
    geom_vline(xintercept = log(lasso_model$lambda.min), linetype = "dashed", color = primary_color) +
    geom_vline(xintercept = log(lasso_model$lambda.1se), linetype = "dashed", color = secondary_color) +
    scale_x_reverse() +
    labs(x = "Log(Lambda)", y = "Mean Squared Error") +
    theme_minimal() +
    theme_prism(base_size = base_size)

  print(p)

  if (save_plots) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    plot_path <- file.path(save_dir, "lasso_path_visualization.pdf")
    cat("Saving plot to:", plot_path, "\n")
    tryCatch({
      ggsave(filename = plot_path,
             plot = p,
             width = plot_width,
             height = plot_height,
             device = "pdf")
      cat("Plot saved successfully.\n")
    }, error = function(e) {
      cat("Error saving plot: ", e$message, "\n")
    })
  }

  if (inherits(object, 'PrognosiX')) {
    object@feature.result[["lasso_cv_plot"]] <- p
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'feature.result' slot updated.\n")
    return(object)
  }

  return(p)
}


lasso_feature_importance <- function(object,
                                     palette_name = "AsteroidCity1",
                                     base_size = 14,
                                     save_plots = TRUE,
                                     save_dir = here("PrognosiX", "Sel_feature"),
                                     plot_width = 5,
                                     plot_height = 5) {

  cat("Checking the input object...\n")
  if (inherits(object, 'PrognosiX')) {
    lasso_model <- slot(object, "feature.result")[["lasso_results"]][["model"]]
    cat("Input is a 'PrognosiX' object.\n")
  } else if (is.list(object)) {
    lasso_model <- object
    cat("Input is a list object.\n")
  } else {
    stop("Input must be an object of class 'PrognosiX' or a data frame.")
  }

  lasso_coef <- coef(lasso_model, s = "lambda.min")
  lasso_coef_df <- as.data.frame(as.matrix(lasso_coef))

  colnames(lasso_coef_df) <- "Coefficient"
  lasso_coef_df$Feature <- rownames(lasso_coef_df)

  lasso_coef_df <- lasso_coef_df[lasso_coef_df$Coefficient != 0, ]
  lasso_coef_df <- lasso_coef_df[order(abs(lasso_coef_df$Coefficient), decreasing = TRUE), ]

  cat("Creating feature importance plot...\n")
  p <- ggplot(lasso_coef_df, aes(x = Feature, y = Coefficient)) +
    geom_segment(aes(x = Feature, xend = Feature, y = 0, yend = Coefficient), color = "grey50") +
    geom_point(aes(color = Coefficient > 0), size = 3) +
    scale_color_manual(values =wes_palette(palette_name)) +
    coord_flip() +
    theme_bw() +
    labs(title = "LASSO Model Feature Importance",
         y = "Coefficient",  x = "") +
    theme(legend.title = element_blank(),
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12)) +
    theme_prism(base_size = base_size)

  print(p)

  if (save_plots) {
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
    plot_path <- file.path(save_dir, "lasso_feature_importance.pdf")
    cat("Saving plot to:", plot_path, "\n")
    tryCatch({
      ggsave(filename = plot_path,
             plot = p,
             width = plot_width,
             height = plot_height,
             device = "pdf")
      cat("Plot saved successfully.\n")
    }, error = function(e) {
      cat("Error saving plot: ", e$message, "\n")
    })
  }

  if (inherits(object, 'PrognosiX')) {
    object@feature.result[["lasso_feature_plot"]] <- p
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'feature.result' slot updated.\n")

    return(object)
  }

  return(p)
}

