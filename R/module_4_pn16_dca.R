#' Plot DCA Without Save
#'
#' Plots Decision Curve Analysis without saving.
#'
#' @param model_results List of model results.
#' @param palette_name Palette name.
#' @param type Type (unused).
#' @param status_col Status column.
#' @param time_col Time column.
#'
#' @return DCA results list.
#' @export
plot_dca_without_save <- function(model_results,
                                  palette_name ="Royal1",
                                  type = "continuous",
                                  status_col = "status",
                                  time_col = "time") {
  
  
  if (is.null(model_results) || length(model_results) == 0) {
    stop("No model results provided.")
  }
  
  dca_results <- list()
  
  
  for (model_name in names(model_results)) {
    model_data <- model_results[[model_name]]$test_data_risk
    
    
    if (length(unique(model_data[[status_col]])) <= 1 ||
        length(unique(model_data$risk_group)) <= 1) {
      warning(paste("Skipping model:", model_name, "- insufficient levels in status or risk_group"))
      next
    }
    
    model_data[[status_col]] <- as.numeric(model_data[[status_col]])
    model_data$risk_group <- factor(model_data$risk_group)
    
    dca_results[[model_name]] <- rmda::decision_curve(
      formula = as.formula(paste(status_col, "~ risk_score + risk_group")),
      data = model_data,
      family = binomial(link = "logit"),
      thresholds = seq(0, 1, by = 0.01),
      confidence.intervals = 0.95,
      study.design = "cohort"
    )
  }
  
  
  if (length(dca_results) == 0) {
    stop("No valid decision curve results available.")
  }
  
  
  palette_colors <- wesanderson::wes_palette(name = palette_name, type = "discrete")
  
  
  p <- rmda::plot_decision_curve(
    dca_results,
    curve.names = names(dca_results),
    col = palette_colors,
    lwd = 2,
    xlim = c(0, 1),
    ylim = c(0, 1),
    xlab = "Risk Threshold",
    ylab = "Net Benefit",
    legend.position = "bottomleft",
    confidence.intervals = FALSE,
    main = "Decision Curve Analysis"
  )
  
  
  return(dca_results)
}


#' Plot DCA for Best Model
#'
#' Plots DCA for the best model in PrognosiX object.
#'
#' @param object PrognosiX object.
#' @param palette_name Palette.
#' @param status_col Status column.
#' @param time_col Time column.
#'
#' @return Updated PrognosiX object.
#' @export
plot_dca_best_model <- function(object,
                                palette_name ="Royal1",
                                status_col = "status",  
                                time_col = "time") {    
  
  best_model <- object@best.model[["model"]]
  best_model_name <- object@best.model[["model_name"]]
  status_col <- object@status_col
  time_col <- object@time_col
  cat(paste("Plotting DCA for best model:", best_model_name, "\n"))
  
  model_results <- list()
  model_results[[best_model_name]] <- list(
    model = best_model,
    test_data_risk = object@survival.model[[best_model_name]]$results_km$data_risk
  )
  
  dca_re <- plot_dca_without_save(
    model_results = model_results,
    palette_name = palette_name,
    status_col = status_col,  
    time_col = time_col       
  )
  
  cat("DCA plot generated successfully for the best model.\n")
  if (inherits(object, 'PrognosiX')) {  
    object@best.model[["dca_results"]] <- dca_re
    
    cat("Updating 'PrognosiX' object...\n")
    cat("The 'PrognosiX' object has been updated with the following slots:\n")
    cat("- 'best.model' slot updated.\n")
    return(object)
  } else {
    warning("No valid best model found.")
  }
}
