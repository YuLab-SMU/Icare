#' Create Forest Plot
#'
#' @param hr_results HR results.
#' @param plot_title Plot title.
#' @param save_plot Save plot.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param palette_name Palette.
#' @param base_size Base size.
#' @param hr_limit HR limit.
#' @param ci_range_limit CI range limit.
#' @import ggplot2 dplyr tidyr ggtext patchwork here
#' @export
create_forest_plot <- function(hr_results,
                               plot_title = "Univariate Analysis: Clinical Factors and Hazard Ratios with 95% Confidence Intervals",
                               save_plot = FALSE,
                               save_dir = here('SurObj', "univariate_analysis"),
                               plot_width = 11,
                               plot_height = 5,
                               palette_name = "AsteroidCity1",
                               base_size = 14,
                               hr_limit = c(0, 3),
                               ci_range_limit = 1000) {
  hr_results$P_value <- round(hr_results$P_value, 3)
  
  
  hr_results <- hr_results %>%
    filter(is.finite(HR) & is.finite(CI_lower) & is.finite(CI_upper)) %>%
    filter(HR >= hr_limit[1] & HR <= hr_limit[2]) %>%
    filter((CI_upper - CI_lower) <= ci_range_limit)
  
  cat(paste(nrow(hr_results), "rows remain after filtering. Rows with extreme or invalid values were excluded.\n"))
  
  tab_header <- c("Clinical factors", "HR (95% CI)", "P-value")
  tab_header_bold <- TRUE
  
  
  tmp_df <- hr_results %>%
    dplyr::select(Variable, HR_95CI, P_value) %>%
    rbind(tab_header, .) %>%
    cbind("clinical" = .[, "Variable"], .)
  
  tmp_df$clinical <- factor(tmp_df$clinical, levels = rev(tmp_df$clinical))
  
  tmp_df <- tmp_df %>%
    pivot_longer(cols = 2:ncol(.), names_to = "x", values_to = "label")
  
  
  tmp_df$label_bold <- sapply(tmp_df$label, function(x) {
    if (x %in% unlist(tab_header)) {
      if (tab_header_bold) {
        paste0("<b>", x, "</b>")
      } else {
        x
      }
    } else {
      x
    }
  }, simplify = T)
  
  p_tab <- ggplot(data = tmp_df, aes(x = x, y = clinical)) +
    geom_tile(color = "white", fill = "white") +
    geom_richtext(aes(label = label_bold), label.color = NA) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line.x = element_line(linewidth = 1),
      axis.line.y = element_blank(),
      text = element_text()
    ) +
    geom_hline(yintercept = c(length(unique(tmp_df$clinical)) - 0.5), linewidth = 1)
  
  tmp_error <- rbind(rep(NA, ncol(hr_results)), hr_results)
  tmp_error$Variable[1] <- tab_header[[1]]
  
  tmp_error$Variable <- factor(tmp_error$Variable, levels = rev(tmp_error$Variable))
  error_bar_height <- 0.2
  ref_line <- 1
  
  p_errorbar <- ggplot(data = tmp_error) +
    geom_point(aes(x = HR, y = Variable, color = HR > 1), shape = "diamond", size = 4) +
    scale_color_manual(values = wes_palette(palette_name)) +
    geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper, y = Variable), linewidth = 1, height = error_bar_height) +
    theme_classic(base_size = base_size) +
    theme(
      title = element_blank(),
      axis.line.x = element_line(linewidth = 1),
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
    ) +
    expand_limits(x = hr_limit) +
    coord_cartesian(xlim = hr_limit) +
    geom_vline(xintercept = ref_line, color = "grey50")
  
  p <- p_tab + p_errorbar +
    plot_annotation(title = plot_title,
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  
  print(p)
  
  if (save_plot) {
    pdf_path <- file.path(save_dir, "forest_plot.pdf")
    ggsave(filename = pdf_path, plot = p, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved at: ", pdf_path, "\n")
  }
}

#' Forest Plot Univariate Analysis
#'
#' @param object SurObj.
#' @param plot_title Title.
#' @param time_col Time col.
#' @param status_col Status col.
#' @param var_col Var col.
#' @param palette_name Palette.
#' @param save_plot Save plot.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param use_subgroup_data Use subgroup.
#' @param hr_limit HR limit.
#' @param result_type Result type.
#' @export
forest_plot_univariate_analysis <- function(object,
                                            plot_title = "Univariate Analysis: Clinical Factors and Hazard Ratios with 95% Confidence Intervals",
                                            time_col = "time",
                                            status_col = "status",
                                            var_col = NULL,
                                            palette_name = "AsteroidCity1",
                                            save_plot = TRUE,
                                            save_dir = here('SurObj', "univariate_analysis"),
                                            plot_width = 14,
                                            plot_height = 7,
                                            base_size = 14,
                                            use_subgroup_data = FALSE,
                                            hr_limit = c(0.1, 3),
                                            result_type = "all") {  # Added result_type parameter
  
  if (inherits(object, 'SurObj')) {
    all_results <- slot(object, "univariate.analysis")[["all_univariate_results"]]
    
    if (result_type == "all") {
      hr_results <- all_results$all_results
    } else if (result_type == "significant") {
      hr_results <- all_results$significant_results
    } else {
      stop("Invalid result_type. Please choose 'all' or 'significant'.")
    }
    
    if (is.null(hr_results) || nrow(hr_results) == 0) {
      stop("The hr_results in the SurObj object is empty.")
    }
    
    cat("Using SurObj object for univariate analysis data...\n")
    
  } else if (is.data.frame(object)) {
    hr_results <- object
    if (nrow(hr_results) == 0) {
      stop("The provided data frame is empty.")
    }
    cat("Using provided data frame for univariate analysis data...\n")
    
  } else {
    stop("Input must be an object of class 'SurObj' or a data frame.")
  }
  
  p <- create_forest_plot(hr_results,
                          hr_limit=hr_limit,
                          plot_title = plot_title,
                          save_plot = save_plot,
                          save_dir = save_dir,
                          plot_width = plot_width,
                          plot_height = plot_height,
                          base_size = base_size)
  
  if (inherits(object, 'SurObj')) {
    cat("Updating 'SurObj' object with forest plot...\n")
    object@univariate.analysis[["forest_plot"]] <- p
    cat("The 'SurObj' object has been updated with the following slots:\n")
    cat("- 'univariate.analysis' slot updated.\n")
    return(object)
  }
  
  return(p)
}
