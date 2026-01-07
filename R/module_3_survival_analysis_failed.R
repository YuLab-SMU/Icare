#' Survival Analysis Var Plot
#'
#' @param data Data.
#' @param time_col Time col.
#' @param status_col Status col.
#' @param var_col Var col.
#' @param palette_name Palette.
#' @param save_plots Save plots.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param break.x.by Break x by.
#' @param pval.method Pval method.
#' @param conf.int Conf int.
#' @param title Title.
#' @param surv.median.line Median line.
#' @param risk.table Risk table.
#' @param xlab Xlab.
#' @param legend.position Legend pos.
#' @param legend.title Legend title.
#' @import survival survminer ggplot2 gridExtra
#' @export
survival_analysis_var_plot  <- function(data,
                                       time_col = "time",
                                       status_col = "status",
                                       var_col = NULL,
                                       palette_name = "Dark2",
                                       save_plots = TRUE,
                                       save_dir = here('SurObj', "univariate_analysis"),
                                       plot_width = 6,
                                       plot_height = 8,
                                       base_size = 14,
                                       break.x.by = NULL,
                                       pval.method = TRUE,
                                       conf.int = FALSE,
                                       title = NULL,
                                       surv.median.line = "hv",
                                       risk.table = TRUE,
                                       xlab = "Follow up time",
                                       legend.position = c(0.8, 0.75),
                                       legend.title = NULL) {
  # 参数验证和数据处理
  cat("Starting survival analysis...\n")
  
  # 检查必要列是否存在
  required_cols <- c(var_col, time_col, status_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing:", paste(missing_cols, collapse = ", ")))
  }

  
  if (nrow(data) == 0) {
    stop("Filtered data has no valid rows after removing NA values.")
  }
  
  cat("Data validation successful. Number of valid rows:", nrow(data), "\n")
  data <- data %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[status_col]]),
      !is.na(.data[[var_col]])
    ) %>%
    mutate(
      !!status_col := as.numeric(.data[[status_col]]),
      !!var_col := factor(.data[[var_col]])  # 强制转换为因子
    )
  if (nlevels(data[[var_col]]) < 2) {
  stop("分组变量需至少包含2个水平！")
}
  # 转换状态列为数值型
  data[[status_col]] <- as.numeric(data[[status_col]], levels = c(0, 1))
  
  # 拟合生存模型

  fit <- survfit(as.formula(paste("Surv(", time_col, ",", status_col, ") ~", var_col)), data = data)
  
  
  cat("Survival model fitted successfully.\n")
  
  # 自动确定时间间隔（如果没有提供）
  if (is.null(break.x.by)) {
    max_time <- max(data[[time_col]], na.rm = TRUE)
    break.x.by <- case_when(
      max_time > 365 ~ 365,
      max_time > 100 ~ 50,
      max_time > 50 ~ 10,
      TRUE ~ 5
    )
    cat("break.x.by automatically set to:", break.x.by, "\n")
  }
  
  # 设置默认标题和图例标题
  if (is.null(title)) {
    title <- sprintf("Kaplan-Meier Survival Curve by %s", var_col)
  }
  
  if (is.null(legend.title)) {
    legend.title <- var_col
  }

  # 生成生存曲线图
  km_plot <- ggsurvplot(
    fit,
    data = data,
    conf.int = conf.int,
    pval = TRUE,
    pval.method = pval.method,
    title = title,
    surv.median.line = surv.median.line,
    risk.table = risk.table,
    xlab = xlab,
    legend = legend.position,
    legend.title = legend.title,
    legend.labs = unique(data[[var_col]]),
    break.x.by = break.x.by,
    palette = palette_name,
    risk.table.height = 0.25
  )
  
  cat("Base survival curve plotted.\n")
  
  # 增强绘图主题
  surv_plot <- km_plot$plot +
    theme_prism(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = base_size),
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.8),
      axis.text.y = element_text(size = base_size * 0.8),
      axis.title = element_text(size = base_size),
      legend.position = legend.position,
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size * 0.8)
    )
  
  # 增强风险表主题
  risk_table <- km_plot$table + 
    theme_prism(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = base_size * 0.8),
      axis.text.y = element_text(size = base_size * 0.8),
      axis.title = element_text(size = base_size)
    )
  
  # 组合图形
  combined_plot <- arrangeGrob(
    surv_plot, 
    risk_table,
    ncol = 1,
    heights = c(2, 0.7)
  )
  
  cat("Combined plot created.\n")
  print(combined_plot)
  # 保存图形
  if (save_plots) {
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
      cat("Created directory:", save_dir, "\n")
    }
    
    # 保存组合图
    ggsave(
      filename = file.path(save_dir, paste0("survival_combined_by_", var_col, ".pdf")),
      plot = combined_plot,
      width = plot_width,
      height = plot_height,
      device = "pdf"
    )
    
    # 单独保存生存曲线和风险表
    ggsave(
      filename = file.path(save_dir, paste0("survival_curve_by_", var_col, ".pdf")),
      plot = surv_plot,
      width = plot_width,
      height = plot_height,
      device = "pdf"
    )
    
    ggsave(
      filename = file.path(save_dir, paste0("risk_table_by_", var_col, ".pdf")),
      plot = risk_table,
      width = plot_width,
      height = plot_height * 0.5,
      device = "pdf"
    )
    
    cat("Plots saved to:", save_dir, "\n")
  }
  
  cat("Survival analysis completed successfully.\n")
  
  return(list(
    km_plot = km_plot,
    surv_plot = surv_plot,
    risk_table = risk_table,
    combined_plot = combined_plot
  ))
}



#' Plot Var Kaplan Meier
#'
#' @param object Object.
#' @param time_col Time col.
#' @param status_col Status col.
#' @param var_col Var col.
#' @param palette_name Palette.
#' @param save_plots Save plots.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param break.x.by Break x by.
#' @param pval.method Pval method.
#' @param conf.int Conf int.
#' @param title Title.
#' @param surv.median.line Median line.
#' @param risk.table Risk table.
#' @param xlab Xlab.
#' @param legend.position Legend pos.
#' @param legend.title Legend title.
#' @param use_subgroup_data Use subgroup.
#' @export
plot_var_kaplan_meier <- function(object,
                                  time_col = "time",
                                  status_col = "status",
                                  var_col = NULL,
                                  palette_name = "Dark2",
                                  save_plots = TRUE,
                                  save_dir = here('SurObj', "univariate_analysis"),
                                  plot_width = 5,
                                  plot_height = 5,
                                  base_size = 14,
                                  break.x.by = NULL,
                                  pval.method = TRUE,
                                  conf.int = FALSE,
                                  title = NULL,
                                  surv.median.line = "hv",
                                  risk.table = TRUE,
                                  xlab = "Follow up time",
                                  legend.position = c(0.8, 0.75),
                                  legend.title = NULL,
                                  use_subgroup_data = FALSE) {
  
  # 数据源选择
  if (inherits(object, 'SurObj')) {
    data <- if (use_subgroup_data) {
      cat("Using subgroup analysis data...\n")
      slot(object, "sub.data")
    } else {
      cat("Using original survival data...\n")
      slot(object, "survival.data")
    }
    
    if (is.null(data) || nrow(data) == 0) {
      stop("The selected data in the SurObj object is empty.")
    }
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using the provided data...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'SurObj' or a data frame.")
  }
  
  # 变量检查
  if (is.null(var_col) || !(var_col %in% colnames(data))) {
    stop("The specified var_col is either NULL or not found in the dataset.")
  }
  
  cat("Starting Kaplan-Meier plotting for variable:", var_col, "\n")
  
  # 转换数据为数值型
  data <- convert_to_numeric(data)
  
  # 执行生存分析
  analysis_results <- survival_analysis_var_plot(
    data = data,
    time_col = "time",
    status_col = "status",
    var_col = var_col,
    palette_name = palette_name,
    save_plots =save_plots,
    save_dir =save_dir,
    plot_width = plot_width,
    plot_height = plot_height,
    base_size = base_size,
    break.x.by = break.x.by,
    pval.method = pval.method,
    conf.int = conf.int,
    title = title,
    surv.median.line =surv.median.line,
    risk.table = risk.table,
    xlab = xlab,
    legend.position = legend.position,
    legend.title = legend.title
  )
  
  # 更新SurObj对象（如果输入是SurObj）
  if (inherits(object, 'SurObj')) {
    cat("Updating 'SurObj' object with analysis results...\n")
    result_col <- paste0("km_results_", var_col)
    object@univariate.analysis[["km_results"]][[result_col]] <- analysis_results
    
    cat("The 'SurObj' object has been updated with the following:\n")
    cat("- Added Kaplan-Meier results for variable:", var_col, "\n")
    return(object)
  }
  
  cat("Kaplan-Meier analysis completed for variable:", var_col, "\n")
  return(analysis_results)
}




#' Plot Survival Time Distribution
#'
#' @param data Data.
#' @param time_col Time col.
#' @param var_col Var col.
#' @param palette_name Palette.
#' @param binwidth Binwidth.
#' @param save_plots Save plots.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @export
plot_survival_time_distribution <- function(data,
                                            time_col,
                                            var_col,
                                            palette_name = "AsteroidCity1",
                                            binwidth = 10,
                                            save_plots = TRUE,
                                            save_dir = here('SurObj', "univariate_analysis"),
                                            plot_width = 5,
                                            plot_height = 5,
                                            base_size = 14) {
  data[[var_col]] <- as.factor(data[[var_col]])
  
  mean_df <- data %>%
    group_by(!!sym(var_col)) %>%
    summarise(mean_time = mean(!!sym(time_col), na.rm = TRUE),
              se_time = sd(!!sym(time_col), na.rm = TRUE) / sqrt(n())) %>%
    ungroup()
  
  label_df <- data %>%
    group_by(!!sym(var_col)) %>%
    summarise(label = paste0(round(mean(!!sym(time_col), na.rm = TRUE), 1), "±",
                             round(sd(!!sym(time_col), na.rm = TRUE) / sqrt(n()), 1), " days")) %>%
    ungroup() %>%
    mutate(y = seq(12, 12 - 5 * (n() - 1), length.out = n()))
  
  colors <- wes_palette(palette_name, type = "discrete")
  p <- ggplot(data, aes(x = !!sym(time_col), fill = !!sym(var_col))) +
    geom_histogram(binwidth = binwidth, color = "black", alpha = 0.7) +
    geom_vline(data = mean_df, aes(xintercept = mean_time), color = colors[4], linetype = "solid", size = 1) +
    geom_vline(data = mean_df, aes(xintercept = mean_time - se_time), color = colors[2], linetype = "dashed", size = 0.6) +
    geom_vline(data = mean_df, aes(xintercept = mean_time + se_time), color = colors[2], linetype = "dashed", size = 0.6) +
    scale_fill_manual(values = wes_palette(palette_name)) +
    labs(title = paste("Distribution of Survival Time by", var_col),
         x = "Survival Time (days)",
         y = "Frequency") +
    facet_grid(paste0(var_col, " ~ .")) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.title.y = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none",
      strip.background = element_rect(fill = "white"),
      panel.grid = element_line(linetype = "dotted", color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line.y.left = element_line(color = "black", size = 0.5),
      axis.line.y.right = element_line(color = "black", size = 0.5))
  
  cat("Plotting survival time distribution...\n")
  print(p)
  
  if (save_plots) {
    
    ggsave(file.path(save_dir, paste0("survival_distribution_by_", var_col,".pdf")),
           plot = p, width = plot_width, height = plot_height,device = "pdf")
    cat("Plot saved to:", file.path(save_dir, paste0("survival_distribution_by_", var_col,".pdf"), "\n"))
  }
  
  return(p)
}



#' Plot Var Survival Time
#'
#' @param object Object.
#' @param time_col Time col.
#' @param status_col Status col.
#' @param var_col Var col.
#' @param binwidth Binwidth.
#' @param palette_name Palette.
#' @param save_plots Save plots.
#' @param save_dir Save dir.
#' @param plot_width Width.
#' @param plot_height Height.
#' @param base_size Base size.
#' @param use_subgroup_data Use subgroup.
#' @export
plot_var_survival_time <- function(object,
                                   time_col = "time",
                                   status_col = "status",
                                   var_col=NULL,
                                   binwidth = 10,
                                   palette_name = "AsteroidCity1",
                                   save_plots = TRUE,
                                   save_dir = here('SurObj', "univariate_analysis"),
                                   plot_width = 5,
                                   plot_height = 5,
                                   base_size = 14,
                                   use_subgroup_data = FALSE) {
  
  if (inherits(object, 'SurObj')) {
    if (use_subgroup_data) {
      data <- slot(object, "sub.data")
      cat("Using subgroup analysis data...\n")
    } else {
      data <- slot(object, "survival.data")
      cat("Using original survival data...\n")
    }
    
    if (is.null(data) || nrow(data) == 0) {
      stop("The survival.data in the SurObj object is empty.")
    }
  } else if (is.data.frame(object)) {
    cat("Input is a data frame. Using the provided data...\n")
    data <- object
  } else {
    stop("Input must be an object of class 'SurObj' or a data frame.")
  }
  data <- convert_to_numeric(data)
  
  p <- plot_survival_time_distribution(data=data,
                                       time_col=time_col,
                                       var_col=var_col,
                                       palette_name=palette_name,
                                       binwidth=binwidth,
                                       save_plots=save_plots,
                                       save_dir=save_dir,
                                       plot_width=plot_width,
                                       plot_height=plot_height,
                                       base_size=base_size)
  result_col<- paste0("distribution_time_",var_col)
  if (inherits(object, 'SurObj')) {
    cat("Updating 'SurObj' object...\n")
    
    object@univariate.analysis[["distribution_time"]][[result_col]]<-  p
    cat("The 'SurObj' object has been updated with the following slots:\n")
    cat("- 'univariate.analysis' slot updated.\n")
    return(object)
  }
  
  cat("Plotting function execution completed.\n")
  return(p)
}
