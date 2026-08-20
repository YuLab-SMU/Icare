## ============================================================
##  Publication-quality visualization functions
##  Organized by module:
##     1. Stat       -- distribution, correlation, PCA, DEG
##     2. Train_Model -- ROC, confusion matrix, feature importance, calibration
##     3. Subtyping  -- dim-reduction, cluster heatmap, silhouette, alluvial
##     4. PrognosiX  -- KM, forest, time-ROC, RCS, nomogram, calibration, DCA, risk
##
##  Every function:
##    - accepts an S4 object OR raw data frames
##    - returns a ggplot2 / grob object (further customisable)
##    - saves to PDF / PNG / SVG when save_plot = TRUE
##    - uses theme_prism + wesanderson / RColorBrewer palettes by default
## ============================================================

## -- shared helpers ------------------------------------------------------------

.pub_theme <- function(base_size = 13) {
  ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                            size  = base_size + 1),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, colour = "grey40"),
      axis.title    = ggplot2::element_text(face  = "bold"),
      legend.title  = ggplot2::element_text(face  = "bold"),
      strip.text    = ggplot2::element_text(face  = "bold")
    )
}

.get_palette <- function(palette_name, n) {
  tryCatch(
    as.character(wesanderson::wes_palette(n = n, name = palette_name,
                                          type = if (n > 5) "continuous" else "discrete")),
    error = function(e)
      RColorBrewer::brewer.pal(max(3L, n), palette_name)[seq_len(n)]
  )
}

.save_plot <- function(p, dir, filename, width, height, format = "pdf") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, paste0(tools::file_path_sans_ext(filename), ".", format))
  ggplot2::ggsave(path, plot = p, width = width, height = height, dpi = 300)
  cat("Plot saved to:", path, "\n")
  invisible(path)
}

#' Get Output Directory for Visualization Functions
#' 
#' Internal helper to get the appropriate output directory based on global config
#' @keywords internal
.get_viz_output_dir <- function(type = c("Stat", "Model", "Subtyping", "PrognosiX")) {
  type <- match.arg(type)
  root <- get_output_root()
  if (is.null(root)) {
    viz_dir <- here::here("Figures", type)
  } else {
    viz_dir <- file.path(root, "Figures", type)
  }
  # Create directory on-demand
  if (!dir.exists(viz_dir)) {
    dir.create(viz_dir, recursive = TRUE, showWarnings = FALSE)
  }
  return(viz_dir)
}

## =============================================================================
##   1  STAT  -- distribution - correlation - PCA - DEG
## =============================================================================

#' Grouped Distribution Plot with Violin and Boxplot
#'
#' This function generates a faceted grouped distribution plot for selected numeric features.
#' It visualizes the distribution of each feature across groups using violin plots overlaid
#' with boxplots, and optionally adds statistical significance labels (p-values).
#'
#' The function converts the grouping column to character internally to ensure that
#' custom color names (provided via `group_colors`) can match the group levels correctly.
#'
#' @param object A `Stat` object (from the package) or a data frame containing the data.
#'   If a `Stat` object is provided, the function will extract the clean data and
#'   the default group column from the object.
#' @param features A character vector of feature (column) names to plot. If `NULL`
#'   (default), the first 12 numeric columns (excluding the group column) are used.
#' @param group_col A character string specifying the column name that defines the
#'   grouping variable. If `NULL` (default) and `object` is a `Stat` object, it will
#'   be automatically taken from `object@group_col`. Otherwise, it must be provided.
#' @param test A character string specifying the statistical test to use for pairwise
#'   comparisons. Passed to \code{ggpubr::stat_compare_means()}. Default is `"wilcox.test"`.
#' @param palette_name A character string indicating the color palette to use when
#'   `group_colors` is `NULL`. Supports RColorBrewer palettes (e.g., "Set1") or
#'   wesanderson palettes (e.g., "Zissou1"). Default is `"Royal1"`.
#' @param group_colors An optional named character vector specifying custom colors
#'   for each group. The names must exactly match the group levels (case-sensitive).
#'   If provided, this overrides `palette_name`. Default is `NULL`.
#' @param ncol Number of columns in the faceted plot. Default is 3.
#' @param base_size Base font size for the plot theme. Default is 13.
#' @param save_plot Logical; if `TRUE`, the plot is saved to disk. Default is `FALSE`.
#' @param save_dir Directory path where the plot will be saved. If `NULL` and
#'   `save_plot` is `TRUE`, a default directory is used.
#' @param width Width of the saved plot (in inches). Default is 10.
#' @param height Height of the saved plot (in inches). Default is 6.
#' @param format Output format for the saved plot (`"pdf"` or `"png"`). Default is `"pdf"`.
#'
#' @return A `ggplot` object. The plot is also saved if `save_plot = TRUE`.
#' @export
#'
#' @examples
#' \dontrun{
#' # Using a Stat object with custom colors
#' p <- PlotGroupedDistribution(
#'   object = stat_obj,
#'   features = c("GeneA", "GeneB"),
#'   group_colors = c("0" = "#404040", "1" = "#ca0020")
#' )
#' print(p)
#'
#' # Using a data frame
#' df <- data.frame(
#'   Group = rep(c("A", "B"), each = 20),
#'   Feature1 = rnorm(40),
#'   Feature2 = rnorm(40, mean = 2)
#' )
#' p2 <- PlotGroupedDistribution(df, group_col = "Group")
#' }
PlotGroupedDistribution <- function(object,
                                    features     = NULL,
                                    group_col    = NULL,
                                    test         = "wilcox.test",
                                    palette_name = "Royal1",
                                    group_colors = NULL,
                                    ncol         = 3,
                                    base_size    = 13,
                                    save_plot    = FALSE,
                                    save_dir     = NULL,
                                    width        = 10,
                                    height       = 6,
                                    format       = "pdf") {
  cat("Generating grouped distribution plot...\n")
  if (is.null(save_dir) && save_plot) {
    save_dir <- .get_viz_output_dir("Stat")
  }
  
  # ---- Extract data and group column ----
  if (inherits(object, "Stat")) {
    df <- object@clean.data
    if (is.null(group_col)) group_col <- object@group_col
  } else {
    df <- as.data.frame(object)
    if (is.null(group_col)) {
      stop("'group_col' must be provided when 'object' is not a Stat object.")
    }
  }
  if (!group_col %in% colnames(df)) {
    stop("Group column '", group_col, "' not found in data.")
  }
  
  # ---- Convert group column to character for reliable matching ----
  df[[group_col]] <- as.character(df[[group_col]])
  
  # ---- Select numeric features ----
  num_cols <- names(df)[sapply(df, is.numeric)]
  num_cols <- setdiff(num_cols, group_col)
  if (is.null(features)) features <- head(num_cols, 12)
  features <- intersect(features, num_cols)
  if (length(features) == 0) {
    stop("No valid numeric features found.")
  }
  
  # ---- Determine group levels and colors ----
  groups <- unique(df[[group_col]])
  n_groups <- length(groups)
  
  if (!is.null(group_colors)) {
    # Check if named vector and match groups
    if (!is.null(names(group_colors))) {
      matched <- intersect(names(group_colors), groups)
      if (length(matched) == n_groups) {
        cols <- group_colors[groups]   # reorder to match group order
      } else {
        warning("Names in 'group_colors' do not match all group levels. Falling back to palette.")
        cols <- .get_palette(palette_name, n_groups)
      }
    } else {
      # Unnamed: assume order matches group levels
      if (length(group_colors) == n_groups) {
        cols <- group_colors
        names(cols) <- groups
        message("Using 'group_colors' in the order of group levels: ",
                paste(groups, collapse = ", "))
      } else {
        warning("Length of 'group_colors' does not match number of groups. Falling back to palette.")
        cols <- .get_palette(palette_name, n_groups)
      }
    }
  } else {
    cols <- .get_palette(palette_name, n_groups)
  }
  names(cols) <- groups   # ensure names
  
  # ---- Reshape data to long format ----
  long_df <- tidyr::pivot_longer(
    df[, c(group_col, features), drop = FALSE],
    cols = -dplyr::all_of(group_col),
    names_to = "Feature",
    values_to = "Value"
  )
  long_df[[group_col]] <- as.factor(long_df[[group_col]])
  
  # ---- Build the plot ----
  p <- ggplot2::ggplot(
    long_df,
    ggplot2::aes(
      x = .data[[group_col]],
      y = Value,
      fill = .data[[group_col]]
    )
  ) +
    ggplot2::geom_violin(alpha = 0.5, colour = NA) +
    ggplot2::geom_boxplot(
      width = 0.2,
      outlier.shape = 21,
      outlier.size = 1.2,
      colour = "grey30"
    ) +
    ggpubr::stat_compare_means(
      method = test,
      label = "p.signif",
      comparisons = utils::combn(levels(long_df[[group_col]]), 2, simplify = FALSE)
    ) +
    ggplot2::facet_wrap(~ Feature, scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = group_col, y = "Value", fill = group_col) +
    .pub_theme(base_size)   # internal theme function
  
  # ---- Save plot if requested ----
  if (save_plot) {
    .save_plot(p, save_dir, "grouped_distribution", width, height, format)   # internal helper
  }
  
  return(p)
}


#' Publication-ready correlation heatmap using the corrplot package (enhanced + stable)
#'
#' Generate a highly customisable correlation heatmap leveraging the `corrplot`
#' package. Supports multiple matrix shapes, visualisation methods, clustering
#' orders, significance masking, and rectangle drawing.
#'
#' @param object A Stat object or data frame.
#' @param features Features to include (NULL = all numeric, max 40).
#' @param method Correlation method: `"spearman"` (default) or `"pearson"`.
#' @param vis_method Visualisation method for the tiles:
#'   `"circle"`, `"square"`, `"ellipse"`, `"number"`, `"shade"`, `"color"`, `"pie"`.
#'   Default `"color"`.
#' @param matrix_type One of `"full"`, `"lower"`, `"upper"`. Default `"lower"`.
#' @param order Order of variables:
#'   `"original"`, `"AOE"`, `"FPC"`, `"hclust"`, `"alphabet"`. Default `"hclust"`.
#' @param hclust_method Agglomeration method for `order = "hclust"`:
#'   `"ward"`, `"ward.D"`, `"ward.D2"`, `"single"`, `"complete"`, `"average"`,
#'   `"mcquitty"`, `"median"`, `"centroid"`. Default `"complete"`.
#' @param addrect Integer, number of rectangles to draw around clusters (requires
#'   `order = "hclust"`). Default `NULL` (no rectangles).
#' @param rect.col Colour of rectangle borders. Default `"navy"`.
#' @param rect.lwd Line width of rectangle borders. Default `2`.
#' @param p_thresh P-value threshold for masking non-significant pairs.
#'   Only used when `insig != "n"`.
#' @param insig Character, how to handle insignificant correlations:
#'   `"blank"` (hide tile), `"p-value"` (print p-value), `"pch"` (print symbol),
#'   `"label_sig"` (print significance stars), `"n"` (do nothing). Default `"blank"`.
#' @param pch Symbol to use when `insig = "pch"`. Default `4` (x).
#' @param sig.level Significance level(s) for `insig = "label_sig"`. Can be a vector,
#'   e.g. `c(0.001, 0.01, 0.05)`. Default `0.05`.
#' @param palette_name RColorBrewer diverging palette name (used when `color_scheme = "brew"`).
#'   Default `"RdYlBu"`.
#' @param color_scheme Colour scheme: `"viridis"`, `"brew"`, `"gradient"`, `"scientific"`.
#'   Default `"scientific"`.
#' @param tl.cex Size of variable labels. Default `0.8`.
#' @param tl.srt Rotation angle of variable labels (degrees). Default `45`.
#' @param tl.col Text color for variable labels. Default `"black"`.
#' @param cl.cex Size of colour legend labels. Default `0.8`.
#' @param cl.pos Position of colour legend: `"r"` (right), `"b"` (bottom), `"n"` (none).
#'   Default `"b"` (bottom).
#' @param diag Logical, whether to show diagonal tiles. Default `FALSE`.
#' @param save_plot Logical. Save the plot? Default `FALSE`.
#' @param save_dir Output directory. When `NULL` and `save_plot = TRUE`,
#'   the default Stat figure directory is used.
#' @param width,height Plot dimensions (inches). Default `8` x `7`.
#' @param format File format (`"pdf"`, `"png"`, `"svg"`). Default `"pdf"`.
#' @param min_obs Minimum number of observations for pairwise correlations. Default `3`.
#'
#' @return Invisibly returns the correlation matrix used for plotting.
#'   The plot is drawn on the current graphics device.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Lower triangle with hierarchical clustering and rectangles
#' # (Color legend now at bottom by default)
#' PlotCorrelationHeatmap(stat_obj_test, matrix_type = "lower", order = "hclust",
#' addrect = 3, color_scheme = "scientific")
#'
#' # Upper triangle with ellipse method and significance stars
#' # Black text labels for better readability
#' PlotCorrelationHeatmap(stat_obj_test, matrix_type = "lower", vis_method = "ellipse",
#' insig = "label_sig", sig.level = c(0.001, 0.01, 0.05))
#'
#' # Custom text color if needed
#' PlotCorrelationHeatmap(stat_obj_test, matrix_type = "upper",tl.col = "navy", cl.pos = "r")
#' PlotCorrelationHeatmap(stat_obj_test, matrix_type = "full",tl.col = "navy", cl.pos = "r")
#' }
PlotCorrelationHeatmap <- function(object,
                                   features      = NULL,
                                   method        = "spearman",
                                   vis_method    = c("color", "circle", "square", "ellipse",
                                                     "number", "shade", "pie"),
                                   matrix_type   = c("lower", "full", "upper"),
                                   order         = c("hclust", "original", "AOE", "FPC", "alphabet"),
                                   hclust_method = "complete",
                                   addrect       = NULL,
                                   rect.col      = "navy",
                                   rect.lwd      = 2,
                                   p_thresh      = 0.05,
                                   insig         = c("blank", "p-value", "pch", "label_sig", "n"),
                                   pch           = 4,
                                   sig.level     = 0.05,
                                   palette_name  = "RdYlBu",
                                   color_scheme  = c("scientific", "viridis", "brew", "gradient"),
                                   tl.cex        = 0.8,
                                   tl.srt        = 45,
                                   tl.col        = "black",
                                   cl.cex        = 0.8,
                                   cl.pos        = c("b", "r", "n"),
                                   diag          = FALSE,
                                   save_plot     = FALSE,
                                   save_dir      = NULL,
                                   width         = 8,
                                   height        = 7,
                                   format        = "pdf",
                                   min_obs       = 3) {
  
  cat("Generating correlation heatmap (corrplot backend)...\n")
  vis_method    <- match.arg(vis_method)
  matrix_type   <- match.arg(matrix_type)
  order         <- match.arg(order)
  insig         <- match.arg(insig)
  color_scheme  <- match.arg(color_scheme)
  cl.pos        <- match.arg(cl.pos)
  
  # ------ 1. Data extraction and validation ------
  tryCatch({
    df <- if (inherits(object, "Stat")) object@clean.data else as.data.frame(object)
    if (!is.data.frame(df)) stop("Input must be a Stat object or data frame.")
    
    cat("Original data dimensions:", nrow(df), "x", ncol(df), "\n")
    
    numeric_cols <- names(df)[sapply(df, is.numeric)]
    if (length(numeric_cols) == 0) {
      stop("No numeric columns found in data.")
    }
    
    if (is.null(features)) {
      features <- head(numeric_cols, 40)
    }
    features <- intersect(features, numeric_cols)
    if (length(features) == 0) {
      stop("No valid numeric features found.")
    }
    
    cat("Selected features:", length(features), "\n")
    
    # ------ 1a. Handle missing values ------
    df_numeric <- df[, features, drop = FALSE]
    
    # Remove rows that are entirely NA
    all_na_rows <- apply(df_numeric, 1, function(x) all(is.na(x)))
    if (any(all_na_rows)) {
      cat("Warning: Removing", sum(all_na_rows), "rows with all missing values.\n")
      df_numeric <- df_numeric[!all_na_rows, ]
    }
    
    # Remove columns with excessive missing values (>50%)
    missing_pct <- colSums(is.na(df_numeric)) / nrow(df_numeric)
    high_missing <- missing_pct > 0.5
    if (any(high_missing)) {
      cat("Warning: Removing columns with >50% missing values:",
          paste(names(df_numeric)[high_missing], collapse = ", "), "\n")
      df_numeric <- df_numeric[, !high_missing, drop = FALSE]
      features <- colnames(df_numeric)
    }
    
    if (nrow(df_numeric) == 0 || ncol(df_numeric) == 0) {
      stop("No valid data remaining after removing missing values.")
    }
    
    cat("Data dimensions after NA handling:", nrow(df_numeric), "x", ncol(df_numeric), "\n")
    
    # ------ 1b. Remove constant columns ------
    zero_var <- apply(df_numeric, 2, function(x) {
      var(x, na.rm = TRUE) == 0 || all(is.na(x))
    })
    
    if (any(zero_var)) {
      cat("Warning: Removing columns with zero variance:",
          paste(names(df_numeric)[zero_var], collapse = ", "), "\n")
      df_numeric <- df_numeric[, !zero_var, drop = FALSE]
      features <- colnames(df_numeric)
    }
    
    if (ncol(df_numeric) < 2) {
      stop("At least 2 numeric features with variance are required for correlation analysis.")
    }
    
    # ------ 1c. Compute correlation matrix ------
    corr_mat <- tryCatch({
      cor(df_numeric, method = method, use = "pairwise.complete.obs")
    }, error = function(e) {
      # Fallback: use complete cases only
      cat("Warning: Using complete cases only due to:", conditionMessage(e), "\n")
      complete_idx <- complete.cases(df_numeric)
      if (sum(complete_idx) < min_obs) {
        stop("Insufficient complete cases for correlation computation (need at least ", min_obs, ")")
      }
      cor(df_numeric[complete_idx, ], method = method)
    })
    
    # Validate correlation matrix
    if (any(is.na(corr_mat))) {
      stop("Correlation matrix contains NaN values. Check data quality and variance.")
    }
    if (any(is.infinite(corr_mat))) {
      stop("Correlation matrix contains infinite values.")
    }
    
    cat("Correlation matrix computed successfully (",
        nrow(corr_mat), "x", ncol(corr_mat), ")\n")
    
  }, error = function(e) {
    stop("Error in data preparation: ", conditionMessage(e), call. = FALSE)
  })
  
  # ------ 2. Prepare p-value matrix if needed ------
  p_mat <- NULL
  if (insig != "n") {
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      warning("Package 'Hmisc' is required for significance masking. Install with: install.packages('Hmisc')")
      insig <- "n"
    } else {
      tryCatch({
        # Use the cleaned numeric data for p-value calculation
        corr_full <- Hmisc::rcorr(as.matrix(df_numeric), type = method)
        p_mat <- corr_full$P
        
        # Ensure p_mat aligns with corr_mat
        if (!all(rownames(p_mat) == rownames(corr_mat))) {
          p_mat <- p_mat[rownames(corr_mat), colnames(corr_mat)]
        }
        
        cat("P-value matrix computed successfully.\n")
        
        # Issue warning about order constraints
        if (order != "original") {
          cat("Note: Using order='hclust' with p-value masking. ",
              "P-values are aligned based on clustering order.\n")
        }
      }, error = function(e) {
        cat("Warning: Could not compute p-value matrix:", conditionMessage(e), "\n")
        cat("Proceeding without significance masking.\n")
        insig <<- "n"
      })
    }
  }
  
  # ------ 3. Colour scheme mapping ------
  col_vector <- tryCatch({
    switch(color_scheme,
           viridis    = {
             if (!requireNamespace("viridisLite", quietly = TRUE)) {
               warning("Package 'viridisLite' not found. Using default palette.")
               grDevices::colorRampPalette(c("#440154", "#31688e", "#35b779", "#fde724"))(200)
             } else {
               viridisLite::viridis(200, option = "C", direction = -1)
             }
           },
           brew       = {
             if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
               warning("Package 'RColorBrewer' not found. Using default palette.")
               grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200)
             } else {
               corrplot::COL2(palette_name, 200)
             }
           },
           gradient   = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
           scientific = grDevices::colorRampPalette(c("#3B528B", "#5A6BB5", "#8A8DC9", "#F7F7F7",
                                                      "#EAA582", "#D6604D", "#B2182B"))(200)
    )
  }, error = function(e) {
    cat("Warning: Color scheme error, using default palette:", conditionMessage(e), "\n")
    grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200)
  })
  
  # ------ 4. Save plot if requested ------
  if (save_plot) {
    tryCatch({
      if (is.null(save_dir)) {
        save_dir <- getwd()  # Use current working directory as fallback
        cat("No save_dir specified. Using current directory:", save_dir, "\n")
      }
      
      if (!dir.exists(save_dir)) {
        dir.create(save_dir, recursive = TRUE)
        cat("Created output directory:", save_dir, "\n")
      }
      
      plot_file <- file.path(save_dir, paste0("correlation_heatmap_", 
                                              format(Sys.time(), "%Y%m%d_%H%M%S"), 
                                              ".", format))
      
      # Close any open devices
      if (!is.null(grDevices::dev.list())) {
        try(grDevices::dev.off(), silent = TRUE)
      }
      
      # Open appropriate device
      if (format == "pdf") {
        grDevices::pdf(plot_file, width = width, height = height)
      } else if (format == "png") {
        grDevices::png(plot_file, width = width, height = height, units = "in", res = 300)
      } else if (format == "svg") {
        grDevices::svg(plot_file, width = width, height = height)
      } else {
        stop("Unsupported format. Use 'pdf', 'png', or 'svg'.")
      }
      
      on.exit({
        if (!is.null(grDevices::dev.list())) {
          try(grDevices::dev.off(), silent = TRUE)
        }
      }, add = TRUE)
      
    }, error = function(e) {
      stop("Error setting up plot output: ", conditionMessage(e), call. = FALSE)
    })
  }
  
  # ------ 5. Build and execute corrplot call ------
  tryCatch({
    # Base arguments with defaults
    args <- list(
      corr         = corr_mat,
      method       = vis_method,
      type         = matrix_type,
      order        = if (order == "original") "original" else order,
      col          = col_vector,
      col.lim      = c(-1, 1),
      tl.cex       = tl.cex,
      tl.srt       = tl.srt,
      tl.col       = tl.col,
      cl.cex       = cl.cex,
      cl.pos       = cl.pos,
      diag         = diag,
      addgrid.col  = "grey90",
      mar          = c(0, 0, 0, 0),
      na.label     = "?"     # Label NA cells
    )
    
    # Add hierarchical clustering method if needed
    if (order %in% c("hclust", "AOE", "FPC")) {
      args$hclust.method <- hclust_method
    }
    
    # Add p-value and significance arguments (only if p_mat is available)
    if (!is.null(p_mat) && insig != "n") {
      args$p.mat <- p_mat
      args$sig.level <- p_thresh
      args$insig <- insig
      if (insig == "pch") args$pch <- pch
      if (insig == "label_sig") args$sig.level <- sig.level
    }
    
    # Add rectangle arguments if requested (only valid with hclust order)
    if (!is.null(addrect) && addrect > 0 && order == "hclust") {
      args$addrect <- addrect
      args$rect.col <- rect.col
      args$rect.lwd <- rect.lwd
    }
    
    # Execute the plot
    cat("Drawing correlation heatmap...\n")
    do.call(corrplot::corrplot, args)
    
    if (save_plot) {
      cat("[OK] Heatmap saved to:", plot_file, "\n")
    }
    
  }, error = function(e) {
    stop("Error generating corrplot: ", conditionMessage(e), call. = FALSE)
  })
  
  cat("[OK] Correlation heatmap completed successfully.\n")
  invisible(corr_mat)
}
#' PCA Plot with Custom Group Colors
#'
#' This function performs Principal Component Analysis (PCA) on numeric columns
#' of a dataset and generates a scatter plot of the first two principal components.
#' Points can be colored and shaped by categorical variables. Custom colors
#' for the grouping variable can be specified via `group_colors`.
#'
#' @param object A `Stat` object (from the package) or a data frame.
#' @param color_by Name of the column in `info.data` (or the data frame) used
#'   to color the points. If NULL and `object` is a `Stat` object, the
#'   `group_col` slot is used.
#' @param shape_by Optional column name for shaping points. Default is NULL.
#' @param pcs Integer vector of length 2 indicating which principal components
#'   to plot. Default is c(1, 2).
#' @param ellipse Logical; if TRUE, draws a 95% confidence ellipse for each group.
#'   Default is TRUE.
#' @param label_points Logical; if TRUE, labels each point with its row name.
#'   Default is FALSE.
#' @param palette_name Name of the palette to use when `group_colors` is not
#'   provided. Supports RColorBrewer (e.g., "Set1") or wesanderson (e.g.,
#'   "Darjeeling1"). Default is "Darjeeling1".
#' @param base_size Base font size for the plot theme. Default is 13.
#' @param save_plot Logical; if TRUE, saves the plot to disk. Default is FALSE.
#' @param save_dir Directory to save the plot. If NULL and save_plot is TRUE,
#'   a default directory is used.
#' @param width Width of the saved plot (inches). Default is 7.
#' @param height Height of the saved plot (inches). Default is 6.
#' @param format Output format ("pdf" or "png"). Default is "pdf".
#' @param group_colors An optional named character vector specifying custom
#'   colors for the groups. Names must match the levels of the `color_by`
#'   variable. If provided, it overrides `palette_name`. Default is NULL.
#'
#' @return A `ggplot` object. The plot is also saved if `save_plot = TRUE`.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- PlotPCA(stat_obj, color_by = "group",
#'              group_colors = c("0" = "#404040", "1" = "#ca0020"))
#' }
PlotPCA <- function(object,
                    color_by      = NULL,
                    shape_by      = NULL,
                    pcs           = c(1, 2),
                    ellipse       = TRUE,
                    label_points  = FALSE,
                    palette_name  = "Darjeeling1",
                    base_size     = 13,
                    save_plot     = FALSE,
                    save_dir      = NULL,
                    width         = 7,
                    height        = 6,
                    format        = "pdf",
                    group_colors  = NULL) {   # new parameter
  
  cat("Generating PCA plot...\n")
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  
  # ---- Extract data ----
  if (inherits(object, "Stat")) {
    mat      <- object@clean.data[, sapply(object@clean.data, is.numeric)]
    info     <- object@info.data
    if (is.null(color_by)) color_by <- object@group_col
  } else {
    mat  <- as.data.frame(object)[, sapply(as.data.frame(object), is.numeric)]
    info <- data.frame(row.names = rownames(mat))
  }
  
  # ---- Remove zero-variance columns ----
  keep_cols <- apply(mat, 2, stats::var, na.rm = TRUE) > 0
  if (any(!keep_cols)) {
    removed <- colnames(mat)[!keep_cols]
    message("Removed ", length(removed), " zero-variance column(s): ",
            paste(removed, collapse = ", "))
    mat <- mat[, keep_cols, drop = FALSE]
  }
  if (ncol(mat) < 2) stop("After removing constant columns, fewer than 2 variables remain.")
  
  # ---- PCA ----
  pca_res  <- stats::prcomp(mat, scale. = TRUE, center = TRUE)
  var_exp  <- round(summary(pca_res)$importance[2, ] * 100, 1)
  scores   <- as.data.frame(pca_res$x[, pcs])
  colnames(scores) <- paste0("PC", pcs)
  
  # ---- Prepare color column ----
  if (nrow(info) == nrow(scores) && color_by %in% colnames(info)) {
    scores[[color_by]] <- as.factor(info[rownames(scores), color_by])
  } else {
    scores[[color_by]] <- factor("all")
  }
  
  # ---- Shape column ----
  if (!is.null(shape_by) && shape_by %in% colnames(info)) {
    scores[[shape_by]] <- as.factor(info[rownames(scores), shape_by])
  }
  
  # ---- Determine colors ----
  group_levels <- levels(scores[[color_by]])
  n_grp <- length(group_levels)
  
  if (!is.null(group_colors)) {
    # Check if group_colors is named and matches all groups
    if (!is.null(names(group_colors))) {
      matched <- intersect(names(group_colors), group_levels)
      if (length(matched) == n_grp) {
        cols <- group_colors[group_levels]   # reorder to match group order
      } else {
        warning("Names in 'group_colors' do not match all group levels. Falling back to palette.")
        cols <- .get_palette(palette_name, n_grp)
      }
    } else {
      # Unnamed: assume order matches group levels
      if (length(group_colors) == n_grp) {
        cols <- group_colors
        names(cols) <- group_levels
        message("Using 'group_colors' in the order of group levels: ",
                paste(group_levels, collapse = ", "))
      } else {
        warning("Length of 'group_colors' does not match number of groups. Falling back to palette.")
        cols <- .get_palette(palette_name, n_grp)
      }
    }
  } else {
    cols <- .get_palette(palette_name, n_grp)
  }
  names(cols) <- group_levels
  
  # ---- Axis labels ----
  xlab <- paste0("PC", pcs[1], " (", var_exp[pcs[1]], "%)")
  ylab <- paste0("PC", pcs[2], " (", var_exp[pcs[2]], "%)")
  
  # ---- Aesthetics ----
  if (!is.null(shape_by) && shape_by %in% colnames(scores)) {
    aes_map <- ggplot2::aes(.data[[paste0("PC", pcs[1])]],
                            .data[[paste0("PC", pcs[2])]],
                            colour = .data[[color_by]],
                            shape  = .data[[shape_by]])
  } else {
    aes_map <- ggplot2::aes(.data[[paste0("PC", pcs[1])]],
                            .data[[paste0("PC", pcs[2])]],
                            colour = .data[[color_by]])
  }
  
  # ---- Build plot ----
  p <- ggplot2::ggplot(scores, aes_map) +
    ggplot2::geom_point(size = 2.5, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(title = "PCA", x = xlab, y = ylab, colour = color_by) +
    .pub_theme(base_size)
  
  if (ellipse) p <- p + ggplot2::stat_ellipse(level = 0.95, linewidth = 0.7, linetype = 2)
  if (label_points) p <- p + ggrepel::geom_text_repel(ggplot2::aes(label = rownames(scores)),
                                                      size = 2.5, max.overlaps = 20)
  if (save_plot) .save_plot(p, save_dir, "PCA", width, height, format)
  return(p)
}

#' Feature selection plot: AUC vs. -log10(p-value) (enhanced)
#'
#' Calculates per-feature ROC AUC and plots it against -log10(p-value).
#' Points that satisfy both \code{auc_thresh} and \code{p_thresh} are highlighted.
#'
#' @param deg_df       Data frame with columns \code{feature}, \code{logFC}, \code{p.adjust} / \code{p_value}.
#' @param mat_test     Expression matrix containing the grouping column.
#' @param group_col    Name of the grouping column in \code{mat_test}.
#' @param auc_thresh   AUC threshold for selection. Default 0.55.
#' @param p_thresh     P-value threshold. Default 0.05.
#' @param selected_fill Fill colour for selected features.
#' @param removed_fill  Fill colour for removed features.
#' @param arrow        Logical. Draw arrows from labels to points.
#' @param base_size    Base font size.
#' @param text_repel_size    Annotation text font size.
#' @param save_plot    Logical. Save the plot? Default \code{FALSE}.
#' @param save_dir     Output directory. When \code{NULL} and \code{save_plot=TRUE}
#'   the default Stat figure directory is used.
#' @param width,height Plot dimensions (inches).
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' last_sig <- ExtractLastTestSig(stat_obj)
#' last_sig$feature <- last_sig$id
#' PlotAUCPval(last_sig, stat_obj_test@clean.data, group_col = "SWAB", 
#' save_plot = FALSE,p_thresh = 0.1, auc_thresh = 0.5)
#' }
PlotAUCPval <- function(deg_df, mat_test,
                        group_col      = "group",
                        auc_thresh     = 0.55,
                        p_thresh       = 0.05,
                        selected_fill  = "#e84118",
                        removed_fill   = "#7f8fa6",
                        arrow          = TRUE,
                        base_size      = 13,
                        text_repel_size=3,
                        save_plot      = FALSE,
                        save_dir       = NULL,
                        width          = 7,
                        height         = 6,
                        format         = "pdf") {
  cat("Generating AUC-P value selection plot...\n")
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  df <- as.data.frame(deg_df)
  if (!"feature" %in% colnames(df)) {
    possible_cols <- c("id", "gene", "symbol", "rownames")
    for (col in possible_cols) {
      if (col %in% colnames(df)) {
        df$feature <- df[[col]]
        break
      }
    }
    if (!"feature" %in% colnames(df)) {
      df$feature <- rownames(df)   
    }
  }
  p_col <- if ("p.adjust" %in% colnames(df)) "p.adjust" else "p_value"
  df$neg_log10p <- -log10(df[[p_col]] + 1e-300)
  
  mat <- as.data.frame(mat_test)
  if (!group_col %in% colnames(mat)) stop("Group column not found in mat_test.")
  mat[[group_col]] <- as.factor(mat[[group_col]])
  lev <- levels(mat[[group_col]])
  
  auc_vec <- sapply(df$feature, function(feat) {
    if (!feat %in% colnames(mat)) return(NA_real_)
    tem_dat <- mat[, c(feat, group_col)]
    colnames(tem_dat)[1] <- "feature"
    colnames(tem_dat)[2] <- "group"
    tem_dat$feature <- as.numeric(tem_dat$feature)
    roc_obj <- tryCatch(
      pROC::roc(group ~ feature, data = tem_dat, levels = lev, direction = "auto", quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(roc_obj)) return(NA_real_)
    round(as.numeric(pROC::auc(roc_obj)), 3)
  })
  
  df$auc <- auc_vec
  if (any(is.na(df$auc))) cat("Features with missing AUC:", sum(is.na(df$auc)), "\n")
  
  df$group <- ifelse(!is.na(df$auc) & df$auc >= auc_thresh & df[[p_col]] < p_thresh,
                     "selected", "removed")
  selected_df <- df[df$group == "selected", , drop = FALSE]
  n_sel <- nrow(selected_df)
  cat("Selected features:", n_sel, "\n")
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = auc, y = neg_log10p)) +
    ggplot2::geom_vline(xintercept = auc_thresh, linetype = 2,
                        colour = "black", linewidth = 0.8, alpha = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = 2,
                        colour = "black", linewidth = 0.8, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(fill = group, size = neg_log10p),
                        colour = "black", stroke = 0.3, shape = 21, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(selected = selected_fill, removed = removed_fill),
                               name = NULL) +
    ggplot2::scale_size(range = c(1, 10), guide = "none") +
    ggrepel::geom_text_repel(
      data = selected_df,
      ggplot2::aes(label = feature),
      size = text_repel_size, colour = "black", fontface = "bold.italic",
      min.segment.length = 0,
      arrow = if (arrow) ggplot2::arrow(length = ggplot2::unit(0.01, "npc")) else NULL,
      show.legend = FALSE, max.overlaps = 200,
      segment.curvature = -0.05, segment.square = FALSE, segment.inflect = TRUE,
      force = 0.1, nudge_x = 0.02, direction = "y", hjust = 0
    ) +
    ggplot2::annotate("text",
                      x = max(df$auc, na.rm = TRUE) * 0.95,
                      y = max(df$neg_log10p, na.rm = TRUE) * 0.95,
                      label = paste0("Selected: ", n_sel),
                      colour = selected_fill, fontface = "bold", size = 4, hjust = 1) +
    .pub_theme(base_size) +
    ggplot2::labs(title = "Feature Selection by AUC & Significance",
                  x = "ROC AUC", y = expression(-log[10]~"(P-value)")) +
    ggplot2::theme(legend.position = "right")
  
  if (save_plot) .save_plot(p, save_dir, "PlotAUCPval", width, height, format)
  return(p)
}

#' Boxplot of top differential features with significance stars (enhanced)
#'
#' @param deg_results DEG data frame with columns 'id', 'change', 'logFC'.
#' @param expr_data   Expression matrix (data frame with features + group column).
#' @param group_col   Column name for the group variable.
#' @param top_n       Number of top features to show.
#' @param palette_fill Colour vector of length 2 for groups.
#' @param base_size   Base font size.
#' @param save_plot   Logical. Save the plot? Default \code{FALSE}.
#' @param save_dir    Output directory. When \code{NULL} and \code{save_plot=TRUE}
#'   the default Stat figure directory is used.
#' @param width,height Plot dimensions.
#' @param format      Output format.
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' stat_obj <- stat_var_feature(stat_obj_test)
#' last_sig <- ExtractLastTestSig(stat_obj)
#' PlotDegBoxplot(last_sig, stat_obj@clean.data, group_col = "SWAB", save_plot = FALSE)
#' }
PlotDegBoxplot <- function(deg_results,
                           expr_data,
                           group_col   = "group",
                           top_n       = 5,
                           palette_fill = c("#999999", "#ef8a62"),
                           base_size   = 12,
                           save_plot   = FALSE,
                           save_dir    = NULL,
                           width       = 7,
                           height      = 5,
                           format      = "pdf") {
  cat("Generating DEG boxplot...\n")
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("Stat")
  
  sig_df <- deg_results[deg_results$change != "Stable", ]
  sig_df <- sig_df[order(-abs(sig_df$logFC)), ]
  n_up   <- sum(sig_df$change == "Up")
  n_down <- sum(sig_df$change == "Down")
  top_up   <- head(sig_df[sig_df$change == "Up", ], min(ceiling(top_n/2), n_up))
  top_down <- head(sig_df[sig_df$change == "Down", ], min(floor(top_n/2), n_down))
  selected <- rbind(top_up, top_down)
  features <- selected$id
  
  box_df <- expr_data[, c(group_col, features), drop = FALSE]
  box_df[[group_col]] <- as.factor(box_df[[group_col]])
  long_df <- tidyr::pivot_longer(box_df, cols = -dplyr::all_of(group_col),
                                 names_to = "feature", values_to = "value")
  
  p_vals <- long_df %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(p_value = wilcox.test(value ~ .data[[group_col]])$p.value,
                     .groups = "drop") %>%
    dplyr::mutate(significance = dplyr::case_when(
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE           ~ ""
    ))
  
  long_df$log_val <- log10(long_df$value + 1)
  y_pos <- long_df %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(y_position = max(log_val, na.rm = TRUE) + 0.2, .groups = "drop")
  p_vals <- dplyr::left_join(p_vals, y_pos, by = "feature")
  long_df <- dplyr::left_join(long_df, p_vals, by = "feature")
  
  long_df$feature <- factor(long_df$feature,
                            levels = p_vals$feature[order(p_vals$y_position, decreasing = TRUE)])
  
  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = feature, y = log_val, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.size = 0.5, outlier.colour = "grey50") +
    ggplot2::geom_text(
      data = dplyr::distinct(long_df, feature, .keep_all = TRUE),
      ggplot2::aes(x = feature, y = y_position, label = significance),
      size = 5, colour = "black", inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = palette_fill) +
    ggplot2::labs(title = paste0("Top ", top_n, " DEGs"),
                  x = NULL, y = expression(log[10]~(Expression + 1)), fill = group_col) +
    .pub_theme(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,vjust = 1))
  print(p)
  if (save_plot) .save_plot(p, save_dir, "deg_boxplot", width, height, format)
  return(p)
}

#' Feature Heatmap with Annotations and Marker Highlighting
#'
#' Creates a heatmap of selected features (e.g., genes, proteins) with optional
#' sample annotations (top/bottom) and feature annotations (left/right). Supports
#' highlighting specific marker features using annotation marks on the left side.
#'
#' @param object A data frame (samples x features), numeric matrix, or a `Stat` object
#'   (which must contain a `clean.data` component with numeric columns).  
#'   Rows = samples, columns = features.
#' @param features Character vector of feature names to display.  
#'   If `NULL`, the top 50 features with highest variance are selected.
#' @param group_col (Legacy) Character string; column name in the `Stat` object's
#'   info data used for sample grouping. This will be added to the top annotation.
#' @param clinical_data Data frame with sample annotations (rows = samples).  
#'   Column names are variable names, row names must match sample names in `object`.
#'   Added as top annotation.
#' @param left_annotation,right_annotation,top_annotation,bottom_annotation  
#'   Data frames for row/column annotations.  
#'   - Row annotations (left/right): row names must match feature names (colnames of `object`).  
#'   - Column annotations (top/bottom): row names must match sample names (rownames of `object`).
#' @param left_columns,right_columns,top_columns,bottom_columns  
#'   Character vectors specifying which columns to use from the corresponding
#'   annotation data frames. If `NULL`, all columns are used.
#' @param marker_features Character vector of feature names to highlight on the left side
#'   using `anno_mark`. If `NULL`, no marker annotation is added.
#' @param marker_cex Font size for marker labels. If `NULL`, automatically scaled
#'   based on the number of markers (range 0.5-1.2).
#' @param scale_features Logical. Should each feature (column) be Z-score scaled
#'   across samples? Default `TRUE`.
#' @param exclude_constant Logical. If `TRUE`, features with zero variance
#'   (constant values) are removed before analysis. Default `TRUE`.
#' @param cluster_rows Logical. Cluster rows (features). Default `TRUE`.
#' @param cluster_cols Logical. Cluster columns (samples). Default `TRUE`.
#' @param show_rownames Logical. Show row (feature) names. Default `FALSE`
#'   (often hidden when using `marker_features`).
#' @param show_colnames Logical. Show column (sample) names. Default `FALSE`.
#' @param palette_name RColorBrewer diverging palette name, used when `color_palette = NULL`.
#'   Default `"RdYlBu"`.
#' @param color_palette Character vector of 4 colours:  
#'   `c(NA_colour, low_colour, mid_colour, high_colour)`.  
#'   Default: `c("grey90", "#7B1FA2", "#000000", "#FFEB3B")` (purple-black-yellow).
#' @param color_range Numeric vector of length 3 (low, mid, high) for colour mapping.
#'   If `NULL`, the 5\% and 95\% quantiles of the data are used, with mid = 0.
#' @param na_col Colour for missing values (NA). Default `"grey90"`.
#' @param ann_colors Optional named list of colours for annotations.  
#'   For discrete variables: a named vector (levels -> colours).  
#'   For continuous variables: a `colorRamp2` function.  
#'   Example: `list(Sex = c(M = "blue", F = "red"), Age = colorRamp2(c(20,50,80), c("blue","white","red")))`.
#' @param ann_continuous_palette Character; RColorBrewer palette for continuous annotation
#'   variables when `ann_colors` is not provided. Default `"RdYlBu"`.
#' @param base_fontsize Base font size for legends and annotation labels.
#' @param save_plot Logical. Save the plot to a file? Default `FALSE`.
#' @param save_dir Output directory. If `NULL` and `save_plot = TRUE`, files are written
#'   to `"./plots"`.
#' @param width,height Dimensions of the saved plot in inches.
#' @param format File format: `"pdf"` or `"png"`.
#' @param ... Additional arguments passed to `ComplexHeatmap::Heatmap`.
#'
#' @return Invisibly, a `ComplexHeatmap::Heatmap` object. The heatmap is drawn
#'   (or saved) as a side effect.
#'
#' @importFrom ComplexHeatmap Heatmap rowAnnotation anno_mark max_text_width ht_opt
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid gpar unit
#' @export
#' 
#' @examples
#' \dontrun{
#' stat_obj_test@info.data=stat_obj_test@clean.data
#' PlotFeatureHeatmap(stat_obj_test, clinical_data = stat_obj@info.data, save_plot = FALSE)
#' # Example with a simple data frame
#' set.seed(123)
#' expr <- data.frame(
#'   Gene1 = rnorm(20), Gene2 = rnorm(20), Gene3 = rnorm(20),
#'   Gene4 = rnorm(20), Gene5 = rnorm(20)
#' )
#' rownames(expr) <- paste0("Sample", 1:20)
#'
#' # Clinical annotation
#' clin <- data.frame(
#'   Age = sample(20:80, 20, replace = TRUE),
#'   Sex = sample(c("M","F"), 20, replace = TRUE),
#'   row.names = rownames(expr)
#' )
#'
#' # Basic heatmap
#' PlotFeatureHeatmap(expr, clinical_data = clin)
#'
#' # Highlight specific markers
#' markers <- c("Gene1", "Gene5")
#' PlotFeatureHeatmap(expr, clinical_data = clin, marker_features = markers)
#'
#' # Custom annotation colours (like your example)
#' my_colors <- list(
#'   Sex = c(M = "#fc8d62", F = "#8da0cb"),
#'   Age = circlize::colorRamp2(c(20, 50, 80), c("#ef8a62", "#ffffff", "#999999"))
#' )
#' PlotFeatureHeatmap(expr, clinical_data = clin, ann_colors = my_colors)
#'
#' # Using a Stat object (assuming it exists)
#' # PlotFeatureHeatmap(stat_obj, clinical_data = stat_obj@clean.data[, c("GENDER","SWAB")])
#' }
PlotFeatureHeatmap <- function(object,
                               features          = NULL,
                               group_col         = NULL,
                               clinical_data     = NULL,
                               left_annotation   = NULL,
                               right_annotation  = NULL,
                               top_annotation    = NULL,
                               bottom_annotation = NULL,
                               left_columns      = NULL,
                               right_columns     = NULL,
                               top_columns       = NULL,
                               bottom_columns    = NULL,
                               marker_features   = NULL,
                               marker_cex        = NULL,
                               scale_features    = TRUE,
                               exclude_constant  = TRUE,
                               cluster_rows      = TRUE,
                               cluster_cols      = TRUE,
                               show_rownames     = FALSE,
                               show_colnames     = FALSE,
                               palette_name      = "RdYlBu",
                               color_palette     =  c("lightgray", "#7B1FA2", "black", "gold"),
                               color_range       = NULL,
                               na_col            = "grey90",
                               ann_colors        = NULL,
                               ann_continuous_palette = "RdYlBu",
                               base_fontsize     = 10,
                               save_plot         = FALSE,
                               save_dir          = NULL,
                               width             = 9,
                               height            = 8,
                               format            = "pdf",
                               ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    stop("Package 'ComplexHeatmap' is required.")
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package 'circlize' is required.")
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    stop("Package 'RColorBrewer' is required.")
  
  tryCatch({
    ComplexHeatmap::ht_opt$message <- FALSE
    ComplexHeatmap::ht_opt$verbose <- FALSE
  }, error = function(e) NULL)
  
  cat("Generating feature heatmap...\n")
  
  if (save_plot && is.null(save_dir)) {
    save_dir <- "./plots"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  }
  
  # ------------------------- 1. Extract numeric matrix -------------------------
  if (inherits(object, "Stat")) {
    if (isS4(object)) {
      raw <- object@clean.data
      info <- object@info.data
      if (is.null(group_col) && !is.null(object@group_col)) group_col <- object@group_col
    } else {
      raw <- object$clean.data
      info <- object$info.data
      if (is.null(group_col) && !is.null(object$group_col)) group_col <- object$group_col
    }
    if (is.null(raw)) stop("No 'clean.data' found in Stat object.")

    num_idx <- sapply(raw, is.numeric)
    if (sum(num_idx) == 0) stop("No numeric columns in clean.data.")
    mat <- as.matrix(raw[, num_idx, drop = FALSE])
    if (is.null(info)) info <- data.frame()
  } else {
    raw <- as.data.frame(object)
    num_idx <- sapply(raw, is.numeric)
    if (sum(num_idx) == 0) stop("No numeric columns found in object.")
    mat <- as.matrix(raw[, num_idx, drop = FALSE])
    info <- data.frame()
  }
  
  if (is.null(rownames(mat))) rownames(mat) <- as.character(seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("Var", seq_len(ncol(mat)))
  
  # ------------------------- 2. Filter constant features -------------------------
  if (exclude_constant) {
    feat_var <- apply(mat, 2, var, na.rm = TRUE)
    const_feats <- names(feat_var[feat_var == 0 | is.na(feat_var)])
    if (length(const_feats) > 0) {
      warning(paste("Excluded", length(const_feats), "constant features:",
                    paste(head(const_feats, 5), collapse = ", ")))
      mat <- mat[, !colnames(mat) %in% const_feats, drop = FALSE]
    }
  }
  
  # ------------------------- 3. Select top features -------------------------
  if (is.null(features)) {
    feat_var <- apply(mat, 2, var, na.rm = TRUE)
    features <- names(sort(feat_var, decreasing = TRUE))[1:min(50, length(feat_var))]
  }
  features <- intersect(features, colnames(mat))
  if (length(features) == 0) stop("No valid features found.")
  mat <- mat[, features, drop = FALSE]
  
  # ------------------------- 4. Scale features (columns) -------------------------
  if (scale_features) {
    col_mean <- colMeans(mat, na.rm = TRUE)
    col_sd   <- apply(mat, 2, sd, na.rm = TRUE)
    col_sd[col_sd == 0] <- 1
    mat <- sweep(mat, 2, col_mean, "-")
    mat <- sweep(mat, 2, col_sd, "/")
    mat[is.nan(mat)] <- 0
  }
  
  # ------------------------- 5. Helper: generate colors for annotation -------------------------
  generate_anno_colors <- function(df, columns, provided_colors = ann_colors) {
    if (is.null(df)) return(NULL)
    if (!is.data.frame(df)) return(NULL)
    if (!is.null(columns)) {
      columns <- intersect(columns, colnames(df))
      if (length(columns) == 0) return(NULL)
      df <- df[, columns, drop = FALSE]
    }
    col_map <- list()
    for (col in colnames(df)) {
      if (!is.null(provided_colors) && col %in% names(provided_colors)) {
        col_map[[col]] <- provided_colors[[col]]
        next
      }
      if (is.factor(df[[col]]) || is.character(df[[col]])) {
        vals <- unique(df[[col]])
        n <- length(vals)
        palettes <- c("Set3", "Paired", "Set2")
        colors_all <- unlist(lapply(palettes, function(p) {
          RColorBrewer::brewer.pal(min(12, RColorBrewer::brewer.pal.info[p, "maxcolors"]), p)
        }))
        if (n > length(colors_all)) colors_all <- c(colors_all, rainbow(n))
        cols <- colors_all[1:n]
        names(cols) <- as.character(vals)
        col_map[[col]] <- cols
      } else if (is.numeric(df[[col]])) {
        rng <- range(df[[col]], na.rm = TRUE)
        if (diff(rng) == 0) rng <- c(rng[1] - 1, rng[1] + 1)
        col_map[[col]] <- circlize::colorRamp2(
          breaks = c(rng[1], mean(rng), rng[2]),
          colors = RColorBrewer::brewer.pal(3, ann_continuous_palette)
        )
      }
    }
    return(col_map)
  }
  
  make_annotation <- function(df, columns, which_dim, orientation, col_map) {
    if (is.null(df)) return(NULL)
    if (!is.data.frame(df)) {
      warning("Annotation must be a data frame.")
      return(NULL)
    }

    if (which_dim == "feature") {
      ids <- colnames(mat)
      df_ids <- rownames(df)
    } else {
      ids <- rownames(mat)
      df_ids <- rownames(df)
    }
    common <- intersect(df_ids, ids)
    if (length(common) == 0) {
      warning("No overlapping IDs between annotation and ", which_dim, "s.")
      return(NULL)
    }
    df <- df[common, , drop = FALSE]
    df <- df[ids, , drop = FALSE]  
    if (!is.null(columns)) {
      columns <- intersect(columns, colnames(df))
      if (length(columns) == 0) return(NULL)
      df <- df[, columns, drop = FALSE]
    }

    for (col in colnames(df)) {
      if (is.character(df[[col]])) df[[col]] <- as.factor(df[[col]])
    }
    ComplexHeatmap::HeatmapAnnotation(
      df = df,
      col = col_map,
      which = orientation,
      annotation_name_gp = grid::gpar(fontsize = base_fontsize),
      annotation_legend_param = list(
        title_gp = grid::gpar(fontsize = base_fontsize),
        labels_gp = grid::gpar(fontsize = base_fontsize - 1)
      ),
      show_annotation_name = TRUE
    )
  }
  
  if (!is.null(group_col) && nrow(info) > 0 && group_col %in% colnames(info)) {
    legacy_ann <- data.frame(Group = info[[group_col]], row.names = rownames(info))
    clinical_data <- if (is.null(clinical_data)) legacy_ann else cbind(clinical_data, legacy_ann)
  }
  if (!is.null(clinical_data)) {
    top_annotation <- if (is.null(top_annotation)) clinical_data else cbind(top_annotation, clinical_data)
  }
  
  left_colors   <- generate_anno_colors(left_annotation, left_columns)
  right_colors  <- generate_anno_colors(right_annotation, right_columns)
  top_colors    <- generate_anno_colors(top_annotation, top_columns)
  bottom_colors <- generate_anno_colors(bottom_annotation, bottom_columns)
  
  left_ha   <- make_annotation(left_annotation,  left_columns,  "feature", "row", left_colors)
  right_ha  <- make_annotation(right_annotation, right_columns, "feature", "row", right_colors)
  top_ha    <- make_annotation(top_annotation,   top_columns,   "sample",  "column", top_colors)
  bottom_ha <- make_annotation(bottom_annotation, bottom_columns, "sample", "column", bottom_colors)
  
  # Marker annotation (left)
  if (!is.null(marker_features)) {
    hm_mat <- t(mat)
    marker_pos <- which(rownames(hm_mat) %in% marker_features)
    if (length(marker_pos) > 0) {
      marker_labels <- rownames(hm_mat)[marker_pos]
      if (is.null(marker_cex)) {
        n_markers <- length(marker_pos)
        marker_cex <- max(0.5, min(1.2, 30 * 0.65 / max(n_markers, 25)))
      }
      marker_ha <- ComplexHeatmap::rowAnnotation(
        link = ComplexHeatmap::anno_mark(
          at = marker_pos,
          labels = marker_labels,
          labels_gp = grid::gpar(cex = marker_cex, fontface = "italic"),
          link_gp = grid::gpar(col = "black"),
          side = "left",
          extend = c(0, 0.3),
          link_width = grid::unit(0.5, "in")
        ),
        width = grid::unit(0.5, "in") + ComplexHeatmap::max_text_width(marker_labels,
                                                                       gp = grid::gpar(cex = marker_cex))
      )
      left_ha <- if (is.null(left_ha)) marker_ha else c(left_ha, marker_ha)
    } else {
      warning("None of the marker_features found in the selected feature set.")
    }
  }
  
  # ------------------------- 7. Heatmap body colors -------------------------
  hm_mat <- t(mat)
  mat_vals <- as.vector(hm_mat)
  if (is.null(color_range)) {
    q05 <- quantile(mat_vals, 0.05, na.rm = TRUE)
    q95 <- quantile(mat_vals, 0.95, na.rm = TRUE)
    if (q05 == q95) {
      q05 <- q05 - 1e-6
      q95 <- q95 + 1e-6
    }
    color_range <- c(q05, 0, q95)
  }
  if (length(color_palette) != 4) {
    warning("color_palette should be length 4 (NA, low, mid, high). Using default.")
    color_palette <- c("grey90", "#7B1FA2", "#000000", "#FFEB3B")
  }
  col_fun <- circlize::colorRamp2(breaks = color_range, colors = color_palette[2:4])
  
  # ------------------------- 8. Build main heatmap -------------------------
  ht <- ComplexHeatmap::Heatmap(
    matrix = hm_mat,
    name = if (scale_features) "Z-score" else "Expression",
    col = col_fun,
    na_col = color_palette[1],
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_rownames,
    show_column_names = show_colnames,
    row_names_gp = grid::gpar(fontsize = base_fontsize - 1),
    column_names_gp = grid::gpar(fontsize = base_fontsize - 1),
    heatmap_legend_param = list(
      title = if (scale_features) "Z-score" else "Expression",
      title_gp = grid::gpar(fontsize = base_fontsize),
      labels_gp = grid::gpar(fontsize = base_fontsize - 1),
      direction = "horizontal",
      legend_width = grid::unit(1.5, "in"),
      title_position = "topcenter"
    ),
    left_annotation = left_ha,
    right_annotation = right_ha,
    top_annotation = top_ha,
    bottom_annotation = bottom_ha,
    use_raster = nrow(hm_mat) * ncol(hm_mat) > 200000,
    ...
  )
  
  # ------------------------- 9. Draw and save -------------------------
  if (!save_plot) {
    ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
  } else {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    path <- file.path(save_dir, paste0("feature_heatmap.", format))
    if (format == "pdf") {
      grDevices::pdf(path, width = width, height = height)
      ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
      grDevices::dev.off()
    } else if (format == "png") {
      grDevices::png(path, width = width, height = height, units = "in", res = 150)
      ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
      grDevices::dev.off()
    } else {
      stop("Unsupported format. Use 'pdf' or 'png'.")
    }
    cat("Heatmap saved to:", path, "\n")
  }
  
  invisible(ht)
}


## =============================================================================
##   2  TRAIN_MODEL  -- ROC - confusion matrix - feature importance - calibration
## =============================================================================

#' Plot ROC Curves for Multiple Models
#'
#' Generates a comparison plot of ROC curves across all trained models.
#' Handles both original and converted factor levels.
#'
#' @param object        A Train_Model object with trained models.
#' @param test_data     Test data (optional). If NULL, uses object@split.data$test.
#' @param palette_name  Color palette name (e.g., "Darjeeling1").
#' @param show_ci       Logical. Show confidence intervals (currently unused).
#' @param base_size     ggplot2 base font size.
#' @param save_plot     Logical. Save the plot to file?
#' @param save_dir      Directory to save plot. If NULL, auto-detects.
#' @param width         Plot width in inches.
#' @param height        Plot height in inches.
#' @param format        File format ("pdf", "png", "jpg").
#' @return A ggplot2 object showing ROC curves.
#' @export
#' @examples
#' \dontrun{
#' model_obj <- ModelTrainAnalysis(
#' object       = model_obj_test,
#' methods      = c("glm", "rf", "gbm"),
#' control      = list(method = "repeatedcv", number = 5, repeats = 1),
#' save_plots   = TRUE,
#' save_dir     = ".",
#' seed         = 123
#' )
#' PlotMultiROC(model_obj, save_plot = FALSE)
#' }
PlotMultiROC <- function(object,
                         test_data    = NULL,
                         palette_name = "Darjeeling1",
                         show_ci      = FALSE,
                         base_size    = 13,
                         save_plot    = FALSE,
                         save_dir     = NULL,
                         width        = 7,
                         height       = 6,
                         format       = "pdf") {
  
  if (is.null(save_dir)) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (!inherits(object, "Train_Model")) {
    stop("'object' must be a Train_Model.")
  }
  
  if (length(object@train.models) == 0) {
    stop("No trained models found in object@train.models.")
  }
  
  td <- if (!is.null(test_data)) {
    test_data
  } else {
    object@split.data$test
  }
  
  if (is.null(td)) {
    stop("Provide 'test_data' or populate split.data$test.")
  }
  
  group_col <- as.character(object@group_col)
  truth <- factor(td[[group_col]])
  
  truth_levels <- levels(truth)
  
  if (length(truth_levels) != 2) {
    stop(sprintf("Expected 2 factor levels, got %d. This function is for binary classification.",
                 length(truth_levels)))
  }
  
  pos_level <- truth_levels[2]
  
  cat(sprintf("Truth factor levels: %s, %s\n", truth_levels[1], truth_levels[2]))
  cat(sprintf("Using '%s' as positive class\n\n", pos_level))
  
  roc_list <- lapply(names(object@train.models), function(nm) {
    cat(sprintf("Processing model: %s\n", nm))
    
    model <- object@train.models[[nm]]
    
    probs <- tryCatch({
      prob_matrix <- stats::predict(model, newdata = td, type = "prob")
      if (pos_level %in% colnames(prob_matrix)) {
        prob_matrix[, pos_level]
      } else {
        cat(sprintf("  Warning: Column '%s' not found. Using 2nd column instead.\n", pos_level))
        prob_matrix[, 2]
      }
    }, error = function(e) {
      cat(sprintf("  Note: Using raw predictions instead of probabilities\n"))
      as.numeric(stats::predict(model, newdata = td, type = "raw") == pos_level)
    })
    if (length(probs) != nrow(td)) {
      stop(sprintf("Prediction length mismatch for model %s", nm))
    }
    
    if (any(is.na(probs))) {
      warning(sprintf("Model %s produced NA predictions", nm))
      probs[is.na(probs)] <- 0.5
    }
    roc_obj <- pROC::roc(truth, probs, 
                         levels = levels(truth), 
                         quiet = TRUE)
    
    auc_val <- round(as.numeric(pROC::auc(roc_obj)), 3)
    coords <- pROC::coords(roc_obj, "all", 
                           ret = c("specificity", "sensitivity"))
    
    cat(sprintf("  [OK] AUC = %.3f\n", auc_val))
    
    data.frame(
      Model       = paste0(nm, " (AUC=", auc_val, ")"),
      Specificity = coords$specificity,
      Sensitivity = coords$sensitivity,
      stringsAsFactors = FALSE
    )
  })
  
  roc_df <- dplyr::bind_rows(roc_list)
  
  if (nrow(roc_df) == 0) {
    stop("No ROC data generated. Check your models and test data.")
  }
  cols <- .get_palette(palette_name, length(object@train.models))
  p <- ggplot2::ggplot(roc_df,
                       ggplot2::aes(x = 1 - Specificity, y = Sensitivity, colour = Model)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50", linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "ROC Curves -- Model Comparison",
      x = "1 - Specificity (FPR)",
      y = "Sensitivity (TPR)",
      colour = NULL
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = c(0.72, 0.22))
  if (save_plot) {
    .safe_dir(save_dir)
    .save_plot(p, save_dir, "multi_ROC", width, height, format)
    cat(sprintf("\n[OK] Plot saved to: %s\n", save_dir))
  }
  
  return(p)
}

#' Plot Confusion Matrix for a Single Model
#'
#' Generates a confusion matrix heatmap for a trained model.
#' Supports both direct class predictions (type = "raw") and probability-based
#' predictions with user-defined threshold. Additionally calculates positive and
#' negative predictive values, optionally adjusted for a given population prevalence.
#'
#' @param object      A Train_Model object with trained models.
#' @param model_name  Name of the model to plot. If NULL, uses first model.
#' @param test_data   Test data (optional). If NULL, uses object@split.data$test.
#' @param threshold   Numeric threshold for probability predictions (0-1).
#'                    If provided, predictions are obtained via type = "prob"
#'                    and binarised using this threshold.
#'                    If NULL (default), type = "raw" is used.
#' @param prob_class  Character string indicating which class is considered
#'                    positive when using threshold. If NULL, the second level
#'                    of the true labels is used as positive class.
#' @param prevalence  Optional population prevalence (0-1) of the positive class.
#'                    If provided, the function calculates adjusted positive and
#'                    negative predictive values using Bayes' theorem.
#' @param palette     Vector of 2 colors (low, high) for the gradient.
#' @param base_size   ggplot2 base font size.
#' @param save_plot   Logical. Save the plot to file?
#' @param save_dir    Directory to save plot. If NULL, auto-detects.
#' @param width       Plot width in inches.
#' @param height      Plot height in inches.
#' @param format      File format ("pdf", "png", "jpg").
#'
#' @return A ggplot2 object showing the confusion matrix.
#'
#' @export
#' @examples
#' \dontrun{
#' model_obj <- ModelTrainAnalysis(
#' object       = model_obj_test,
#' methods      = c("glm", "rf", "gbm"),
#' control      = list(method = "repeatedcv", number = 5, repeats = 1),
#' save_plots   = TRUE,
#' save_dir     = ".",
#' seed         = 123
#' )
#' PlotConfusionMatrix(model_obj, model_name = names(model_obj@train.models)[1], save_plot = FALSE)
#' }
PlotConfusionMatrix <- function(object,
                                model_name = NULL,
                                test_data  = NULL,
                                threshold  = NULL,
                                prob_class = NULL,
                                prevalence = NULL,
                                palette    = c("#d8b365", "#5ab4ac"),
                                base_size  = 13,
                                save_plot  = FALSE,
                                save_dir   = NULL,
                                width      = 5,
                                height     = 4.5,
                                format     = "pdf") {
  
  # ---------------------------------------------------------------------------
  # Input Validation
  # ---------------------------------------------------------------------------
  if (is.null(save_dir)) {
    save_dir <- .get_viz_output_dir("Model")
  }
  
  if (!inherits(object, "Train_Model")) {
    stop("'object' must be a Train_Model.")
  }
  
  if (length(object@train.models) == 0) {
    stop("No trained models found.")
  }
  
  if (!is.null(prevalence) && (prevalence <= 0 || prevalence >= 1)) {
    warning("Prevalence should be between 0 and 1 (exclusive). Adjusted PPV/NPV may be unreliable.")
  }
  
  # ---------------------------------------------------------------------------
  # Select Model
  # ---------------------------------------------------------------------------
  nm <- if (is.null(model_name)) {
    cat(sprintf("No model specified. Using first model: %s\n", 
                names(object@train.models)[1]))
    names(object@train.models)[1]
  } else {
    model_name
  }
  
  if (!nm %in% names(object@train.models)) {
    stop(sprintf("Model '%s' not found. Available models: %s",
                 nm, paste(names(object@train.models), collapse = ", ")))
  }
  
  model <- object@train.models[[nm]]
  
  # ---------------------------------------------------------------------------
  # Prepare Test Data
  # ---------------------------------------------------------------------------
  td <- if (!is.null(test_data)) {
    test_data
  } else {
    object@split.scale.data$test
  }
  
  if (is.null(td)) {
    stop("Provide 'test_data' or populate split.data$test.")
  }
  
  # ---------------------------------------------------------------------------
  # Extract True Labels
  # ---------------------------------------------------------------------------
  group_col <- as.character(object@group_col)
  
  # Ensure true labels are factor
  truth_raw <- td[[group_col]]
  if (!is.factor(truth_raw)) {
    truth_raw <- factor(truth_raw)
  }
  
  truth_levels <- levels(truth_raw)
  
  if (length(truth_levels) != 2) {
    stop(sprintf("Expected 2 factor levels, got %d. This function is for binary classification.",
                 length(truth_levels)))
  }
  
  truth <- factor(truth_raw, levels = truth_levels)
  
  cat(sprintf("True labels - Levels: [%s]\n", paste(truth_levels, collapse = ", ")))
  
  # Determine positive class for threshold-based predictions
  if (is.null(prob_class)) {
    pos_class <- truth_levels[2]   # second level as positive by convention
    neg_class <- truth_levels[1]
  } else {
    if (!prob_class %in% truth_levels) {
      stop(sprintf("prob_class = '%s' not found in true label levels: %s",
                   prob_class, paste(truth_levels, collapse = ", ")))
    }
    pos_class <- prob_class
    neg_class <- setdiff(truth_levels, prob_class)[1]
  }
  
  # ---------------------------------------------------------------------------
  # Get Predictions (Threshold vs. Raw)
  # ---------------------------------------------------------------------------
  if (!is.null(threshold)) {
    # --- Use probability predictions and threshold ---
    cat(sprintf("Using probability predictions with threshold = %.3f\n", threshold))
    cat(sprintf("Positive class: %s | Negative class: %s\n", pos_class, neg_class))
    
    prob_pred <- tryCatch({
      stats::predict(model, newdata = td, type = "prob")
    }, error = function(e) {
      # Fallback for glm / similar (type = "response")
      stats::predict(model, newdata = td, type = "response")
    })
    
    # Convert to vector of probabilities for the positive class
    if (is.matrix(prob_pred) || is.data.frame(prob_pred)) {
      # Multi-column output (e.g., randomForest, rpart, lda)
      col_names <- colnames(prob_pred)
      if (ncol(prob_pred) == 2) {
        # Try to locate positive class column
        if (pos_class %in% col_names) {
          prob_pos <- prob_pred[, pos_class]
        } else {
          # Assume second column is positive (common default)
          cat("Assuming second column corresponds to positive class.\n")
          prob_pos <- prob_pred[, 2]
        }
      } else {
        stop("Probability output has more than 2 columns. Please specify 'prob_class'.")
      }
    } else if (is.numeric(prob_pred)) {
      # Single numeric vector (glm, etc.) - assumed to be probability of positive class
      prob_pos <- prob_pred
    } else {
      stop("Unable to extract positive class probabilities from model predictions.")
    }
    
    # Binarise using threshold
    pred_char <- ifelse(prob_pos > threshold, pos_class, neg_class)
    pred_raw <- factor(pred_char, levels = truth_levels)
    
    # For consistency with raw path, store as character for later mapping
    pred_char <- as.character(pred_raw)
    pred_levels <- levels(pred_raw)
    
    cat(sprintf("Predicted labels (thresholded) - Levels: [%s]\n", 
                paste(pred_levels, collapse = ", ")))
    
  } else {
    # --- Original raw prediction path ---
    pred_raw <- tryCatch({
      stats::predict(model, newdata = td, type = "raw")
    }, error = function(e) {
      cat(sprintf("Error getting predictions: %s\n", e$message))
      stop("Could not obtain predictions from model. Consider using 'threshold' for probability-based predictions.")
    })
    
    pred_char <- as.character(pred_raw)
    pred_levels <- unique(pred_char)
    cat(sprintf("Predicted labels (raw) - Unique values: [%s]\n", 
                paste(pred_levels, collapse = ", ")))
  }
  
  # ---------------------------------------------------------------------------
  # Handle Factor Level Mismatch (only needed for raw predictions; threshold
  # path already yields factor with correct levels, but we keep the logic
  # for robustness)
  # ---------------------------------------------------------------------------
  if (!identical(sort(unique(pred_char)), sort(truth_levels))) {
    cat("INFO: Prediction levels differ from truth levels. Attempting mapping...\n")
    
    # If predictions use different names (e.g., "0"/"1" vs "X0"/"X1")
    # try to map them based on order
    if (length(unique(pred_char)) == length(truth_levels)) {
      pred_sorted <- sort(unique(pred_char))
      truth_sorted <- sort(truth_levels)
      
      mapping <- setNames(truth_sorted, pred_sorted)
      cat(sprintf("Mapping: %s -> %s\n", 
                  paste(names(mapping), collapse=", "),
                  paste(mapping, collapse=", ")))
      
      pred_char <- mapping[as.character(pred_char)]
    }
  }
  
  # Convert to factor with truth's levels
  pred <- factor(pred_char, levels = truth_levels)
  
  if (length(pred) != length(truth)) {
    stop(sprintf("Prediction length (%d) != truth length (%d)",
                 length(pred), length(truth)))
  }
  
  cat(sprintf("Final predicted levels: [%s]\n\n", paste(levels(pred), collapse = ", ")))
  
  # ---------------------------------------------------------------------------
  # Create Confusion Matrix
  # ---------------------------------------------------------------------------
  cm <- table(Actual    = truth,
              Predicted = pred)
  
  cm_df <- as.data.frame(cm)
  colnames(cm_df) <- c("Actual", "Predicted", "Freq")
  
  cm_df <- cm_df %>%
    dplyr::group_by(Actual) %>%
    dplyr::mutate(Pct = round(Freq / sum(Freq) * 100, 1)) %>%
    dplyr::ungroup()
  
  # ---------------------------------------------------------------------------
  # Calculate Classification Metrics
  # ---------------------------------------------------------------------------
  # Extract counts using the already defined positive/negative classes
  tn <- cm[neg_class, neg_class]
  fp <- cm[neg_class, pos_class]
  fn <- cm[pos_class, neg_class]
  tp <- cm[pos_class, pos_class]
  
  total <- tp + tn + fp + fn
  accuracy <- (tp + tn) / total
  sensitivity <- tp / (tp + fn)   # also recall, TPR
  specificity <- tn / (tn + fp)   # TNR
  precision <- tp / (tp + fp)     # PPV
  npv <- tn / (tn + fn)           # Negative Predictive Value
  
  # Handle division by zero
  if (is.nan(sensitivity)) sensitivity <- 0
  if (is.nan(specificity)) specificity <- 0
  if (is.nan(precision)) precision <- 0
  if (is.nan(npv)) npv <- 0
  
  # Test set prevalence (observed)
  obs_prevalence <- (tp + fn) / total
  
  cat("----------------------------------------\n")
  cat("Performance Metrics (test set)\n")
  cat("----------------------------------------\n")
  cat(sprintf("TP = %d | TN = %d | FP = %d | FN = %d\n", tp, tn, fp, fn))
  cat(sprintf("Total N = %d\n", total))
  cat(sprintf("Accuracy     = %.3f (%.1f%%)\n", accuracy, accuracy * 100))
  cat(sprintf("Sensitivity  = %.3f (%.1f%%)   [Recall / TPR]\n", sensitivity, sensitivity * 100))
  cat(sprintf("Specificity  = %.3f (%.1f%%)   [TNR]\n", specificity, specificity * 100))
  cat(sprintf("PPV (Precision) = %.3f (%.1f%%)\n", precision, precision * 100))
  cat(sprintf("NPV             = %.3f (%.1f%%)\n", npv, npv * 100))
  cat(sprintf("Prevalence (test) = %.3f (%.1f%%)\n", obs_prevalence, obs_prevalence * 100))
  
  # ---------------------------------------------------------------------------
  # Adjusted predictive values using given population prevalence
  # ---------------------------------------------------------------------------
  if (!is.null(prevalence) && prevalence > 0 && prevalence < 1) {
    # Bayes' theorem adjustments
    # P(D+|T+) = (sens * prev) / (sens * prev + (1-spec) * (1-prev))
    adj_ppv <- (sensitivity * prevalence) / 
      (sensitivity * prevalence + (1 - specificity) * (1 - prevalence))
    # P(D-|T-) = (spec * (1-prev)) / ((1-sens) * prev + spec * (1-prev))
    adj_npv <- (specificity * (1 - prevalence)) / 
      ((1 - sensitivity) * prevalence + specificity * (1 - prevalence))
    
    if (is.nan(adj_ppv)) adj_ppv <- 0
    if (is.nan(adj_npv)) adj_npv <- 0
    
    cat("\n----------------------------------------\n")
    cat(sprintf("Adjusted for population prevalence = %.3f (%.1f%%)\n", prevalence, prevalence * 100))
    cat("----------------------------------------\n")
    cat(sprintf("Adjusted PPV (true positive prediction) = %.3f (%.1f%%)\n", adj_ppv, adj_ppv * 100))
    cat(sprintf("Adjusted NPV (true negative prediction) = %.3f (%.1f%%)\n", adj_npv, adj_npv * 100))
  }
  cat("\n")
  
  # ---------------------------------------------------------------------------
  # Plot Confusion Matrix
  # ---------------------------------------------------------------------------
  p <- ggplot2::ggplot(cm_df,
                       ggplot2::aes(x = Actual, y = Predicted, fill = Freq)) +
    ggplot2::geom_tile(colour = "white", linewidth = 1) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(Freq, "\n(", Pct, "%)")),
      size = 4,
      fontface = "bold",
      colour = "black"
    ) +
    ggplot2::scale_fill_gradient(
      low = palette[1],
      high = palette[2],
      name = "Count"
    ) +
    ggplot2::scale_x_discrete(limits = truth_levels) +
    ggplot2::scale_y_discrete(limits = rev(truth_levels)) +
    ggplot2::labs(
      title = paste("Confusion Matrix --", nm),
      subtitle = sprintf("Accuracy: %.1f%% | Sensitivity: %.1f%% | Specificity: %.1f%% | PPV: %.1f%% | NPV: %.1f%%",
                         accuracy * 100, sensitivity * 100, specificity * 100,
                         precision * 100, npv * 100),
      x = "Actual",
      y = "Predicted"
    ) +
    .pub_theme(base_size) +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = base_size),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, colour = "grey50")
    )
  
  # ---------------------------------------------------------------------------
  # Save Plot
  # ---------------------------------------------------------------------------
  if (save_plot) {
    .safe_dir(save_dir)
    filename <- paste0("confusion_", nm, 
                       ifelse(!is.null(threshold), paste0("_thr", gsub("\\.", "", sprintf("%.2f", threshold))), ""),
                       ifelse(!is.null(prevalence), paste0("_prev", gsub("\\.", "", sprintf("%.2f", prevalence))), ""),
                       ".", format)
    filepath <- file.path(save_dir, filename)
    ggplot2::ggsave(filepath, p, width = width, height = height, device = format)
    cat(sprintf("[OK] Plot saved to: %s\n", filepath))
  }
  
  return(p)
}

# Helper functions (.get_viz_output_dir, .safe_dir, .pub_theme) remain unchanged.
#' Feature importance plot
#'
#' Bar chart of variable importance scores across trained models. Each model
#' produces one facet panel.
#'
#' @param object      A \code{Train_Model} object.
#' @param top_n       Show top N features per model. Default \code{20}.
#' @param palette_name Colour palette.
#' @param base_size   Base font size.
#' @param save_plot   Logical.
#' @param save_dir    Output directory.
#' @param width,height Inches.
#' @param format      File format.
#' @returns A \code{ggplot} object.
#' @export
#' @export
#' @examples
#' \dontrun{
#' model_obj <- ModelTrainAnalysis(
#' object       = train_obj_test,
#' methods      = c("glm", "rf", "gbm"),
#' control      = list(method = "repeatedcv", number = 5, repeats = 1),
#' save_plots   = TRUE,
#' save_dir     = ".",
#' seed         = 123
#' )
#' PlotFeatureImportance(model_obj, top_n = 10, save_plot = FALSE)
#' }
PlotFeatureImportance <- function(object,
                                  top_n        = 20,
                                  palette_name = "Zissou1",
                                  base_size    = 12,
                                  save_plot    = FALSE,
                                  save_dir = NULL, 
                                  width        = 10,
                                  height       = 6,
                                  format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Model")

  if (!inherits(object, "Train_Model"))
    stop("'object' must be a Train_Model.")
  if (length(object@train.models) == 0)
    stop("No trained models found.")

  imp_list <- lapply(names(object@train.models), function(nm) {
    model <- object@train.models[[nm]]
    imp   <- tryCatch(
      caret::varImp(model, scale = TRUE)$importance,
      error = function(e) NULL
    )
    if (is.null(imp)) return(NULL)
    imp_df <- data.frame(Feature    = rownames(imp),
                         Importance = imp[, 1],
                         Model      = nm,
                         row.names  = NULL)
    head(imp_df[order(-imp_df$Importance), ], top_n)
  })

  imp_df <- dplyr::bind_rows(imp_list[!sapply(imp_list, is.null)])
  if (nrow(imp_df) == 0) stop("Could not extract variable importance.")

  imp_df <- imp_df %>%
    dplyr::group_by(Model) %>%
    dplyr::mutate(Feature = factor(Feature, levels = rev(unique(Feature)))) %>%
    dplyr::ungroup()

  cols <- .get_palette(palette_name, length(unique(imp_df$Model)))

  p <- ggplot2::ggplot(imp_df,
    ggplot2::aes(Importance, Feature, fill = Model)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::facet_wrap(~ Model, scales = "free") +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(title = paste0("Top ", top_n, " Feature Importance"),
                  x = "Importance (scaled)", y = NULL) +
    .pub_theme(base_size)

  if (save_plot)
    .save_plot(p, save_dir, "feature_importance", width, height, format)
  return(p)
}


#' Calibration Curve and Prediction Distribution Plot
#'
#' Evaluates the calibration of a binary classifier by comparing predicted
#' probabilities with observed outcomes. Produces a calibration curve
#' (observed proportion vs. mean predicted probability) and a histogram
#' of predicted probabilities. The smooth curve can be customized
#' (linear, LOESS, or none) and binning automatically adapts to sample
#' size and probability distribution.
#'
#' @param object A \code{Train_Model} object containing trained models.
#' @param model_name Character. Name of the model to evaluate. If \code{NULL},
#'   uses the first model in \code{object@train.models}.
#' @param test_data Data frame. Test data used for evaluation. If \code{NULL},
#'   uses \code{object@split.data$test}.
#' @param n_bins Integer or \code{NULL}. Number of bins for calibration curve.
#'   If \code{NULL}, automatically chosen based on sample size (range 4–10).
#' @param bin_method Character. One of \code{"auto"}, \code{"equal_width"},
#'   or \code{"quantile"}. If \code{"auto"}, switches to quantile binning
#'   when predicted probabilities are highly concentrated (spread < 0.15 or
#'   unique values < 10), otherwise uses equal-width bins.
#' @param hist_colors Character vector of length 2 for histogram fill colors.
#'   Default \code{c("#f0f0f0", "#1b9e77")}.
#' @param curve_colors Character vector of length 2: point color and smooth
#'   line color. Default \code{c("#1b9e77", "#c51b8a")}.
#' @param base_size Numeric. Base font size for \code{ggprism::theme_prism()}.
#'   Default 14.
#' @param combine_plots Logical. If \code{TRUE}, returns a patchwork composite
#'   of histogram + calibration curve. Default \code{TRUE}.
#' @param save_plot Logical. Save the plot to a PDF file. Default \code{FALSE}.
#' @param save_dir Character. Directory to save the plot. Required if
#'   \code{save_plot = TRUE}.
#' @param se Logical. Display standard error ribbon for the smooth curve.
#'   Default \code{FALSE}.
#' @param width Numeric. Width of saved plot in inches. Default 10.
#' @param height Numeric. Height of saved plot in inches. Default 5.
#' @param format Character. File format for saved plot ("pdf", "png", etc.).
#'   Default "pdf".
#' @param show_stats_on_plot Logical. Annotate the calibration curve with
#'   Brier score, calibration intercept, slope, and Eavg. Default \code{TRUE}.
#' @param smooth_method Character. Method for the smooth curve:
#'   \code{"lm"} (linear), \code{"loess"} (local polynomial), or \code{"none"}.
#'   Default \code{"lm"}.
#'
#' @return If \code{combine_plots = TRUE}, invisibly returns a \code{patchwork}
#'   object. Otherwise, returns a list with \code{histogram} and
#'   \code{calibration} ggplot objects.
#'
#' @examples
#' \dontrun{
#' library(caret)
#' data("mtcars")
#' mtcars$am <- as.factor(mtcars$am)
#' model <- CreateModelObject(data = mtcars, group_col = "am")
#' set.seed(123)
#' idx <- sample(1:nrow(mtcars), 20)
#' model@split.data <- list(training = mtcars[idx, ], testing = mtcars[-idx, ])
#' trained <- ModelTrainAnalysis(model, methods = c("glm", "rf"),
#'                               control = list(method = "cv", number = 3),
#'                               save_plots = FALSE)
#' PlotCalibration(trained, model_name = "glm", test_data = mtcars[-idx, ],
#'                 n_bins = 8, smooth_method = "lm", save_plot = FALSE)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_abline geom_point
#'   geom_smooth geom_text scale_fill_manual scale_size_continuous
#'   coord_equal xlim ylim labs theme element_blank element_text
#' @importFrom ggprism theme_prism
#' @importFrom dplyr group_by summarise n
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom stats predict coef sd quantile
#' @export
PlotCalibration <- function(object,
                            model_name = NULL,
                            test_data = NULL,
                            n_bins = NULL,
                            bin_method = c("auto", "equal_width", "quantile"),
                            hist_colors = c("#f0f0f0", "#1b9e77"),
                            curve_colors = c("#1b9e77", "#c51b8a"),
                            base_size = 14,
                            combine_plots = TRUE,
                            save_plot = FALSE,
                            save_dir = NULL,
                            se = FALSE,
                            width = 10,
                            height = 5,
                            format = "pdf",
                            show_stats_on_plot = TRUE,
                            smooth_method = c("lm", "loess", "none")) {
  
  smooth_method <- match.arg(smooth_method)
  bin_method <- match.arg(bin_method)
  
  # ---- 1. Parameter Extraction ----
  if (!inherits(object, "Train_Model"))
    stop("'object' must be a Train_Model.")
  
  if (is.null(model_name))
    model_name <- names(object@train.models)[1]
  model <- object@train.models[[model_name]]
  
  if (is.null(test_data))
    test_data <- object@split.data$test
  
  group_col <- as.character(object@group_col)
  if (!group_col %in% colnames(test_data))
    stop("group_col '", group_col, "' not found in test_data.")
  
  # ---- 2. Truth Label Extraction & Encoding ----
  truth_factor <- test_data[[group_col]]
  if (!is.factor(truth_factor))
    truth_factor <- as.factor(truth_factor)
  
  levels_true <- levels(truth_factor)
  cat("===== True Class Levels =====\n")
  print(levels_true)
  
  if (length(levels_true) != 2)
    warning("Grouping variable is not binary; calibration curve may not be applicable.")
  
  pos_level <- levels_true[2]
  truth_numeric <- as.integer(truth_factor) - 1L  # 0/1 encoding
  cat("\nPositive Level (Assumed levels[2]):", pos_level, "\n")
  cat("Positive Samples:", sum(truth_numeric == 1), "/", length(truth_numeric), "\n")
  
  # ---- 3. Predicted Probability Extraction ----
  pred_prob <- tryCatch({
    prob_mat <- stats::predict(model, newdata = test_data, type = "prob")
    if (is.matrix(prob_mat) || is.data.frame(prob_mat)) {
      if (pos_level %in% colnames(prob_mat)) {
        prob_mat[, pos_level]
      } else {
        warning("Prediction matrix column names do not include positive level; using second column.")
        prob_mat[, 2]
      }
    } else if (is.numeric(prob_mat)) {
      prob_mat
    } else {
      stop("Unable to parse return type of predict()")
    }
  }, error = function(e) {
    stop("Prediction failed: ", e$message)
  })
  
  # ---- 4. Basic Statistics ----
  cat("\n===== Predicted Probability Distribution =====\n")
  print(summary(pred_prob))
  cat("Standard Deviation:", sd(pred_prob), "\n")
  
  brier <- mean((truth_numeric - pred_prob)^2)
  cat("\nBrier Score (Lower is better, range 0-0.25):", round(brier, 4), "\n")
  
  # Calibration Slope & Intercept (Logit Calibration)
  cal_df <- data.frame(truth = truth_numeric, prob = pred_prob)
  cal_df$prob_clip <- pmax(pmin(cal_df$prob, 1 - 1e-6), 1e-6)
  cal_glm <- suppressWarnings(
    glm(truth ~ log(prob_clip/(1 - prob_clip)), family = binomial(), data = cal_df)
  )
  intercept <- stats::coef(cal_glm)[1]
  slope     <- stats::coef(cal_glm)[2]
  
  cat("\n===== Calibration Parameters (Logit) =====\n")
  cat("Intercept:", round(intercept, 3), "(Ideal: 0)\n")
  cat("Slope:    ", round(slope, 3), "(Ideal: 1)\n")
  
  if (abs(intercept) > 0.5) cat("Warning: Large intercept suggests global probability shift.\n")
  if (abs(slope - 1) > 0.2) cat("Warning: Slope deviation suggests poor calibration at extremes.\n")
  
  # ---- 5. Intelligent Binning ----
  if (is.null(n_bins)) {
    n <- length(pred_prob)
    base_bins <- min(10, max(4, floor(n / 15)))
    
    if (bin_method == "auto") {
      n_unique <- length(unique(round(pred_prob, 2)))
      spread <- diff(quantile(pred_prob, probs = c(0.1, 0.9), na.rm = TRUE))
      
      if (spread < 0.15 || n_unique < 10) {
        bin_method <- "quantile"
        n_bins <- min(6, max(3, base_bins - 1))
        cat(sprintf("Auto-switched to quantile binning (spread=%.2f, n_unique=%d), n_bins=%d\n",
                    spread, n_unique, n_bins))
      } else {
        bin_method <- "equal_width"
        n_bins <- base_bins
        cat(sprintf("Auto-selected equal-width binning, n_bins=%d\n", n_bins))
      }
    } else {
      n_bins <- base_bins
    }
  }
  
  # Create bins
  if (bin_method == "quantile") {
    breaks <- unique(quantile(pred_prob, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
    if (min(breaks) > 0) breaks <- c(0, breaks)
    if (max(breaks) < 1) breaks <- c(breaks, 1)
    breaks <- unique(breaks)
  } else {
    breaks <- seq(0, 1, length.out = n_bins + 1)
  }
  cal_df$bin <- cut(pred_prob, breaks = breaks, include.lowest = TRUE)
  
  cal_sum <- cal_df %>%
    dplyr::group_by(bin) %>%
    dplyr::summarise(
      mean_pred = mean(prob),
      obs_rate  = mean(truth),
      n         = dplyr::n(),
      .groups   = "drop"
    )
  
  cat("\n===== Binning Calibration Summary (n_bins =", n_bins, ") =====\n")
  print(cal_sum)
  
  e_avg <- mean(abs(cal_sum$obs_rate - cal_sum$mean_pred), na.rm = TRUE)
  cat("\nMean Absolute Calibration Error (Eavg):", round(e_avg, 4), "\n")
  
  # ---- 6. Visualization ----
  hist_df <- data.frame(prob = pred_prob,
                        truth = factor(truth_numeric, levels = c(0,1), labels = levels_true))
  
  p1 <- ggplot2::ggplot(hist_df, ggplot2::aes(x = prob, fill = truth)) +
    ggplot2::geom_histogram(alpha = 0.6, bins = 30, position = "stack",
                            color = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = hist_colors, name = "True Class") +
    ggplot2::labs(title = paste0("Probability Distribution: ", model_name),
                  x = "Predicted Probability", y = "Count") +
    ggprism::theme_prism(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom")
  
  p2 <- ggplot2::ggplot(cal_sum, ggplot2::aes(x = mean_pred, y = obs_rate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, alpha = 0.6,
                         linetype = "dashed", color = "#d9d9d9", linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(size = n), color = curve_colors[1], alpha = 0.8)
  
  if (smooth_method != "none") {
    p2 <- p2 + ggplot2::geom_smooth(method = smooth_method, se = se,
                                    color = curve_colors[2], linewidth = 1.5,
                                    linetype = "dotted")
  }
  
  p2 <- p2 +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "n") +
    ggplot2::coord_equal() +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::labs(title = paste0("Calibration Curve: ", model_name),
                  subtitle = paste0("Bins = ", n_bins, " (", bin_method, ")"),
                  x = "Mean Predicted Probability", y = "Observed Proportion") +
    ggprism::theme_prism(base_size = base_size)
  
  if (show_stats_on_plot) {
    stats_text <- paste0(
      "Brier: ", round(brier, 4), "\n",
      "Intercept: ", round(intercept, 3), "\n",
      "Slope: ", round(slope, 3), "\n",
      "Eavg: ", round(e_avg, 4)
    )
    p2 <- p2 + ggplot2::annotate("text", x = 0.05, y = 0.95, label = stats_text,
                                 hjust = 0, vjust = 1, size = base_size * 0.3,
                                 family = "mono", fontface = "bold")
  }
  
  # ---- 7. Output Management ----
  if (combine_plots) {
    if (!requireNamespace("patchwork", quietly = TRUE))
      stop("Please install 'patchwork' package.")
    combined <- p1 + p2 + patchwork::plot_layout(widths = c(1, 1))
    
    if (save_plot) {
      if (is.null(save_dir)) save_dir <- getwd()
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      out_path <- file.path(save_dir, paste0("Calibration_", model_name, ".", format))
      ggplot2::ggsave(out_path, combined, width = width, height = height)
      cat("\nPlot saved to:", out_path, "\n")
    }
    
    print(combined)
    return(invisible(combined))
  } else {
    return(list(histogram = p1, calibration = p2))
  }
}




## =============================================================================
##   3  SUBTYPING  -- dim-reduction - cluster heatmap - silhouette - alluvial
## =============================================================================

#' t-SNE / UMAP scatter coloured by cluster and features
#'
#' Draws 2-D embedding points coloured by cluster label plus optional
#' metadata overlay. Returns a \code{patchwork} composite when
#' \code{overlay_features} is supplied.
#'
#' @param object           A \code{Subtyping} object.
#' @param reduction        \code{"tsne"} or \code{"umap"}.
#' @param color_by         \code{"cluster"} or a column name in
#'   \code{info.data}.
#' @param overlay_features Character vector of features to colour additional
#'   panels.
#' @param point_size       Point size. Default \code{1.8}.
#' @param palette_name     Palette.
#' @param base_size        Base font size.
#' @param save_plot        Logical.
#' @param save_dir         Output directory.
#' @param width,height     Inches.
#' @param format           File format.
#' @returns A \code{ggplot} or \code{patchwork} object.
#' @export
#' @examples
#' \dontrun{
#' subtype_obj_test=Sub_tsne_analyse(subtype_obj_test)
#' subtype_obj_test=Sub_umap_analyse(subtype_obj_test)
#' PlotDimReduction(subtype_obj_test,reduction='tsne',color_by='SWAB')
#' PlotDimReduction(subtype_obj_test,reduction='umap',color_by='SWAB')
#' }

PlotDimReduction <- function(object,
                             reduction        = c("tsne", "umap"),
                             color_by         = "cluster",
                             overlay_features = NULL,
                             point_size       = 1.8,
                             palette_name     = "Darjeeling1",
                             base_size        = 13,
                             save_plot        = FALSE,
                             save_dir = NULL,  # default: auto-detect from config
                             width            = 6,
                             height           = 5,
                             format           = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("Subtyping")

  if (!inherits(object, "Subtyping")) stop("'object' must be a Subtyping object.")
  reduction <- match.arg(reduction)

  vis <- object@visualization.results
  key <- paste0(reduction, ".df")
  if (is.null(vis) || (!key %in% names(vis) && !reduction %in% names(vis)))
    stop("Run Sub_tsne_analyse() or Sub_umap_analyse() first.")
  
  if (key %in% names(vis)) {
    coords_df <- vis[[key]]
  } else {
    coords_df <- vis[[reduction]]$coords
    if (is.null(coords_df))
      coords_df <- vis[[reduction]][, 1:2]
  }
  
  if(ncol(coords_df) >= 2) {
    colnames(coords_df)[1:2] <- c("Dim1", "Dim2")
  }

  cdata <- object@clustered.data
  if (!is.null(cdata) && nrow(cdata) > 0) {
    if ("cluster" %in% colnames(cdata)) {
      coords_df$cluster <- as.factor(cdata[rownames(coords_df), "cluster"])
    } else if ("group" %in% colnames(cdata)) {
      coords_df$cluster <- as.factor(cdata[rownames(coords_df), "group"])
    }
  }
  info <- object@info.data
  if (nrow(info) > 0)
    coords_df <- cbind(coords_df, info[rownames(coords_df), , drop = FALSE])

  ax <- toupper(reduction)
  .one_panel <- function(clr_col, is_factor = TRUE, title_extra = "") {
    if (is_factor) {
      n_lv <- length(unique(coords_df[[clr_col]]))
      cols  <- .get_palette(palette_name, n_lv)
      p <- ggplot2::ggplot(coords_df,
        ggplot2::aes(Dim1, Dim2, colour = as.factor(.data[[clr_col]]))) +
        ggplot2::scale_colour_manual(values = cols,
                                     name   = clr_col)
    } else {
      p <- ggplot2::ggplot(coords_df,
        ggplot2::aes(Dim1, Dim2,
                     colour = as.numeric(.data[[clr_col]]))) +
        ggplot2::scale_colour_viridis_c(name = clr_col, option = "D")
    }
    p + ggplot2::geom_point(size = point_size, alpha = 0.8) +
      ggplot2::labs(title = paste0(ax, title_extra),
                    x = paste(ax, "1"), y = paste(ax, "2")) +
      .pub_theme(base_size)
  }

  main_p <- .one_panel(color_by,
                        is_factor = !color_by %in% names(coords_df)[
                          sapply(coords_df, is.numeric)])

  if (!is.null(overlay_features) && length(overlay_features) > 0) {
    feat_panels <- lapply(overlay_features, function(ft) {
      if (!ft %in% colnames(coords_df)) {
        cd <- object@clean.data
        if (ft %in% colnames(cd))
          coords_df[[ft]] <<- cd[rownames(coords_df), ft]
        else return(NULL)
      }
      .one_panel(ft, is_factor = FALSE,
                  title_extra = paste0(" - ", ft))
    })
    feat_panels <- feat_panels[!sapply(feat_panels, is.null)]
    all_panels  <- c(list(main_p), feat_panels)
    out_p       <- patchwork::wrap_plots(all_panels)
  } else {
    out_p <- main_p
  }

  if (save_plot)
    .save_plot(out_p, save_dir, paste0(reduction, "_", color_by),
               width * max(1, 1 + length(overlay_features)), height, format)
  return(out_p)
}



#' Plot Cluster Heatmap for Clinlabomics
#'
#' @description
#' Generates an enhanced cluster heatmap designed for clinical laboratory omics data. 
#' Features automated Z-score thresholding, custom log transformations, split-group ordering, 
#' and multi-layered top annotations for clinical metadata.
#'
#' @param object An S4 object, a list, or a matrix containing data. If an S4 object or a list, 
#'   it must contain slots or elements named \code{clean.data} (data matrix) and \code{info.data} (metadata frame).
#' @param deg_df A data frame containing Differential Expression Analysis results. Must include 
#'   columns: \code{id}, \code{p.adjust}, \code{logFC}, and \code{target_group}.
#' @param group_by Character. The column name in metadata used for primary patient clustering. 
#'   Default is \code{"cluster_lpa"}.
#' @param annotation_cols Character vector. Optional column names from the metadata to display 
#'   as tracking tracks (e.g., \code{c("Age", "Gender", "Stage")}) above the heatmap.
#' @param top_n Integer. Maximum number of top features to slice per cluster based on LogFC. Default is 3.
#' @param p_cutoff Numeric. Adjusted p-value significance threshold for feature selection. Default is 0.1.
#' @param logfc_cutoff Numeric. Log fold-change threshold for feature selection. Default is 0.
#' @param log_transform Logical. If \code{TRUE}, applies a \code{log10(x + 1)} transformation to the matrix before scaling.
#' @param custom_levels Character vector. Explicitly sets the factor level order of the clusters.
#' @param state_palette Named character vector. Colors assigned to the primary cluster levels.
#' @param annotation_palette A named list of color vectors or color mapping functions corresponding 
#'   to columns in \code{annotation_cols}.
#' @param heatmap_palette Character vector of 3 colors representing the low, middle, and high values of the Z-score.
#' @param z_limit Numeric. Fixed cutoff limit for standard deviations. If \code{NULL}, automatically uses the 99.5th percentile.
#' @param save_path Character. Output file path for the generated PDF plot.
#' @param show_gene_names Logical. If \code{TRUE}, utilizes \code{anno_mark} to display feature labels cleanly on the right.
#' @param raster_quality Numeric. Quality multiplier for rasterization. Default is 5.
#'
#' @return Invisibly returns a structured \code{\link[ComplexHeatmap:Heatmap-class]{Heatmap}}object.
#' @export
#' 
#' @importFrom methods slot
#' @importFrom gtools mixedsort
#' @importFrom dplyr mutate group_by arrange desc slice_head ungroup
#' @importFrom scales hue_pal
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit grid.rect
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_mark draw
#' @examples
#' \dontrun{
#' # 1. Prepare dummy omics matrix
#' set.seed(42)
#' mock_mat <- matrix(matrix(rnorm(200, mean = 2)), ncol = 10)
#' colnames(mock_mat) <- c("wbc", "lymphocytes", "neutrophils", "crp", "ggt", paste0("Lab", 6:10))
#' rownames(mock_mat) <- paste0("Patient_", 1:20)
#' 
#' # 2. Prepare metadata with clinical tracks
#' mock_meta <- data.frame(
#'   PatientID = rownames(mock_mat),
#'   cluster_lpa = rep(c("Cluster_1", "Cluster_2"), each = 10),
#'   Age = sample(30:80, 20, replace = TRUE),
#'   Gender = sample(c("Male", "Female"), 20, replace = TRUE),
#'   row.names = rownames(mock_mat)
#' )
#' 
#' # 3. Prepare dummy DEG frame
#' mock_deg <- data.frame(
#'   id = colnames(mock_mat),
#'   p.adjust = runif(10, 0, 0.02),
#'   logFC = c(runif(5, 1, 3), runif(5, -3, -1)),
#'   target_group = rep(c("Cluster_1", "Cluster_2"), each = 5)
#' )
#' 
#' # 4. Generate Heatmap with clinical info tracking
#' PlotClusterHeatmap(
#'   object = list(clean.data = mock_mat, info.data = mock_meta),
#'   deg_df = mock_deg,
#'   annotation_cols = c("Age", "Gender"),
#'   annotation_palette = list(
#'     Gender = c("Male" = "#4682B4", "Female" = "#FF69B4"),
#'     Age = circlize::colorRamp2(c(30, 80), c("white", "#E65100"))
#'   ),
#'   save_path = "./Clinlabomics_Heatmap_Output.pdf"
#' )
#' }
PlotClusterHeatmap <- function(
    object,
    deg_df,
    group_by           = "cluster_lpa",
    annotation_cols    = NULL,  
    top_n              = 3,
    p_cutoff           = 0.1,
    logfc_cutoff       = 0,
    log_transform      = FALSE,  
    custom_levels      = NULL,
    state_palette      = NULL,
    annotation_palette = NULL, 
    heatmap_palette    = c("#2166AC", "black", "#FFEC00"),
    z_limit            = NULL,   
    save_path          = "./Cluster_Heatmap_Final.pdf",
    show_gene_names    = TRUE,
    raster_quality     = 5
) {
  # Safe Package Checks 
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop("Package 'ComplexHeatmap' required.")
  if (!requireNamespace("circlize", quietly = TRUE))       stop("Package 'circlize' required.")
  if (!requireNamespace("dplyr", quietly = TRUE))          stop("Package 'dplyr' required.")
  if (!requireNamespace("gtools", quietly = TRUE))         stop("Package 'gtools' required.")
  
  # [1] Data Parsing -------------------------------------------------------
  if (isS4(object)) {
    mat_df  <- methods::slot(object, "clean.data")
    meta_df <- methods::slot(object, "info.data")
  } else if (is.list(object) && !is.data.frame(object)) {
    mat_df  <- object$clean.data
    meta_df <- object$info.data
  } else {
    mat_df  <- object
    meta_df <- object
  }
  
  # [2] Optional Log Transformation -----------------------------------------
  if (log_transform) {
    message("Applying log10(x + 1) transformation...")
    mat_df <- log10(mat_df + 1)
  }
  
  # [3] Group Ordering -----------------------------------------------------
  if (!group_by %in% colnames(meta_df)) {
    # Previously, a missing group_by column silently produced a ZERO-LENGTH
    # vector here (as.character(NULL) == character(0), not a column of NAs),
    # which then propagated through group-level computation, row ordering,
    # and z-score scaling, only surfacing 5-6 steps later as an unrelated
    # error deep inside plotting code (scales::hue_pal(): "Must request at
    # least one colour from a hue palette."), with no indication that the
    # real cause was a missing/misspelled column. Fail clearly here instead.
    cluster_like <- grep("^cluster_", colnames(meta_df), value = TRUE)
    stop(
      "group_by = \"", group_by, "\" is not a column of this object's ",
      "info.data (it has ", ncol(meta_df), " columns). ",
      if (length(cluster_like) > 0) {
        paste0("Available clustering column(s): ", paste(cluster_like, collapse = ", "), ". ")
      } else {
        paste0("No columns starting with \"cluster_\" were found -- make sure you've ",
               "run a clustering step (e.g. Sub_nmf_assign_subtypes(), ",
               "Sub_kmeans_with_optimal_k(), Sub_lpa_with_optimal_k()) before calling ",
               "PlotClusterHeatmap(). ")
      },
      "Pass the correct column via `group_by =`.",
      call. = FALSE
    )
  }
  raw_group <- as.character(meta_df[[group_by]])
  if (all(is.na(raw_group))) {
    stop("group_by = \"", group_by, "\" exists but is entirely NA/missing -- ",
         "no samples have been assigned to a cluster under this column yet.",
         call. = FALSE)
  }
  lvls <- if (!is.null(custom_levels)) custom_levels else gtools::mixedsort(unique(raw_group))
  
  # [4] Feature Filtering & Slicing ----------------------------------------
  available_features <- colnames(mat_df)
  filtered_deg <- deg_df %>%
    dplyr::filter(id %in% available_features) %>% 
    dplyr::filter(p.adjust < p_cutoff, logFC > logfc_cutoff) %>%
    dplyr::mutate(target_group = factor(as.character(target_group), levels = lvls)) %>%
    dplyr::group_by(target_group) %>%
    dplyr::arrange(dplyr::desc(logFC), .by_group = TRUE) %>% 
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(target_group)
  
  if (nrow(filtered_deg) == 0) stop("No matching clinical features found with the current thresholds.")
  
  # [5] Matrix Computations & Z-score Scaling -------------------------------
  sample_order_idx <- order(factor(raw_group, levels = lvls))
  ordered_meta     <- meta_df[sample_order_idx, , drop = FALSE]
  ordered_group_tags <- factor(raw_group[sample_order_idx], levels = lvls)
  
  plot_mat_raw <- as.matrix(mat_df[rownames(ordered_meta), filtered_deg$id])
  plot_mat_z   <- scale(plot_mat_raw)
  
  # [6] Adaptive Z-limit Calculation ---------------------------------------
  if (is.null(z_limit)) {
    z_limit <- quantile(abs(plot_mat_z), 0.995, na.rm = TRUE)
    message(paste0("Auto-calculated z_limit (99.5th percentile): ", round(z_limit, 2)))
  }
  
  plot_mat <- t(plot_mat_z)
  plot_mat[plot_mat > z_limit]  <- z_limit
  plot_mat[plot_mat < -z_limit] <- -z_limit
  
  # [7] Visual Legends & Dynamic Palette Allocation -----------------------
  col_fun <- circlize::colorRamp2(c(-z_limit, 0, z_limit), heatmap_palette)
  
  if (is.null(state_palette)) {
    state_colors <- scales::hue_pal()(length(lvls))
  } else {
    state_colors <- state_palette
  }
  names(state_colors) <- lvls
  
  # Initialize the color map registry list for annotations
  color_registry_list <- list(Cluster = state_colors)
  if (!is.null(annotation_palette)) {
    color_registry_list <- c(color_registry_list, annotation_palette)
  }
  
  # Isolate metadata columns specified for track representation
  top_anno_df <- ordered_meta[, c(group_by, annotation_cols), drop = FALSE]
  colnames(top_anno_df)[1] <- "Cluster"
  
  common_lgd_param <- list(
    direction = "horizontal", nrow = 1, title_position = "topleft",
    title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    labels_gp = grid::gpar(fontsize = 9)
  )
  
  # Construct Top Clinical Tracks Layout
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    df = top_anno_df,
    col = color_registry_list,
    simple_anno_size = grid::unit(5, "mm"),
    annotation_legend_param = color_registry_list, 
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 11, fontface = "bold"),
    annotation_name_side = "left"
  )
  
  # [8] Construct Heatmap Body ----------------------------------------------
  ht <- ComplexHeatmap::Heatmap(
    plot_mat,
    name = "Z-score",
    col = col_fun,
    heatmap_legend_param = common_lgd_param,
    column_title = paste0("Clusters (", group_by, ")"),
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    column_split = ordered_group_tags,
    row_split    = filtered_deg$target_group, 
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    cluster_column_slices = FALSE,
    cluster_row_slices    = FALSE,
    cluster_columns = TRUE, 
    cluster_rows    = FALSE,
    show_column_names = FALSE,
    show_row_names    = FALSE, 
    column_gap = grid::unit(1.5, "mm"),
    row_gap    = grid::unit(1.5, "mm"),
    border = TRUE,
    top_annotation = top_annotation,
    right_annotation = if (show_gene_names) {
      ComplexHeatmap::rowAnnotation(
        mark = ComplexHeatmap::anno_mark(
          at = 1:nrow(plot_mat), 
          labels = filtered_deg$id,
          labels_gp = grid::gpar(fontsize = 12, fontface = "bold.italic")
        )
      )
    } else { NULL },
    use_raster = TRUE, 
    raster_quality = raster_quality,
    layer_fun = function(j, i, x, y, w, h, fill) {
      grid::grid.rect(gp = grid::gpar(lwd = 0.5, fill = "transparent", col = "white"))
    }
  )
  
  # [9] Compilation and File Export -----------------------------------------
  if (!is.null(save_path)) {
    if (!dir.exists(dirname(save_path))) dir.create(dirname(save_path), recursive = TRUE)
    pdf(save_path, width = 12, height = 9)
    ComplexHeatmap::draw(
      ht, 
      heatmap_legend_side = "bottom", 
      annotation_legend_side = "bottom", 
      legend_gap = grid::unit(15, "mm"), 
      merge_legend = TRUE
    )
    dev.off()
    message(paste("Success: PlotClusterHeatmap generated at:", save_path))
  }
  
  return(invisible(ht))
}




#' Plot Silhouette Score for Subtyping Results
#'
#' @description
#' Calculates and visualizes silhouette widths to evaluate clustering cohesion.
#' Automatically aligns cluster labels from metadata if missing in the numeric matrix.
#'
#' @param object A \code{Subtyping} S4 object.
#' @param group_by Character. The column name in \code{info.data} containing cluster labels. 
#'   Default is \code{"cluster_lpa"}.
#' @param dist_method Character. Distance metric for \code{stats::dist()}. Default is \code{"euclidean"}.
#' @param palette_name Character. Palette name (e.g., from Wes Anderson or RColorBrewer).
#' @param base_size Numeric. Base font size for the plot.
#' @param save_plot Logical. Whether to save the output as a file.
#' @param save_dir Character. Directory path to save the plot.
#' @param width,height,format Plot export settings.
#'
#' @return A \code{ggplot} object.
#' @export
PlotSilhouette <- function(
    object,
    group_by     = "cluster_lpa",
    dist_method  = "euclidean",
    palette_name = "Darjeeling1",
    base_size    = 13,
    save_plot    = FALSE,
    save_dir     = NULL,
    width        = 7,
    height       = 5,
    format       = "pdf"
) {
  # 1. Validation ----------------------------------------------------------
  if (!inherits(object, "Subtyping")) {
    stop("Input 'object' must be a Subtyping S4 object.")
  }
  
  # 2. Data Preparation ----------------------------------------------------
  # Extract numeric data and metadata
  mat  <- as.matrix(object@clean.data)
  info <- object@info.data
  
  if (!(group_by %in% colnames(info))) {
    stop(sprintf("Column '%s' not found in info.data. Available: %s", 
                 group_by, paste(colnames(info), collapse = ", ")))
  }
  
  # Ensure cluster labels are factor/numeric
  cl_labels <- info[[group_by]]
  if (any(is.na(cl_labels))) {
    stop("Cluster labels contain missing values (NA). Please check clustering results.")
  }
  
  # 3. Silhouette Calculation ----------------------------------------------
  # Calculate distance matrix
  d_mat   <- stats::dist(mat, method = dist_method)
  cl_int  <- as.integer(as.factor(cl_labels))
  
  sil     <- cluster::silhouette(cl_int, d_mat)
  sil_df  <- as.data.frame(sil[, ])
  
  # Reconstruct group names for the plot
  lvl_names <- levels(as.factor(cl_labels))
  sil_df$cluster_name <- lvl_names[sil_df$cluster]
  
  # Sort for visualization
  sil_df  <- sil_df[order(sil_df$cluster, sil_df$sil_width), ]
  sil_df$order <- seq_len(nrow(sil_df))
  avg_sil <- round(mean(sil_df$sil_width), 3)
  
  # 4. Visualization -------------------------------------------------------
  cols <- .get_palette(palette_name, length(lvl_names))
  
  p <- ggplot2::ggplot(sil_df, ggplot2::aes(x = order, y = sil_width, fill = cluster_name)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::geom_hline(yintercept = avg_sil, linetype = "dashed", color = "red") +
    ggplot2::annotate("text", x = 5, y = avg_sil + 0.05, 
                      label = paste0("Average Width: ", avg_sil), 
                      hjust = 0, fontface = "bold") +
    ggplot2::scale_fill_manual(values = cols, name = "Cluster") +
    ggplot2::labs(
      title = paste("Silhouette Analysis:", group_by),
      subtitle = paste("Distance Method:", dist_method),
      x = "Samples", 
      y = "Silhouette Width"
    ) +
    .pub_theme(base_size)
  
  # 5. Export --------------------------------------------------------------
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- "./"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    
    file_name <- paste0("Silhouette_", group_by, ".", format)
    dest_path <- file.path(save_dir, file_name)
    
    ggplot2::ggsave(dest_path, plot = p, width = width, height = height)
    message("Silhouette plot saved to: ", dest_path)
  }
  
  return(p)
}


#' Plot Multi-Group Alluvial Diagram
#'
#' @description
#' Visualizes sample transitions across multiple clustering methods using 
#' the 'Lodes' format for ggalluvial. 
#'
#' @param object A \code{Subtyping} S4 object or a data frame.
#' @param cols_list Character vector. Columns to use as axes (e.g., c("Kmeans", "LPA", "NMF")).
#' @param palette_name Character. Color palette name.
#' @param base_size Numeric. Base font size.
#' @param save_plot Logical. Save to file?
#' @param save_dir Character. Output directory.
#' @param width,height,format Export dimensions.
#'
#' @export
PlotMultiAlluvial <- function(
    object,
    cols_list    = c("cluster_kmeans", "cluster_lpa", "cluster_nmf"),
    palette_name = "Darjeeling1",
    base_size    = 13,
    save_plot    = FALSE,
    save_dir     = NULL,
    width        = 10,
    height       = 7,
    format       = "pdf"
) {
  # 1. Dependencies --------------------------------------------------------
  requireNamespace("ggalluvial", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  
  # 2. Data Extraction -----------------------------------------------------
  df_meta <- if (inherits(object, "Subtyping")) {
    as.data.frame(methods::slot(object, "info.data"))
  } else {
    as.data.frame(object)
  }
  
  if (!all(cols_list %in% colnames(df_meta))) {
    stop("Some columns in 'cols_list' were not found in the metadata.")
  }
  
  # 3. Transform to Long Format (Lodes) ------------------------------------
  # This is the "ggalluvial" way to handle multiple axes
  df_meta$ID <- seq_len(nrow(df_meta))
  
  plot_df <- df_meta %>%
    dplyr::select(ID, dplyr::all_of(cols_list)) %>%
    # Convert to lodes form: creates 'stratum', 'alluvium', and 'x' (axis name)
    ggalluvial::to_lodes_form(key = "Method", axes = cols_list) %>%
    dplyr::mutate(stratum = as.factor(stratum),
                  Method = factor(Method, levels = cols_list))
  
  # Identify initial clusters for coloring (to track flows from the first method)
  first_method_map <- df_meta %>% 
    dplyr::select(ID, fill_group = !!cols_list[1])
  
  plot_df <- plot_df %>%
    dplyr::left_join(first_method_map, by = "ID") %>%
    dplyr::mutate(fill_group = as.factor(fill_group))
  
  # 4. Visualization -------------------------------------------------------
  cols <- .get_palette(palette_name, nlevels(plot_df$fill_group))
  
  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = Method, 
                                    stratum = stratum, 
                                    alluvium = ID, 
                                    y = 1)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill = fill_group), 
                              width = 1/12, alpha = 0.7) +
    ggalluvial::geom_stratum(width = 1/12, fill = "grey90", color = "white") +
    ggplot2::geom_text(stat = ggalluvial::StatStratum, 
                       ggplot2::aes(label = ggplot2::after_stat(stratum)),
                       size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = cols, name = paste("Initial:", cols_list[1])) +
    ggplot2::labs(title = "Multi-Method Cluster Consistency",
                  subtitle = paste("Tracking flow from:", cols_list[1]),
                  x = "Clustering Method",
                  y = "Sample Count") +
    .pub_theme(base_size)
  
  # 5. Export --------------------------------------------------------------
  if (save_plot) {
    if (is.null(save_dir)) save_dir <- "./"
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    dest <- file.path(save_dir, paste0("multi_alluvial.", format))
    ggplot2::ggsave(dest, plot = p, width = width, height = height)
    message("Multi-group alluvial plot saved to: ", dest)
  }
  
  return(p)
}


## =============================================================================
##   4  PROGNOSIS  -- KM - forest - time-ROC - RCS - nomogram - calibration - DCA - risk
## =============================================================================

#' Kaplan-Meier Survival Curve with Clinical Group Support
#'
#' Draws a KM curve stratified by a grouping variable. If a \code{PrognosiX}
#' object is provided, the grouping variable can be from either \code{survival.data}
#' or \code{info.data}.
#'
#' @param object A \code{PrognosiX} object or a data frame with \code{time} and
#'   \code{status} columns.
#' @param group_col Character; name of the grouping column. If \code{NULL},
#'   a single overall curve is drawn.
#' @param time_unit Label for x‑axis (e.g., "days", "months").
#' @param conf_int Logical; show confidence bands.
#' @param pairwise_p Logical; add pairwise log‑rank p‑values.
#' @param palette_name Color palette name.
#' @param base_size Base font size.
#' @param save_plot Logical; save to PDF.
#' @param save_dir Output directory.
#' @param width,height Plot dimensions.
#' @param format File format.
#' @return A \code{ggsurvplot} object.
#' @export
PlotKaplanMeier <- function(object,
                            group_col    = NULL,
                            time_unit    = "months",
                            conf_int     = FALSE,
                            pairwise_p   = FALSE,
                            palette_name = "Darjeeling1",
                            base_size    = 13,
                            save_plot    = FALSE,
                            save_dir = NULL,
                            width        = 7,
                            height       = 8,
                            format       = "pdf") {
  
  if (is.null(save_dir) && save_plot) save_dir <- .get_viz_output_dir("PrognosiX")
  
  # ---- Extract data ----
  if (inherits(object, "PrognosiX")) {
    surv_df <- as.data.frame(object@survival.data)
    info_df <- as.data.frame(object@info.data)
    # If group_col is provided, merge it from info.data
    if (!is.null(group_col) && group_col %in% colnames(info_df)) {
      common <- intersect(rownames(surv_df), rownames(info_df))
      if (length(common) == 0) stop("No matching rows between survival.data and info.data.")
      surv_df <- surv_df[common, , drop = FALSE]
      info_df <- info_df[common, , drop = FALSE]
      surv_df[[group_col]] <- info_df[[group_col]]
    }
  } else {
    surv_df <- as.data.frame(object)
    if (!all(c("time", "status") %in% colnames(surv_df)))
      stop("Data frame must contain 'time' and 'status' columns.")
  }
  
  # ---- Build formula ----
  if (is.null(group_col)) {
    form_str <- "Surv(time, status) ~ 1"
  } else {
    if (!group_col %in% colnames(surv_df)) stop("Group column not found.")
    form_str <- paste0("Surv(time, status) ~ ", group_col)
  }
  fit <- survival::survfit(as.formula(form_str), data = surv_df)
  
  n_grp <- if (is.null(group_col)) 1L else length(unique(surv_df[[group_col]]))
  cols <- .get_palette(palette_name, n_grp)
  
  # ---- Plot ----
  km <- survminer::ggsurvplot(
    fit,
    data             = surv_df,
    palette          = cols,
    conf.int         = conf_int,
    pval             = !is.null(group_col),
    pval.method      = !is.null(group_col),
    risk.table       = TRUE,
    risk.table.col   = "strata",
    surv.median.line = "hv",
    xlab             = paste0("Time (", time_unit, ")"),
    ylab             = "Overall Survival Probability",
    legend.title     = if (!is.null(group_col)) group_col else "",
    ggtheme          = .pub_theme(base_size),
    tables.theme     = ggplot2::theme_minimal(base_size = base_size - 2) +
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold"))
  )
  
  if (pairwise_p && !is.null(group_col)) {
    pw <- survminer::pairwise_survdiff(as.formula(form_str), data = surv_df)
    km$plot <- km$plot + ggplot2::labs(caption = paste(
      "Pairwise log-rank p:\n",
      paste(capture.output(print(pw$p.value, digits = 3)), collapse = "\n")
    ))
  }
  
  if (save_plot) .save_plot(km$plot, save_dir, paste0("KM_", if (!is.null(group_col)) group_col else "overall"), width, height, format)
  return(km)
}


#' Forest plot -- univariate / multivariate Cox HR
#'
#' Renders a publication-style forest plot with HR, 95 \% CI, and p-value
#' columns alongside the graphical HR display.
#'
#' @param object       A \code{PrognosiX} object or a data frame with columns
#'   \code{Variable}, \code{HR}, \code{CI_lower}, \code{CI_upper},
#'   \code{P_value}, \code{HR_95CI}.
#' @param analysis     \code{"univariate"} (default) or \code{"multivariate"}.
#'   Selects the result stored in \code{univariate.analysis}.
#' @param hr_limit     x-axis range for HR. Default \code{c(0.1, 10)}.
#' @param log_scale    Use log10 x-axis. Default \code{TRUE}.
#' @param palette_name Palette (2 colours: sig / NS).
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotForestPlot <- function(object,

                           analysis     = "univariate",
                           hr_limit     = c(0.1, 10),
                           log_scale    = TRUE,
                           palette_name = "Royal1",
                           base_size    = 12,
                           save_plot    = FALSE,
                           save_dir = NULL,  # default: auto-detect from config
                           width        = 10,
                           height       = 7,
                           format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")

  if (inherits(object, "PrognosiX")) {
    hr_df <- object@univariate.analysis$hr_results
    if (is.null(hr_df) && analysis == "multivariate")
      hr_df <- object@best.model$hr_results
  } else {
    hr_df <- as.data.frame(object)
  }

  req_cols <- c("Variable", "HR", "CI_lower", "CI_upper", "P_value")
  if (!all(req_cols %in% colnames(hr_df)))
    stop("hr_df must contain: ", paste(req_cols, collapse = ", "))

  hr_df <- hr_df %>%
    dplyr::filter(is.finite(HR) & is.finite(CI_lower) & is.finite(CI_upper)) %>%
    dplyr::filter(HR >= hr_limit[1] & HR <= hr_limit[2]) %>%
    dplyr::mutate(
      Sig      = P_value < 0.05,
      P_label  = ifelse(P_value < 0.001, "<0.001",
                        as.character(round(P_value, 3))),
      HR_label = if ("HR_95CI" %in% colnames(.)) HR_95CI
                 else paste0(round(HR, 2), " (",
                             round(CI_lower, 2), "-",
                             round(CI_upper, 2), ")"),
      Variable = factor(Variable, levels = rev(Variable))
    )

  cols <- .get_palette(palette_name, 2)

  p <- ggplot2::ggplot(hr_df,
    ggplot2::aes(y = Variable, x = HR,
                 xmin = CI_lower, xmax = CI_upper,
                 colour = Sig)) +
    ggplot2::geom_vline(xintercept = 1,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::geom_errorbarh(height = 0.25, linewidth = 0.8) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "grey60", `TRUE` = cols[1]),
      labels = c(`FALSE` = "p >= 0.05", `TRUE` = "p < 0.05"),
      name   = NULL) +
    ggplot2::geom_text(ggplot2::aes(x = hr_limit[2] * 1.1, label = HR_label),
                       hjust = 0, size = 3, colour = "black") +
    ggplot2::geom_text(ggplot2::aes(x = hr_limit[2] * 1.5, label = P_label),
                       hjust = 0, size = 3, colour = "black") +
    {if (log_scale) ggplot2::scale_x_log10(limits = c(hr_limit[1],
                                                        hr_limit[2] * 2))
     else ggplot2::scale_x_continuous(limits = c(hr_limit[1],
                                                   hr_limit[2] * 2))} +
    ggplot2::labs(title = paste(tools::toTitleCase(analysis), "Cox Forest Plot"),
                  x = "Hazard Ratio (95% CI)", y = NULL) +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = "top")

  if (save_plot)
    .save_plot(p, save_dir, paste0("forest_", analysis), width, height, format)
  return(p)
}


#' Time-dependent ROC curves
#'
#' Draws AUC-annotated ROC curves for multiple time-points using
#' \pkg{timeROC}.
#'
#' @param object       A \code{PrognosiX} object or a survival data frame.
#' @param marker_col   Column to use as the marker / risk score.
#' @param time_points  Numeric vector of time-points. When \code{NULL} the
#'   25th, 50th, 75th percentiles of \code{time} are used.
#' @param time_unit    Unit label. Default \code{"months"}.
#' @param palette_name Palette.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotTimeROC <- function(object,

                        marker_col   = "risk_score",
                        time_points  = NULL,
                        time_unit    = "months",
                        palette_name = "Darjeeling1",
                        base_size    = 13,
                        save_plot    = FALSE,
                        save_dir = NULL,  # default: auto-detect from config
                        width        = 6,
                        height       = 5.5,
                        format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")


  surv_df <- if (inherits(object, "PrognosiX")) object@survival.data
              else as.data.frame(object)
  if (!marker_col %in% colnames(surv_df))
    stop("Column '", marker_col, "' not found.")

  time   <- surv_df[["time"]]
  status <- as.numeric(surv_df[["status"]])
  marker <- surv_df[[marker_col]]

  tps <- if (is.null(time_points)) {
    round(stats::quantile(time, c(0.25, 0.5, 0.75)), 1)
  } else {
    time_points
  }

  roc_res <- timeROC::timeROC(T = time, delta = status,
                               marker = marker, cause = 1,
                               weighting = "marginal", times = tps)

  auc_vals <- round(roc_res$AUC, 3)
  cols     <- .get_palette(palette_name, length(tps))

  roc_list <- lapply(seq_along(tps), function(i) {
    fp <- roc_res$FP[, i]
    tp <- roc_res$TP[, i]
    data.frame(FPR       = fp,
               TPR       = tp,
               TimePoint = paste0("t=", tps[i], " ", time_unit,
                                  " (AUC=", auc_vals[i], ")"))
  })
  roc_df <- dplyr::bind_rows(roc_list)

  p <- ggplot2::ggplot(roc_df,
    ggplot2::aes(FPR, TPR, colour = TimePoint)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0,
                         linetype = "dashed", colour = "grey50") +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::coord_equal() +
    ggplot2::labs(title    = "Time-Dependent ROC",
                  subtitle = paste("Marker:", marker_col),
                  x = "False Positive Rate",
                  y = "True Positive Rate") +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = c(0.68, 0.22))

  if (save_plot)
    .save_plot(p, save_dir, "time_ROC", width, height, format)
  return(p)
}


#' Restricted Cubic Spline (RCS) dose-response plot
#'
#' Fits an RCS model (linear, logistic, or Cox) and plots the non-linear
#' association with 95 \% CI ribbon.
#'
#' @param object       A \code{PrognosiX} object or data frame.
#' @param x_col        Continuous predictor column.
#' @param y_col        Outcome column (event status for Cox).
#' @param method       \code{"linear"}, \code{"logistic"}, or \code{"cox"}.
#' @param knots        Number of RCS knots. Default \code{4}.
#' @param adjust_vars  Additional adjustment covariates.
#' @param ref_point    Reference value for OR/HR = 1. Default: median.
#' @param palette_name Single colour palette.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotRCS <- function(object,

                    x_col        = NULL,
                    y_col        = "status",
                    method       = c("cox", "logistic", "linear"),
                    knots        = 4,
                    adjust_vars  = NULL,
                    ref_point    = NULL,
                    palette_name = "Royal1",
                    base_size    = 13,
                    save_plot    = FALSE,
                    save_dir = NULL,  # default: auto-detect from config
                    width        = 6,
                    height       = 5,
                    format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")

  method  <- match.arg(method)
  surv_df <- if (inherits(object, "PrognosiX")) object@survival.data
              else as.data.frame(object)
  if (is.null(x_col)) stop("Specify 'x_col'.")

  ref <- if (is.null(ref_point)) stats::median(surv_df[[x_col]], na.rm = TRUE)
          else ref_point

  dd <- rms::datadist(surv_df)
  options(datadist = "dd")

  fit <- if (method == "cox") {
    rms::cph(as.formula(paste0("Surv(time, status) ~ rms::rcs(", x_col, ",", knots, ")",
                               if (!is.null(adjust_vars)) paste0("+", paste(adjust_vars, collapse = "+")))),
              data = surv_df, x = TRUE, y = TRUE)
  } else if (method == "logistic") {
    rms::lrm(as.formula(paste0(y_col, " ~ rms::rcs(", x_col, ",", knots, ")",
                               if (!is.null(adjust_vars)) paste0("+", paste(adjust_vars, collapse = "+")))),
              data = surv_df, x = TRUE, y = TRUE)
  } else {
    rms::ols(as.formula(paste0(y_col, " ~ rms::rcs(", x_col, ",", knots, ")",
                               if (!is.null(adjust_vars)) paste0("+", paste(adjust_vars, collapse = "+")))),
              data = surv_df, x = TRUE, y = TRUE)
  }

  pred <- rms::Predict(fit, !!rlang::sym(x_col), ref.zero = TRUE, fun = exp)
  pred_df <- as.data.frame(pred)

  col <- .get_palette(palette_name, 1)
  y_label <- switch(method, cox = "Hazard Ratio", logistic = "Odds Ratio", "beta")

  hist_df <- data.frame(x = surv_df[[x_col]])

  p <- ggplot2::ggplot(pred_df,
    ggplot2::aes(x = .data[[x_col]], y = yhat)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = paste0(col, "30"), colour = NA) +
    ggplot2::geom_line(colour = col, linewidth = 1.2) +
    ggplot2::geom_rug(data = hist_df, ggplot2::aes(x = x, y = NULL),
                      sides = "b", alpha = 0.3, colour = "grey30") +
    ggplot2::labs(title    = paste("RCS --", x_col, "vs", y_label),
                  subtitle = paste0(knots, " knots, ref = ", round(ref, 2)),
                  x = x_col, y = y_label) +
    .pub_theme(base_size)

  if (save_plot)
    .save_plot(p, save_dir, paste0("RCS_", x_col), width, height, format)
  return(p)
}


#' Risk score distribution plot
#'
#' Draws a dot plot of per-sample risk scores sorted by value, coloured by
#' high / low risk group, with a KM curve and time-ROC in a composite layout.
#'
#' @param object       A \code{PrognosiX} object or a data frame with
#'   \code{risk_score}, \code{time}, \code{status}.
#' @param score_col    Column name for the risk score. Default
#'   \code{"risk_score"}.
#' @param cutoff_method Method to determine high / low boundary:
#'   \code{"median"} (default), \code{"mean"}, or a numeric value.
#' @param palette_name Palette (2 colours: high / low).
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{patchwork} composite.
#' @export
PlotRiskScore <- function(object,

                          score_col      = "risk_score",
                          cutoff_method  = "median",
                          palette_name   = "Royal1",
                          base_size      = 12,
                          save_plot      = FALSE,
                          save_dir = NULL,  # default: auto-detect from config
                          width          = 10,
                          height         = 4,
                          format         = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")


  surv_df <- if (inherits(object, "PrognosiX")) object@survival.data
              else as.data.frame(object)
  if (!score_col %in% colnames(surv_df))
    stop("Column '", score_col, "' not found.")

  cutoff <- if (is.numeric(cutoff_method)) cutoff_method
             else if (cutoff_method == "mean") mean(surv_df[[score_col]], na.rm = TRUE)
             else stats::median(surv_df[[score_col]], na.rm = TRUE)

  surv_df <- surv_df[order(surv_df[[score_col]]), ]
  surv_df$Risk  <- ifelse(surv_df[[score_col]] >= cutoff, "High", "Low")
  surv_df$Index <- seq_len(nrow(surv_df))

  cols <- .get_palette(palette_name, 2)
  names(cols) <- c("Low", "High")

  # panel A - dot plot
  p_score <- ggplot2::ggplot(surv_df,
    ggplot2::aes(Index, .data[[score_col]], colour = Risk)) +
    ggplot2::geom_point(size = 0.7, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = which(surv_df$Risk == "High")[1] - 0.5,
                        linetype = "dashed", colour = "grey40") +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::labs(title = "Risk Score Distribution",
                  x = "Patients (ranked)", y = "Risk Score",
                  colour = "Risk Group") +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = "right")

  # panel B - survival status tiles
  p_status <- ggplot2::ggplot(surv_df,
    ggplot2::aes(Index, 1,
                 colour = as.factor(as.numeric(surv_df[["status"]])))) +
    ggplot2::geom_point(shape = 3, size = 0.6) +
    ggplot2::scale_colour_manual(values = c("0" = "grey70", "1" = "black"),
                                  labels = c("0" = "Censored", "1" = "Event"),
                                  name   = NULL) +
    ggplot2::labs(x = NULL, y = "Event") +
    ggplot2::scale_y_continuous(breaks = NULL) +
    .pub_theme(base_size - 2) +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  out_p <- patchwork::wrap_plots(p_score, p_status, ncol = 1,
                                  heights = c(4, 1))

  if (save_plot)
    .save_plot(out_p, save_dir, "risk_score", width, height, format)
  return(out_p)
}


#' Decision Curve Analysis (DCA) plot
#'
#' Plots net benefit vs. threshold probability for one or more models
#' alongside the "treat all" and "treat none" reference lines.
#'
#' @param object       A \code{PrognosiX} object, or a data frame with outcome
#'   and predictor columns.
#' @param predictors   Character vector of predictor / risk score column names.
#' @param outcome_col  Binary outcome column. Default \code{"status"}.
#' @param thresholds   Sequence of decision thresholds. Default
#'   \code{seq(0, 0.5, by = 0.01)}.
#' @param palette_name Palette.
#' @param base_size    Base font size.
#' @param save_plot    Logical.
#' @param save_dir     Output directory.
#' @param width,height Inches.
#' @param format       File format.
#' @returns A \code{ggplot} object.
#' @export
PlotDCA <- function(object,

                    predictors   = "risk_score",
                    outcome_col  = "status",
                    thresholds   = seq(0, 0.5, by = 0.01),
                    palette_name = "Darjeeling1",
                    base_size    = 13,
                    save_plot    = FALSE,
                    save_dir = NULL,  # default: auto-detect from config
                    width        = 7,
                    height       = 5.5,
                    format       = "pdf") {
  if (is.null(save_dir)) save_dir <- .get_viz_output_dir("PrognosiX")

  df <- if (inherits(object, "PrognosiX")) object@survival.data
         else as.data.frame(object)
  df[[outcome_col]] <- as.numeric(df[[outcome_col]])
  n <- nrow(df)
  prev <- mean(df[[outcome_col]], na.rm = TRUE)

  dca_list <- lapply(thresholds, function(pt) {
    treat_all_nb <- prev - (1 - prev) * pt / (1 - pt)
    rows <- lapply(predictors, function(pr) {
      if (!pr %in% colnames(df)) return(NULL)
      test_pos   <- df[[pr]] >= pt
      true_pos   <- sum(test_pos & df[[outcome_col]] == 1, na.rm = TRUE)
      false_pos  <- sum(test_pos & df[[outcome_col]] == 0, na.rm = TRUE)
      net_benefit <- (true_pos / n) - (false_pos / n) * (pt / (1 - pt))
      data.frame(Threshold = pt, Model = pr, NetBenefit = net_benefit)
    })
    ref_row <- data.frame(Threshold  = pt,
                          Model      = "Treat All",
                          NetBenefit = treat_all_nb)
    dplyr::bind_rows(c(rows, list(ref_row)))
  })

  dca_df <- dplyr::bind_rows(dca_list)
  dca_df <- dplyr::bind_rows(
    dca_df,
    data.frame(Threshold  = thresholds,
               Model      = "Treat None",
               NetBenefit = 0)
  )

  model_levels <- c(predictors, "Treat All", "Treat None")
  dca_df$Model <- factor(dca_df$Model, levels = model_levels)

  n_cols <- length(model_levels)
  cols   <- .get_palette(palette_name, max(3, length(predictors)))
  cols   <- c(cols[seq_along(predictors)], "grey40", "grey70")
  names(cols) <- model_levels

  ltys <- c(rep("solid", length(predictors)), "dashed", "dotted")
  names(ltys) <- model_levels

  p <- ggplot2::ggplot(dca_df,
    ggplot2::aes(Threshold, NetBenefit, colour = Model, linetype = Model)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(values = cols)  +
    ggplot2::scale_linetype_manual(values = ltys) +
    ggplot2::coord_cartesian(ylim = c(-0.05, max(dca_df$NetBenefit, na.rm = TRUE) + 0.05)) +
    ggplot2::labs(title    = "Decision Curve Analysis",
                  x = "Threshold Probability",
                  y = "Net Benefit",
                  colour   = NULL,
                  linetype = NULL) +
    .pub_theme(base_size) +
    ggplot2::theme(legend.position = c(0.75, 0.75))

  if (save_plot)
    .save_plot(p, save_dir, "DCA", width, height, format)
  return(p)
}

#' Plot Group Mean Heatmap
#'
#' Creates a heatmap of group-wise mean expression values for selected features.
#' The function filters features based on differential expression results,
#' computes mean expression per group, optionally Z-scores by row, and visualizes
#' using ComplexHeatmap with auto-annotation.
#'
#' @param object An S4 object (e.g., \code{Subtyping} or \code{Stat}) or a list
#'   containing \code{clean.data} (expression matrix) and \code{info.data}
#'   (metadata). The \code{clean.data} slot/data should have samples as rows and
#'   features as columns.
#' @param deg_df A data frame containing differential expression analysis results.
#'   Must include columns: \code{id} (feature names), \code{p.adjust} (adjusted
#'   p-values), \code{logFC} (log fold change), and \code{target_group}
#'   (group labels).
#' @param group_by Character string specifying the column name in \code{info.data}
#'   used for grouping samples. Default is \code{"cluster_lpa"}.
#' @param top_n Integer. Number of top features to select per group based on
#'   \code{logFC}. Default is \code{5}.
#' @param p_cutoff Numeric. Adjusted p-value threshold for feature selection.
#'   Default is \code{0.05}.
#' @param logfc_cutoff Numeric. Log fold-change threshold for feature selection.
#'   Only features with absolute logFC > \code{logfc_cutoff} are retained.
#'   Default is \code{0}.
#' @param custom_levels Character vector. Optional custom ordering of group
#'   levels. If \code{NULL}, groups are ordered alphabetically using
#'   \code{gtools::mixedsort}. Default is \code{NULL}.
#' @param z_score_type Character string specifying whether to Z-score rows.
#'   Currently only \code{"row"} is supported. Default is \code{"row"}.
#' @param heatmap_palette Character vector of length 3 specifying colors for
#'   the low, middle, and high values in the heatmap. Default is
#'   \code{c("#2166AC", "white", "#B2182B")}.
#' @param save_path Character string specifying the file path to save the PDF.
#'   If \code{NULL}, the plot is not saved. Default is
#'   \code{"./Group_Mean_Heatmap.pdf"}.
#'
#' @return A \code{ComplexHeatmap::Heatmap} object (invisibly). The heatmap
#'   is also drawn on the current graphics device and saved to PDF if
#'   \code{save_path} is provided.
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom gtools mixedsort
#' @importFrom scales hue_pal
#' @importFrom grid gpar
#' @importFrom methods slot
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'sub_obj' is a Subtyping object and 'deg_results' is a data frame
#' PlotGroupMeanHeatmap(
#'   object = sub_obj,
#'   deg_df = deg_results,
#'   group_by = "cluster_lpa",
#'   top_n = 10,
#'   p_cutoff = 0.01
#' )
#' }
PlotGroupMeanHeatmap <- function(
    object,
    deg_df,
    group_by           = "cluster_lpa",
    top_n              = 5,
    p_cutoff           = 0.05,
    logfc_cutoff       = 0,
    custom_levels      = NULL,
    z_score_type       = "row",
    heatmap_palette    = c("#2166AC", "white", "#B2182B"),
    save_path          = "./Group_Mean_Heatmap.pdf"
) {
  
  # [1] Data Extraction
  if (isS4(object)) {
    mat_df  <- methods::slot(object, "clean.data")
    meta_df <- methods::slot(object, "info.data")
  } else {
    mat_df  <- object$clean.data
    meta_df <- object$info.data
  }
  
  raw_group <- as.character(meta_df[[group_by]])
  lvls <- if (!is.null(custom_levels)) custom_levels else gtools::mixedsort(unique(raw_group))
  
  # [2] Feature Selection
  filtered_deg <- deg_df %>%
    dplyr::filter(id %in% colnames(mat_df)) %>% 
    dplyr::filter(p.adjust < p_cutoff, logFC > logfc_cutoff) %>%
    dplyr::group_by(target_group) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(match(target_group, lvls))
  
  selected_genes <- unique(filtered_deg$id)
  calc_mat <- as.matrix(mat_df[, selected_genes, drop = FALSE])
  
  # [3] Group Mean Calculation
  plot_mat <- sapply(lvls, function(g) {
    colMeans(calc_mat[which(raw_group == g), , drop = FALSE], na.rm = TRUE)
  })
  
  if (z_score_type == "row") {
    plot_mat <- t(scale(t(plot_mat)))
  }
  
  # [4] Auto-Annotations
  group_colors <- setNames(hue_pal()(length(lvls)), lvls)
  
  ha = HeatmapAnnotation(
    Group = lvls,
    col = list(Group = group_colors),
    show_legend = TRUE,
    annotation_name_side = "left"
  )
  
  # [5] Heatmap Construction
  col_fun <- colorRamp2(c(min(plot_mat), 0, max(plot_mat)), heatmap_palette)
  
  ht <- Heatmap(
    plot_mat,
    name = "Mean Exp",
    col = col_fun,
    top_annotation = ha,           
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = factor(filtered_deg$target_group[match(rownames(plot_mat), filtered_deg$id)], levels = lvls),
    column_names_side = "top",     
    row_names_side = "left",
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 1)
  )
  
  # [6] Save
  if (!is.null(save_path)) {
    pdf(save_path, width = 7, height = 8)
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
    dev.off()
  }
  return(ht)
}

