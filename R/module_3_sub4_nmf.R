#' Generate NMF Plot
#'
#' @param estimate NMF estimate.
#' @param save_dir Directory.
#' @importFrom NMF plot
#' @export
generate_nmf_plot <- function(estimate,
                              save_dir = here('Subtyping',"cluster_results","nmf_results" )) {
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  nmf_plot_path <- file.path(save_dir, "nmf_plot.pdf")
  cat("Saving NMF plot to: ", nmf_plot_path, "\n")
  
  p <- plot(estimate)
  
  ggsave(nmf_plot_path, plot = p, width = 8, height = 6)
  
  cat("NMF plot saved successfully using ggsave.\n")
}

#' Generate Consensus Map Plot
#'
#' @param estimate NMF estimate.
#' @param save_dir Directory.
#' @importFrom NMF consensusmap
#' @importFrom grDevices pdf dev.off
#' @export
generate_consensus_map_plot <- function(estimate,
                                        save_dir = here('Subtyping',"cluster_results","nmf_results" )) {
  
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }
  
  consensus_plot_path <- file.path(save_dir, "consensus_map.pdf")
  cat("Saving Consensus Map plot to: ", consensus_plot_path, "\n")
  
  pdf(file = consensus_plot_path)
  consensusmap(estimate,
               annRow = NA,
               annCol = NA,
               main = "Consensus matrix",
               info = FALSE)
  cat("Consensus Map plot saved successfully.\n")
  dev.off()
}

#' NMF Estimate
#'
#' @param data Data matrix.
#' @param rank_range Rank range.
#' @param seed Seed.
#' @param nrun nrun.
#' @param method Method.
#' @param save_dir Directory.
#' @importFrom NMF nmf
#' @export
nmf_estimate <- function(data,
                         rank_range = 2:4,
                         seed = 8891,
                         nrun = 10,
                         method = "brunet",
                         save_dir = here('Subtyping',"cluster_results","nmf_results" )) {
  
  if (any(data < 0)) {
    warning("The data contains negative values, which have been replaced with 0.")
    data[data < 0] <- 0
  }
  
  cat("Starting NMF decomposition...\n")
  estimate <- nmf(data,
                  rank = rank_range,
                  method = method,
                  nrun = nrun,
                  seed = seed)
  cat("NMF decomposition completed.\n")
  
  generate_nmf_plot(estimate=estimate, save_dir=save_dir)
  
  generate_consensus_map_plot(estimate=estimate, save_dir=save_dir)
  
  return(estimate)
}

#' Sub NMF Estimate
#'
#' @param object Subtyping object.
#' @param rank_range Rank range.
#' @param seed Seed.
#' @param nrun nrun.
#' @param method Method.
#' @param save_dir Directory.
#' @export
Sub_nmf_estimate <- function(object,
                             rank_range = 2:4,
                             seed = 8891,
                             nrun = 10,
                             method = "brunet",
                             save_dir = here('Subtyping',"cluster_results","nmf_results" )) {
  
  if (inherits(object, "Subtyping")) {
    data <- slot(object, "clean.data")
  } else if (is.data.frame(object)) {
    data <- object
  } else {
    stop("Input must be an object of class 'Subtyping' or a data frame")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  estimate <- nmf_estimate(data=data,
                           rank_range=rank_range,
                           seed=seed,
                           nrun=nrun,
                           method=method,
                           save_dir=save_dir)
  
  if (inherits(object, "Subtyping")) {
    object@cluster.results[["nmf.result"]] <- list(estimate =estimate)
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'cluster.results' slot updated.\n")
    
    return(object)
  } else {
    return(estimate)
  }
}

#' NMF Best Rank Analysis
#'
#' @param data Data.
#' @param estimate Estimate.
#' @param rank_range Rank range.
#' @param nrun nrun.
#' @param seed Seed.
#' @param method Method.
#' @param palette_name2 Palette.
#' @param palette_name1 Palette.
#' @param save_dir Directory.
#' @importFrom NMF basismap
#' @export
nmf_best_rank_analysis <- function(data,
                                   estimate,
                                   rank_range = 2:4,
                                   nrun = 10,
                                   seed = 8891,
                                   method = "brunet",
                                   palette_name2="Zissou1",
                                   palette_name1 = "Set3",
                                   save_dir = here('Subtyping',"cluster_results","nmf_results" )) {
  
  if (any(data < 0)) {
    warning("The data contains negative values, which have been replaced with 0.")
    data[data < 0] <- 0
  }
  
  cophenetic_scores <- estimate[["measures"]][["cophenetic"]]
  best_rank_index <- which.max(cophenetic_scores)
  best_rank <- rank_range[best_rank_index]
  cat("Best rank is: ", best_rank, "\n")
  
  estimate_best <- nmf(data, rank = best_rank, method = method, nrun = nrun, seed = seed)
  
  group_colors <- RColorBrewer::brewer.pal(best_rank+1, palette_name1)
  annColors_group <- setNames(group_colors, as.character(1:best_rank))
  
  save_path <- file.path(save_dir, "basismap_plot.pdf")
  
  cat("Saving the basis map as PDF...\n")
  
  pdf(file = save_path)
  p3 <- basismap(estimate_best,
                 cexCol = 1.5,
                 cexRow = 1.5,
                 annColors = list(group = annColors_group),
                 color = wes_palette(palette_name2, 100, type = "continuous")
  )
  dev.off()
  
  cat("Basis map has been saved as PDF: ", save_path, "\n")
  
  return(list(
    best_estimate = estimate_best,
    best_rank = best_rank
  ))
}

#' Sub Best Rank Analysis
#'
#' @param object Subtyping object.
#' @param nrun nrun.
#' @param seed Seed.
#' @param method Method.
#' @param palette_name2 Palette.
#' @param palette_name1 Palette.
#' @param save_dir Directory.
#' @export
Sub_best_rank_analysis <- function(object,
                                   nrun = 10,
                                   seed = 8891,
                                   method = "brunet",
                                   palette_name2 = "Zissou1",
                                   palette_name1 = "Set3",
                                   save_dir = here('Subtyping',"cluster_results","nmf_results" )) {
  
  if (inherits(object, "Subtyping")) {
    data <- slot(object, "clean.data")
    estimate <- slot(object, "cluster.results")[["nmf.result"]][["estimate"]]
  } else if (is.list(object)) {
    data <- object$data
    estimate <- object$estimate
  } else {
    stop("Input must be an object of class 'Subtyping' or a list")
  }
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found in the input")
  }
  
  re1 <- nmf_best_rank_analysis(estimate = estimate,
                                data = data,
                                nrun = nrun,
                                seed = seed,
                                method = method,
                                palette_name2 = palette_name2,
                                palette_name1 = palette_name1,
                                save_dir = save_dir)
  
  if (inherits(object, "Subtyping")) {
    object@cluster.results[["nmf.result"]][["best_estimate"]]<-re1$best_estimate
    object@cluster.results[["nmf.result"]][["best_rank"]]<-re1$best_rank
    object@Optimal.cluster<-re1$best_rank
    cat("Updating 'Subtyping' object...\n")
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'cluster.results' slot updated.\n")
    cat("- 'Optimal.cluster' slot updated.\n")
    
    return(object)
  } else {
    return(re1)
  }
}

#' Assign Groups NMF
#'
#' @param best_estimate Best estimate.
#' @param data Data.
#' @export
assign_groups_nmf <- function(best_estimate,
                              data) {
  
  cat("Extracting the W matrix from the NMF estimate...\n")
  W_matrix <- best_estimate@fit@W
  
  cat("Assigning samples to groups based on the highest weight factor...\n")
  group_assignment <- apply(W_matrix, 1, which.max)
  
  data$group <- factor(group_assignment)
  
  cat("Outputting group counts:\n")
  print(table(data$group))
  
  return(data)
}

#' Sub Assign Groups
#'
#' @param object Subtyping object.
#' @export
Sub_assign_groups <- function(object) {
  
  cat("Checking input type...\n")
  if (inherits(object, "Subtyping")) {
    cat("Detected 'Subtyping' object...\n")
    data <- slot(object, "clean.data")
    best_estimate <- slot(object, "cluster.results")[["nmf.result"]][["best_estimate"]]
  } else if (is.list(object)) {
    cat("Detected a list input...\n")
    data <- object$data
    best_estimate <- object$estimate
  } else {
    stop("Input must be an object of class 'Subtyping' or a list containing data and estimation results")
  }
  
  
  if (is.null(data) || nrow(data) == 0) {
    stop("No valid data found")
  }
  
  
  clustered_data <- assign_groups_nmf(best_estimate = best_estimate, data = data)
  
  if (inherits(object, "Subtyping")) {
    cat("Updating 'Subtyping' object with clustered data...\n")
    object@clustered.data <- clustered_data
    cat("The 'Subtyping' object has been updated with the following slots:\n")
    cat("- 'clustered.data' slot updated.\n")
    
    return(object)
  } else {
    cat("Returning clustered data as a data frame...\n")
    return(clustered_data)
  }
}
