#file           :clustering_maps.R
#description    :Perform clustering on raster data, visualize dendrogram, and generate maps
#author         :Giovanni Bonafè | ARPA FVG
#created        :2025-03-24
#last update    :2025-04-15
#version        :1.1
#dependencies   :terra, ncdf4, dendextend, RColorBrewer, classInt, tidyr, readr
#notes          :Uses Jaccard dissimilarity and Ward clustering
#                module load r-rgdal geos netcdf-fortran/4.5.2/gcc
#============================================================================

# Define command-line options
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-s", "--specie"), type = "character", default = "O3",
              help = "Pollutant to analyze [default: %default]"),
  make_option(c("-p", "--pattern"), type = "character", 
              default = "data/processed/{specie}_*_unbiased_scenario.tif",
              help = "File pattern for input rasters [default: '%default']"),
  make_option(c("-T", "--topo_file"), type = "character", 
              default = "data/fairmode-wg5-exercise-202501/topography.nc",
              help = "Topography NetCDF file [default: '%default']"),
  make_option(c("-t", "--thresholds"), type = "character",
              help = "Optional thresholds for Jaccard similarity (comma-separated)"),
  make_option(c("-c", "--clusters"), type = "integer", default = 10,
              help = "Maximum number of clusters [default: %default]")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))
specie <- opt$specie
thresholds <- if(!is.null(opt$thresholds)){
  as.numeric(strsplit(opt$thresholds, ",")[[1]])
} else {
  NULL
}
file_pattern <- opt$pattern
nc_topo_file <- opt$topo_file
max_nclu <- opt$clusters

# Script for Jaccard distance
source("R/distance_scenarios.R")

# Load packages
library(terra)
library(ncdf4)
library(dendextend)
library(RColorBrewer)
library(glue)
library(classInt)
library(tidyr)
library(readr)
library(dplyr)
library(futile.logger)

# Function to load topography from NetCDF
topography_from_nc <- function(nc_file, varname = "topography", grid) {
  nc <- nc_open(nc_file)
  topo_vals <- ncvar_get(nc, varname)
  nc_close(nc)
  topo_rast <- rast(t(topo_vals))
  ext(topo_rast) <- ext(grid)
  return(topo_rast)
}

# Function to load raster stack from files
load_rasters <- function(files) {
  rasters <- lapply(files, rast)
  names(rasters) <- sapply(basename(files), function(x) strsplit(x, "_")[[1]][2])
  return(rast(rasters))
}

# Function to compute global quantile thresholds
compute_thresholds <- function(r_stack, probs) {
  thrs <- signif(quantile(values(r_stack, mat=FALSE), probs = probs, na.rm = TRUE), 4)
  flog.info(glue("Thresholds: {paste(thrs,collapse=', ')} ({paste(names(thrs), collapse=', ')})"))
  return(thrs)
}

# Function to find the optimal number of clusters
find_optimal_clusters <- function(dend, max_k = 10) {
  heights <- rev(sort(get_branches_heights(dend)))
  gaps <- diff(heights)
  max_k <- min(max_k, length(gaps))
  optimal_k <- which.max(gaps[1:max_k]) + 1
  return(optimal_k)
}

# Function to perform clustering analysis on species distribution data
# 
# Args:
#   specie: Character, name of the species to analyze
#   file_pattern: Character, pattern to find species distribution files (glue syntax)
#   nc_topo_file: Character, path to NetCDF file with topography data
#   thresholds: Numeric vector, optional thresholds for Jaccard similarity calculation
#
# Returns:
#   A list containing clustering results and intermediate objects:
#     - dend: dendrogram object from hierarchical clustering
#     - r_stack: original raster stack of species distribution
#     - r_stack_land: raster stack masked to land areas only
#     - clu_idx: cluster assignments ordered by dendrogram
#     - nclu: optimal number of clusters found
#     - thresholds: thresholds used for Jaccard similarity

perform_clustering <- function(
    specie, 
    file_pattern = "data/processed/{specie}_*_unbiased_scenario.tif", 
    nc_topo_file = "data/fairmode-wg5-exercise-202501/topography.nc",
    thresholds = NULL) {
  
  flog.info(glue("Processing specie {specie}..."))

    # Load and prepare raster data ----
  # Find all files matching the pattern for the given species
  files <- Sys.glob(glue(file_pattern))
  
  # Load raster files into a stack
  r_stack <- load_rasters(files)
  
  # Load topography data and prepare land mask ----
  # Extract topography data matching the raster extent
  topo <- topography_from_nc(nc_topo_file, grid = ext(r_stack))
  crs(topo) <- crs(r_stack)  # Ensure consistent coordinate reference system
  
  # Create land mask (values > 1 indicate land)
  land_mask <- topo > 1
  
  # Apply land mask to the raster stack (remove ocean areas)
  r_stack_land <- mask(x = r_stack, mask = land_mask, maskvalues = 0)
  
  # Calculate thresholds if not provided ----
  # Compute quantile thresholds for Jaccard similarity calculation
  if(is.null(thresholds)) {
    thresholds <- compute_thresholds(r_stack_land, probs = c(0.5, 0.75, 0.95))
  }
  nthr <- length(thresholds)
  
  # Calculate similarity matrix ----
  # Compute average Jaccard similarity across all thresholds
  flog.info("Calculating dissimilarity matrix")
  d <- Reduce("+", lapply(1:nthr, function(i) {
    jaccard_matrix(r_stack_land, thresholds[i])
  })) / nthr
  
  # Hierarchical clustering ----
  flog.info("Performing hierachical clustering")
  clustering <- hclust(d, method = "ward.D2")  # Ward's method for compact clusters
  dend <- as.dendrogram(clustering)  # Convert to dendrogram object
  
  # Determine optimal number of clusters ----
  nclu <- find_optimal_clusters(dend, max_k = max_nclu)
  
  # Cut dendrogram to get cluster assignments ----
  clu_idx <- cutree(dend, k = nclu)
  
  # Process cluster assignments ----
  # Create output table with ordered cluster indices
  out <- tibble(name = names(clu_idx), clu_idx = unname(clu_idx)) %>%
    # Join with dendrogram order information
    left_join(tibble(name = labels(dend), dend_order = length(labels(dend)):1)) %>%
    # Sort by dendrogram order (from top to bottom)
    arrange(dend_order) %>%
    # Convert cluster IDs to factors ordered by first appearance in dendrogram
    mutate(clu_idx_ordered = factor(clu_idx, levels = unique(clu_idx))) %>%
    # Convert to numeric IDs (1 = first cluster in dendrogram, etc.)
    mutate(dend_idx = as.integer(clu_idx_ordered)) %>%
    # Remove temporary column
    select(-clu_idx_ordered)
  
  # Return results ----
  return(list(
    dend = dend,               # Dendrogram object
    r_stack = r_stack[[out$name]],# Original raster stack, reordered according to dendogram
    r_stack_land = r_stack_land[[out$name]], # Land-only raster stack
    clu_idx = out$dend_idx,    # Cluster IDs ordered by dendrogram
    raster_name = out$name,    # Name of the raster
    nclu = nclu,               # Optimal number of clusters
    thresholds = thresholds    # Thresholds used for similarity
  ))
}

# Function to find medoid
calculate_jaccard_medoid <- function(r_stack, thresholds, cluster_members) {
  if (length(cluster_members) == 1) {
    return(cluster_members[1])
  }
  
  cluster_stack <- r_stack[[cluster_members]]
  
  # Calculate Jaccard dissimilarity matrix
  d_list <- lapply(thresholds, function(thr) {
    jaccard_matrix(cluster_stack, thr)
  })
  
  # Ensure we're working with matrices and average them
  d_matrices <- lapply(d_list, function(x) {
    if (inherits(x, "dist")) as.matrix(x) else x
  })
  
  d <- Reduce("+", d_matrices) / length(thresholds)
  
  # Validate matrix structure
  if (!is.matrix(d) || nrow(d) != length(cluster_members)) {
    stop("Invalid dissimilarity matrix structure")
  }
  
  # Find medoid
  avg_distances <- rowMeans(d)
  medoid_idx <- which.min(avg_distances)
  return(cluster_members[medoid_idx])
}

# Function to plot and save dendrogram
plot_dendrogram <- function(dend, nclu, specie) {
  colors <- brewer.pal(min(nclu, 8), "Dark2")
  dend <- dend %>% 
    color_labels(k = nclu, col = colors) %>% 
    color_branches(k = nclu, col = colors, groupLabels= (nclu:1)) %>%
    set("labels_cex", 0.6) %>% 
    set("branches_lwd", 2)
  flog.info(glue("Plotting dendrogram in file dendrogram_{specie}.pdf"))
  pdf(glue("dendrogram_{specie}.pdf"), width = 8, height = 6)
  par(mar = c(2, 2, 2, 8))
  plot(dend, horiz = TRUE, axes = FALSE, 
       xlim = c(max(get_branches_heights(dend)) * 1.2, min(get_branches_heights(dend))))
  dev.off()
}

# Function to compute cluster statistics
compute_cluster_stats <- function(r_stack, clu_idx, thresholds) {
  unique_clusters <- unique(clu_idx)
  mean_rasters <- list()
  sd_rasters <- list()
  medoids <- list()  # To store most representative members
  
  for (k in unique_clusters) {
    flog.info(glue("Processing cluster {k}..."))
    layers_in_cluster <- which(clu_idx == k)
    member_names <- names(r_stack)[layers_in_cluster]
    
    # Identify medoid using Jaccard distance
    medoids[[k]] <- calculate_jaccard_medoid(r_stack, thresholds, member_names)
    
    # Calculate mean and standard deviation rasters
    if (length(layers_in_cluster) == 1) {
      mean_r <- r_stack[[layers_in_cluster]]
      sd_r <- mean_r * NA  # SD is undefined for single-member clusters
    } else {
      mean_r <- app(r_stack[[layers_in_cluster]], mean, na.rm = TRUE)
      sd_r <- app(r_stack[[layers_in_cluster]], sd, na.rm = TRUE)
    }
    
    mean_rasters[[paste0("Cluster_", k)]] <- mean_r
    sd_rasters[[paste0("Cluster_", k)]] <- sd_r
  }
  
  return(list(
    mean = mean_rasters, 
    sd = sd_rasters,
    medoids = medoids  # Most representative members
  ))
}

# Function to generate cluster maps
generate_cluster_maps <- function(cluster_stats, specie, n_clas = 7) {
  cluster_mean <- rast(cluster_stats$mean)
  cluster_sd <- rast(cluster_stats$sd)
  pp <- c(0, 0.01, (1:(n_clas-3))/(n_clas-2), 0.99, 1)
  mean_bins <- unique(round(signif(quantile(values(cluster_mean), pp), 2), 1))
  sd_bins <- unique(round(signif(quantile(values(cluster_sd), pp, na.rm = TRUE), 2), 2))
  mean_pal <- rev(brewer.pal(n_clas, "Spectral"))  
  sd_pal <- brewer.pal(n_clas, "PuOr")
  flog.info(glue("Plotting maps in file cluster_maps_{specie}.pdf:"))
  pdf(glue("cluster_maps_{specie}.pdf"), width = 8, height = 6)
  for (i in 1:length(cluster_stats$mean)) {
    flog.info(glue("cluster {i}"))
    k <- names(cluster_stats$mean)[i]
    par(mfrow = c(1, 2), mar = c(2, 2, 3, 4))
    plot(cluster_mean[[i]], main = paste("Mean -", k), col = mean_pal, breaks = mean_bins, legend = "topleft")
    plot(cluster_sd[[i]], main = paste("Std. Dev. -", k), col = sd_pal, breaks = sd_bins, legend = "topleft")
  }
  dev.off()
}

write_clustering <- function(clustering_results, cluster_stats, specie) {
  # Create medoid information data frame
  medoid_df <- data.frame(
    cluster_index = 1:length(cluster_stats$medoids),
    medoid_member = unlist(cluster_stats$medoids),
    stringsAsFactors = FALSE
  )
  
  # Prepare main output table
  out <- data.frame(
    member = clustering_results$raster_name,
    cluster_index = clustering_results$clu_idx, 
    method = "Jaccard distance, ward.D2",
    thresholds = paste(clustering_results$thresholds, collapse="|"),
    row.names = NULL) %>%
    separate(col = member, sep = "\\.", 
             into = c("sequence","calibration","correction","spatialization"), 
             fill = "right", remove = FALSE) %>%
    left_join(medoid_df, by = "cluster_index") %>%
    mutate(is_medoid = ifelse(member == medoid_member, "yes", "no"))
  
  # Write output file
  flog.info(glue("Saving clustering results to clustering_{specie}.csv"))
  write_csv(out, file = glue("clustering_{specie}.csv"))
  
  # Log medoid information
  flog.info("Most representative members (medoids) for each cluster:")
  for (k in 1:length(cluster_stats$medoids)) {
    flog.info(glue("Cluster {k}: {cluster_stats$medoids[[k]]}"))
  }
}

# execution
clustering_results <- perform_clustering(
  specie = specie, 
  thresholds = thresholds,
  file_pattern = file_pattern,
  nc_topo_file = nc_topo_file
)
plot_dendrogram(
  clustering_results$dend, 
  specie = specie, 
  nclu = clustering_results$nclu)
cluster_stats <- compute_cluster_stats(
  clustering_results$r_stack, 
  clustering_results$clu_idx,
  clustering_results$thresholds)
write_clustering(
  clustering_results, 
  cluster_stats, 
  specie)
generate_cluster_maps(
  cluster_stats, 
  specie = specie)
