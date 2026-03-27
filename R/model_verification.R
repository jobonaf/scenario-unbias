# =============================================================================
# title          :model_verification
# description    :Verification and selection of air quality models
#                 based on global metrics, percentiles, and error distribution analysis.
# author         :Giovanni Bonafe'
# date           :20250912
# version        :0.4
# notes          :Requires model raster files and reference raster data
# R_version      :3.5.2
# =============================================================================

# Define command-line options
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("glue"))
option_list <- list(
  make_option(c("-p", "--pollutant"), type = "character", default = "PM25",
              help = "Pollutant to analyze [default: %default]"),
  make_option(c("-r", "--reference"), type = "character", 
              default = "data/esercizio-dominio-italiano/Scenario_Reference_Gridded/SCEN_SC_A_c_{pollutant}_YEARLY.nc",
              help = "Reference data path [default: %default]"),
  make_option(c("-m", "--model_pattern"), type = "character", 
              default = "data/processed_italian/{pollutant}_*_unbiased_scenario.tif",
              help = "File pattern for model rasters [default: %default]"),
  make_option(c("-z", "--zones"), type = "character", 
              default = "data/models_verification/homogeneous_zones.tif",
              help = "Homogeneous zones raster path [default: %default]"),
  make_option(c("-b", "--boundary"), type = "character", 
              default = "/atlas/arpa/bonafeg/data/geo/LimitiAmministrativi/Italy_WGS84_LatLong/ITA_adm0.shp",
              help = "Italy boundary shapefile [default: %default]"),
  make_option(c("-R", "--regions"), type = "character", 
              default = "/atlas/arpa/bonafeg/data/geo/LimitiAmministrativi/Italy_WGS84_LatLong/ITA_adm1.shp",
              help = "Italy regions shapefile [default: %default]"),
  make_option(c("-o", "--output"), type = "character", 
              default = "data/models_verification/italian_exercise",
              help = "Output directory [default: %default]")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "modules required: r-rgdal geos netcdf-fortran/4.5.2/gcc"))

# Interpolate paths with pollutant name
pollutant <- opt$pollutant
reference_path <- glue(opt$reference)
model_pattern <- glue(opt$model_pattern)
zones_path <- opt$zones
italy_boundary_path <- opt$boundary
italy_regions_path <- opt$regions
output_dir <- opt$output

# Load required libraries
library(terra)    # For raster operations
library(sf)       
library(futile.logger) # For logging
library(dplyr)    # For data manipulation
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(forcats)

# =============================================================================
# Source external functions
# =============================================================================
source("R/read_netcdf_as_raster.R")

# =============================================================================
# Initialize logging
# =============================================================================
flog.appender(appender.console())  # Write to standard output
flog.threshold(INFO)
flog.info("Starting model verification workflow")
flog.info("Using parameters: %s", paste(commandArgs(trailingOnly = TRUE), collapse = " "))

# =============================================================================
# Phase 0: Data Preparation
# =============================================================================
flog.info("Phase 0: Data Preparation")

# Pollutant extended
poll_ext <- case_match(
  pollutant,
  "PM25" ~ "PM2.5",
  "NO2" ~ "nitrogen dioxide",
  "O3" ~ "ozone",
  .default = pollutant
)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load reference data using the custom read_nc function
flog.info(glue("Loading reference data for {pollutant}..."))
reference <- read_nc(reference_path) # Assuming read_nc returns a SpatRaster

# Load Italy boundary and regions from shapefiles
flog.info("Loading Italy boundary...")
italy_boundary <- st_read(italy_boundary_path, quiet = TRUE) %>% 
  vect()

flog.info("Loading Italian regions...")
italy_regions <- st_read(italy_regions_path, quiet = TRUE) %>% 
  vect()

# Create abbreviated region names
region_abbreviations <- abbreviate(italy_regions$NAME_1, minlength = 6)
flog.info(glue("Region abbreviations: {paste(region_abbreviations, collapse = ', ')}"))

# Load homogeneous zones mask
flog.info("Loading homogeneous zones mask...")
zones_mask <- rast(zones_path)

# Prima del loop principale, crea un raster delle regioni con la stessa risoluzione e estensione
flog.info("Creating regions raster...")
regions_raster <- rast(ext = ext(zones_mask), resolution = res(zones_mask), crs = crs(zones_mask))
italy_regions_raster <- rasterize(italy_regions, regions_raster, field = "ID_1")
italy_regions_raster <- crop(italy_regions_raster, italy_boundary)
italy_regions_raster <- mask(italy_regions_raster, italy_boundary)

# Define zone abbreviations
zone_abbreviations <- c(
  "1" = "HiMount",
  "2" = "HillAdr",  # Hills and Adriatic coast
  "3" = "Thyrren", # Tyrrhenian coast and Salento
  "4" = "PoUrban",  # Po Valley, Rome and Naples
  "5" = "LoMount"  # Prealps and Apennines
)

flog.info(glue("Zone abbreviations: {paste(zone_abbreviations, collapse = ', ')}"))

# Crop all spatial data to Italy boundary
flog.info("Cropping data to Italy boundary...")
reference <- crop(reference, italy_boundary)
reference <- mask(reference, italy_boundary)
zones_mask <- crop(zones_mask, italy_boundary)
zones_mask <- mask(zones_mask, italy_boundary)

# Get list of model files (TIF format)
model_files <- Sys.glob(model_pattern)
model_names <- gsub(glue("^{pollutant}_|_unbiased_scenario$"), "", 
                    tools::file_path_sans_ext(basename(model_files)))

flog.info(glue("Found {length(model_files)} model files for {pollutant}"))

# =============================================================================
# Helper function: Check for valid raster data
# =============================================================================
check_raster_validity <- function(model_path, boundary) {
  flog.info(glue("Checking file {basename(model_path)}"))
  tryCatch({
    model_rast <- rast(model_path)
    
    if (is.null(model_rast) || nlyr(model_rast) == 0) {
      return(FALSE)
    }
    
    model_rast <- crop(model_rast, boundary)
    model_rast <- mask(model_rast, boundary)
    
    valid_vals <- values(model_rast)
    valid_count <- sum(!is.na(valid_vals))
    
    return(valid_count > 0)
    
  }, error = function(e) {
    flog.error(glue("Error with file {basename(model_path)}: {e$message}"))
    return(FALSE)
  })
}
# Check all models for validity
valid_models <- sapply(model_files, check_raster_validity, boundary = italy_boundary)
invalid_models <- model_files[!valid_models]

if (length(invalid_models) > 0) {
  flog.warn(glue("Found {length(invalid_models)} invalid models (only NA/NaN values):"))
  for (model in invalid_models) {
    flog.warn(glue("  - {basename(model)}"))
  }
}

# Keep only valid models
model_files <- model_files[valid_models]
model_names <- model_names[valid_models]
flog.info(glue("Proceeding with {length(model_files)} valid models"))

# =============================================================================
# Helper function: Create selection plot
# =============================================================================
create_selection_plot <- function(data, phase_name, output_dir, pollutant, nrows = 1, 
                                  metrics_to_plot, facet_scales = "free_x", plot_median = TRUE) {
  flog.info(glue("Creating selection plot for {phase_name}"))
  
  # Select only the requested metrics
  valid_metrics <- metrics_to_plot[metrics_to_plot %in% names(data)]
  if (length(valid_metrics) == 0) {
    flog.err("No valid metrics found for plotting. Using default metrics.")
  }
  
  # Prepare data for ggplot
  plot_data <- data %>%
    select(model, all_of(valid_metrics)) %>%
    pivot_longer(cols = -model, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = valid_metrics))
  
  # Order models based on the first metric in the list
  primary_metric <- valid_metrics[1]
  model_order <- data %>%
    arrange(desc(.data[[primary_metric]])) %>%
    pull(model)
  
  plot_data$model <- factor(plot_data$model, levels = model_order)
  
  # Calculate PDF dimensions based on number of models
  n_models <- length(unique(plot_data$model))
  n_metrics <- length(valid_metrics)
  
  # Dynamic dimensions
  pdf_width <- 1 + (n_metrics/nrows * 1.7)  # Width based on number of metrics
  pdf_height <- 1 + (n_models*nrows * 0.2)   # Height based on number of models
  
  # Limit minimum and maximum dimensions
  pdf_width <- max(4, min(pdf_width, 12))
  pdf_height <- max(4, min(pdf_height, 12))
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = value, y = model)) +
    geom_point(size = 2, color = "steelblue")  +
    facet_wrap(~ metric, scales = facet_scales, nrow = nrows) +
    labs(title = glue("Pollutant: {poll_ext}"),
         x = "Metric Value", 
         y = "Models") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, hjust = 0),
          plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          strip.text = element_text(size = 9))
  if(plot_median) p <- p +
    geom_vline(data = plot_data %>% group_by(metric) %>% 
                 summarize(median_val = median(value, na.rm = TRUE)),
               aes(xintercept = median_val), color = "red", linetype = "dashed", linewidth = 0.5)
  
  # Save plot as PDF
  output_file <- glue("{output_dir}/selection_plot_{phase_name}_{pollutant}.pdf")
  ggsave(output_file, p, width = pdf_width, height = pdf_height, units = "in", dpi = 300)
  
  flog.info(glue("Selection plot saved: {output_file} (Dimensions: {pdf_width}x{pdf_height} inches)"))
  flog.info(glue("Metrics plotted: {paste(valid_metrics, collapse = ', ')}"))
  
  return(p)
}

# =============================================================================
# Phase 1: Global Screening (Standard Metrics)
# =============================================================================
flog.info("Phase 1: Global Screening")

calculate_global_metrics <- function(model_path, reference_rast, boundary) {
  # Load and prepare model data
  model_rast <- rast(model_path)
  model_rast <- crop(model_rast, boundary)
  model_rast <- mask(model_rast, boundary)
  
  # Get valid values
  model_vals <- values(model_rast)
  ref_vals <- values(reference_rast)
  valid_idx <- !is.na(model_vals) & !is.na(ref_vals)
  
  model_vals <- model_vals[valid_idx]
  ref_vals <- ref_vals[valid_idx]
  
  # Calculate metrics
  error <- model_vals - ref_vals
  n <- length(ref_vals)
  
  metrics <- list(
    ME = mean(error),
    MAE = mean(abs(error)),
    RMSE = sqrt(mean(error^2)),
    R2 = 1 - (sum((ref_vals - model_vals)^2) / sum((ref_vals - mean(ref_vals))^2)),
    IOA = 1 - (sum((ref_vals - model_vals)^2) / 
                 sum((abs(model_vals - mean(ref_vals)) + abs(ref_vals - mean(ref_vals)))^2)),
    correlation = cor(model_vals, ref_vals, method = "pearson"),
    n_pixels = n
  )
  
  return(metrics)
}

# Calculate metrics for all models
global_results <- list()
for (i in seq_along(model_files)) {
  flog.info(glue("Calculating global scores for model {i}/{length(model_files)}: {model_names[i]}"))
  global_results[[model_names[i]]] <- calculate_global_metrics(
    model_files[i], reference, italy_boundary
  )
}

# Convert to dataframe and rank by IOA
global_df <- do.call(rbind, lapply(global_results, as.data.frame))
global_df$model <- rownames(global_df)
rownames(global_df) <- NULL
global_df %>%
  mutate(model = ordered(model)) %>%
  mutate(model = fct_reorder(model, IOA)) -> global_df

# Select models better than median for key metrics INCLUDING R2
median_ioa <- median(global_df$IOA, na.rm = TRUE)
median_rmse <- median(global_df$RMSE, na.rm = TRUE)
median_cor <- median(global_df$correlation, na.rm = TRUE)
median_r2 <- median(global_df$R2, na.rm = TRUE)  # Added R2 median

selected_phase1 <- global_df %>%
  filter(IOA > median_ioa,
         RMSE < median_rmse,
         correlation > median_cor,
         R2 > median_r2)  %>%
  droplevels()

flog.info(glue("Selected {nrow(selected_phase1)} models from Phase 1 (better than median)"))

# Create selection plot
create_selection_plot(global_df, "Phase1_GlobalScreening", output_dir, pollutant,
                      metrics_to_plot = c("IOA", "MAE", "RMSE", "R2", "correlation"))

# =============================================================================
# Phase 2: Percentile Screening (Jaccard)
# =============================================================================
flog.info("Phase 2: Percentile Screening")

calculate_jaccard <- function(model_path, reference_rast, boundary, percentiles = c(0.5, 0.75, 0.95)) {
  # Load and prepare model data
  model_rast <- rast(model_path)
  model_rast <- crop(model_rast, boundary)
  model_rast <- mask(model_rast, boundary)
  
  # Calculate thresholds
  ref_vals <- values(reference_rast)
  ref_vals <- ref_vals[!is.na(ref_vals)]
  thresholds <- quantile(ref_vals, probs = percentiles)
  
  # Calculate Jaccard for each percentile
  jaccard_scores <- numeric(length(percentiles))
  names(jaccard_scores) <- paste0("P", percentiles * 100)
  
  for (j in seq_along(percentiles)) {
    # Create binary masks
    ref_binary <- reference_rast > thresholds[j]
    model_binary <- model_rast > thresholds[j]
    
    # Calculate intersection and union
    intersection <- sum(values(ref_binary) & values(model_binary), na.rm = TRUE)
    union <- sum(values(ref_binary) | values(model_binary), na.rm = TRUE)
    
    # Jaccard similarity
    jaccard_scores[j] <- intersection / union
  }
  
  return(jaccard_scores)
}

# Calculate Jaccard for selected models
jaccard_results <- list()
for (i in seq_along(selected_phase1$model)) {
  model_name <- as.character(selected_phase1$model[i])
  model_idx <- which(model_names == model_name)
  flog.info(glue("Calculating Jaccard for model {i}/{nrow(selected_phase1)}: {model_name}"))
  jaccard_results[[model_name]] <- calculate_jaccard(
    model_files[model_idx], reference, italy_boundary
  )
}

# Add Jaccard scores to results
jaccard_df <- do.call(rbind, jaccard_results)
selected_phase1 <- cbind(selected_phase1, jaccard_df)

# Select models better than median for ALL Jaccard percentiles
median_p50 <- median(selected_phase1$P50, na.rm = TRUE)
median_p75 <- median(selected_phase1$P75, na.rm = TRUE)
median_p95 <- median(selected_phase1$P95, na.rm = TRUE)

selected_phase2 <- selected_phase1 %>%
  filter(P50 > median_p50,
         P75 > median_p75,
         P95 > median_p95) %>%
  droplevels()

flog.info(glue("Selected {nrow(selected_phase2)} models from Phase 2 (better than median for all Jaccard percentiles)"))

# Create selection plot
create_selection_plot(selected_phase1, "Phase2_PercentileScreening", output_dir, pollutant,
                      metrics_to_plot = c("P95", "P75", "P50"))

# =============================================================================
# Phase 3: Error Distribution Analysis (Boxplots by Zone/Region)
# =============================================================================
flog.info("Phase 3: Error Distribution Analysis")

# Function to extract error values by zone/region
extract_errors_by_area <- function(model_path, reference_rast, area_mask, area_type = "zone", 
                                   area_names, model_name = NULL) {
  # Load and prepare model data
  model_rast <- rast(model_path)
  model_rast <- crop(model_rast, area_mask)
  model_rast <- mask(model_rast, area_mask)
  
  # Calculate error raster
  error_rast <- model_rast - reference_rast
  
  # Get unique area IDs
  area_ids <- unique(values(area_mask))
  area_ids <- area_ids[!is.na(area_ids)]
  
  # Initialize list for results
  error_list <- list()
  
  # Extract errors for each area
  for (area_id in area_ids) {
    # Create mask for current area
    area_mask_current <- area_mask == area_id
    
    # Extract error values for this area
    error_vals <- values(mask(error_rast, area_mask_current, maskvalues = FALSE))
    error_vals <- error_vals[!is.na(error_vals)]
    
    if (length(error_vals) > 0) {
      area_name <- area_names[as.character(area_id)]
      
      error_list[[area_name]] <- data.frame(
        model = if(!is.null(model_name)) model_name else basename(model_path),
        area = unname(area_name),
        area_type = area_type,
        error = error_vals,
        stringsAsFactors = F
      )
    }
  }
  
  # Combine all results
  error_df <- bind_rows(error_list)
  
  return(error_df)
}

# Extract errors for all selected models by homogeneous zones
zone_errors_list <- list()
region_errors_list <- list()

for (i in seq_along(selected_phase2$model)) {
  model_name <- selected_phase2$model[i]
  model_idx <- which(model_names == model_name)
  flog.info(glue("Extracting errors for model {i}/{nrow(selected_phase2)}: {model_name}"))
  
  # Extract errors by homogeneous zones
  zone_errors <- extract_errors_by_area(
    model_files[model_idx], reference, zones_mask, "zone", 
    zone_abbreviations, model_name
  )
  zone_errors_list[[model_name]] <- zone_errors
  
  # Extract errors by Italian regions
  reg_abbr <- region_abbreviations
  names(reg_abbr) <- 1:length(reg_abbr)
  region_errors <- extract_errors_by_area(
    model_files[model_idx], reference, italy_regions_raster, "region", 
    reg_abbr, model_name
  )
  region_errors_list[[model_name]] <- region_errors
}

# Combine all error data
zone_errors_all <- bind_rows(zone_errors_list)
region_errors_all <- bind_rows(region_errors_list)

# Create boxplots for zones
flog.info("Creating boxplots for homogeneous zones")
ll <- quantile(zone_errors_all$error, probs=c(0.01,0.99), na.rm=T)
zone_boxplot <- ggplot(zone_errors_all %>%
                         mutate(area = ordered(area)) %>%
                         mutate(area = fct_reorder(area, error)), 
                       aes(x = area, y = error, group = area)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~ model, ncol = 3) +
  scale_y_continuous(limits = ll) +
  labs(title = "Error Distribution by Homogeneous Zone",
       subtitle = glue("Pollutant: {poll_ext}"),
       x = "Zone", 
       y = "Error (Model - Reference)") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "none")

# Save zone boxplot
zone_output_file <- glue("{output_dir}/error_boxplot_zones_{pollutant}.pdf")
ggsave(zone_output_file, zone_boxplot, width = 10, height = 8, units = "in", dpi = 300)
flog.info(glue("Zone boxplot saved: {zone_output_file}"))

# Create boxplots for regions
flog.info("Creating boxplots for regions")
region_boxplot <- ggplot(region_errors_all %>%
                           mutate(area = ordered(area)) %>%
                           mutate(area = fct_reorder(area, error)), 
                         aes(x = area, y = error, fill = area, group = area)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~ model, ncol = 4) +
  scale_fill_manual(values = rep(brewer.pal(8, "Set2"), length.out = n_distinct(region_errors_all$area))) +
  scale_y_continuous(limits = ll) +
  labs(title = "Error Distribution by Region",
       subtitle = glue("Pollutant: {poll_ext}"),
       x = "Region", 
       y = "Error (Model - Reference)") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "none")

# Save region boxplot
region_output_file <- glue("{output_dir}/error_boxplot_regions_{pollutant}.pdf")
ggsave(region_output_file, region_boxplot, width = 12, height = 8, units = "in", dpi = 300)
flog.info(glue("Region boxplot saved: {region_output_file}"))

# Phase 3 is descriptive only - no model selection
selected_phase3 <- selected_phase2
flog.info(glue("Phase 3 completed (error distribution analysis). {nrow(selected_phase3)} models proceed to next phase"))

# =============================================================================
# Phase 4: Detailed Analysis of Finalists
# =============================================================================
flog.info("Phase 4: Detailed Analysis")

# Select top models for detailed analysis (max 10)
n_finalists <- min(10, nrow(selected_phase3))
finalists <- selected_phase3 %>%
  arrange(desc(IOA), desc(P95)) %>%
  slice_head(n = n_finalists)

flog.info(glue("Selected {nrow(finalists)} finalists for detailed analysis"))

generate_model_maps <- function(model_files, model_names, reference_rast, boundary, output_dir, pollutant, n_clas = 7) {
  flog.info(glue("Generating comprehensive PDF report for {pollutant}"))
  
  # Define output PDF path
  pdf_file <- glue("{output_dir}/model_comparison_{pollutant}.pdf")
  pdf(pdf_file, width = 10, height = 8)
  
  # Calculate consistent color breaks for concentrations and errors across all models
  all_model_vals <- unlist(lapply(model_files, function(f) {
    model_rast <- rast(f)
    model_rast <- crop(model_rast, boundary)
    model_rast <- mask(model_rast, boundary)
    vals <- values(model_rast)
    vals[!is.na(vals)]
  }))
  
  ref_vals <- values(reference_rast)
  ref_vals <- ref_vals[!is.na(ref_vals)]
  
  # Calculate quantile-based breaks for concentrations
  conc_pp <- c(0, 0.01, (1:(n_clas-3))/(n_clas-2), 0.99, 1)
  conc_bins <- unique(round(signif(quantile(c(all_model_vals, ref_vals), conc_pp, na.rm = TRUE), 2), 1))
  n_clas <- length(conc_bins)-1
  
  # Calculate symmetric breaks for errors (centered on zero) using absolute values
  all_errors <- unlist(lapply(model_files, function(f) {
    model_rast <- rast(f)
    model_rast <- crop(model_rast, boundary)
    model_rast <- mask(model_rast, boundary)
    error_vals <- values(model_rast - reference_rast)
    error_vals[!is.na(error_vals)]
  }))
  
  # Use absolute maximum error to create symmetric bins centered on zero
  n_clas_err <- 6
  err_pp <- c((0:(n_clas_err-1)/2)/((n_clas_err+1)/2),0.99,1)
  error_bins <- unique(round(signif(quantile(abs(all_errors), err_pp, na.rm = TRUE), 2), 1))
  error_bins <- setdiff(unique(sort(c(-error_bins,error_bins))), 0)
  max_err <- ceiling(max(abs(all_errors)))
  n_clas_err <- length(error_bins)-1
  error_bins[1] <- -max_err
  error_bins[length(error_bins)] <- max_err
  
  # Define color palettes
  conc_pal <- rev(brewer.pal(n_clas, "Spectral"))  # For concentrations
  error_pal <- rev(brewer.pal(n_clas_err, "RdBu"))     # For errors (diverging)
  
  # Plot reference map first
  par(mfrow = c(1, 1))
  plot(reference_rast, main = glue("Reference: {poll_ext}"), 
       col = conc_pal, breaks = conc_bins, 
       legend = "topright",
       mar = c(3, 1, 2, 1))
  mtext(text = bquote(.(glue("Concentration range: {round(min(values(reference_rast), na.rm = TRUE), 1)}",
                             " to {round(max(values(reference_rast), na.rm = TRUE), 1)}"))~mu*g/m^3), 
        side = 1, line = 4, cex = 0.8)
  plot(boundary, add = TRUE)
  
  # Plot each model and its error
  for (i in seq_along(model_files)) {
    flog.info(glue("Plotting model {i}/{length(model_files)}: {model_names[i]}"))
    
    # Load and prepare model data
    model_rast <- rast(model_files[i])
    model_rast <- crop(model_rast, boundary)
    model_rast <- mask(model_rast, boundary)
    
    error_rast <- model_rast - reference_rast
    
    # Create two-panel plot for each model
    par(mfrow = c(1, 2))
    
    # Model concentration map
    plot(model_rast, main = glue("Model: {model_names[i]}"), 
         col = conc_pal, breaks = conc_bins,
         mar = c(2, 1, 2, 1), 
         legend = "topright")
    mtext(text = bquote(.(glue("Concentration range: {round(min(values(model_rast), na.rm = TRUE), 1)}",
                               " to {round(max(values(model_rast), na.rm = TRUE), 1)}"))~mu*g/m^3), 
          side = 1, line = 4, cex = 0.8)
    plot(boundary, add = TRUE)
    
    # Error map with symmetric bins centered on zero
    plot(error_rast, main = glue("Error: {model_names[i]} - Reference"), 
         col = error_pal, breaks = error_bins,
         mar = c(2, 1, 2, 1), 
         legend = "topright")
    mtext(text = bquote(.(glue("Error range: {round(min(values(error_rast), na.rm = TRUE), 1)}",
                               " to {round(max(values(error_rast), na.rm = TRUE), 1)}"))~mu*g/m^3), 
          side = 1, line = 4, cex = 0.8)
    plot(boundary, add = TRUE)
  }
  
  dev.off()
  flog.info(glue("PDF report saved: {pdf_file}"))
}

# Generate reports for all finalists
model_idx <- match(finalists$model, model_names)
generate_model_maps(model_files[model_idx], model_names[model_idx], 
                    reference, italy_boundary, output_dir, pollutant)

# =============================================================================
# Phase 5: Scorecard and Final Selection
# =============================================================================
flog.info("Phase 5: Scorecard and Final Selection")

# Create comprehensive scorecard with ALL statistics
scorecard <- finalists %>%
  select(model, IOA, RMSE, ME, MAE, correlation, R2, P50, P75, P95, n_pixels)

# Save scorecard
write.csv(scorecard, glue("{output_dir}/model_scorecard_{pollutant}.csv"), row.names = FALSE)

# Save all results (including invalid models info)
results <- list(
  global_results = global_results,
  jaccard_results = jaccard_results,
  zone_errors = zone_errors_all,
  region_errors = region_errors_all,
  finalists = finalists,
  scorecard = scorecard,
  invalid_models = invalid_models
)

save(results, file = glue("{output_dir}/verification_results_{pollutant}.RData"))

flog.info(glue("Workflow completed for {pollutant}. Results saved in {output_dir}"))
flog.info("Please review visual reports and fill strengths/weaknesses in the scorecard")

# =============================================================================
# End of script
# =============================================================================