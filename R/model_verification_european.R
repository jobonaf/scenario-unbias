# =============================================================================
# title          :model_verification_european
# description    :Calculation of air quality model skill scores for European domain
#                 Simplified version - computes only global metrics and saves to CSV
# author         :Giovanni Bonafe'
# date           :20250203
# version        :1.1
# notes          :Requires model raster files and reference raster data;
#                 requires modules r-rgdal geos netcdf-fortran/4.5.2/gcc
# R_version      :3.5.2
# =============================================================================

# Define command-line options
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("glue"))
option_list <- list(
  make_option(c("-p", "--pollutant"), type = "character", default = "NO2",
              help = "Pollutant to analyze (NO2, O3, PM25) [default: %default]"),
  make_option(c("-r", "--reference_dir"), type = "character", 
              default = "data/fairmode-wg5-exercise-202501/YEARLY/Scenario_Reference_Gridded/",
              help = "Reference data directory [default: %default]"),
  make_option(c("-m", "--model_dir"), type = "character", 
              default = "data/processed_fairmode",
              help = "Model rasters directory [default: %default]"),
  make_option(c("-T", "--topo_file"), type = "character", 
              default = "data/fairmode-wg5-exercise-202501/topography.nc",
              help = "Topography NetCDF file for land mask [default: %default]"),
  make_option(c("-o", "--output_dir"), type = "character", 
              default = "data/models_verification/fairmode_exercise",
              help = "Output directory [default: %default]"),
  make_option(c("-L", "--use_land_mask"), type = "logical", 
              default = TRUE,
              help = "Apply land mask (topography > 1) [default: %default]")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Calculate skill scores for European domain models"))

# Extract parameters
pollutant <- opt$pollutant
reference_dir <- opt$reference_dir
model_dir <- opt$model_dir
topo_file <- opt$topo_file
output_dir <- opt$output_dir
use_land_mask <- opt$use_land_mask

# Load required libraries
library(terra)           # For raster operations
library(futile.logger)   # For logging
library(dplyr)           # For data manipulation

# =============================================================================
# Source external functions
# =============================================================================
source("R/read_netcdf_as_raster.R")
source("R/topography_from_nc.R")

# =============================================================================
# Initialize logging
# =============================================================================
flog.appender(appender.console())  # Write to standard output
flog.threshold(INFO)
flog.info("Starting European domain skill scores calculation")
flog.info("Pollutant: %s", pollutant)
flog.info("Using parameters: %s", paste(commandArgs(trailingOnly = TRUE), collapse = " "))

# =============================================================================
# Data Preparation
# =============================================================================
flog.info("Loading reference data...")

# Handle special naming for PM25
pollutant_ref <- ifelse(pollutant == "PM25", "PM25_rh50", pollutant)
reference_file <- glue("{reference_dir}/SCEN_REF_{pollutant_ref}_YEARLY.nc")

# Check if reference file exists
if (!file.exists(reference_file)) {
  flog.error(glue("Reference file not found: {reference_file}"))
  stop(glue("Reference file not found: {reference_file}"))
}

# Load reference data
flog.info(glue("Loading reference: {reference_file}"))
reference <- read_nc(reference_file)

# Get model files
model_pattern <- glue("{model_dir}/{pollutant}_*_scenario.tif")
model_files <- Sys.glob(model_pattern)

if (length(model_files) == 0) {
  flog.error(glue("No model files found matching pattern: {model_pattern}"))
  stop(glue("No model files found matching pattern: {model_pattern}"))
}

# Extract model names - remove pollutant prefix and both "_unbiased_scenario" and "_scenario" suffixes
model_names <- basename(model_files)
model_names <- tools::file_path_sans_ext(model_names)
model_names <- gsub(glue("^{pollutant}_"), "", model_names)
model_names <- gsub("_unbiased_scenario$", "", model_names)
model_names <- gsub("_scenario$", "", model_names)

flog.info(glue("Found {length(model_files)} model files for {pollutant}"))

# Load topography and create land mask (optional)
land_mask <- NULL
if (use_land_mask) {
  flog.info("Creating land mask from topography...")
  
  if (!file.exists(topo_file)) {
    flog.warn(glue("Topography file not found: {topo_file}. Proceeding without land mask."))
    use_land_mask <- FALSE
  } else {
    # Load topography
    topo <- topography_from_nc(topo_file, grid = ext(reference))
    crs(topo) <- crs(reference)
    
    # Create land mask (values > 1 indicate land)
    land_mask <- topo > 1
    
    # Apply land mask to reference
    reference <- mask(x = reference, mask = land_mask, maskvalues = 0)
    
    flog.info("Land mask applied successfully")
  }
}

# =============================================================================
# Helper function: Check for valid raster data
# =============================================================================
check_raster_validity <- function(model_path, reference_rast, land_mask = NULL) {
  flog.info(glue("Checking file {basename(model_path)}"))
  tryCatch({
    model_rast <- rast(model_path)
    
    if (is.null(model_rast) || nlyr(model_rast) == 0) {
      return(FALSE)
    }
    
    # Match extent and resolution with reference
    model_rast <- resample(model_rast, reference_rast, method = "bilinear")
    
    # Apply land mask if provided
    if (!is.null(land_mask)) {
      model_rast <- mask(x = model_rast, mask = land_mask, maskvalues = 0)
    }
    
    valid_vals <- values(model_rast)
    valid_count <- sum(!is.na(valid_vals))
    
    return(valid_count > 0)
    
  }, error = function(e) {
    flog.error(glue("Error with file {basename(model_path)}: {e$message}"))
    return(FALSE)
  })
}

# =============================================================================
# Helper function: Calculate global metrics
# =============================================================================
calculate_global_metrics <- function(model_path, reference_rast, land_mask = NULL) {
  # Load and prepare model data
  model_rast <- rast(model_path)
  
  # Match extent and resolution with reference
  model_rast <- resample(model_rast, reference_rast, method = "bilinear")
  
  # Apply land mask if provided
  if (!is.null(land_mask)) {
    model_rast <- mask(x = model_rast, mask = land_mask, maskvalues = 0)
  }
  
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

# =============================================================================
# Check validity of all models
# =============================================================================
flog.info("Checking validity of model files...")
valid_models <- sapply(model_files, check_raster_validity, 
                       reference_rast = reference, 
                       land_mask = land_mask)
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
# Calculate skill scores for all valid models
# =============================================================================
flog.info("Calculating skill scores...")

skill_scores <- list()
for (i in seq_along(model_files)) {
  flog.info(glue("Processing model {i}/{length(model_files)}: {model_names[i]}"))
  
  skill_scores[[model_names[i]]] <- calculate_global_metrics(
    model_files[i], 
    reference, 
    land_mask
  )
}

# =============================================================================
# Convert to dataframe and save
# =============================================================================
flog.info("Preparing output...")

# Convert to dataframe
skill_df <- do.call(rbind, lapply(skill_scores, as.data.frame))
skill_df$model <- rownames(skill_df)
rownames(skill_df) <- NULL

# Reorder columns
skill_df <- skill_df %>%
  select(model, ME, MAE, RMSE, R2, IOA, correlation, n_pixels)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save to CSV
output_file <- glue("{output_dir}/skill_scores_{pollutant}.csv")
write.csv(skill_df, output_file, row.names = FALSE)

flog.info(glue("Skill scores saved to: {output_file}"))
flog.info(glue("Processed {nrow(skill_df)} models successfully"))

# Print summary statistics
flog.info("Summary statistics:")
flog.info(glue("  Mean RMSE: {round(mean(skill_df$RMSE), 3)}"))
flog.info(glue("  Mean MAE: {round(mean(skill_df$MAE), 3)}"))
flog.info(glue("  Mean correlation: {round(mean(skill_df$correlation), 3)}"))
flog.info(glue("  Mean IOA: {round(mean(skill_df$IOA), 3)}"))
flog.info(glue("  Mean R2: {round(mean(skill_df$R2), 3)}"))

flog.info("European domain skill scores calculation completed successfully!")

# =============================================================================
# End of script
# =============================================================================