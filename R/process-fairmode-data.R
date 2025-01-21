
# Define command-line arguments
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-p", "--pollutants"), type = "character", default = "NO2,O3,PM25",
              help = "Comma-separated list of pollutants to process [default: %default]"),
  make_option(c("-o", "--output_dir"), type = "character", default = "data/processed",
              help = "Output directory for processed data [default: %default]"),
  make_option(c("-u", "--unbias_sequences"), type = "character", default = "SCA,CSA,CAS,CA",
              help = "Comma-separated list of unbias sequences [default: %default]"),
  make_option(c("-c", "--calibration_methods"), type = "character", default = "Point,Grid,Cell",
              help = "Comma-separated list of calibration methods [default: %default]"),
  make_option(c("-a", "--correction_algorithms"), type = "character", default = "Add,Mult,Lin",
              help = "Comma-separated list of correction algorithms [default: %default]"),
  make_option(c("-s", "--spatialization_methods"), type = "character", default = "tps,idw,ok,ked",
              help = "Comma-separated list of spatialization methods [default: %default]")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))
pollutants <- strsplit(opt$pollutants, ",")[[1]]
output_dir <- opt$output_dir
unbias_sequences <- strsplit(opt$unbias_sequences, ",")[[1]]
calibration_methods <- strsplit(opt$calibration_methods, ",")[[1]]
correction_algorithms <- strsplit(opt$correction_algorithms, ",")[[1]]
spatialization_methods <- strsplit(opt$spatialization_methods, ",")[[1]]

# Load necessary libraries
library(dplyr)
library(terra)
library(glue)
library(futile.logger)

# scripts
source("R/read-fairmode-data.R")
source("R/unbias-aq-scenario.R")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Function to process data for a single combination
process_combination <- function(pollutant, output_dir, unbias_sequence, 
                                calibration_method, correction_algorithm, 
                                spatialization_method) {
  flog.info("Processing pollutant: %s with combination: %s.%s.%s.%s", 
            pollutant, unbias_sequence, calibration_method, correction_algorithm, spatialization_method)
  
  # Check if the combination is valid
  if (!is_valid_combination(unbias_sequence, calibration_method)) {
    flog.warn("Invalid combination for pollutant %s: %s.%s.%s.%s", 
              pollutant, unbias_sequence, calibration_method, 
              correction_algorithm, spatialization_method)
    return(NULL)
  }
  
  # Read and process data
  flog.info("Reading data for pollutant: %s", pollutant)
  data_list <- read_data(pollutant)
  
  flog.info("Applying unbiasing function for pollutant: %s", pollutant)
  unbias_result <- process_data(
    observed_data = data_list$observed_data,
    base_case = data_list$base_case,
    scenario = data_list$scenario,
    unbias_sequence = unbias_sequence,
    calibration_method = calibration_method,
    correction_algorithm = correction_algorithm,
    spatialization_method = spatialization_method
  )
  
  # Save processed data
  fileout <- glue(
    "{output_dir}/{pollutant}_",
    "{unbias_sequence}.{calibration_method}.{correction_algorithm}.{spatialization_method}",
    "_unbiased_scenario.tif")
  flog.info("Saving processed raster to file: %s", fileout)
  writeRaster(unbias_result, filename = fileout, overwrite = TRUE)
  
  flog.info("Completed processing for pollutant: %s", pollutant)
}

# Function to validate combination
is_valid_combination <- function(unbias_sequence, calibration_method) {
  if (unbias_sequence == "SCA" && calibration_method == "Point") {
    return(FALSE)
  }
  if (unbias_sequence %in% c("CSA", "CAS", "CA") && calibration_method != "Point") {
    return(FALSE)
  }
  return(TRUE)
}

# Main script
flog.info("Starting data processing...")
for (pollutant in pollutants) {
  for (unbias_sequence in unbias_sequences) {
    for (calibration_method in calibration_methods) {
      for (correction_algorithm in correction_algorithms) {
        for (spatialization_method in spatialization_methods) {
          process_combination(
            pollutant = pollutant,
            output_dir = output_dir,
            unbias_sequence = unbias_sequence,
            calibration_method = calibration_method,
            correction_algorithm = correction_algorithm,
            spatialization_method = spatialization_method
          )
        }
      }
    }
  }
}
flog.info("Data processing completed. Results saved to: %s", output_dir)
