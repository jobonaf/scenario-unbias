suppressPackageStartupMessages(library("optparse"))

# Define command-line options
option_list <- list(
  make_option(c("-p", "--pollutants"), type = "character", default = "NO2,O3,PM25",
              help = "Comma-separated list of pollutants to process [default: %default]"),
  make_option(c("-o", "--output_dir"), type = "character", default = "data/processed",
              help = "Output directory for processed data [default: %default]"),
  make_option(c("-u", "--unbias_sequences"), type = "character", default = "SCA,CSA,CAS,CA",
              help = "Comma-separated list of unbias sequences [default: %default]"),
  make_option(c("-c", "--calibration_methods"), type = "character", default = "All,Each,Grid,Cell",
              help = "Comma-separated list of calibration methods [default: %default]"),
  make_option(c("-a", "--correction_algorithms"), type = "character", default = "Add,Mult,Lin",
              help = "Comma-separated list of correction algorithms [default: %default]"),
  make_option(c("-s", "--spatialization_methods"), type = "character", default = "tps,idw,ok,ked",
              help = "Comma-separated list of spatialization methods [default: %default]")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Convert comma-separated string inputs into lists
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

# Load external scripts containing necessary functions
source("R/read-fairmode-data.R")
source("R/unbias-aq-scenario.R")

# Create output directory if it does not exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Function to check if a given combination of unbias sequence and calibration method is valid
is_valid_combination <- function(unbias_sequence, calibration_method, correction_algorithm) {
  if (!correction_algorithm %in% c("Add", "Mult") && calibration_method %in% c("Each", "Cell")) {
    return(FALSE)
  }
  if (unbias_sequence == "SCA" && calibration_method %in% c("Each", "All")) {
    return(FALSE)
  }
  if (unbias_sequence %in% c("CAS", "CA") && !calibration_method %in% c("Each", "All")) {
    return(FALSE)
  }
  if (unbias_sequence == "CSA" && calibration_method != "Each") {
    return(FALSE)
  }
  return(TRUE)
}

# Function to process a specific combination of parameters
process_combination <- function(pollutant, output_dir, unbias_sequence, 
                                calibration_method, correction_algorithm, 
                                spatialization_method = NULL) {
  flog.info("Processing pollutant: %s with combination: %s.%s.%s%s", 
            pollutant, unbias_sequence, calibration_method, correction_algorithm,
            ifelse(is.null(spatialization_method), "", paste0(".", spatialization_method)))
  
  # Check if the combination is valid
  if (!is_valid_combination(unbias_sequence, calibration_method, correction_algorithm)) {
    flog.warn("Invalid combination for pollutant %s: %s.%s.%s%s", 
              pollutant, unbias_sequence, calibration_method, 
              correction_algorithm, ifelse(is.null(spatialization_method), "", paste0(".", spatialization_method)))
    return(NULL)
  }
  
  # Read input data for the pollutant
  flog.info("Reading data for pollutant: %s", pollutant)
  data_list <- read_data(pollutant)
  
  # Apply the unbiasing function
  flog.info("Applying unbiasing function for pollutant: %s", pollutant)
  unbias_result <- process_data(
    observed_data = data_list$observed_data,
    base_case = data_list$base_case,
    scenario = data_list$scenario,
    unbias_sequence = unbias_sequence,
    calibration_method = calibration_method,
    correction_algorithm = correction_algorithm,
    spatialization_method = ifelse(unbias_sequence == "CA", NULL, spatialization_method)
  )
  
  # Define the output file name based on parameters
  fileout <- if (unbias_sequence == "CA") {
    glue("{output_dir}/{pollutant}_{unbias_sequence}.{calibration_method}.{correction_algorithm}")
  } else {
    glue("{output_dir}/{pollutant}_{unbias_sequence}.{calibration_method}.{correction_algorithm}.{spatialization_method}")
  }
  
  # Determine the type of output to save based on the data structure
  if (inherits(unbias_result, "SpatRaster")) {
    fileout <- paste0(fileout, "_unbiased_scenario.tif")
    flog.info("Saving processed raster to file: %s", fileout)
    writeRaster(unbias_result, filename = fileout, overwrite = TRUE)
  } else if (is.data.frame(unbias_result)) {
    fileout <- paste0(fileout, "_unbiased_scenario.csv")
    flog.info("Saving processed data frame to file: %s", fileout)
    write.csv(unbias_result, fileout, row.names = FALSE)
  } else {
    stop("Unexpected data type for 'unbias_result'. Expected SpatRaster or data.frame.")
  }
}

# Main processing loop
flog.info("Starting data processing...")
for (pollutant in pollutants) {
  for (unbias_sequence in unbias_sequences) {
    for (calibration_method in calibration_methods) {
      for (correction_algorithm in correction_algorithms) {
        
        if (unbias_sequence == "CA") {
          # Call process_combination only once without iterating over spatialization_methods
          process_combination(
            pollutant = pollutant,
            output_dir = output_dir,
            unbias_sequence = unbias_sequence,
            calibration_method = calibration_method,
            correction_algorithm = correction_algorithm
          )
        } else {
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
}
flog.info("Data processing completed. Results saved to: %s", output_dir)
