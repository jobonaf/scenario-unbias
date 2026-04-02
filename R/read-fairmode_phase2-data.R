# Load necessary packages
library(terra)
library(readr)
library(dplyr)
library(glue)
library(ncdf4)
library(futile.logger)

# Function to read gridded (.nc) and observed (.csv) data for a specific parameter and scenario
read_data <- function(parameter,
                      scenario_year = NULL,  # Can be "2015", "2022", "2023", or "2024"
                      data_path = "/u/arpa/bonafeg/src/scenario-unbias/data/fairmode-wg5-exercise-202602/YEARLY") {
  
  # Helper function to read NetCDF file as SpatRaster
  read_nc <- function(file) {
    # Check if file exists
    if (!file.exists(file)) {
      stop("File not found: ", file)
    }
    
    # Open the NetCDF file
    nc_data <- nc_open(file)
    
    # Extract variable data (assuming it's the first variable)
    var_data <- ncvar_get(nc_data, nc_data$var[[1]]$name)
    
    # Get the latitude and longitude coordinates
    lon <- ncvar_get(nc_data, "lon")  # Longitude
    lat <- ncvar_get(nc_data, "lat")  # Latitude
    
    # Compute resolution (assuming uniform grid spacing)
    res_x <- mean(diff(lon))  # Resolution in X direction
    res_y <- mean(diff(lat))  # Resolution in Y direction
    
    # Compute new extent by shifting from center-based to corner-based convention
    xmin_new <- min(lon) - res_x / 2
    xmax_new <- max(lon) + res_x / 2
    ymin_new <- min(lat) - res_y / 2
    ymax_new <- max(lat) + res_y / 2
    
    # Convert data to SpatRaster (flipping Y direction if necessary)
    raster_data <- rast(t(var_data)[ncol(var_data):1,], crs="EPSG:4326")
    ext(raster_data) <- c(xmin_new, xmax_new, ymin_new, ymax_new)
    
    # Close the NetCDF file
    nc_close(nc_data)
    
    return(raster_data)
  }
  
  # For fairmode_phase2, scenario_year must be provided
  if (is.null(scenario_year)) {
    stop("scenario_year must be specified for fairmode_phase2 exercise")
  }
  
  # Define file paths for each type of data
  obs_code <- case_when(
    parameter == "NO2" ~ "NO2",
    parameter == "O3" ~ "O3",
    parameter == "PM25" ~ "PM25"
  )
  
  # Load base case (2015) gridded data
  base_case_gridded_file <- glue("{data_path}/BaseCase_2015_Gridded/EMEP_yearly_2015.nc")
  base_case <- read_nc(base_case_gridded_file)
  
  # Load observed data from 2015
  observed_data_file <- glue("{data_path}/BaseCase_2015_Points/yearly_{obs_code}_2015.csv")
  observed_data <- read_csv(observed_data_file, show_col_types = FALSE) %>%
    rename(x = Longitude, y = Latitude, value = Average)
  
  # Get model bounding box (use base_case grid)
  bb <- ext(base_case)
  
  # Filter observed stations within bounding box
  observed_data <- observed_data %>%
    filter(
      x >= bb$xmin,
      x <= bb$xmax,
      y >= bb$ymin,
      y <= bb$ymax
    )
  
  # Load scenario gridded data based on scenario_year
  if (scenario_year == "2015") {
    # If scenario_year is 2015, use the base case gridded data as scenario
    scenario <- base_case
    flog.info("Using base case (2015) gridded data as scenario data")
  } else {
    # Load scenario gridded data for the specified year
    scenario_gridded_file <- glue("{data_path}/Scenario_{scenario_year}_Gridded/EMEP_yearly_{scenario_year}.nc")
    scenario <- read_nc(scenario_gridded_file)
    flog.info("Loaded scenario gridded data for year: %s", scenario_year)
  }
  
  # Return the data as a list with the SAME structure as fairmode and italian exercises
  # This ensures complete backward compatibility with the calling script
  return(list(
    observed_data = observed_data,
    base_case = base_case,      # Always the 2015 base case
    scenario = scenario         # The target grid (could be 2015 or scenario year)
  ))
}