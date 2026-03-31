# Load necessary packages
library(terra)
library(readr)
library(dplyr)
library(glue)
library(ncdf4)

# Function to read gridded (.nc) and observed (.csv) data for a specific parameter
read_data <- function(parameter,
                      data_path="data/fairmode-wg5-exercise-202501/YEARLY") {
  
  # Define file paths for each type of data
  obs_code <- case_when(
    parameter == "NO2" ~ "ug_NO2",
    parameter == "O3" ~ "ppb_O3",
    parameter == "PM25" ~ "ug_PM25_rh50"
  )
  mod_code <- case_when(
    parameter == "PM25" ~ "PM25_rh50",
    .default = parameter
  )
  base_case_gridded  <- glue("{data_path}/BaseCase_Perturbed_Gridded/BaseCase_PERT_{mod_code}_YEARLY.nc")
  scenario_gridded   <- glue("{data_path}/Scenario_Perturbed_Gridded/SCEN_PERT_{mod_code}_YEARLY.nc")
  observed_data_file <- glue("{data_path}/BaseCase_Reference_Points/yearly_SURF_{obs_code}.csv")
  
  # Function to read a NetCDF file as SpatRaster, adjusting extent if needed
  read_nc <- function(file) {
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
  
  # Read the gridded data
  base_case <- read_nc(base_case_gridded)
  scenario <- read_nc(scenario_gridded)
  
  # Read the observed data and rename columns to x, y, and value
  observed_data <- read_csv(observed_data_file, show_col_types = FALSE) %>%
    rename(x = 1, y = 2, value = 3)  # Assuming x, y, value are in the first three columns
  
  # Return the data as a list
  return(list(observed_data = observed_data, base_case = base_case, scenario = scenario))
}

