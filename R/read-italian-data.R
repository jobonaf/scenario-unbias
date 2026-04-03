# Load necessary packages
library(terra)
library(readr)
library(dplyr)
library(glue)
library(ncdf4)

# Function to read gridded (.nc) and observed (.csv) data for a specific parameter
# Modified for the new file structure with different naming conventions
read_data <- function(parameter,
                      data_path = "data/esercizio-dominio-italiano",
                      preproc_obs = c("exclude_industrial", "exclude_traffic", 
                                      "closest_to_center")) {
  
  # Define file name components for each type of data
  # Note: The new structure uses different naming patterns
  obs_file_suffix <- case_when(
    parameter == "NO2" ~ "NO2_BASE_extract_height.csv",
    parameter == "O3" ~ "O3_BASE_extract_height.csv",
    parameter == "PM10" ~ "PM10_BASE_extract_height.csv",
    parameter == "PM25" ~ "PM25_BASE_extract_height.csv"
  )
  
  # Define file patterns for gridded data
  base_case_perturbed <- glue("{data_path}/BaseCase_Perturbed_Gridded/SCEN_SC_B_c_{parameter}_YEARLY.nc")
  scenario_perturbed <- glue("{data_path}/Scenario_Perturbed_Gridded/SCEN_SC_AB_c_{parameter}_YEARLY.nc")
  base_case_reference <- glue("{data_path}/BaseCase_Reference_Gridded/SCEN_BASE_c_{parameter}_YEARLY.nc")
  scenario_reference <- glue("{data_path}/Scenario_Reference_Gridded/SCEN_SC_A_c_{parameter}_YEARLY.nc")
  observed_data_file <- glue("{data_path}/BaseCase_Reference_Points/{obs_file_suffix}")
  
  # Function to read a NetCDF file as SpatRaster, adjusting extent if needed
  source("R/read_netcdf_as_raster.R")
  
  # Read the gridded data
  base_case <- read_nc(base_case_perturbed)
  scenario <- read_nc(scenario_perturbed)
  base_case_reference <- read_nc(base_case_reference)
  scenario_reference <- read_nc(scenario_reference)
  
  # Read the observed data and rename columns to x, y, and value
  observed_data <- read_csv(observed_data_file, show_col_types = FALSE) 
  
  # Exclude stations by type
  if("exclude_industrial" %in% preproc_obs) observed_data <- observed_data %>% filter(AirQualityStationType != "industrial")
  if("exclude_traffic" %in% preproc_obs)    observed_data <- observed_data %>% filter(AirQualityStationType != "traffic")
  if("exclude_background" %in% preproc_obs) observed_data <- observed_data %>% filter(AirQualityStationType != "background")
  observed_data <- observed_data %>%
    transmute(x = Longitude, y = Latitude, value = Value_sampled)  
  
  # Keep only one point for each cell, the closest to the cell center
  if("closest_to_center" %in% preproc_obs) {
    observed_data <- observed_data %>%
      mutate(cell = cellFromXY(base_case, cbind(x, y))) %>%
      group_by(cell) %>%
      slice_min(order_by = sqrt((x - xyFromCell(base_case, cell)[,1])^2 +
                                  (y - xyFromCell(base_case, cell)[,2])^2)) %>%
      ungroup() %>%
      select(-cell)
  }
  
  # Move each point to the cell center (suitable if data are "synthetic" observations)
  if("move_to_center" %in% preproc_obs) {
    observed_data <- observed_data %>%
      mutate(cell = cellFromXY(base_case, cbind(x, y)),
             x = xyFromCell(base_case, cell)[,1],
             y = xyFromCell(base_case, cell)[,2]) %>%
      select(-cell)
  }
  
  # Return the data as a list
  return(list(observed_data = observed_data, 
              base_case = base_case, 
              scenario = scenario, 
              base_case_reference = base_case_reference, 
              scenario_reference = scenario_reference))
}

# Example usage:
# no2_data <- read_data("NO2")
# pm25_data <- read_data("PM25")