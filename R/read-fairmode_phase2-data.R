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
                      data_path = "/u/arpa/bonafeg/src/scenario-unbias/data/fairmode-wg5-exercise-202602/YEARLY",
                      preproc_obs = c("exclude_industrial", "exclude_traffic", 
                                      "closest_to_center")) {
  
  # Helper function to read NetCDF file as SpatRaster
  read_nc <- function(file, varname) {
    if (!file.exists(file)) {
      stop("File not found: ", file)
    }
    
    nc_data <- nc_open(file)
    
    # Check variable exists
    if (!(varname %in% names(nc_data$var))) {
      stop("Variable ", varname, " not found in file: ", file)
    }
    
    var_data <- ncvar_get(nc_data, varname)
    
    lon <- ncvar_get(nc_data, "lon")
    lat <- ncvar_get(nc_data, "lat")
    
    res_x <- mean(diff(lon))
    res_y <- mean(diff(lat))
    
    xmin_new <- min(lon) - res_x / 2
    xmax_new <- max(lon) + res_x / 2
    ymin_new <- min(lat) - res_y / 2
    ymax_new <- max(lat) + res_y / 2
    
    raster_data <- rast(t(var_data)[ncol(var_data):1,], crs="EPSG:4326")
    ext(raster_data) <- c(xmin_new, xmax_new, ymin_new, ymax_new)
    
    nc_close(nc_data)
    
    return(raster_data)
  }
  
  # For fairmode_phase2, scenario_year must be provided
  if (is.null(scenario_year)) {
    stop("scenario_year must be specified for fairmode_phase2 exercise")
  }
  
  # Define codes for parameters
  obs_code <- case_when(
    parameter == "NO2" ~ "NO2",
    parameter == "O3" ~ "O3",
    parameter == "PM25" ~ "PM25"
  )
  mod_code <- case_when(
    parameter == "NO2" ~ "SURF_ug_NO2",
    parameter == "O3" ~ "SURF_ug_O3",
    parameter == "PM25" ~ "SURF_ug_PM25_rh50"
  )
  
  # Load base case (2015) gridded data
  base_case_gridded_file <- glue("{data_path}/BaseCase_2015_Gridded/EMEP_yearly_2015.nc")
  base_case <- read_nc(base_case_gridded_file, mod_code)
  
  # Load observed data from 2015
  observed_data_file <- glue("{data_path}/BaseCase_2015_Points/yearly_{obs_code}_2015.csv")
  observed_data <- read_csv(observed_data_file, show_col_types = FALSE) %>%
    rename(x = Longitude, y = Latitude, value = Average)
  
  # Exclude stations by type
  if("exclude_industrial" %in% preproc_obs) observed_data <- observed_data %>% filter(Type != "Industrial")
  if("exclude_traffic" %in% preproc_obs)    observed_data <- observed_data %>% filter(Type != "Traffic")
  if("exclude_background" %in% preproc_obs) observed_data <- observed_data %>% filter(Type != "Background")
  
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
    scenario <- read_nc(scenario_gridded_file, mod_code)
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