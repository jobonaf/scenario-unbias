# Load required packages
library(terra)
library(ncdf4)
library(glue)
library(futile.logger)
library(stringr)
library(dplyr)

# Function to convert a TIFF to NetCDF with grid rounding
tiff_to_netcdf <- function(input_file, output_file = NULL, 
                           decimals = 5, overwrite = FALSE,
                           varname = "conc", unit = "ug/m3") {
  # If output_file is not specified, use the input filename with .nc extension
  if (is.null(output_file)) {
    output_file <- sub("\\.tif$", ".nc", input_file, ignore.case = TRUE)
  }
  
  # Check if the output file already exists
  if (file.exists(output_file) && !overwrite) {
    stop("Output file already exists. Use overwrite = TRUE to replace it.")
  }
  
  # Load the raster file
  raster_data <- rast(input_file)
  
  # Round the grid extent to the specified number of decimal places
  xmin_new <- round(ext(raster_data)$xmin, decimals)
  xmax_new <- round(ext(raster_data)$xmax, decimals)
  ymin_new <- round(ext(raster_data)$ymin, decimals)
  ymax_new <- round(ext(raster_data)$ymax, decimals)
  
  # Update raster extent
  ext(raster_data) <- ext(xmin_new, xmax_new, ymin_new, ymax_new)
  
  # Get raster dimensions
  nx <- ncol(raster_data)
  ny <- nrow(raster_data)
  
  # NetCDF refers to cell center, while GeoTIFF to cell border
  xmin_nc <- xmin_new + xres(raster_data) * 0.5
  xmax_nc <- xmax_new - xres(raster_data) * 0.5
  ymin_nc <- ymin_new + yres(raster_data) * 0.5
  ymax_nc <- ymax_new - yres(raster_data) * 0.5
  
  # Generate coordinate values
  lon_values <- round(seq(xmin_nc, xmax_nc, length.out = nx), digits = decimals)
  lat_values <- round(seq(ymin_nc, ymax_nc, length.out = ny), digits = decimals)

  # Define NetCDF dimensions
  lon_dim <- ncdim_def(name = "lon", units = "degrees_east", vals = lon_values, longname = "Longitude")
  lat_dim <- ncdim_def(name = "lat", units = "degrees_north", vals = lat_values, longname = "Latitude")
  
  # Define the NetCDF variable
  var_def <- ncvar_def(name = varname, units = unit, dim = list(lon_dim, lat_dim), 
                       prec = "double", missval = NA_real_)
  
  # Create the NetCDF file
  nc <- nc_create(output_file, vars = list(var_def))
  
  # Write coordinate values
  ncvar_put(nc, "lon", lon_values)
  ncvar_put(nc, "lat", lat_values)
  
  # Convert to matrix preserving spatial structure
  corrected_matrix <- as.matrix(raster_data, wide = TRUE)
  
  # Transpose to match NetCDF convention (lon, lat)
  corrected_matrix <- t(corrected_matrix)
  
  # Reverse rows to ensure latitude is from south to north
  corrected_matrix <- corrected_matrix[, ny:1]
  
  # Write raster data
  ncvar_put(nc, varname, corrected_matrix)
  
  # Close the NetCDF file
  nc_close(nc)
  
  cat("NetCDF file saved:", output_file, "\n")
}

# Example usage:
# tiff_to_netcdf("input.tif", decimals = 5, overwrite = TRUE)

exercise <- "phase2"

if(exercise=="phase1") {
  #---------------------
  # FAIRMODE WG5 Phase 1
  # Naming convention:
  # <scenario>_<group>_<species>_<info>_CORR_YEARLY.nc
  # group:   ITAWG
  # species: O3, PM25, NO2
  # info:    <sequence>.<calibration>.<correction>.<spatialization>
  scen_in <- "unbiased_basecase"
  scen_out <- "BaseCase"
  group <- "ITAWG"
  for (f in Sys.glob(glue("data/processed_fairmode/*_{scen_in}.tif"))) {
    specie <- strsplit(basename(f), "_")[[1]][1]
    varname <- case_match(specie, "O3"~"SURF_ppb_O3", "PM25"~"SURF_ug_PM25_rh50", "NO2"~"SURF_ug_NO2")
    unit <- case_match(specie, "O3"~"ppb", "PM25"~"ug/m3", "NO2"~"ug/m3")
    info <- strsplit(basename(f), "_")[[1]][2]
    infos <- strsplit(info,"\\.")[[1]]
    if(length(infos)==4) infos[4] <- str_to_upper(infos[4])
    fileout <- glue("data/fairmode-wg5-exercise-output/{scen_out}_{group}_{specie}_{paste(infos,collapse='.')}_CORR_YEARLY.nc")
    tiff_to_netcdf(input_file = f, output_file = fileout, overwrite = T, varname = varname, unit = unit)
    flog.info(glue("Written file {fileout}"))
  }
}

if(exercise=="phase2") {
  #---------------------
  # FAIRMODE WG5 Phase 2
  # Naming convention:
  # Scen_<year>_<group>_<species>_<info>_CORR_YEARLY.nc
  # group:   ITAWG
  # species: O3, PM25, NO2
  # info:    <sequence>.<calibration>.<correction>.<spatialization>
  scen_in <- "unbiased_scenario"
  dir_in <- "data/processed_phase2"
  dir_out <- "data/fairmode-wg5-exercise-phase2-output"
  group <- "ITAWG"
  for (f in Sys.glob(glue("{dir_in}/*_{scen_in}_????.tif"))) {
    specie <- strsplit(basename(f), "_")[[1]][1]
    varname <- specie
    unit <- case_match(specie, "O3"~"ug/m3", "PM25"~"ug/m3", "NO2"~"ug/m3")
    info <- strsplit(basename(f), "_")[[1]][2]
    infos <- strsplit(info,"\\.")[[1]]
    year <- tools::file_path_sans_ext(strsplit(basename(f), "_")[[1]][5])
    fileout <- glue("{dir_out}/Scen_{year}_{group}_{specie}_{paste(infos,collapse='.')}_CORR_YEARLY.nc")
    tiff_to_netcdf(input_file = f, output_file = fileout, overwrite = T, varname = varname, unit = unit)
    flog.info(glue("Written file {fileout}"))
  }
}
