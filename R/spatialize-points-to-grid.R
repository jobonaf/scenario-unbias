library(tidyverse)
library(terra)
library(fields)
library(gstat)
library(futile.logger)
library(doParallel)
library(rlang)
library(glue)

# Configura logging
flog.threshold(INFO)
flog.appender(appender.console())

# Spatialization function (updated with detailed logging)
spatialize <- function(points_data, grid_data, 
                       method = c("tps", "idw", "ok", "ked"), 
                       params = list(), ...) {
  
  # Validate input
  if (!inherits(grid_data, "SpatRaster")) {
    flog.error("'grid_data' must be a SpatRaster object.")
    stop("'grid_data' must be a SpatRaster object.")
  }
  if (!is.data.frame(points_data) || !all(c("x", "y", "value") %in% colnames(points_data))) {
    flog.error("'points_data' must be a data frame with 'x', 'y' and 'value' columns.")
    stop("'points_data' must be a data frame with 'x', 'y' and 'value' columns.")
  }
  
  # Prepare output raster (same grid as grid_data)
  empty_grid <- rast(grid_data)
  result <- empty_grid
  
  # Log input data details
  flog.info(glue("Starting spatialization with method '{method}'"))
  flog.info(glue("Points data: {nrow(points_data)} observations"))
  flog.info(glue("Grid dimensions: {nrow(grid_data)} x {ncol(grid_data)} cells"))
  flog.info(glue("Grid extent: {as.character(ext(grid_data))}"))
  flog.info(glue("Grid CRS: {crs(grid_data)}"))
  
  # Debug data validation
  flog.debug(glue("NA values in points_data$value: {sum(is.na(points_data$value))}"))
  flog.debug(glue("Value range: [{min(points_data$value, na.rm = TRUE)}, {max(points_data$value, na.rm = TRUE)}]"))
  flog.debug(glue("X coordinate range: [{min(points_data$x)}, {max(points_data$x)}]"))
  flog.debug(glue("Y coordinate range: [{min(points_data$y)}, {max(points_data$y)}]"))
  
  # Check spatial overlap
  points_ext <- ext(c(range(points_data$x), range(points_data$y)))
  grid_ext <- ext(grid_data)
  overlap_check <- points_ext[1] >= grid_ext[1] & points_ext[2] <= grid_ext[2] &
    points_ext[3] >= grid_ext[3] & points_ext[4] <= grid_ext[4]
  
  if (!overlap_check) {
    flog.warn(glue("Points extent ({as.character(points_ext)}) does not fully overlap with grid extent ({as.character(grid_ext)})"))
  }
  
  # Select spatialization method
  flog.info(glue("Spatializing with method '{method}'"))
  
  if (method == "tps") {
    # Thin Plate Spline spatialization
    flog.debug("Fitting Thin Plate Spline model")
    tps_model <- Tps(points_data[, c("x", "y")], points_data$value, ...)
    result <- interpolate(empty_grid, tps_model)
    
  } else if (method == "idw") {
    # Inverse Distance Weighted spatialization
    idp  <- params$idp %||% 2
    nmax <- params$nmax %||% Inf
    flog.debug(glue("Performing IDW with idp = {idp}, nmax = {nmax}"))
    
    gstat_model <- gstat(NULL, id = "var", formula = var ~ 1, locations = ~x + y, 
                         data = data.frame(points_data, var = points_data$value), 
                         nmax = nmax, set = list(idp = idp), ...)
    result <- interpolate(empty_grid, gstat_model, index = 1)
    
  } else if (method == "ok") {
    # Ordinary Kriging spatialization with detailed debugging
    flog.info("Starting Ordinary Kriging process")
    
    # Prepare data for gstat
    gstat_data <- data.frame(x = points_data$x, y = points_data$y, var = points_data$value)
    
    # Empirical variogram
    flog.debug("Calculating empirical variogram")
    sample_variogram <- tryCatch({
      v <- variogram(var ~ 1, ~x + y, data = gstat_data)
      flog.debug(glue("Empirical variogram calculated with {nrow(v)} bins"))
      flog.debug(glue("Variogram distance range: [{min(v$dist)}, {max(v$dist)}]"))
      v
    }, error = function(e) {
      flog.error(glue("Error calculating variogram: {e$message}"))
      stop(e)
    })
    
    # Variogram model parameters with sensible defaults
    v_psill  <- params$psill %||% NA
    v_model  <- params$model %||% "Sph"
    v_range  <- params$range %||% NA
    v_nugget <- params$nugget %||% NA
    
    # Set reasonable defaults if NA
    if (is.na(v_psill)) v_psill <- 0.8 * var(points_data$value, na.rm = TRUE)
    if (is.na(v_range)) v_range <- max(sample_variogram$dist) / 2
    if (is.na(v_nugget)) v_nugget <- 0.2 * var(points_data$value, na.rm = TRUE)
    
    flog.debug(glue("Variogram parameters - psill: {v_psill}, range: {v_range}, nugget: {v_nugget}, model: {v_model}"))
    
    variogram_model <- vgm(psill = v_psill, model = v_model, range = v_range, nugget = v_nugget, ...)
    
    # Fit variogram
    flog.debug("Fitting variogram model")
    fit_variogram <- tryCatch({
      fit <- fit.variogram(sample_variogram, variogram_model)
      flog.info(glue("Variogram fit successful - psill: {fit$psill[2]}, range: {fit$range[2]}, nugget: {fit$psill[1]}"))
      fit
    }, error = function(e) {
      flog.warn(glue("Variogram fit failed: {e$message}. Using initial model."))
      variogram_model
    })
    
    # Create gstat model
    flog.debug("Creating gstat model")
    gstat_model <- gstat(NULL, "var", var ~ 1, 
                         data = gstat_data, 
                         locations = ~x + y, model = fit_variogram)
    
    # Perform interpolation
    flog.debug("Starting kriging interpolation")
    num_cores <- max(1, detectCores() - 1)
    flog.debug(glue("Using {num_cores} cores for interpolation"))
    
    result <- tryCatch({
      terra::interpolate(grid_data, gstat_model, index = 1, cores = num_cores, 
                         cpkgs = c("terra", "gstat"))
    }, error = function(e) {
      flog.warn(glue("Parallel interpolation failed: {e$message}. Trying single core."))
      terra::interpolate(grid_data, gstat_model, index = 1, cores = 1)
    })
    
  } else if (method == "ked") {
    # Kriging with External Drift spatialization
    flog.info("Starting Kriging with External Drift")
    
    # Extract drift values
    drift_values <- as.data.frame(terra::extract(grid_data, points_data[, c("x", "y")])[, -1])
    names(drift_values) <- names(grid_data)
    
    if (ncol(drift_values) == 0) {
      flog.error("No drift variables extracted from grid")
      stop("No drift variables available for KED")
    }
    
    flog.debug(glue("Extracted {ncol(drift_values)} drift variables: {paste(names(drift_values), collapse = ', ')}"))
    
    # Prepare data
    ked_data <- data.frame(points_data, var = points_data$value, drift_values)
    model_formula <- as.formula(paste0("var ~ ", paste(names(drift_values), collapse = " + ")))
    
    flog.debug(glue("KED formula: {deparse(model_formula)}"))
    
    # Variogram calculation
    sample_variogram <- tryCatch({
      variogram(model_formula, ~x + y, data = ked_data)
    }, error = function(e) {
      flog.error(glue("Error calculating KED variogram: {e$message}"))
      stop(e)
    })
    
    # Variogram model parameters
    v_psill  <- params$psill %||% NA
    v_model  <- params$model %||% "Sph"
    v_range  <- params$range %||% NA
    v_nugget <- params$nugget %||% NA
    
    if (is.na(v_psill)) v_psill <- 0.8 * var(points_data$value, na.rm = TRUE)
    if (is.na(v_range)) v_range <- max(sample_variogram$dist) / 2
    if (is.na(v_nugget)) v_nugget <- 0.2 * var(points_data$value, na.rm = TRUE)
    
    variogram_model <- vgm(psill = v_psill, model = v_model, range = v_range, nugget = v_nugget, ...)
    
    # Fit variogram
    fit_variogram <- tryCatch({
      fit.variogram(sample_variogram, variogram_model)
    }, error = function(e) {
      flog.warn(glue("KED variogram fit failed: {e$message}. Using initial model."))
      variogram_model
    })
    
    # Create gstat model and interpolate
    gstat_model <- gstat(NULL, "var", model_formula, 
                         data = ked_data, 
                         locations = ~x + y, model = fit_variogram)
    
    num_cores <- max(1, detectCores() - 1)
    result <- tryCatch({
      terra::interpolate(grid_data, gstat_model, index = 1, cores = num_cores, 
                         cpkgs = c("terra", "gstat"))
    }, error = function(e) {
      flog.warn(glue("KED parallel interpolation failed: {e$message}. Trying single core."))
      terra::interpolate(grid_data, gstat_model, index = 1, cores = 1)
    })
    
  } else {
    flog.error(glue("Unsupported spatialization method: {method}"))
    stop("Unsupported spatialization method.")
  }
  
  # Apply mask to maintain original grid shape
  result <- mask(result, grid_data)
  
  # Check output
  result_values <- values(result)
  valid_cells <- sum(!is.na(result_values))
  total_cells <- ncell(result)
  if (inherits(result_values, "array")) {
    result_values <- as.numeric(result_values)
  }
  
  if (valid_cells == 0) {
    flog.warn(glue("Spatialization output is entirely NA: {valid_cells}/{total_cells} valid cells"))
    flog.debug(glue("Result extent: {as.character(ext(result))}"))
    flog.debug(glue("Result CRS: {crs(result)}"))
  } else {
    vals_clean <- na.omit(result_values)
    flog.info(glue("Raster value summary: Min={min(vals_clean)}, Mean={mean(vals_clean)}, Max={max(vals_clean)}"))
    flog.info(glue("Spatialization produced {valid_cells}/{total_cells} valid cells ({round(100 * valid_cells / total_cells, 1)}%)"))
  }
  
  return(result)
}

# Helper function for variogram debugging
debug_variogram <- function(points_data) {
  flog.info("=== VARIOGRAM DEBUGGING ===")
  
  gstat_data <- data.frame(x = points_data$x, y = points_data$y, var = points_data$value)
  
  sample_variogram <- tryCatch({
    v <- variogram(var ~ 1, ~x + y, data = gstat_data)
    flog.info(glue("Empirical variogram with {nrow(v)} distance bins"))
    flog.info(glue("Distance range: {min(v$dist)} to {max(v$dist)}"))
    flog.info(glue("Gamma range: {min(v$gamma)} to {max(v$gamma)}"))
    print(v)
    v
  }, error = function(e) {
    flog.error(glue("Variogram calculation failed: {e$message}"))
    NULL
  })
  
  return(sample_variogram)
}

# Test function with simple data
test_spatialization <- function() {
  flog.info("=== TESTING WITH SIMPLE DATA ===")
  
  test_points <- data.frame(
    x = c(0, 1, 0, 1),
    y = c(0, 0, 1, 1), 
    value = c(10, 20, 15, 25)
  )
  
  test_grid <- rast(ext(0, 1, 0, 1), res = 0.1)
  
  flog.info("Testing OK with simple data...")
  test_result <- spatialize(test_points, test_grid, method = "ok")
  valid_cells <- sum(!is.na(values(test_result)))
  flog.info(glue("Test result: {valid_cells} valid cells"))
  
  return(test_result)
}