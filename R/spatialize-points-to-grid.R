library(tidyverse)
library(terra)
library(fields)
library(gstat)
library(sf)

# Spatialization function (updated)
spatialize <- function(points_data, grid_data, 
                              method = c("tps", "idw", "ok", "ked"), 
                              params = list(), ...) {
  # Validate input
  if (!inherits(grid_data, "SpatRaster")) stop("'grid_data' must be a SpatRaster object.")
  if (!is.data.frame(points_data) || !all(c("x", "y") %in% colnames(points_data))) stop("'points_data' must be a data frame with 'x' and 'y' columns.")
  
  # Prepare output raster (same grid as grid_data)
  empty_grid <- rast(grid_data)
  result <- empty_grid  # Start with an empty raster grid
  
  # Select spatialization method
  if (method == "tps") {
    # Thin Plate Spline spatialization
    tps_model <- Tps(points_data[, c("x", "y")], points_data$values, ...)
    result <- interpolate(empty_grid, tps_model)
    
  } else if (method == "idw") {
    # Inverse Distance Weighted spatialization
    idp  <- params$idp %||% 2     # Default power parameter for IDW
    nmax <- params$nmax %||% Inf  # Default nmax
    gstat_model <- gstat(NULL, id = "var", formula = var ~ 1, locations = ~x + y, 
                         data = data.frame(points_data, var = points_data$values), 
                         nmax = nmax, 
                         set = list(idp = idp), ...)
    result <- interpolate(empty_grid, gstat_model, index = 1)
    
  } else if (method == "ok") {
    # Ordinary Kriging spatialization
    v_psill  <- params$psill  %||% NA    # Default sill of the variogram model component
    v_model  <- params$model  %||% "Sph" # Default model type
    v_range  <- params$range  %||% NA    # Default range parameter of the variogram model component
    v_nugget <- params$nugget %||% NA    # Default nugget component of the variogram
    sample_variogram <- variogram(var ~ 1, ~x + y, 
                                  data = data.frame(points_data, var = points_data$values))
    variogram_model <- vgm(psill = v_psill, model = v_model, range = v_range, nugget = v_nugget, ...)
    fit_variogram <- fit.variogram(sample_variogram, variogram_model)
    gstat_model <- gstat(NULL, "var", var ~ 1, 
                         data = data.frame(points_data, var = points_data$values), 
                         locations = ~x + y, model = fit_variogram)
    result <- interpolate(grid_data, gstat_model, index = 1)
    
  } else if (method == "ked") {
    # Kriging with External Drift spatialization
    v_psill  <- params$psill  %||% NA    # Default sill of the variogram model component
    v_model  <- params$model  %||% "Sph" # Default model type
    v_range  <- params$range  %||% NA    # Default range parameter of the variogram model component
    v_nugget <- params$nugget %||% NA    # Default nugget component of the variogram
    drift_values <- as.data.frame(extract(grid_data, points_data[, c("x", "y")])[, -1])
    names(drift_values) <- names(grid_data)
    model_formula <- as.formula(paste0("var ~ ", paste(names(drift_values), collapse = " + ")))
    sample_variogram <- variogram(model_formula, ~x + y, 
                                  data = data.frame(points_data, var = points_data$values, drift_values))
    variogram_model <- vgm(psill = v_psill, model = v_model, range = v_range, nugget = v_nugget, ...)
    fit_variogram <- fit.variogram(sample_variogram, variogram_model)
    gstat_model <- gstat(NULL, "var", model_formula, 
                         data = data.frame(points_data, var = points_data$values, drift_values), 
                         locations = ~x + y, model = fit_variogram)
    result <- interpolate(grid_data, gstat_model, index = 1)
    
  } else {
    stop("Unsupported spatialization method.")
  }
  
  # Apply mask to maintain original grid shape
  result <- mask(result, grid_data)
  
  return(result)
}
