library(terra)

## SPLIT Point INTO PtGlo e PtLoc
################
calibrate <- function(obs, mod, 
                      calibration_method = c("Point", "Grid", "Cell", "Neigh"), 
                      correction_algorithm = c("Add", "Mult", "Lin"),
                      neigh_radius = 3) {
  
  # Match arguments
  calibration_method <- match.arg(calibration_method)
  correction_algorithm <- match.arg(correction_algorithm)
  
  # Validate inputs
  if (!inherits(mod, "SpatRaster")) {
    stop("'mod' must be a SpatRaster object.")
  }
  if (!inherits(obs, c("SpatRaster", "data.frame"))) {
    stop("'obs' must be either a SpatRaster or a data.frame.")
  }
  if (calibration_method %in% c("Grid", "Cell", "Neigh") && !inherits(obs, "SpatRaster")) {
    stop(paste0("For '", calibration_method, "' calibration, 'obs' must be a SpatRaster."))
  }
  if (calibration_method %in% c("Grid", "Cell", "Neigh") && !terra::compareGeom(obs, mod)) {
    stop("'obs' and 'mod' must have the same geometry for Grid, Cell, or Neigh calibration.")
  }
  
  # Calculate obs_values and mod_values
  if (calibration_method == "Point") {
    # Extract model values at observation points
    obs_values <- obs$value  # Assuming 'value' column contains observed values
    mod_values <- terra::extract(mod, obs[, c("x", "y")], xy = FALSE, ID=FALSE)[[1]]
  } else if (calibration_method == "Grid") {
    # Use raster values directly for grid-based calibration
    obs_values <- terra::values(obs)
    mod_values <- terra::values(mod)
  } else if (calibration_method == "Cell") {
    # Use cell-wise values for calibration
    obs_values <- obs
    mod_values <- mod
  } else if (calibration_method == "Neigh") {
    # Use neighborhood-based values for calibration
    weights <- function(dist) 1 / pmax(dist^2, 0.5^2)  # Avoid too small weights
    w <- terra::focalWeight(mod, neigh_radius, type = "circle")
    w <- weights(w)  # Apply inverse-square weighting with threshold
    
    obs_values <- terra::focal(obs, w, fun = weighted.mean, na.rm = TRUE)
    mod_values <- terra::focal(mod, w, fun = weighted.mean, na.rm = TRUE)
  }
  
  # Calculate correction coefficients
  if (correction_algorithm == "Add") {
    if (inherits(obs_values, "SpatRaster")) {
      coefficient <- obs_values - mod_values
    } else {
      coefficient <- mean(obs_values, na.rm = TRUE) - mean(mod_values, na.rm = TRUE)
    }
  } else if (correction_algorithm == "Mult") {
    if (inherits(obs_values, "SpatRaster")) {
      coefficient <- obs_values / mod_values
    } else {
      coefficient <- mean(obs_values, na.rm = TRUE) / mean(mod_values, na.rm = TRUE)
    }
  } else if (correction_algorithm == "Lin") {
    fit <- lm(obs_values ~ mod_values, na.action = na.exclude, 
              data = data.frame(obs_values, mod_values))
    coefficient <- list(intercept = coef(fit)[1], slope = coef(fit)[2])
  }
  
  return(coefficient)
}
