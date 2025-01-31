library(terra)

calibrate <- function(obs, mod, 
                      calibration_method = c("Each", "All", "Grid", "Cell", "Neigh"), 
                      correction_algorithm = c("Add", "Mult", "Lin"),
                      neigh_radius = 3) {
  
  # Match arguments
  calibration_method <- match.arg(calibration_method)
  correction_algorithm <- match.arg(correction_algorithm)
  
  # Validate inputs
  if (!inherits(mod, "SpatRaster")) {
    stop("'mod' must be a SpatRaster object.")
  }
  if (calibration_method %in% c("Each", "All") && (!inherits(obs, "data.frame") || 
                                                   !all(c("x", "y", "value") %in% colnames(obs)))) {
    stop("For 'Each' and 'All' calibration, 'obs' must be a data.frame with 'x', 'y', and 'value' columns.")
  }
  if (calibration_method %in% c("Grid", "Cell", "Neigh") && !inherits(obs, "SpatRaster")) {
    stop(paste0("For '", calibration_method, "' calibration, 'obs' must be a SpatRaster."))
  }
  if (calibration_method %in% c("Grid", "Cell", "Neigh") && !terra::compareGeom(obs, mod)) {
    stop("'obs' and 'mod' must have the same geometry for Grid, Cell, or Neigh calibration.")
  }
  
  # Calculate obs_values and mod_values
  if (calibration_method %in% c("Each", "All")) {
    obs_values <- obs$value
    mod_values <- terra::extract(mod, obs[, c("x", "y")], xy = FALSE, ID = FALSE)[[1]]
  } else if (calibration_method == "Grid") {
    obs_values <- terra::values(obs)
    mod_values <- terra::values(mod)
  } else if (calibration_method == "Cell") {
    obs_values <- obs
    mod_values <- mod
  } else if (calibration_method == "Neigh") {
    weights <- function(dist) 1 / pmax(dist^2, 0.5^2)  # Avoid too small weights
    w <- terra::focalWeight(mod, neigh_radius, type = "circle")
    w <- weights(w)  # Apply inverse-square weighting with threshold
    
    obs_values <- terra::focal(obs, w, fun = weighted.mean, na.rm = TRUE)
    mod_values <- terra::focal(mod, w, fun = weighted.mean, na.rm = TRUE)
  }
  
  # Calculate correction coefficients
  if (calibration_method == "All") {
    if (correction_algorithm == "Add") {
      coefficient <- mean(obs_values, na.rm = TRUE) - mean(mod_values, na.rm = TRUE)
    } else if (correction_algorithm == "Mult") {
      coefficient <- mean(obs_values, na.rm = TRUE) / mean(mod_values, na.rm = TRUE)
    } else if (correction_algorithm == "Lin") {
      fit <- lm(obs_values ~ mod_values, na.action = na.exclude, 
                data = data.frame(obs_values, mod_values))
      coefficient <- list(intercept = coef(fit)[1], slope = coef(fit)[2])
    }
  } else if (calibration_method == "Each") {
    coefficient <- data.frame(x = obs$x, y = obs$y)
    if (correction_algorithm == "Add") {
      coefficient$value <- obs_values - mod_values
    } else if (correction_algorithm == "Mult") {
      coefficient$value <- obs_values / mod_values
    } else if (correction_algorithm == "Lin") {
      stop("'Lin' correction is not supported for 'Each' calibration.")
    }
  } else if (calibration_method %in% c("Grid", "Cell", "Neigh")) {
    if (correction_algorithm == "Add") {
      coefficient <- obs_values - mod_values
    } else if (correction_algorithm == "Mult") {
      coefficient <- obs_values / mod_values
    } else if (correction_algorithm == "Lin") {
      fit <- lm(obs_values ~ mod_values, na.action = na.exclude, 
                data = data.frame(obs_values, mod_values))
      coefficient <- list(intercept = coef(fit)[1], slope = coef(fit)[2])
    }
  }
  
  return(coefficient)
}
