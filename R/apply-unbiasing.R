apply_correction <- function(scenario, coefficients, 
                             calibration_method, 
                             correction_algorithm = c("Add", "Mult", "Lin")) {
  # Validate correction algorithm
  correction_algorithm <- match.arg(correction_algorithm)
  
  # Determine scenario type
  if (inherits(scenario, "SpatRaster")) {
    scenario_type <- "raster"
  } else if (is.data.frame(scenario) && all(c("x", "y", "value") %in% colnames(scenario))) {
    scenario_type <- "data.frame"
  } else {
    stop("'scenario' must be either a SpatRaster object or a data.frame with columns 'x', 'y', 'value'.")
  }
  
  # Apply correction based on calibration method and algorithm
  if (calibration_method %in% c("Grid", "All")) {
    if (!is.numeric(coefficients) && !is.list(coefficients)) {
      stop("For 'Grid' or 'All', coefficients must be numeric or a list (for 'Lin').")
    }
    
    if (correction_algorithm == "Add") {
      corrected_scenario <- scenario + coefficients
    } else if (correction_algorithm == "Mult") {
      corrected_scenario <- scenario * coefficients
    } else if (correction_algorithm == "Lin") {
      corrected_scenario <- coefficients$intercept + scenario * coefficients$slope
    }
    
  } else if (calibration_method == "Each") {
    if (!is.data.frame(coefficients) || !all(c("x", "y", "value") %in% colnames(coefficients))) {
      stop("For 'Each', coefficients must be a data.frame with columns 'x', 'y', 'value'.")
    }
    
    scenario_values <- terra::extract(scenario, coefficients[, c("x", "y")], xy = FALSE, ID = FALSE)[[1]]
    
    if (correction_algorithm == "Add") {
      corrected_values <- scenario_values + coefficients$value
    } else if (correction_algorithm == "Mult") {
      corrected_values <- scenario_values * coefficients$value
    } else if (correction_algorithm == "Lin") {
      stop("Linear correction is not supported for 'Each'.")
    }
    
    corrected_scenario <- data.frame(x = coefficients$x, y = coefficients$y, value = corrected_values)
    
  } else if (calibration_method %in% c("Cell", "Neigh")) {
    if (!inherits(coefficients, "SpatRaster")) {
      stop("For 'Cell' or 'Neigh', coefficients must be a SpatRaster.")
    }
    
    if (correction_algorithm == "Add") {
      corrected_scenario <- scenario + coefficients
    } else if (correction_algorithm == "Mult") {
      corrected_scenario <- scenario * coefficients
    } else if (correction_algorithm == "Lin") {
      corrected_scenario <- coefficients$intercept + scenario * coefficients$slope
    }
  } else {
    stop("Unsupported calibration method.")
  }
  
  return(corrected_scenario)
}
