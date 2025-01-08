apply_correction <- function(scenario, coefficients, 
                             correction_algorithm = c("Add", "Mult", "Lin")) {
  # Validate correction algorithm
  correction_algorithm <- match.arg(correction_algorithm)
  
  # Check if 'scenario' is a SpatRaster or a data.frame
  if (inherits(scenario, "SpatRaster")) {
    # SpatRaster scenario
    scenario_type <- "raster"
  } else if (is.data.frame(scenario) && all(c("x", "y", "value") %in% colnames(scenario))) {
    # data.frame scenario
    scenario_type <- "data.frame"
  } else {
    stop("'scenario' must be either a SpatRaster object or a data.frame with columns 'x', 'y', 'value'.")
  }
  
  # Apply correction based on the algorithm
  if (correction_algorithm == "Add") {
    # Add coefficients (global or local)
    if (!is.numeric(coefficients) && !inherits(coefficients, "SpatRaster") && !is.numeric(coefficients[1])) {
      stop("For 'Add', coefficients must be numeric or a SpatRaster.")
    }
    
    if (scenario_type == "raster") {
      corrected_scenario <- scenario + coefficients
    } else if (scenario_type == "data.frame") {
      scenario$value <- scenario$value + coefficients
      corrected_scenario <- scenario
    }
    
  } else if (correction_algorithm == "Mult") {
    # Multiply by coefficients (global or local)
    if (!is.numeric(coefficients) && !inherits(coefficients, "SpatRaster") && !is.numeric(coefficients[1])) {
      stop("For 'Mult', coefficients must be numeric or a SpatRaster.")
    }
    
    if (scenario_type == "raster") {
      corrected_scenario <- scenario * coefficients
    } else if (scenario_type == "data.frame") {
      scenario$value <- scenario$value * coefficients
      corrected_scenario <- scenario
    }
    
  } else if (correction_algorithm == "Lin") {
    # Apply linear model (intercept and slope)
    if (!is.list(coefficients) || !all(c("intercept", "slope") %in% names(coefficients))) {
      stop("For 'Lin', coefficients must be a list with 'intercept' and 'slope'.")
    }
    
    if (scenario_type == "raster") {
      corrected_scenario <- coefficients$intercept + scenario * coefficients$slope
    } else if (scenario_type == "data.frame") {
      scenario$value <- coefficients$intercept + scenario$value * coefficients$slope
      corrected_scenario <- scenario
    }
    
  } else {
    stop("Unsupported correction algorithm.")
  }
  
  return(corrected_scenario)
}
