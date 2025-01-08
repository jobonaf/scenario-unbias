apply_correction <- function(scenario, coefficients, 
                              correction_algorithm = c("Add", "Mult", "Lin")) {
  # Validate correction algorithm
  correction_algorithm <- match.arg(correction_algorithm)
  
  # Check inputs
  if (!inherits(scenario, "SpatRaster")) {
    stop("'scenario' must be a SpatRaster object.")
  }
  
  # Apply correction based on the algorithm
  if (correction_algorithm == "Add") {
    # Add coefficients (global or local)
    if (!is.numeric(coefficients) || !inherits(coefficients, "SpatRaster")) {
      stop("For 'Add', coefficients must be numeric or a SpatRaster.")
    }
    corrected_scenario <- scenario + coefficients
    
  } else if (correction_algorithm == "Mult") {
    # Multiply by coefficients (global or local)
    if (!is.numeric(coefficients) || !inherits(coefficients, "SpatRaster")) {
      stop("For 'Mult', coefficients must be numeric or a SpatRaster.")
    }
    corrected_scenario <- scenario * coefficients
    
  } else if (correction_algorithm == "Lin") {
    # Apply linear model (intercept and slope)
    if (!is.list(coefficients) || !all(c("intercept", "slope") %in% names(coefficients))) {
      stop("For 'Lin', coefficients must be a list with 'intercept' and 'slope'.")
    }
    corrected_scenario <- coefficients$intercept + scenario * coefficients$slope
    
  } else {
    stop("Unsupported correction algorithm.")
  }
  
  return(corrected_scenario)
}
