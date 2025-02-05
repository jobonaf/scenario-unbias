library(futile.logger)

apply_correction <- function(scenario, coefficients, correction_algorithm = c("Add", "Mult", "Lin")) {
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
  
  # Determine coefficient type
  if (inherits(coefficients, "SpatRaster")) {
    coeff_type <- "raster"
  } else if (is.data.frame(coefficients) && all(c("x", "y", "value") %in% colnames(coefficients))) {
    coeff_type <- "data.frame"
  } else if (is.list(coefficients) && all(c("intercept", "slope") %in% names(coefficients))) {
    coeff_type <- "lin"
  } else if (is.numeric(coefficients)) {
    coeff_type <- "numeric"
  } else {
    stop("Unsupported coefficients format.")
  }
  flog.info("Applying correction algorithm '%s' on scenario of type '%s' with coefficients of type '%s'", 
            correction_algorithm, scenario_type, coeff_type)
  
  # Apply correction based on scenario type and coefficient type
  if (scenario_type == "raster") {
    if (coeff_type == "numeric") {
      if (correction_algorithm == "Add") {
        corrected_scenario <- scenario + coefficients
      } else if (correction_algorithm == "Mult") {
        corrected_scenario <- scenario * coefficients
      } else {
        stop("Linear correction requires a list of coefficients with 'intercept' and 'slope'.")
      }
    } else if (coeff_type == "lin") {
      corrected_scenario <- coefficients$intercept + scenario * coefficients$slope
    } else if (coeff_type == "raster") {
      if (correction_algorithm == "Add") {
        corrected_scenario <- scenario + coefficients
      } else if (correction_algorithm == "Mult") {
        corrected_scenario <- scenario * coefficients
      } else if (correction_algorithm == "Lin") {
        if (nlyr(coefficients) < 2) {
          stop("For 'Lin' correction, raster coefficients must have two layers: intercept and slope.")
        }
        corrected_scenario <- coefficients[[1]] + scenario * coefficients[[2]]
      }
    } else if (coeff_type == "data.frame") {
      stop("Point-wise coefficients (data.frame) cannot be directly applied to a raster scenario.")
    }
    names(corrected_scenario) <- "value"
  } else if (scenario_type == "data.frame") {
    if (coeff_type == "numeric") {
      if (correction_algorithm == "Add") {
        corrected_scenario <- scenario
        corrected_scenario$value <- scenario$value + coefficients
      } else if (correction_algorithm == "Mult") {
        corrected_scenario <- scenario
        corrected_scenario$value <- scenario$value * coefficients
      } else {
        stop("Linear correction requires a list of coefficients with 'intercept' and 'slope'.")
      }
    } else if (coeff_type == "lin") {
      corrected_scenario <- scenario
      corrected_scenario$value <- coefficients$intercept + scenario$value * coefficients$slope
    } else if (coeff_type == "data.frame") {
      if (!all(scenario[, c("x", "y")] == coefficients[, c("x", "y")])) {
        stop("Scenario and coefficients must have matching x, y locations.")
      }
      if (correction_algorithm == "Add") {
        corrected_scenario <- scenario
        corrected_scenario$value <- scenario$value + coefficients$value
      } else if (correction_algorithm == "Mult") {
        corrected_scenario <- scenario
        corrected_scenario$value <- scenario$value * coefficients$value
      } else {
        stop("Linear correction is not supported for point-wise coefficients.")
      }
    } else if (coeff_type == "raster") {
      stop("Raster coefficients cannot be applied directly to a data.frame scenario.")
    }
  }
  
  return(corrected_scenario)
}
