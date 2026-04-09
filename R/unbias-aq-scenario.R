# External functions
source("R/spatialize-points-to-grid.R")
source("R/calibrate-unbias-coefficients.R")
source("R/apply-unbiasing.R")

# Function to check if a given combination of unbias sequence, calibration method,
# and correction algorithm is valid according to the classification scheme described in
# https://doi.org/10.5281/zenodo.15188017
is_valid_combination <- function(unbias_sequence, calibration_method, correction_algorithm) {
  
  # SCA sequences spatialize first, so calibration must operate on gridded data
  # (Grid, Cell, or Neigh); point-based strategies (Each, All) are not permitted
  if (unbias_sequence == "SCA" && calibration_method %in% c("Each", "All")) {
    return(FALSE)
  }
  
  # CAS, CA, and CSA sequences calibrate at observation points, so only point-based
  # strategies (Each, All) are valid; grid-based strategies are not permitted
  if (unbias_sequence %in% c("CAS", "CA", "CSA") && !calibration_method %in% c("Each", "All")) {
    return(FALSE)
  }
  
  # Complex adjustment algorithms (Lin, Quant) require pooled data to fit their parameters
  # and cannot be calibrated at a single point or cell; only All or Grid are appropriate.
  # Simple algorithms (Add, Mult) are compatible with any calibration strategy
  if (!correction_algorithm %in% c("Add", "Mult") && calibration_method %in% c("Each", "Cell")) {
    return(FALSE)
  }
  
  # CSA with All calibration produces a scalar coefficient that cannot be spatialized;
  # use CA with All instead
  if (unbias_sequence == "CSA" && calibration_method == "All") {
    return(FALSE)
  }
  
  return(TRUE)
}

# Main process function
process_data <- function(observed_data, base_case, scenario, 
                         unbias_sequence       = c("SCA", "CSA", "CAS", "CA"), 
                         calibration_method    = c("All", "Each", "Grid", "Cell", "Neigh"), 
                         correction_algorithm  = c("Add", "Mult", "Lin"), 
                         spatialization_method = c("tps", "idw", "ok", "ked", "scm")) {
  
  # Validate inputs
  unbias_sequence       <- match.arg(unbias_sequence)
  correction_algorithm  <- match.arg(correction_algorithm)
  calibration_method    <- match.arg(calibration_method)
  if (unbias_sequence != "CA") {
    spatialization_method <- match.arg(spatialization_method)
  }

  # Check if the combination is valid
  if (!is_valid_combination(unbias_sequence, calibration_method, correction_algorithm)) {
    flog.error("Invalid combination for pollutant %s: %s.%s.%s%s", 
              pollutant, unbias_sequence, calibration_method, 
              correction_algorithm, 
              ifelse(is.null(spatialization_method), "", paste0(".", spatialization_method))
    )
    stop(1)
  }
  
  # Execute based on the chosen unbias_sequence
  if (unbias_sequence == "SCA") {
    # Spatialize the observed data (scattered points), calibrate coefficients, then apply correction
    spatialized_data <- spatialize(observed_data, scenario, spatialization_method)  # Spatialize observed data
    calibrated_coefficients <- calibrate(
      obs = spatialized_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )
    # Apply correction
    corrected_data <- apply_correction(
      scenario = scenario, 
      coefficients = calibrated_coefficients, 
      correction_algorithm = correction_algorithm)  
    return(corrected_data)
    
  } else if (unbias_sequence == "CSA") {
    # Calibrate coefficients (using observed data), spatialize the coefficients, then apply correction
    calibrated_coefficients <- calibrate(
      obs = observed_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )
    spatialized_coefficients <- spatialize(calibrated_coefficients, scenario, spatialization_method)  # Spatialize coefficients
    # Apply correction
    corrected_data <- apply_correction(
      scenario = scenario, 
      coefficients = spatialized_coefficients, 
      correction_algorithm = correction_algorithm)  
    return(corrected_data)
    
  } else if (unbias_sequence == "CAS") {
    # Calibrate coefficients (using observed data), apply correction, then spatialize the corrected data
    
    # Extract values from the 'scenario' based on the coordinates in 'observed_data'
    scenario_values <- terra::extract(scenario, observed_data[, c("x", "y")], xy = FALSE, ID=FALSE)
    scenario_sparse <- data.frame(observed_data[, c("x", "y")], value = unname(scenario_values))
    
    # Apply correction to the sparse scenario
    calibrated_coefficients <- calibrate(
      obs = observed_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )

    # Apply correction
    corrected_sparse <- apply_correction(
      scenario = scenario_sparse, 
      coefficients = calibrated_coefficients, 
      correction_algorithm = correction_algorithm)  
    
    # Spatialize the corrected sparse data
    spatialized_data <- spatialize(corrected_sparse, scenario, spatialization_method)  # Spatialize corrected data
    return(spatialized_data)
    
  } else if (unbias_sequence == "CA") {
    # Calibrate coefficients (using observed data) and apply correction
    calibrated_coefficients <- calibrate(
      obs = observed_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )

    # If 'calibrated_coefficients' is a data.frame object,
    # extract values from scenario before applying correction
    if (inherits(calibrated_coefficients, "data.frame")) {
      scenario_sparse <- terra::extract(scenario, calibrated_coefficients[, c("x", "y")], xy = FALSE, ID=FALSE)
      scenario <- data.frame(calibrated_coefficients[, c("x", "y")], value = unname(scenario_sparse))
    }

    # Apply correction
    corrected_data <- apply_correction(
      scenario = scenario, 
      coefficients = calibrated_coefficients, 
      correction_algorithm = correction_algorithm)  
    return(corrected_data)
  }
}
