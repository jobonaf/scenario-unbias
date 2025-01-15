# Main process function
process_data <- function(observed_data, base_case, scenario, 
                         unbias_sequence       = c("SCA", "CSA", "CAS", "CA"), 
                         calibration_method    = c("Point", "Grid", "Cell", "Neigh"), 
                         correction_algorithm  = c("Add", "Mult", "Lin"), 
                         spatialization_method = c("tps", "idw", "ok", "ked", "scm")) {
  
  # Validate inputs
  unbias_sequence       <- match.arg(unbias_sequence)
  correction_algorithm  <- match.arg(correction_algorithm)
  calibration_method    <- match.arg(calibration_method)
  spatialization_method <- match.arg(spatialization_method)
  
  # Check calibration method restrictions for each unbias_sequence
  if (unbias_sequence == "SCA" && calibration_method == "Point") {
    stop("For sequence 'SCA', calibration method 'Point' is not allowed.")
  }
  if (unbias_sequence %in% c("CSA", "CAS", "CA") && calibration_method != "Point") {
    stop("For sequences 'CSA', 'CAS', and 'CA', calibration method must be 'Point'.")
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
    corrected_data <- apply_correction(scenario, calibrated_coefficients, correction_algorithm)  # Apply correction
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
    corrected_data <- apply_correction(scenario, spatialized_coefficients, correction_algorithm)  # Apply correction
    return(corrected_data)
    
  } else if (unbias_sequence == "CAS") {
    # Calibrate coefficients (using observed data), apply correction, then spatialize the corrected data
    
    # Extract values from the 'scenario' based on the coordinates in 'observed_data'
    scenario_sparse <- terra::extract(scenario, observed_data[, c("x", "y")], xy = TRUE)
    
    # Apply correction to the sparse scenario
    calibrated_coefficients <- calibrate(
      obs = observed_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )
    corrected_sparse <- apply_correction(scenario_sparse, calibrated_coefficients, correction_algorithm)  # Apply correction
    
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
    corrected_data <- apply_correction(scenario, calibrated_coefficients, correction_algorithm)  # Apply correction
    return(corrected_data)
  }
}
