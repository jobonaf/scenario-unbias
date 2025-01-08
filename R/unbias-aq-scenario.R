# Main process function
process_data <- function(observed_data, base_case, scenario, 
                         unbias_sequence      = c("ICA", "CIA", "CAI", "CA"), 
                         calibration_method   = c("Point", "Grid", "Cell", "Neigh"), 
                         correction_algorithm = c("Add", "Mult", "Lin"), 
                         interpolation_method = c("tps", "idw", "ok", "ked", "scm")) {
  
  # Validate inputs
  unbias_sequence      <- match.arg(unbias_sequence)
  correction_algorithm <- match.arg(correction_algorithm)
  calibration_method   <- match.arg(calibration_method)
  interpolation_method <- match.arg(interpolation_method)
  
  # Check calibration method restrictions for each unbias_sequence
  if (unbias_sequence == "ICA" && calibration_method == "Point") {
    stop("For sequence 'ICA', calibration method 'Point' is not allowed.")
  }
  if (unbias_sequence %in% c("CIA", "CAI", "CA") && calibration_method != "Point") {
    stop("For sequences 'CIA', 'CAI', and 'CA', calibration method must be 'Point'.")
  }
  
  # Execute based on the chosen unbias_sequence
  if (unbias_sequence == "ICA") {
    # Interpolate the observed data (scattered points), calibrate coefficients, then apply correction
    interpolated_data <- do_interpolation(observed_data, scenario, interpolation_method)  # Interpolate observed data
    calibrated_coefficients <- calibrate(
      obs = interpolated_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )
    corrected_data <- apply_correction(scenario, calibrated_coefficients, correction_algorithm)  # Apply correction
    return(corrected_data)
    
  } else if (unbias_sequence == "CIA") {
    # Calibrate coefficients (using observed data), interpolate the coefficients, then apply correction
    calibrated_coefficients <- calibrate(
      obs = observed_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )
    interpolated_coefficients <- do_interpolation(calibrated_coefficients, scenario, interpolation_method)  # Interpolate coefficients
    corrected_data <- apply_correction(scenario, interpolated_coefficients, correction_algorithm)  # Apply correction
    return(corrected_data)
    
  } else if (unbias_sequence == "CAI") {
    # Calibrate coefficients (using observed data), apply correction, then interpolate the corrected data
    calibrated_coefficients <- calibrate(
      obs = observed_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )
    corrected_sparse <- apply_correction(scenario, calibrated_coefficients, correction_algorithm)  # Apply correction
    interpolated_data <- do_interpolation(corrected_sparse, scenario, interpolation_method)  # Interpolate corrected data
    return(interpolated_data)
    
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
