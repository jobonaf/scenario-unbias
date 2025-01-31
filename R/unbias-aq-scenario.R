# External functions
source("R/spatialize-points-to-grid.R")
source("R/calibrate-unbias-coefficients.R")
source("R/apply-unbiasing.R")

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
  spatialization_method <- match.arg(spatialization_method)
  
  # Check calibration method restrictions for each unbias_sequence
  if (unbias_sequence == "SCA" && !calibration_method %in% c("Grid", "Cell", "Neigh")) {
    stop("For sequence 'SCA', calibration method must be one of 'Grid', 'Cell', or 'Neigh'.")
  }
  if (unbias_sequence %in% c("CSA", "CAS", "CA") && !calibration_method %in% c("All", "Each")) {
    stop("For sequences 'CSA', 'CAS', and 'CA', calibration method must be 'All' or 'Each'.")
  }
  # Check compatibility between calibration method and correction algorithm
  if (calibration_method %in% c("Each", "Cell") && !correction_algorithm %in% c("Add", "Mult")) {
    stop("Calibration methods 'Each' and 'Cell' are only compatible with correction algorithms 'Add' and 'Mult'.")
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
      correction_algorithm = correction_algorithm, 
      calibration_method = calibration_method)  
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
      coefficients = calibrated_coefficients, 
      correction_algorithm = correction_algorithm, 
      calibration_method = calibration_method)  
    return(corrected_data)
    
  } else if (unbias_sequence == "CAS") {
    # Calibrate coefficients (using observed data), apply correction, then spatialize the corrected data
    
    # Extract values from the 'scenario' based on the coordinates in 'observed_data'
    scenario_values <- terra::extract(scenario, observed_data[, c("x", "y")], xy = TRUE)
    scenario_sparse <- data.frame(observed_data[, c("x", "y")], value = scenario_values)
    
    # Apply correction to the sparse scenario
    calibrated_coefficients <- calibrate(
      obs = observed_data, 
      mod = base_case, 
      calibration_method = calibration_method, 
      correction_algorithm = correction_algorithm
    )
    # Apply correction
    corrected_sparse <- apply_correction(
      scenario = scenario, 
      coefficients = calibrated_coefficients, 
      correction_algorithm = correction_algorithm, 
      calibration_method = calibration_method)  
    
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
    # Apply correction
    corrected_data <- apply_correction(
      scenario = scenario, 
      coefficients = calibrated_coefficients, 
      correction_algorithm = correction_algorithm, 
      calibration_method = calibration_method)  
    return(corrected_data)
  }
}
