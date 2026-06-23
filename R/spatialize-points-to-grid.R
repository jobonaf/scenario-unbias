# =============================================================================
# title          :spatialize-points-to-grid
# description    :Spatial interpolation of point data to a regular grid using
#                 TPS, IDW, Ordinary Kriging (OK) and Kriging with External
#                 Drift (KED). Includes robust variogram fitting and automatic
#                 reprojection to planar coordinates (EPSG:3035) for geographic
#                 CRS to ensure numerical stability.
# author         :Giovanni Bonafe'
# date           :20250420
# version        :1.7
# notes          :Requires tidyverse, terra, fields, gstat, futile.logger,
#                 parallel, rlang, glue.
# R_version      :3.5.2
# =============================================================================

library(tidyverse)
library(terra)
library(fields)
library(gstat)
library(futile.logger)
library(parallel)
library(rlang)
library(glue)

# Set logging
flog.threshold(INFO)      # Adjust as needed (DEBUG for more details)
flog.appender(appender.console())

# ──────────────────────────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────────────────────────

# Safely detect number of available cores for parallel processing
safe_cores <- function(requested = NULL) {
  cores <- tryCatch({
    parallel::detectCores()
  }, error = function(e) NA_integer_)
  
  if (is.na(cores) || cores < 2) {
    flog.warn("Parallelization not available, using single core")
    return(1L)
  }
  
  cores <- max(1L, cores - 1L)
  
  ok <- tryCatch({
    cl <- parallel::makeCluster(1L)
    parallel::stopCluster(cl)
    TRUE
  }, error = function(e) FALSE)
  
  if (!ok) {
    flog.warn("Cluster creation failed, disabling parallelization")
    return(1L)
  }
  
  return(cores)
}

# Check whether a fitted variogram is degenerate (implausible parameters)
is_degenerate_variogram <- function(vgm_fit, max_range = NULL) {
  if (is.null(vgm_fit) || nrow(vgm_fit) < 2) return(TRUE)
  
  nugget    <- vgm_fit$psill[1]
  sill_comp <- vgm_fit$psill[2]
  sill      <- sum(vgm_fit$psill)
  range     <- vgm_fit$range[2]
  
  if (is.na(range)    || range     <= 0)        return(TRUE)
  if (is.na(sill)     || sill      <= 0)        return(TRUE)
  if (is.na(sill_comp)|| sill_comp <= 0)        return(TRUE)
  if (nugget >= 0.99 * sill)                    return(TRUE)
  
  # If max_range is provided, check against it
  if (!is.null(max_range) && range > max_range) return(TRUE)
  
  return(FALSE)
}

# ──────────────────────────────────────────────────────────────────────────────
# Variogram fitting (shared by OK and KED)
# ──────────────────────────────────────────────────────────────────────────────

# Robust variogram fitting with ordered strategies.
# Returns a fitted variogram model or a numerically stable fallback.
fit_variogram_model <- function(formula, gstat_data, values, params,
                                label    = "variogram",
                                max_range = NULL) {
  
  flog.debug(glue("[{label}] Starting variogram fitting with {nrow(gstat_data)} points"))
  
  # Compute adaptive max_range if not provided
  if (is.null(max_range)) {
    coords <- gstat_data[, c("x", "y")]
    diag   <- sqrt(diff(range(coords$x))^2 + diff(range(coords$y))^2)
    max_range <- diag / 3
    flog.debug(glue("[{label}] Adaptive max_range set to {round(max_range, 4)}"))
  }
  
  # Preliminary variogram (to determine cutoff)
  v0 <- tryCatch(
    variogram(formula, ~x + y, data = gstat_data),
    error = function(e) {
      flog.error(glue("[{label}] Preliminary variogram failed: {e$message}"))
      stop(e)
    }
  )
  
  if (nrow(gstat_data) < 10) {
    flog.error(glue("[{label}] Too few points, returning NULL"))
    return(NULL)
  }
  
  cutoff_val <- quantile(v0$dist, 0.9, na.rm = TRUE)
  width_val  <- cutoff_val / 15
  
  # Empirical variogram (Cressie robust estimator)
  sv <- tryCatch(
    variogram(formula, ~x + y, data = gstat_data,
              cressie = TRUE, cutoff = cutoff_val, width = width_val),
    error = function(e) {
      variogram(formula, ~x + y, data = gstat_data, cressie = TRUE, cutoff = cutoff_val)
    }
  )
  
  data_var <- var(values, na.rm = TRUE)
  
  # Ordered list of fitting strategies (from most to least informative)
  strategies <- list(
    list(psill = params$psill  %||% (0.8 * data_var),
         range = params$range  %||% median(sv$dist, na.rm = TRUE),
         nugget= params$nugget %||% (0.2 * data_var),
         model = params$model  %||% "Sph"),
    list(psill = data_var, range = quantile(sv$dist, 0.3), nugget = 0.1*data_var, model = "Sph"),
    list(psill = data_var, range = quantile(sv$dist, 0.5), nugget = 0.1*data_var, model = "Exp"),
    list(psill = data_var, range = quantile(sv$dist, 0.3), nugget = 0, model = "Sph"),
    list(psill = 0.5*data_var, range = quantile(sv$dist, 0.2), nugget = 0.05*data_var, model = "Mat")
  )
  
  for (i in seq_along(strategies)) {
    s <- strategies[[i]]
    vgm_init <- vgm(psill = s$psill, model = s$model, range = s$range, nugget = s$nugget)
    
    # Constrained fit (method 7: L-BFGS-B with bounds)
    fit <- tryCatch(
      fit.variogram(sv, vgm_init, fit.method = 7, fit.ranges = c(0.001, max_range)),
      error = function(e) NULL
    )
    
    # Fallback to unconstrained fit
    if (is.null(fit)) {
      fit <- tryCatch(
        fit.variogram(sv, vgm_init),
        error = function(e) NULL
      )
    }
    
    if (!is.null(fit) && !is_degenerate_variogram(fit, max_range)) {
      flog.info(glue("[{label}] Strategy {i} succeeded"))
      return(fit)
    }
  }
  
  # If all strategies fail, return a numerically stable exponential model
  flog.warn(glue("[{label}] Using manual stable exponential model"))
  vgm(psill = 0.8 * data_var, 
      model = "Exp", 
      range = max_range / 2,          
      nugget = 0.2 * data_var)
}

# ──────────────────────────────────────────────────────────────────────────────
# Kriging execution
# ──────────────────────────────────────────────────────────────────────────────

# Execute kriging using terra::interpolate with parallel fallback
run_kriging <- function(grid_data, gstat_model, label = "kriging") {
  num_cores <- safe_cores()
  flog.debug(glue("[{label}] Using {num_cores} cores"))
  
  tryCatch(
    terra::interpolate(grid_data, gstat_model, index = 1,
                       cores = num_cores, cpkgs = c("terra", "gstat")),
    error = function(e) {
      flog.warn(glue(
        "[{label}] Parallel interpolation failed: {e$message}. ",
        "Retrying single core."
      ))
      terra::interpolate(grid_data, gstat_model, index = 1, cores = 1L)
    }
  )
}

# Create an all-NA raster with the same geometry as template
na_raster <- function(template) {
  r <- rast(template)
  values(r) <- NA_real_
  r
}

# Clip extreme tails of a raster (0.1% and 99.9% quantiles) for numerical stability
# NOT USED ANYMORE
clip_data <- function(grid_data) {
  vals <- values(grid_data)
  
  q_low  <- quantile(vals, 0.001, na.rm = TRUE)
  q_high <- quantile(vals, 0.999, na.rm = TRUE)
  
  n_low  <- sum(vals < q_low,  na.rm = TRUE)
  n_high <- sum(vals > q_high, na.rm = TRUE)
  
  if (n_low + n_high > 0) {
    flog.debug(glue(
      "Clipping raster tails: ",
      "{n_low} below {signif(q_low,4)}, ",
      "{n_high} above {signif(q_high,4)}"
    ))
    vals[vals < q_low]  <- q_low
    vals[vals > q_high] <- q_high
    values(grid_data) <- vals
  }
  return(grid_data)
}

# ──────────────────────────────────────────────────────────────────────────────
# Main spatialization function
# ──────────────────────────────────────────────────────────────────────────────

spatialize <- function(points_data, grid_data,
                       method = c("tps", "idw", "ok", "ked"),
                       params = list(), ...) {
  
  method <- match.arg(method)
  
  # Input validation
  if (!inherits(grid_data, "SpatRaster")) {
    flog.error("'grid_data' must be a SpatRaster object.")
    stop("'grid_data' must be a SpatRaster object.")
  }
  if (!is.data.frame(points_data) ||
      !all(c("x", "y", "value") %in% colnames(points_data))) {
    flog.error("'points_data' must be a data frame with 'x', 'y' and 'value' columns.")
    stop("'points_data' must be a data frame with 'x', 'y' and 'value' columns.")
  }
  
  # Store original grid geometry for final masking
  grid_data_orig <- grid_data
  crs_orig <- crs(grid_data)
  
  # ----------------------------------------------------------------------------
  # Automatic reprojection to planar CRS (EPSG:3035) if data are in lon/lat
  # ----------------------------------------------------------------------------
  reprojected <- FALSE
  if (is.lonlat(grid_data)) {
    flog.debug("Reprojecting data from lon/lat to EPSG:3035 (LAEA Europe) for stable kriging")
    
    # Reproject raster grid
    grid_data <- project(grid_data, "EPSG:3035", method = "bilinear")
    
    # Reproject point coordinates
    pts_vect <- vect(points_data[, c("x", "y")], geom = c("x", "y"), crs = crs_orig)
    pts_proj <- project(pts_vect, "EPSG:3035")
    coords  <- crds(pts_proj)
    points_data$x <- coords[, 1]
    points_data$y <- coords[, 2]
    
    reprojected <- TRUE
  }
  
  empty_grid <- rast(grid_data)
  
  flog.info(glue("Starting spatialization with method '{method}'"))
  flog.info(glue("Points data: {nrow(points_data)} observations"))
  flog.info(glue("Grid dimensions: {nrow(grid_data)} x {ncol(grid_data)} cells"))
  flog.info(glue("Grid extent: {as.character(ext(grid_data))}"))
  flog.info(glue("Grid CRS: {crs(grid_data)}"))
  
  flog.debug(glue("NA values in points_data$value: {sum(is.na(points_data$value))}"))
  flog.debug(glue("Value range: [{min(points_data$value, na.rm=TRUE)}, {max(points_data$value, na.rm=TRUE)}]"))
  flog.debug(glue("X range: [{min(points_data$x)}, {max(points_data$x)}]"))
  flog.debug(glue("Y range: [{min(points_data$y)}, {max(points_data$y)}]"))
  
  points_ext <- ext(c(range(points_data$x), range(points_data$y)))
  grid_ext   <- ext(grid_data)
  if (!(points_ext[1] >= grid_ext[1] & points_ext[2] <= grid_ext[2] &
        points_ext[3] >= grid_ext[3] & points_ext[4] <= grid_ext[4])) {
    flog.warn(glue(
      "Points extent ({as.character(points_ext)}) does not fully overlap ",
      "with grid extent ({as.character(grid_ext)})"
    ))
  }
  
  # ----------------------------------------------------------------------------
  # Thin Plate Spline
  # ----------------------------------------------------------------------------
  if (method == "tps") {
    flog.debug("Fitting Thin Plate Spline model")
    tps_model <- Tps(points_data[, c("x", "y")], points_data$value, ...)
    result    <- interpolate(empty_grid, tps_model)
    
    # ----------------------------------------------------------------------------
    # Inverse Distance Weighting
    # ----------------------------------------------------------------------------
  } else if (method == "idw") {
    idp  <- params$idp  %||% 2
    nmax <- params$nmax %||% Inf
    flog.debug(glue("IDW: idp={idp}, nmax={nmax}"))
    
    gstat_model <- gstat(
      NULL, id = "var", formula = var ~ 1, locations = ~x + y,
      data = data.frame(points_data, var = points_data$value),
      nmax = nmax, set = list(idp = idp), ...
    )
    result <- interpolate(empty_grid, gstat_model, index = 1)
    
    # ----------------------------------------------------------------------------
    # Ordinary Kriging
    # ----------------------------------------------------------------------------
  } else if (method == "ok") {
    flog.info("Starting Ordinary Kriging")
    
    gstat_data <- data.frame(x   = points_data$x,
                             y   = points_data$y,
                             var = points_data$value)
    
    fit_vgm <- fit_variogram_model(
      formula    = var ~ 1,
      gstat_data = gstat_data,
      values     = points_data$value,
      params     = params,
      label      = "OK"
    )
    
    if (is.null(fit_vgm)) {
      flog.error("[OK] Variogram fitting failed entirely. Returning all-NA raster.")
      return(na_raster(empty_grid))
    }
    
    # Use nmax from params or a sensible default; maxdist is optional
    nmax_ok   <- params$nmax   %||% 200
    maxdist_ok <- params$maxdist %||% 2000000
    
    gstat_model <- gstat(NULL, "var", var ~ 1,
                         data      = gstat_data,
                         locations = ~x + y,
                         model     = fit_vgm,
                         nmax      = nmax_ok,
                         maxdist   = maxdist_ok)
    
    options(gstat.cn_max = 1e6)
    result <- run_kriging(grid_data, gstat_model, label = "OK")
    
    # ----------------------------------------------------------------------------
    # Kriging with External Drift
    # ----------------------------------------------------------------------------
  } else if (method == "ked") {
    flog.info("Starting Kriging with External Drift")
    
    drift_values <- as.data.frame(
      terra::extract(grid_data, points_data[, c("x", "y")])[, -1]
    )
    names(drift_values) <- names(grid_data)
    
    if (ncol(drift_values) == 0) {
      flog.error("No drift variables extracted from grid")
      stop("No drift variables available for KED")
    }
    flog.debug(glue(
      "Drift variables ({ncol(drift_values)}): ",
      "{paste(names(drift_values), collapse=', ')}"
    ))
    
    ked_data      <- data.frame(points_data, var = points_data$value, drift_values)
    model_formula <- as.formula(paste0("var ~ ", paste(names(drift_values), collapse = " + ")))
    flog.debug(glue("KED formula: {deparse(model_formula)}"))
    
    fit_vgm <- fit_variogram_model(
      formula    = model_formula,
      gstat_data = ked_data,
      values     = points_data$value,
      params     = params,
      label      = "KED"
    )
    
    if (is.null(fit_vgm)) {
      flog.error("[KED] Variogram fitting failed entirely. Returning all-NA raster.")
      return(na_raster(empty_grid))
    }
    
    nmax_ked   <- params$nmax   %||% 200
    maxdist_ked <- params$maxdist %||% 2000000
    
    gstat_model <- gstat(NULL, "var", model_formula,
                         data      = ked_data,
                         locations = ~x + y,
                         model     = fit_vgm,
                         nmax      = nmax_ked,
                         maxdist   = maxdist_ked)
    
    options(gstat.cn_max = 1e6)
    result <- run_kriging(grid_data, gstat_model, label = "KED")
  }
  
  # ----------------------------------------------------------------------------
  # Reproject result back to original CRS if needed
  # ----------------------------------------------------------------------------
  if (reprojected) {
    flog.debug("Reprojecting result back to original CRS and aligning to original grid")
    result <- project(result, grid_data_orig, method = "bilinear")
  }
  
  # Apply domain mask using the original grid geometry
  domain_mask        <- rast(grid_data_orig)
  values(domain_mask) <- 1L
  result             <- mask(result, domain_mask)
  
  # Log output statistics
  result_values <- as.numeric(values(result))
  valid_cells   <- sum(!is.na(result_values))
  total_cells   <- ncell(result)
  
  if (valid_cells == 0L) {
    flog.warn(glue(
      "Spatialization output is entirely NA: 0/{total_cells} valid cells | ",
      "method='{method}'"
    ))
  } else {
    vals_clean <- na.omit(result_values)
    flog.info(glue(
      "Output: {valid_cells}/{total_cells} valid cells ",
      "({round(100 * valid_cells / total_cells, 1)}%) | ",
      "min={round(min(vals_clean),4)}, ",
      "mean={round(mean(vals_clean),4)}, ",
      "max={round(max(vals_clean),4)}"
    ))
  }
  
  return(result)
}

# ──────────────────────────────────────────────────────────────────────────────
# Debugging and testing utilities (optional)
# ──────────────────────────────────────────────────────────────────────────────

debug_variogram <- function(points_data, formula = var ~ 1) {
  flog.info("=== VARIOGRAM DEBUGGING ===")
  
  gstat_data <- data.frame(x   = points_data$x,
                           y   = points_data$y,
                           var = points_data$value)
  
  tryCatch({
    v <- variogram(formula, ~x + y, data = gstat_data)
    flog.info(glue(
      "Bins: {nrow(v)} | dist: {min(v$dist)} to {max(v$dist)} | ",
      "gamma: {min(v$gamma)} to {max(v$gamma)}"
    ))
    print(v)
    v
  }, error = function(e) {
    flog.error(glue("Variogram calculation failed: {e$message}"))
    NULL
  })
}

test_spatialization <- function() {
  flog.info("=== TESTING SPATIALIZATION METHODS ===")
  
  test_grid <- rast(ext(-5, 5, 45, 55), res = 0.5)
  values(test_grid) <- 1
  
  set.seed(123)
  n <- 50
  pts <- data.frame(
    x = runif(n, -4, 4),
    y = runif(n, 46, 54)
  )
  pts$value <- with(pts, 2 + 0.1 * x + 0.05 * y + rnorm(n, 0, 0.5))
  
  methods <- c("tps", "idw", "ok", "ked")
  results <- list()
  
  for (m in methods) {
    flog.info(glue("Testing method: {m}"))
    res <- tryCatch(
      spatialize(pts, test_grid, method = m, params = list(idp = 2, nmax = 20)),
      error = function(e) {
        flog.error(glue("Method {m} failed: {e$message}"))
        return(NULL)
      }
    )
    if (!is.null(res)) {
      vc <- sum(!is.na(values(res)))
      flog.info(glue("  Valid cells: {vc} / {ncell(res)} ({round(100*vc/ncell(res),1)}%)"))
    }
    results[[m]] <- res
  }
  
  flog.info("=== TEST COMPLETED ===")
  invisible(results)
}