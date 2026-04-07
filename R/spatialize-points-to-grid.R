library(tidyverse)
library(terra)
library(fields)
library(gstat)
library(futile.logger)
library(parallel)
library(rlang)
library(glue)

# Set logging
flog.threshold(INFO)
flog.appender(appender.console())

# ── Helpers ───────────────────────────────────────────────────────────────────

# Safely detect available cores
safe_cores <- function(requested = NULL) {
  cores <- tryCatch({
    parallel::detectCores()
  }, error = function(e) NA_integer_)
  
  if (is.na(cores) || cores < 2) {
    flog.warn("Parallelization not available, using single core")
    return(1)
  }
  
  cores <- max(1, cores - 1)
  
  ok <- tryCatch({
    cl <- parallel::makeCluster(1)
    parallel::stopCluster(cl)
    TRUE
  }, error = function(e) FALSE)
  
  if (!ok) {
    flog.warn("Cluster creation failed, disabling parallelization")
    return(1)
  }
  
  return(cores)
}

# Check whether a fitted variogram is degenerate
is_degenerate_variogram <- function(vgm_fit) {
  if (is.null(vgm_fit) || nrow(vgm_fit) < 2) return(TRUE)
  
  nugget <- vgm_fit$psill[1]
  sill   <- sum(vgm_fit$psill)
  range  <- vgm_fit$range[2]
  
  if (is.na(range) || range <= 0) return(TRUE)
  if (is.na(sill)  || sill  <= 0) return(TRUE)
  if (nugget >= 0.99 * sill)      return(TRUE)
  
  return(FALSE)
}

# ── Shared variogram fitting ──────────────────────────────────────────────────
#
# formula   : gstat formula (var ~ 1 for OK, var ~ drift1 + ... for KED)
# gstat_data: data frame with x, y, var (and drift columns for KED)
# values    : numeric vector of the target variable (for variance estimates)
# params    : list of optional overrides (psill, range, nugget, model)
# label     : short string used in log messages ("OK" or "KED")
#
# Returns a fitted (or fallback) variogram model ready for gstat().

fit_variogram_model <- function(formula, gstat_data, values, params,
                                label = "variogram") {
  
  # ── Step 1: preliminary variogram to derive a data-driven cutoff ────────────
  v0 <- tryCatch(
    variogram(formula, ~x + y, data = gstat_data),
    error = function(e) {
      flog.error(glue("[{label}] Preliminary variogram failed: {e$message}"))
      stop(e)
    }
  )
  
  cutoff_val <- quantile(v0$dist, 0.8, na.rm = TRUE)
  flog.debug(glue("[{label}] Variogram cutoff (80th pct of distances): {round(cutoff_val, 4)}"))
  
  # ── Step 2: robust empirical variogram ─────────────────────────────────────
  sv <- tryCatch(
    variogram(formula, ~x + y, data = gstat_data,
              cressie = TRUE, cutoff = cutoff_val),
    error = function(e) {
      flog.error(glue("[{label}] Empirical variogram failed: {e$message}"))
      stop(e)
    }
  )
  
  flog.debug(glue(
    "[{label}] Variogram: {nrow(sv)} bins | ",
    "dist [{round(min(sv$dist),3)}, {round(max(sv$dist),3)}] | ",
    "gamma [{round(min(sv$gamma),4)}, {round(max(sv$gamma),4)}]"
  ))
  
  data_var <- var(values, na.rm = TRUE)
  
  # ── Step 3: first fit attempt (user params or data-driven defaults) ─────────
  v_psill  <- params$psill  %||% (0.8 * data_var)
  v_range  <- params$range  %||% median(sv$dist, na.rm = TRUE)
  v_nugget <- params$nugget %||% (0.2 * data_var)
  v_model  <- params$model  %||% "Sph"
  
  flog.debug(glue(
    "[{label}] Initial model: psill={round(v_psill,4)}, ",
    "range={round(v_range,4)}, nugget={round(v_nugget,4)}, model={v_model}"
  ))
  
  vgm_init <- vgm(psill = v_psill, model = v_model,
                  range = v_range, nugget = v_nugget)
  
  fit <- tryCatch(
    fit.variogram(sv, vgm_init),
    error = function(e) {
      flog.warn(glue("[{label}] First fit failed: {e$message}"))
      NULL
    }
  )
  
  if (!is.null(fit) && !is_degenerate_variogram(fit)) {
    flog.info(glue(
      "[{label}] Fit success -> nugget={round(fit$psill[1],4)}, ",
      "psill={round(fit$psill[2],4)}, range={round(fit$range[2],4)}"
    ))
    return(fit)
  }
  
  flog.warn(glue("[{label}] First fit degenerate, trying fallback"))
  
  # ── Step 4: fallback — different range scale (30th percentile) ─────────────
  v_range_fb  <- quantile(sv$dist, 0.3, na.rm = TRUE)
  v_psill_fb  <- data_var
  v_nugget_fb <- 0.1 * data_var
  
  flog.debug(glue(
    "[{label}] Fallback model: psill={round(v_psill_fb,4)}, ",
    "range={round(v_range_fb,4)}, nugget={round(v_nugget_fb,4)}, model=Sph"
  ))
  
  vgm_fb <- vgm(psill = v_psill_fb, model = "Sph",
                range = v_range_fb, nugget = v_nugget_fb)
  
  fit_fb <- tryCatch(
    fit.variogram(sv, vgm_fb),
    error = function(e) {
      flog.warn(glue("[{label}] Fallback fit failed: {e$message}. Using initial fallback model."))
      vgm_fb
    }
  )
  
  flog.info(glue(
    "[{label}] Fallback fit -> nugget={round(fit_fb$psill[1],4)}, ",
    "psill={round(fit_fb$psill[2],4)}, range={round(fit_fb$range[2],4)}"
  ))
  
  return(fit_fb)
}

# ── Kriging interpolation with parallel/single-core fallback ─────────────────

run_kriging <- function(grid_data, gstat_model, label = "kriging") {
  num_cores <- safe_cores()
  flog.debug(glue("[{label}] Using {num_cores} cores"))
  
  tryCatch(
    terra::interpolate(grid_data, gstat_model, index = 1,
                       cores = num_cores, cpkgs = c("terra", "gstat")),
    error = function(e) {
      flog.warn(glue("[{label}] Parallel interpolation failed: {e$message}. Retrying single core."))
      terra::interpolate(grid_data, gstat_model, index = 1, cores = 1)
    }
  )
}

# ── Stabilize raster values (bilateral clipping) ─────────────────
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

# ── Main spatialization function ──────────────────────────────────────────────

spatialize <- function(points_data, grid_data,
                       method = c("tps", "idw", "ok", "ked"),
                       params = list(), ...) {
  
  # Validate input
  if (!inherits(grid_data, "SpatRaster")) {
    flog.error("'grid_data' must be a SpatRaster object.")
    stop("'grid_data' must be a SpatRaster object.")
  }
  if (!is.data.frame(points_data) ||
      !all(c("x", "y", "value") %in% colnames(points_data))) {
    flog.error("'points_data' must be a data frame with 'x', 'y' and 'value' columns.")
    stop("'points_data' must be a data frame with 'x', 'y' and 'value' columns.")
  }
  
  empty_grid <- rast(grid_data)
  
  # Log input data details
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
  
  flog.info(glue("Spatializing with method '{method}'"))
  
  # ── TPS ────────────────────────────────────────────────────────────────────
  if (method == "tps") {
    flog.debug("Fitting Thin Plate Spline model")
    tps_model <- Tps(points_data[, c("x", "y")], points_data$value, ...)
    result <- interpolate(empty_grid, tps_model)
    
    # ── IDW ────────────────────────────────────────────────────────────────────
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
    
    # ── OK ─────────────────────────────────────────────────────────────────────
  } else if (method == "ok") {
    flog.info("Starting Ordinary Kriging")
    
    gstat_data <- data.frame(x = points_data$x, y = points_data$y,
                             var = points_data$value)
    
    fit_vgm <- fit_variogram_model(
      formula    = var ~ 1,
      gstat_data = gstat_data,
      values     = points_data$value,
      params     = params,
      label      = "OK"
    )
    
    gstat_model <- gstat(NULL, "var", var ~ 1,
                         data = gstat_data, locations = ~x + y,
                         model = fit_vgm, nmax = 50)
    
    result <- run_kriging(clip_data(grid_data), gstat_model, label = "OK")
    
    # ── KED ────────────────────────────────────────────────────────────────────
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
    
    gstat_model <- gstat(NULL, "var", model_formula,
                         data = ked_data, locations = ~x + y,
                         model = fit_vgm, nmax = 50)
    
    result <- run_kriging(clip_data(grid_data), gstat_model, label = "KED")
    
  } else {
    flog.error(glue("Unsupported spatialization method: {method}"))
    stop("Unsupported spatialization method.")
  }
  
  # ── Output validation ───────────────────────────────────────────────────────
  result <- mask(result, grid_data)
  
  result_values <- as.numeric(values(result))
  valid_cells   <- sum(!is.na(result_values))
  total_cells   <- ncell(result)
  
  if (valid_cells == 0) {
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

# ── Debugging helpers ─────────────────────────────────────────────────────────

debug_variogram <- function(points_data, formula = var ~ 1) {
  flog.info("=== VARIOGRAM DEBUGGING ===")
  
  gstat_data <- data.frame(x = points_data$x, y = points_data$y,
                           var = points_data$value)
  
  tryCatch({
    v <- variogram(formula, ~x + y, data = gstat_data)
    flog.info(glue("Bins: {nrow(v)} | dist: {min(v$dist)} to {max(v$dist)} | gamma: {min(v$gamma)} to {max(v$gamma)}"))
    print(v)
    v
  }, error = function(e) {
    flog.error(glue("Variogram calculation failed: {e$message}"))
    NULL
  })
}

test_spatialization <- function() {
  flog.info("=== TESTING WITH SIMPLE DATA ===")
  
  test_points <- data.frame(
    x = c(0, 1, 0, 1),
    y = c(0, 0, 1, 1),
    value = c(10, 20, 15, 25)
  )
  test_grid <- rast(ext(0, 1, 0, 1), res = 0.1)
  
  flog.info("Testing OK...")
  result <- spatialize(test_points, test_grid, method = "ok")
  flog.info(glue("Valid cells: {sum(!is.na(values(result)))}"))
  
  return(result)
}