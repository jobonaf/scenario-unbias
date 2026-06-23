# =============================================================================
# file           :generate_cluster_jaccard_maps.R
# description    :Generate Jaccard cluster maps using precomputed CSV clusters
#                 and FAIRMODE scenario read via read_data()
# author         :Giovanni Bonafe'
# created        :2026-02-06
# version        :1.4
# dependencies   :terra, futile.logger, glue, RColorBrewer, readr, dplyr, ncdf4, rnaturalearth, sf
# notes          :No clustering recomputed, European domain, with boundaries
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(futile.logger)
  library(glue)
  library(RColorBrewer)
  library(readr)
  library(dplyr)
  library(ncdf4)
  library(scales)
  library(rnaturalearth)
  library(sf)
})

# =============================================================================
# Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript generate_cluster_jaccard_maps.R <INQUINANTE>")
}

specie <- args[1]

# Format pollutant name for display
specie_display <- ifelse(specie == "PM25", "PM2.5", specie)

# =============================================================================
# Paths (coerenti con gli script esistenti)
# =============================================================================

cluster_csv <- glue("data/clustering/clustering_{specie}.csv")
raster_dir <- "data/processed_fairmode" # BCM unbiased scenarios
output_dir  <- "output/maps-cluster-jaccard"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Logging
# =============================================================================

flog.threshold(INFO)
flog.info(glue("Generating Jaccard cluster maps for {specie}"))

# =============================================================================
# Load Europe boundaries
# =============================================================================

flog.info("Loading Europe boundaries")
europe_sf <- ne_countries(scale = "medium", continent = c("Europe","Africa"), 
                          returnclass = "sf")
# Convert to SpatVector for terra compatibility
europe <- vect(europe_sf)

# =============================================================================
# FAIRMODE data reader (as provided)
# =============================================================================

read_data <- function(parameter,
                      data_path="data/fairmode-wg5-exercise-202501/YEARLY") {
  
  obs_code <- case_when(
    parameter == "NO2" ~ "ug_NO2",
    parameter == "O3"  ~ "ppb_O3",
    parameter == "PM25" ~ "ug_PM25_rh50"
  )
  mod_code <- case_when(
    parameter == "PM25" ~ "PM25_rh50",
    .default = parameter
  )
  
  base_case_gridded  <- glue("{data_path}/BaseCase_Perturbed_Gridded/BaseCase_PERT_{mod_code}_YEARLY.nc")
  scenario_gridded   <- glue("{data_path}/Scenario_Perturbed_Gridded/SCEN_PERT_{mod_code}_YEARLY.nc")
  
  read_nc <- function(file) {
    nc_data <- nc_open(file)
    
    var_data <- ncvar_get(nc_data, nc_data$var[[1]]$name)
    lon <- ncvar_get(nc_data, "lon")
    lat <- ncvar_get(nc_data, "lat")
    
    res_x <- mean(diff(lon))
    res_y <- mean(diff(lat))
    
    ext_vals <- c(
      min(lon) - res_x / 2,
      max(lon) + res_x / 2,
      min(lat) - res_y / 2,
      max(lat) + res_y / 2
    )
    
    r <- rast(t(var_data)[ncol(var_data):1,], crs = "EPSG:4326")
    ext(r) <- ext_vals
    
    nc_close(nc_data)
    r
  }
  
  list(
    base_case = read_nc(base_case_gridded),
    scenario  = read_nc(scenario_gridded)
  )
}

# =============================================================================
# Helper functions
# =============================================================================

read_cluster_csv <- function(file) {
  if (!file.exists(file)) {
    flog.error(glue("Cluster CSV not found: {file}"))
    stop()
  }
  
  df <- read.csv(file, stringsAsFactors = FALSE)
  
  required <- c("member", "cluster_index", "is_medoid")
  if (!all(required %in% names(df))) {
    flog.error(glue("CSV missing required columns: {paste(required, collapse=', ')}"))
    stop()
  }
  
  df$cluster_index <- factor(df$cluster_index)
  df
}

compute_conc_breaks <- function(r_stack, r_ref, n_clas = 7) {
  # Collect all values from stack and reference
  all_vals <- c(values(r_stack, na.rm = TRUE),
                values(r_ref,   na.rm = TRUE))
  
  # Quantile-based breaks to avoid color recycling
  conc_pp <- c(0, 0.01, (1:(n_clas-3))/(n_clas-2), 0.99, 1)
  breaks <- unique(round(signif(quantile(all_vals, conc_pp, na.rm = TRUE), 2), 1))
  
  breaks
}

compute_diff_breaks <- function(all_diffs, n_clas_err = 6) {
  # Use absolute maximum error to create symmetric bins centered on zero
  err_pp <- c((0:(n_clas_err-1)/2)/((n_clas_err+1)/2), 0.99, 1)
  error_bins <- unique(round(signif(quantile(abs(all_diffs), err_pp, na.rm = TRUE), 2), 1))
  error_bins <- setdiff(unique(sort(c(-error_bins, error_bins))), 0)
  max_err <- ceiling(max(abs(all_diffs), na.rm = TRUE))
  
  error_bins[1] <- -max_err
  error_bins[length(error_bins)] <- max_err
  
  error_bins
}

plot_map <- function(r, breaks, cols, title, boundary = NULL) {
  plot(
    r,
    breaks = breaks,
    col = cols,
    axes = FALSE,
    legend = TRUE,
    mar = c(3, 1, 3, 1),
    main = title
  )
  
  # Add Europe boundaries if provided
  if (!is.null(boundary)) {
    plot(boundary, add = TRUE, border = "gray40", lwd = 0.5)
  }
}

# =============================================================================
# Read inputs
# =============================================================================

flog.info("Reading clustering CSV")
clu <- read_cluster_csv(cluster_csv)

flog.info("Loading BCM rasters")
bcm_files <- glue("{raster_dir}/{specie}_{clu$member}_unbiased_scenario.tif")

missing <- bcm_files[!file.exists(bcm_files)]

if (length(missing) > 0) {
  flog.error("Missing BCM rasters:")
  for (f in missing) flog.error(glue("  - {basename(f)}"))
  stop("Cannot proceed: missing BCM raster files")
}

r_stack <- rast(bcm_files)
names(r_stack) <- clu$member

# FAIRMODE scenario (reference)
flog.info("Reading FAIRMODE scenario via read_data()")
fairmode <- read_data(specie)
r_ref <- fairmode$scenario

# =============================================================================
# Compute global color scales
# =============================================================================

flog.info("Computing global color scales")

# Concentration breaks (quantile-based, non-recycling)
n_clas <- 7
conc_breaks <- compute_conc_breaks(r_stack, r_ref, n_clas = n_clas)
n_clas <- length(conc_breaks) - 1
conc_cols <- rev(brewer.pal(n_clas, "Spectral"))

# Standard deviation breaks (global across all clusters)
all_sd_vals <- NULL
for (cl in levels(clu$cluster_index)) {
  idx <- which(clu$cluster_index == cl)
  r_clu <- r_stack[[idx]]
  r_sd <- app(r_clu, sd, na.rm = TRUE)
  all_sd_vals <- c(all_sd_vals, values(r_sd, na.rm = TRUE))
}

sd_pp <- c(0, 0.01, (1:(n_clas-3))/(n_clas-2), 0.99, 1)
sd_breaks <- unique(round(signif(quantile(all_sd_vals, sd_pp, na.rm = TRUE), 2), 1))
n_clas_sd <- length(sd_breaks) - 1
sd_cols <- rev(brewer.pal(n_clas_sd, "YlOrRd"))

# Difference breaks (symmetric, centered on zero)
all_diffs <- NULL
for (cl in levels(clu$cluster_index)) {
  idx <- which(clu$cluster_index == cl)
  medoid_idx <- idx[clu$is_medoid[idx] %in% c("yes", "YES", TRUE)]
  if (length(medoid_idx) == 1) {
    medoid_name <- clu$member[medoid_idx]
    r_medoid <- r_stack[[medoid_name]]
    r_diff <- r_medoid - r_ref
    all_diffs <- c(all_diffs, values(r_diff, na.rm = TRUE))
  }
}

diff_breaks <- compute_diff_breaks(all_diffs, n_clas_err = 5)
n_clas_err <- length(diff_breaks) - 1
diff_cols <- rev(brewer.pal(max(3, n_clas_err), "RdBu"))

# =============================================================================
# Output
# =============================================================================

pdf_file <- glue("{output_dir}/cluster_Jaccard_{specie}.pdf")
pdf(pdf_file, width = 7, height = 6)

# =============================================================================
# Cluster loop
# =============================================================================

for (cl in levels(clu$cluster_index)) {
  
  flog.info(glue("Processing cluster {cl}"))
  
  idx <- which(clu$cluster_index == cl)
  r_clu <- r_stack[[idx]]
  
  r_mean   <- app(r_clu, mean,   na.rm = TRUE)
  r_sd     <- app(r_clu, sd,     na.rm = TRUE)
  r_median <- app(r_clu, median, na.rm = TRUE)
  
  medoid_idx <- idx[clu$is_medoid[idx] %in% c("yes", "YES", TRUE)]
  if (length(medoid_idx) != 1) {
    flog.error(glue("Cluster {cl}: zero or multiple medoids found"))
    stop()
  }
  
  medoid_name <- clu$member[medoid_idx]
  r_medoid <- r_stack[[medoid_name]]
  r_diff   <- r_medoid - r_ref
  
  # Plot maps with pollutant name in title and Europe boundaries
  plot_map(r_mean, conc_breaks, conc_cols,
           glue("{specie_display} - Cluster {cl} - Mean"),
           boundary = europe)
  
  plot_map(r_sd, sd_breaks, sd_cols,
           glue("{specie_display} - Cluster {cl} - Standard deviation"),
           boundary = europe)
  
  plot_map(r_median, conc_breaks, conc_cols,
           glue("{specie_display} - Cluster {cl} - Median"),
           boundary = europe)
  
  plot_map(r_medoid, conc_breaks, conc_cols,
           glue("{specie_display} - Cluster {cl} - Medoid ({medoid_name})"),
           boundary = europe)
  
  plot_map(r_diff, diff_breaks, diff_cols,
           glue("{specie_display} - Cluster {cl} - Medoid minus reference"),
           boundary = europe)
}

dev.off()

flog.info(glue("PDF written: {pdf_file}"))
flog.info("Done.")