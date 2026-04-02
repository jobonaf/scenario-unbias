# make_boxplots.R
#
# Reads all GeoTIFF files produced by the scenario-unbias pipeline
# and extracts global raster statistics.
# Outputs a CSV including also invalid rasters.

library(terra)
library(dplyr)
library(stringr)
library(glue)

# ── Configuration ─────────────────────────────────────────────────────────────

INPUT_DIR  <- "data/processed_phase2"
OUTPUT_CSV <- "output/fairmode_phase2/boxplots_pollutant_year_stats.csv"

# ── Filename parser ───────────────────────────────────────────────────────────

parse_tif_filename <- function(path) {
  fname <- basename(path)
  stem  <- sub("\\.tif$", "", fname)
  
  m <- str_match(
    stem,
    "^([^_]+)_([A-Z]+)\\.([^\\.]+)\\.([^\\.]+)(?:\\.([^_]+))?_unbiased_scenario_(\\d{4})$"
  )
  
  if (is.na(m[1, 1])) {
    warning("Unrecognised filename format: ", fname)
    return(NULL)
  }
  
  data.frame(
    path                  = path,
    filename              = fname,
    pollutant             = m[1, 2],
    unbias_sequence       = m[1, 3],
    calibration_method    = m[1, 4],
    correction_algorithm  = m[1, 5],
    spatialization_method = ifelse(is.na(m[1, 6]), NA_character_, m[1, 6]),
    scenario_year         = as.integer(m[1, 7]),
    stringsAsFactors      = FALSE
  )
}

# ── Collect file metadata ─────────────────────────────────────────────────────

tif_files <- list.files(INPUT_DIR, pattern = "\\.tif$", full.names = TRUE)
cat(sprintf("GeoTIFF files found    : %d\n", length(tif_files)))

file_meta <- bind_rows(lapply(tif_files, parse_tif_filename))
cat(sprintf("Files parsed           : %d\n", nrow(file_meta)))

# ── Define expected combinations (robust) ─────────────────────────────────────

# combinazioni di metodi (indipendenti da pollutant/year)
method_combos <- file_meta %>%
  distinct(unbias_sequence,
           calibration_method,
           correction_algorithm,
           spatialization_method)

# dimensioni complete di pollutant e year
pollutant_year <- file_meta %>%
  distinct(pollutant, scenario_year)

# prodotto cartesiano
expected_grid <- merge(method_combos, pollutant_year)

# ── Identify missing combinations (no TIFF) ───────────────────────────────────

missing_meta <- expected_grid %>%
  left_join(file_meta,
            by = c("pollutant",
                   "unbias_sequence",
                   "calibration_method",
                   "correction_algorithm",
                   "spatialization_method",
                   "scenario_year")) %>%
  filter(is.na(path)) %>%
  mutate(
    path = NA_character_,
    filename = NA_character_
  )

# ── Extract global raster statistics ─────────────────────────────────────────

extract_global_stats <- function(row) {
  
  file_info <- tryCatch(
    file.info(row$path),
    error = function(e) NULL
  )
  
  r <- tryCatch(
    rast(row$path),
    error = function(e) {
      warning("Cannot read raster: ", row$path, " — ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(r)) {
    return(data.frame(
      row,
      valid_raster = FALSE,
      error_type   = "read_error",
      file_size = if (!is.null(file_info)) file_info$size else NA_real_,
      file_mtime = if (!is.null(file_info)) format(file_info$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_,
      stat_min    = NA_real_,
      stat_q25    = NA_real_,
      stat_median = NA_real_,
      stat_mean   = NA_real_,
      stat_q75    = NA_real_,
      stat_max    = NA_real_,
      stat_sd     = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  pixel_values <- tryCatch(
    values(r),
    error = function(e) {
      warning("Cannot extract values: ", row$path, " — ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(pixel_values)) {
    return(data.frame(
      row,
      valid_raster = FALSE,
      error_type   = "values_error",
      file_size = if (!is.null(file_info)) file_info$size else NA_real_,
      file_mtime = if (!is.null(file_info)) format(file_info$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_,
      stat_min    = NA_real_,
      stat_q25    = NA_real_,
      stat_median = NA_real_,
      stat_mean   = NA_real_,
      stat_q75    = NA_real_,
      stat_max    = NA_real_,
      stat_sd     = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  pixel_values <- as.vector(pixel_values)
  pixel_values <- pixel_values[is.finite(pixel_values)]
  
  if (length(pixel_values) == 0) {
    warning("No finite pixel values in: ", row$path)
    
    return(data.frame(
      row,
      valid_raster = FALSE,
      error_type   = "no_finite_values",
      file_size = if (!is.null(file_info)) file_info$size else NA_real_,
      file_mtime = if (!is.null(file_info)) format(file_info$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_,
      stat_min    = NA_real_,
      stat_q25    = NA_real_,
      stat_median = NA_real_,
      stat_mean   = NA_real_,
      stat_q75    = NA_real_,
      stat_max    = NA_real_,
      stat_sd     = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    row,
    valid_raster = TRUE,
    error_type   = NA_character_,
    file_size = if (!is.null(file_info)) file_info$size else NA_real_,
    file_mtime = if (!is.null(file_info)) format(file_info$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_,
    stat_min    = min(pixel_values),
    stat_q25    = quantile(pixel_values, 0.25),
    stat_median = median(pixel_values),
    stat_mean   = mean(pixel_values),
    stat_q75    = quantile(pixel_values, 0.75),
    stat_max    = max(pixel_values),
    stat_sd     = sd(pixel_values),
    stringsAsFactors = FALSE
  )
}

cat("Extracting global raster statistics...\n")

stats_list <- vector("list", nrow(file_meta))
pb <- txtProgressBar(min = 0, max = nrow(file_meta), style = 3)

for (i in seq_len(nrow(file_meta))) {
  stats_list[[i]] <- extract_global_stats(file_meta[i, ])
  setTxtProgressBar(pb, i)
}

close(pb)

stats_df <- bind_rows(stats_list)

cat(sprintf("Rasters processed      : %d\n", nrow(stats_df)))
cat(sprintf("Valid rasters          : %d\n", sum(stats_df$valid_raster, na.rm = TRUE)))
cat(sprintf("Invalid rasters        : %d\n", sum(!stats_df$valid_raster, na.rm = TRUE)))

# ── Add missing TIFF combinations to stats ────────────────────────────────────

missing_stats <- missing_meta %>%
  mutate(
    valid_raster = FALSE,
    error_type   = "missing_tif",
    file_size    = NA_real_,
    file_mtime   = NA_character_,
    stat_min    = NA_real_,
    stat_q25    = NA_real_,
    stat_median = NA_real_,
    stat_mean   = NA_real_,
    stat_q75    = NA_real_,
    stat_max    = NA_real_,
    stat_sd     = NA_real_
  )

stats_df <- bind_rows(stats_df, missing_stats)
cat(sprintf("Missing rasters        : %d\n", sum(is.na(stats_df$filename))))

# ── Save CSV ──────────────────────────────────────────────────────────────────

dir.create(dirname(OUTPUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(stats_df, OUTPUT_CSV, row.names = FALSE)

cat(sprintf("Statistics CSV written : %s\n", OUTPUT_CSV))

stats_df %>%
  filter(!valid_raster) %>%
  group_by(scenario_year, pollutant, 
           method = glue("{unbias_sequence}.{calibration_method}.{correction_algorithm}")) %>%
  summarize(spatialization_method = paste0(spatialization_method, collapse=",")) %>%
  knitr::kable()
