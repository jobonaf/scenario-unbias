# explore_tiff_data.R
#
# Reads all GeoTIFF files produced by the scenario-unbias pipeline
# and extracts global raster statistics.
# Outputs a CSV including also invalid rasters.

library(terra)
library(dplyr)
library(glue)

# ── Configuration ─────────────────────────────────────────────────────────────

INPUT_DIR  <- "data/processed_phase2/"
OUTPUT_CSV <- "output/fairmode_phase2/output_stats.csv"

# ── Filename schemas ──────────────────────────────────────────────────────────
#
# Each schema is a named list with:
#   name    : human-readable label (used in warnings)
#   pattern : regex with named capture groups
#   defaults: named list of fixed values for fields absent from the pattern
#
# Named capture groups map directly to output columns.
# Every schema must produce (or default) all of:
#   pollutant, unbias_sequence, calibration_method, correction_algorithm,
#   spatialization_method, scenario_year, run_type
#
# Schemas are tried in order; the first match wins.

FILENAME_SCHEMAS <- list(
  
  # Schema 1 — with year
  # e.g. NO2_CA.All.Add.IDW_unbiased_scenario_2030.tif
  #      NO2_CA.All.Add_unbiased_basecase_2030.tif
  list(
    name    = "with_year",
    pattern = paste0(
      "^(?<pollutant>[^_]+)",
      "_(?<unbias_sequence>[A-Z]+)",
      "\\.(?<calibration_method>[^\\.]+)",
      "\\.(?<correction_algorithm>[^\\.]+)",
      "(?:\\.(?<spatialization_method>[^_]+))?",
      "_unbiased_(?<run_type>scenario|basecase)",
      "_(?<scenario_year>\\d{4})$"
    ),
    defaults = list()
  ),
  
  # Schema 2 — without year
  # e.g. NO2_CA.All.Add.IDW_unbiased_scenario.tif
  #      NO2_CA.All.Add_unbiased_basecase.tif
  list(
    name    = "without_year",
    pattern = paste0(
      "^(?<pollutant>[^_]+)",
      "_(?<unbias_sequence>[A-Z]+)",
      "\\.(?<calibration_method>[^\\.]+)",
      "\\.(?<correction_algorithm>[^\\.]+)",
      "(?:\\.(?<spatialization_method>[^_]+))?",
      "_unbiased_(?<run_type>scenario|basecase)$"
    ),
    defaults = list(scenario_year = NA_integer_)
  )
  
)

# ── Filename parser ───────────────────────────────────────────────────────────

# All output columns, in canonical order.
CANONICAL_FIELDS <- c(
  "path", "filename",
  "pollutant", "unbias_sequence", "calibration_method",
  "correction_algorithm", "spatialization_method",
  "scenario_year", "run_type"
)

parse_tif_filename <- function(path) {
  fname <- basename(path)
  stem  <- sub("\\.tif$", "", fname)
  
  for (schema in FILENAME_SCHEMAS) {
    rx <- regexpr(schema$pattern, stem, perl = TRUE)
    if (rx == -1L) next                          # no match, try next schema
    
    starts  <- attr(rx, "capture.start")
    lengths <- attr(rx, "capture.length")
    cnames  <- attr(rx, "capture.names")
    
    # Extract each named group; optional groups that didn't participate → NA
    groups <- setNames(
      lapply(seq_along(cnames), function(i) {
        if (starts[i] == -1L) NA_character_
        else substr(stem, starts[i], starts[i] + lengths[i] - 1L)
      }),
      cnames
    )
    groups <- groups[nzchar(cnames)]             # drop any unnamed groups
    
    # Merge with schema defaults (fill fields absent from the pattern)
    row <- c(groups, schema$defaults[!names(schema$defaults) %in% names(groups)])
    
    # Type coercions
    row$scenario_year <-
      if (!is.null(row$scenario_year) && !is.na(row$scenario_year))
        as.integer(row$scenario_year)
    else NA_integer_
    
    # Empty string from a non-participating optional group → NA
    if (is.null(row$spatialization_method) ||
        identical(row$spatialization_method, ""))
      row$spatialization_method <- NA_character_
    
    row$path     <- path
    row$filename <- fname
    
    # Return in canonical column order; pad any still-missing columns
    out <- row[intersect(CANONICAL_FIELDS, names(row))]
    for (col in setdiff(CANONICAL_FIELDS, names(out))) out[[col]] <- NA_character_
    
    return(as.data.frame(out[CANONICAL_FIELDS], stringsAsFactors = FALSE))
  }
  
  warning("Unrecognised filename format: ", fname)
  return(NULL)
}

# ── Collect file metadata ─────────────────────────────────────────────────────

tif_files <- list.files(INPUT_DIR, pattern = "\\.tif$", full.names = TRUE)
cat(sprintf("GeoTIFF files found    : %d\n", length(tif_files)))

file_meta <- bind_rows(lapply(tif_files, parse_tif_filename))
cat(sprintf("Files parsed           : %d\n", nrow(file_meta)))

# ── Define expected combinations (robust) ─────────────────────────────────────

# Unique method combinations (schema-agnostic)
method_combos <- file_meta %>%
  distinct(unbias_sequence, calibration_method,
           correction_algorithm, spatialization_method)

# Unique pollutant × year × run_type combinations
pollutant_year_runtype <- file_meta %>%
  distinct(pollutant, scenario_year, run_type)

# Full Cartesian product
expected_grid <- merge(method_combos, pollutant_year_runtype)

# ── Identify missing combinations (no TIFF) ───────────────────────────────────

missing_meta <- expected_grid %>%
  left_join(file_meta,
            by = c("pollutant", "unbias_sequence", "calibration_method",
                   "correction_algorithm", "spatialization_method",
                   "scenario_year", "run_type")) %>%
  filter(is.na(path)) %>%
  mutate(
    path     = NA_character_,
    filename = NA_character_
  )

# ── Extract global raster statistics ─────────────────────────────────────────

extract_global_stats <- function(row) {
  
  file_info <- tryCatch(file.info(row$path), error = function(e) NULL)
  
  r <- tryCatch(rast(row$path), error = function(e) {
    warning("Cannot read raster: ", row$path, " — ", e$message)
    NULL
  })
  
  make_result <- function(valid, error_type, pixel_values = numeric(0), n_na = 0L) {
    has_values <- length(pixel_values) > 0
    data.frame(
      row,
      valid_raster = valid,
      error_type   = error_type,
      n_na_pixels  = n_na,
      file_size    = if (!is.null(file_info)) file_info$size  else NA_real_,
      file_mtime   = if (!is.null(file_info)) format(file_info$mtime, "%Y-%m-%d %H:%M:%S") else NA_character_,
      stat_min     = if (has_values) min(pixel_values)              else NA_real_,
      stat_q25     = if (has_values) quantile(pixel_values, 0.25)   else NA_real_,
      stat_median  = if (has_values) median(pixel_values)           else NA_real_,
      stat_mean    = if (has_values) mean(pixel_values)             else NA_real_,
      stat_q75     = if (has_values) quantile(pixel_values, 0.75)   else NA_real_,
      stat_max     = if (has_values) max(pixel_values)              else NA_real_,
      stat_sd      = if (has_values) sd(pixel_values)               else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  if (is.null(r)) return(make_result(FALSE, "read_error"))
  
  pixel_values <- tryCatch(as.vector(values(r)), error = function(e) {
    warning("Cannot extract values: ", row$path, " — ", e$message)
    NULL
  })
  
  if (is.null(pixel_values)) return(make_result(FALSE, "values_error"))
  n_na <- sum(is.na(pixel_values) | is.nan(pixel_values) | is.infinite(pixel_values))
  pixel_values <- pixel_values[is.finite(pixel_values)]
  
  pixel_values <- pixel_values[is.finite(pixel_values)]
  
  if (length(pixel_values) == 0) {
    warning("No finite pixel values in: ", row$path)
    return(make_result(FALSE, "no_finite_values", n_na = n_na))
  }
  
  make_result(TRUE,
              error_type = if (n_na > 0) "has_na_pixels" else NA_character_,
              pixel_values,
              n_na = n_na)
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
cat(sprintf("Valid rasters          : %d\n",  sum( stats_df$valid_raster, na.rm = TRUE)))
cat(sprintf("Invalid rasters        : %d\n",  sum(!stats_df$valid_raster, na.rm = TRUE)))

# ── Add missing TIFF combinations to stats ────────────────────────────────────

missing_stats <- missing_meta %>%
  mutate(
    valid_raster = FALSE,
    error_type   = "missing_tif",
    n_na_pixels  = NA_integer_,  
    file_size    = NA_real_,
    file_mtime   = NA_character_,
    stat_min     = NA_real_,
    stat_q25     = NA_real_,
    stat_median  = NA_real_,
    stat_mean    = NA_real_,
    stat_q75     = NA_real_,
    stat_max     = NA_real_,
    stat_sd      = NA_real_
  )

stats_df <- bind_rows(stats_df, missing_stats)
cat(sprintf("Missing rasters        : %d\n", sum(is.na(stats_df$filename))))

# ── Save CSV ──────────────────────────────────────────────────────────────────

dir.create(dirname(OUTPUT_CSV), showWarnings = FALSE, recursive = TRUE)
write.csv(stats_df, OUTPUT_CSV, row.names = FALSE)
cat(sprintf("Statistics CSV written : %s\n", OUTPUT_CSV))

# ── Summary of invalid rasters ────────────────────────────────────────────────

stats_df %>%
  filter(!valid_raster) %>%
  group_by(scenario_year, run_type, pollutant,
           method = glue("{unbias_sequence}.{calibration_method}.{correction_algorithm}")) %>%
  summarize(spatialization_method = paste0(spatialization_method, collapse = ","),
            .groups = "drop") %>%
  knitr::kable()