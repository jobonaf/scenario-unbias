## summarise_skill_scores.R
##
## Reads skill_scores_{pollutant}.csv files, parses the model name into its
## components, then computes summary statistics (mean, median, p10, p25, p75,
## p90) for each skill metric, stratified by every meaningful grouping derived
## from the model name parts.
##
## Input  : {path_in}/skill_scores_{pollutant}.csv
## Output : {path_out}/skill_scores_summary.csv  (long-format)
##
## Model name structure (dot-separated, up to 4 parts):
##   1. sequence           – e.g. CA, CAS
##   2. calibration_scope  – All | Grid | Each | Cell | Neigh
##   3. correction_method  – Add | Lin | Mult
##   4. spatial_algorithm  – idw | ked | <NA> if absent
##
## Compatible with R 3.5.2 (no across(), no group_modify(), no native pipe)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

# ── User configuration ────────────────────────────────────────────────────────

path_in    <- "data/models_verification/fairmode_exercise"  # directory containing skill_scores_*.csv
path_out   <- path_in                                       # directory for output CSV
pollutants <- c("NO2", "PM25", "O3")

# Skill metrics to summarise (all numeric columns except n_pixels)
skill_metrics <- c("ME", "MAE", "RMSE", "R2", "IOA", "correlation")

# Stratification groupings: each element lists the variable(s) to group by.
strata <- list(
  # Single variables
  by_sequence                            = "sequence",
  by_calibration_scope                   = "calibration_scope",
  by_correction_method                   = "correction_method",
  by_spatial_algorithm                   = "spatial_algorithm",
  # Pairwise combinations
  by_sequence_calibration_scope          = c("sequence", "calibration_scope"),
  by_sequence_correction_method          = c("sequence", "correction_method"),
  by_sequence_spatial_algorithm          = c("sequence", "spatial_algorithm"),
  by_calibration_scope_correction_method = c("calibration_scope", "correction_method"),
  by_calibration_scope_spatial_algorithm = c("calibration_scope", "spatial_algorithm"),
  by_correction_method_spatial_algorithm = c("correction_method", "spatial_algorithm")
)

# ── Helpers ───────────────────────────────────────────────────────────────────

#' Parse the 'model' column into its four named components.
#' The four new columns are appended at the end of the data frame.
parse_model_name <- function(df) {
  parts                <- strsplit(df$model, ".", fixed = TRUE)
  df$sequence          <- sapply(parts, function(x) x[1])
  df$calibration_scope <- sapply(parts, function(x) x[2])
  df$correction_method <- sapply(parts, function(x) x[3])
  df$spatial_algorithm <- sapply(parts, function(x) if (length(x) >= 4) x[4] else NA_character_)
  df
}

#' Compute summary statistics for one skill metric from a numeric vector.
summarise_metric <- function(vals, metric_name) {
  data.frame(
    metric = metric_name,
    n      = sum(!is.na(vals)),
    mean   = mean(vals,           na.rm = TRUE),
    median = median(vals,         na.rm = TRUE),
    p10    = quantile(vals, 0.10, na.rm = TRUE, names = FALSE),
    p25    = quantile(vals, 0.25, na.rm = TRUE, names = FALSE),
    p75    = quantile(vals, 0.75, na.rm = TRUE, names = FALSE),
    p90    = quantile(vals, 0.90, na.rm = TRUE, names = FALSE),
    stringsAsFactors = FALSE
  )
}

#' Split df by pollutant + group_vars, compute stats for all metrics, and
#' return a single long data frame with boolean stratum flags.
summarise_stratum <- function(df, group_vars) {
  all_group_vars <- c("pollutant", group_vars)
  
  # One sub-data-frame per unique combination of grouping variables
  groups <- split(df, df[, all_group_vars, drop = FALSE], drop = TRUE)
  
  # Boolean flags indicating which model components define this stratum
  model_cols <- c("sequence", "calibration_scope", "correction_method", "spatial_algorithm")
  flag_df    <- setNames(
    as.data.frame(
      lapply(model_cols, function(col) col %in% group_vars),
      stringsAsFactors = FALSE
    ),
    paste0("by_", model_cols)
  )
  
  result <- lapply(groups, function(g) {
    # Key: first row of grouping columns (identical within the group)
    key          <- g[1L, all_group_vars, drop = FALSE]
    rownames(key) <- NULL
    
    # Stats for every metric stacked into rows
    metric_rows <- do.call(rbind, lapply(skill_metrics, function(m) {
      summarise_metric(g[[m]], m)
    }))
    
    # All four model component columns always present; NA when not in this stratum
    model_df <- setNames(
      as.data.frame(
        lapply(model_cols, function(col) {
          if (col %in% group_vars) g[[col]][1L] else NA_character_
        }),
        stringsAsFactors = FALSE
      ),
      model_cols
    )
    
    # Column order: pollutant | by_* flags | metric stats | model components
    cbind(
      key[, "pollutant", drop = FALSE],
      flag_df[rep(1L, nrow(metric_rows)), , drop = FALSE],
      metric_rows,
      model_df
    )
  })
  
  do.call(rbind, result)
}

# ── Main ──────────────────────────────────────────────────────────────────────

# 1. Read and parse all pollutant files
message("Reading input files...")

raw_list <- lapply(pollutants, function(pol) {
  fpath <- file.path(path_in, paste0("skill_scores_", pol, ".csv"))
  if (!file.exists(fpath)) {
    warning("File not found, skipping: ", fpath)
    return(NULL)
  }
  df           <- read_csv(fpath, col_types = cols())
  df$pollutant <- pol
  df
})

raw <- do.call(rbind, raw_list)
message(sprintf("Read %d rows across %d pollutant(s).",
                nrow(raw), length(unique(raw$pollutant))))

# 2. Parse model name (component columns appended at the end)
parsed <- parse_model_name(raw)

# 3. Compute stratified summaries
message("Computing stratified summaries...")

summary_long <- do.call(rbind, lapply(
  strata,
  function(group_vars) summarise_stratum(parsed, group_vars)
))

rownames(summary_long) <- NULL

# 4. Write output
dir.create(path_out, recursive = TRUE, showWarnings = FALSE)

# Format numeric columns to 4 significant digits before writing
num_cols_long <- c("mean", "median", "p10", "p25", "p75", "p90")
summary_long[, num_cols_long] <- lapply(summary_long[, num_cols_long],
                                        function(x) as.numeric(formatC(x, digits = 4, format = "g")))

out_path <- file.path(path_out, "skill_scores_summary.csv")
write.csv(summary_long, out_path, row.names = FALSE, quote = FALSE)

message("Done. Output written to: ", out_path)
message(sprintf("Rows: %d | Columns: %d", nrow(summary_long), ncol(summary_long)))

# ── Paper table: univariate strata, median ME/RMSE/correlation by pollutant ───

message("Building paper table...")

# Keep only univariate strata (exactly one by_* flag is TRUE)
model_cols <- c("sequence", "calibration_scope", "correction_method", "spatial_algorithm")
flag_cols  <- paste0("by_", model_cols)

n_true <- rowSums(summary_long[, flag_cols])
univar <- summary_long[n_true == 1, ]

# The stratum label is the value of whichever model column is not NA
univar$stratum_label <- apply(univar[, model_cols], 1, function(x) {
  v <- x[!is.na(x)]
  if (length(v) == 1L) v else NA_character_
})

# Keep only the three metrics of interest
univar <- univar[univar$metric %in% c("ME", "RMSE", "correlation"), ]

# One row per (stratum_label, pollutant, metric) — pick median
# (already stored in the 'median' column of summary_long)
univar_sub <- univar[, c("stratum_label", "pollutant", "metric", "n", "median")]

# Pivot to wide: rows = stratum_label, cols = pollutant x metric
# Build manually for R 3.5.2 compatibility
make_wide <- function(df, pol, metrics) {
  sub <- df[df$pollutant == pol, ]
  # n is the same for all metrics within a group — take from first metric
  n_df <- sub[sub$metric == metrics[1], c("stratum_label", "n")]
  names(n_df)[2] <- paste0("n_", pol)
  out <- n_df
  for (m in metrics) {
    col <- sub[sub$metric == m, c("stratum_label", "median")]
    names(col)[2] <- paste0(pol, "_", m)
    out <- merge(out, col, by = "stratum_label", all = TRUE)
  }
  out
}

metrics_ordered <- c("RMSE", "correlation", "ME")
pols            <- c("PM25", "NO2", "O3")

wide <- make_wide(univar_sub, pols[1], metrics_ordered)
for (pol in pols[-1]) {
  wide <- merge(wide, make_wide(univar_sub, pol, metrics_ordered),
                by = "stratum_label", all = TRUE)
}

# Canonical row order: sequence levels, then scope, then method, then algorithm
row_order <- c("SCA", "CSA", "CAS", "CA",
               "All", "Grid", "Each", "Cell", "Neigh",
               "Add", "Mult", "Lin",
               "ked", "ok", "idw", "tps")
wide <- wide[match(row_order, wide$stratum_label), ]
wide <- wide[!is.na(wide$stratum_label), ]

# Use a single n column (same across pollutants for univariate strata)
# Replace three n_* cols with one n taken from PM25 (or whichever is available)
wide$n <- ifelse(!is.na(wide$n_PM25), wide$n_PM25,
                 ifelse(!is.na(wide$n_NO2),  wide$n_NO2,  wide$n_O3))
wide <- wide[, c("stratum_label", "n",
                 paste0("PM25_", metrics_ordered),
                 paste0("NO2_",  metrics_ordered),
                 paste0("O3_",   metrics_ordered))]

# Format: 2 decimals for RMSE and ME, 3 for correlation
rmse_me_cols <- grep("_RMSE$|_ME$",   names(wide), value = TRUE)
corr_cols    <- grep("_correlation$", names(wide), value = TRUE)
wide[, rmse_me_cols] <- lapply(wide[, rmse_me_cols],
                               function(x) as.numeric(formatC(x, digits = 2, format = "f")))
wide[, corr_cols] <- lapply(wide[, corr_cols],
                            function(x) as.numeric(formatC(x, digits = 3, format = "f")))

# Write
table_path <- file.path(path_out, "skill_scores_table.csv")
write.csv(wide, table_path, row.names = FALSE, quote = FALSE)

message("Paper table written to: ", table_path)
message(sprintf("Rows: %d | Columns: %d", nrow(wide), ncol(wide)))