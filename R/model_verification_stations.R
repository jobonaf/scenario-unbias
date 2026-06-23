# =============================================================================
# title          :model_verification_stations
# description    :Skill scores against station observations for 3 years (2022-2024)
#                 Extracts pixel values at station locations, computes metrics
#                 per year then averages across years
# author         :Giovanni Bonafe'
# date           :20260507
# version        :0.1
# R_version      :3.5.2
# =============================================================================

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("glue"))

option_list <- list(
  make_option(c("-p", "--pollutant"), type = "character", default = "NO2",
              help = "Pollutant to analyze (NO2, O3, PM25) [default: %default]"),
  make_option(c("-o", "--obs_dir"), type = "character",
              default = "data/fairmode-wg5-exercise-202602/YEARLY",
              help = "Directory containing observation CSVs [default: %default]"),
  make_option(c("-m", "--model_dir"), type = "character",
              default = "data/processed_phase2",
              help = "Directory containing model TIFFs [default: %default]"),
  make_option(c("-u", "--output_dir"), type = "character",
              default = "data/models_verification/fairmode_exercise_phase2",
              help = "Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Phase 2 skill scores against station observations"))

pollutant  <- opt$pollutant
obs_dir    <- opt$obs_dir
model_dir  <- opt$model_dir
output_dir <- opt$output_dir

library(terra)
library(dplyr)
library(futile.logger)

flog.appender(appender.console())
flog.threshold(INFO)
flog.info("Starting phase 2 skill scores - pollutant: %s", pollutant)

years <- c(2022, 2023, 2024)

# =============================================================================
# Load observations (one CSV per year, same columns)
# =============================================================================
flog.info("Loading observations...")

obs_all <- bind_rows(lapply(years, function(yr) {
  f <- glue("{obs_dir}/Scenario_{yr}_Points/yearly_{pollutant}_{yr}.csv")
  if (!file.exists(f)) {
    flog.warn("Observation file not found: %s", f)
    return(NULL)
  }
  read.csv(f) %>% mutate(year = yr)
}))

flog.info("Loaded %d station-year rows", nrow(obs_all))

# =============================================================================
# Get model files (exclude year 2015)
# =============================================================================
model_files <- Sys.glob(glue("{model_dir}/{pollutant}_*_scenario*.tif"))
model_files <- model_files[!grepl("_2015\\.tif$", model_files)]

if (length(model_files) == 0) {
  flog.error("No model files found in %s", model_dir)
  stop("No model files found")
}

# Parse model name and year from filename
# e.g. NO2_CA.All.Add_unbiased_scenario_2022.tif
parse_tif_name <- function(f) {
  b <- tools::file_path_sans_ext(basename(f))
  b <- gsub(glue("^{pollutant}_"), "", b)
  # year is last 4 digits before extension
  yr  <- as.integer(regmatches(b, regexpr("[0-9]{4}$", b)))
  mdl <- gsub("_unbiased_scenario_[0-9]{4}$|_scenario_[0-9]{4}$", "", b)
  data.frame(file = f, model = mdl, year = yr, stringsAsFactors = FALSE)
}

model_index <- bind_rows(lapply(model_files, parse_tif_name))
flog.info("Found %d model-year combinations", nrow(model_index))

# =============================================================================
# Helper: compute metrics from paired vectors
# =============================================================================
compute_metrics <- function(obs, mod) {
  err  <- mod - obs
  n    <- length(obs)
  denom_ioa <- sum((abs(mod - mean(obs)) + abs(obs - mean(obs)))^2)
  list(
    ME          = mean(err),
    MAE         = mean(abs(err)),
    RMSE        = sqrt(mean(err^2)),
    R2          = 1 - sum((obs - mod)^2) / sum((obs - mean(obs))^2),
    IOA         = 1 - sum((obs - mod)^2) / denom_ioa,
    correlation = cor(mod, obs, method = "pearson"),
    n_stations  = n
  )
}

# =============================================================================
# Extract pixel values at station locations and compute metrics per model-year
# =============================================================================
flog.info("Extracting pixel values and computing metrics...")

results <- list()

for (mdl in unique(model_index$model)) {
  flog.info("Model: %s", mdl)
  yr_metrics <- list()
  
  for (yr in years) {
    row <- model_index[model_index$model == mdl & model_index$year == yr, ]
    if (nrow(row) == 0) {
      flog.warn("  Missing file for year %d", yr)
      next
    }
    
    obs_yr <- obs_all[obs_all$year == yr, ]
    if (nrow(obs_yr) == 0) next
    
    r <- tryCatch(rast(row$file), error = function(e) {
      flog.warn("  Cannot read %s: %s", basename(row$file), e$message)
      NULL
    })
    if (is.null(r)) next
    
    pts      <- vect(obs_yr, geom = c("Longitude", "Latitude"), crs = "EPSG:4326")
    pts      <- project(pts, crs(r))
    extracted <- terra::extract(r, pts)[, 2]
    
    paired <- data.frame(obs = obs_yr$Average, mod = extracted) %>%
      filter(!is.na(obs) & !is.na(mod))
    
    if (nrow(paired) < 3) {
      flog.warn("  Too few paired values for year %d", yr)
      next
    }
    
    yr_metrics[[as.character(yr)]] <- compute_metrics(paired$obs, paired$mod)
  }
  
  if (length(yr_metrics) == 0) next
  
  # Average metrics across years
  metric_names <- names(yr_metrics[[1]])
  avg <- sapply(metric_names, function(m) {
    mean(sapply(yr_metrics, function(x) x[[m]]))
  })
  
  results[[mdl]] <- as.data.frame(t(avg))
}

# =============================================================================
# Assemble and save
# =============================================================================
flog.info("Saving output...")

skill_df <- bind_rows(results, .id = "model") %>%
  select(model, ME, MAE, RMSE, R2, IOA, correlation, n_stations) %>%
  arrange(desc(IOA))

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_file <- glue("{output_dir}/skill_scores_phase2_{pollutant}.csv")
write.csv(skill_df, output_file, row.names = FALSE)

flog.info("Saved %d rows to %s", nrow(skill_df), output_file)
flog.info("Done.")