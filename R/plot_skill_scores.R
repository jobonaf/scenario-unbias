# =============================================================================
# title          :plot_skill_scores
# description    :Create comparison plots of model skill scores from CSV files
# author         :Giovanni Bonafe'
# date           :20250203
# version        :1.1
# notes          :Reads skill scores from CSV and generates multi-panel plots
# R_version      :3.5.2
# =============================================================================

# Define command-line options
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("glue"))

option_list <- list(
  make_option(c("-i", "--input_csv"), type = "character", 
              default = NULL,
              help = "Path to skill scores CSV file [REQUIRED]"),
  make_option(c("-p", "--pollutant"), type = "character", 
              default = NULL,
              help = "Pollutant name for plot title (auto-detected from filename if NULL) [default: %default]"),
  make_option(c("-m", "--metrics"), type = "character", 
              default = "IOA,RMSE,MAE,R2,correlation",
              help = "Comma-separated list of metrics to plot [default: %default]"),
  make_option(c("-r", "--reference_lines"), type = "character", 
              default = NULL,
              help = "Comma-separated reference values for each metric (use 'NA' to skip a metric) [default: %default]"),
  make_option(c("-o", "--output_dir"), type = "character", 
              default = ".",
              help = "Output directory for plots [default: %default]"),
  make_option(c("-n", "--nrows"), type = "integer", 
              default = 1,
              help = "Number of rows for facet layout [default: %default]"),
  make_option(c("-s", "--facet_scales"), type = "character", 
              default = "free_x",
              help = "Facet scales: 'free', 'free_x', 'free_y', 'fixed' [default: %default]"),
  make_option(c("-M", "--plot_median"), type = "logical", 
              default = TRUE,
              help = "Plot median reference line [default: %default]"),
  make_option(c("-w", "--width"), type = "numeric", 
              default = NULL,
              help = "PDF width in inches (NULL for auto) [default: %default]"),
  make_option(c("-H", "--height"), type = "numeric", 
              default = NULL,
              help = "PDF height in inches (NULL for auto) [default: %default]"),
  make_option(c("-f", "--filename"), type = "character", 
              default = NULL,
              help = "Output filename (NULL for auto-generated) [default: %default]")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list,
                               description = "Create skill scores comparison plots"))

# Check required arguments
if (is.null(opt$input_csv)) {
  stop("Input CSV file is required. Use -i or --input_csv to specify the file.")
}

# Extract parameters
input_csv <- opt$input_csv

# Auto-detect pollutant from filename if not provided
if (is.null(opt$pollutant)) {
  # Extract pollutant from filename patterns like "skill_scores_NO2.csv" or "skill_scores_PM25.csv"
  basename_file <- basename(input_csv)
  
  # Try to match common patterns
  if (grepl("_NO2[._]", basename_file, ignore.case = TRUE)) {
    pollutant <- "NO2"
  } else if (grepl("_O3[._]", basename_file, ignore.case = TRUE)) {
    pollutant <- "O3"
  } else if (grepl("_PM25[._]", basename_file, ignore.case = TRUE)) {
    pollutant <- "PM25"
  } else if (grepl("_PM10[._]", basename_file, ignore.case = TRUE)) {
    pollutant <- "PM10"
  } else {
    # Default fallback
    pollutant <- "Unknown"
    warning("Could not auto-detect pollutant from filename. Using 'Unknown'. Use -p to specify.")
  }
  
  message(glue("Auto-detected pollutant: {pollutant}"))
} else {
  pollutant <- opt$pollutant
}

metrics_to_plot <- strsplit(opt$metrics, ",")[[1]]
# Trim whitespace from metrics
metrics_to_plot <- trimws(metrics_to_plot)

# Parse reference lines
reference_values <- NULL
if (!is.null(opt$reference_lines)) {
  ref_strings <- strsplit(opt$reference_lines, ",")[[1]]
  ref_strings <- trimws(ref_strings)
  
  # Convert to numeric, NA for "NA" strings
  reference_values <- sapply(ref_strings, function(x) {
    if (toupper(x) == "NA" || x == "") {
      return(NA)
    } else {
      return(as.numeric(x))
    }
  })
  
  # Check length matches metrics
  if (length(reference_values) != length(metrics_to_plot)) {
    warning(glue("Number of reference values ({length(reference_values)}) does not match number of metrics ({length(metrics_to_plot)}). Reference lines will not be plotted."))
    reference_values <- NULL
  }
}

output_dir <- opt$output_dir
nrows <- opt$nrows
facet_scales <- opt$facet_scales
plot_median <- opt$plot_median
pdf_width <- opt$width
pdf_height <- opt$height
output_filename <- opt$filename

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(futile.logger)

# =============================================================================
# Initialize logging
# =============================================================================
flog.appender(appender.console())
flog.threshold(INFO)
flog.info("Starting skill scores plotting")
flog.info("Input CSV: %s", input_csv)
flog.info("Pollutant: %s", pollutant)
flog.info("Metrics to plot: %s", paste(metrics_to_plot, collapse = ", "))
if (!is.null(reference_values)) {
  flog.info("Reference values: %s", paste(reference_values, collapse = ", "))
}

# =============================================================================
# Load data
# =============================================================================
flog.info("Loading skill scores from CSV...")

if (!file.exists(input_csv)) {
  flog.error(glue("Input file not found: {input_csv}"))
  stop(glue("Input file not found: {input_csv}"))
}

skill_data <- read.csv(input_csv, stringsAsFactors = FALSE)

flog.info(glue("Loaded {nrow(skill_data)} models"))
flog.info(glue("Available columns: {paste(names(skill_data), collapse = ', ')}"))

# =============================================================================
# Pollutant extended names
# =============================================================================
poll_ext <- case_match(
  pollutant,
  "PM25" ~ "PM2.5",
  "PM10" ~ "PM10",
  "NO2" ~ "nitrogen dioxide",
  "O3" ~ "ozone",
  .default = pollutant
)

# =============================================================================
# Function: Create selection plot with multiple panels
# =============================================================================
create_selection_plot <- function(data, pollutant, output_dir, nrows = 1, 
                                  metrics_to_plot, reference_values = NULL,
                                  facet_scales = "free_x", 
                                  plot_median = TRUE, pdf_width = NULL, 
                                  pdf_height = NULL, output_filename = NULL) {
  
  flog.info("Creating skill scores comparison plot...")
  
  # Check if model column exists
  if (!"model" %in% names(data)) {
    flog.error("Data must contain a 'model' column")
    stop("Data must contain a 'model' column")
  }
  
  # Select only the requested metrics
  valid_metrics <- metrics_to_plot[metrics_to_plot %in% names(data)]
  if (length(valid_metrics) == 0) {
    flog.error(glue("None of the requested metrics found in data: {paste(metrics_to_plot, collapse = ', ')}"))
    stop("No valid metrics found for plotting")
  }
  
  if (length(valid_metrics) < length(metrics_to_plot)) {
    missing_metrics <- setdiff(metrics_to_plot, valid_metrics)
    flog.warn(glue("Some metrics not found in data: {paste(missing_metrics, collapse = ', ')}"))
  }
  
  # Prepare data for ggplot
  plot_data <- data %>%
    select(model, all_of(valid_metrics)) %>%
    pivot_longer(cols = -model, names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = valid_metrics))
  
  # Order models based on the first metric in the list
  primary_metric <- valid_metrics[1]
  model_order <- data %>%
    arrange(desc(.data[[primary_metric]])) %>%
    pull(model)
  
  plot_data$model <- factor(plot_data$model, levels = model_order)
  
  # Calculate PDF dimensions dynamically if not provided
  n_models <- length(unique(plot_data$model))
  n_metrics <- length(valid_metrics)
  
  if (is.null(pdf_width)) {
    pdf_width <- 1 + (n_metrics/nrows * 1.7)
    pdf_width <- max(4, min(pdf_width, 12))
  }
  
  if (is.null(pdf_height)) {
    pdf_height <- 1 + (n_models*nrows * 0.2)
    pdf_height <- max(4, min(pdf_height, 12))
  }
  
  flog.info(glue("Plot dimensions: {pdf_width} x {pdf_height} inches"))
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = value, y = model)) +
    geom_point(size = 2, color = "steelblue") +
    facet_wrap(~ metric, scales = facet_scales, nrow = nrows) +
    labs(title = glue("Pollutant: {poll_ext}"),
         x = "metric value", 
         y = "BCMs") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8, hjust = 0),
          plot.title = element_text(hjust = 0.5, size = 14),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          strip.text = element_text(size = 9))
  
  # Add median line if requested
  if (plot_median) {
    median_data <- plot_data %>% 
      group_by(metric) %>% 
      summarize(median_val = median(value, na.rm = TRUE), .groups = "drop")
    
    p <- p +
      geom_vline(data = median_data,
                 aes(xintercept = median_val), 
                 color = "grey60", linetype = "dashed", linewidth = 0.5)
    
    flog.info("Median reference lines added")
  }
  
  # Add custom reference lines if provided
  if (!is.null(reference_values)) {
    # Create data frame for reference lines
    ref_df <- data.frame(
      metric = valid_metrics,
      ref_value = reference_values[match(valid_metrics, metrics_to_plot)],
      stringsAsFactors = FALSE
    )
    
    # Filter out NA values
    ref_df <- ref_df %>% filter(!is.na(ref_value))
    
    if (nrow(ref_df) > 0) {
      ref_df$metric <- factor(ref_df$metric, levels = valid_metrics)
      
      p <- p +
        geom_vline(data = ref_df,
                   aes(xintercept = ref_value), 
                   color = "grey60", linetype = "solid", linewidth = 0.5)
      
      flog.info(glue("Custom reference lines added for {nrow(ref_df)} metrics"))
      for (i in 1:nrow(ref_df)) {
        flog.info(glue("  {ref_df$metric[i]}: {ref_df$ref_value[i]}"))
      }
    }
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    flog.info(glue("Created output directory: {output_dir}"))
  }
  
  # Generate output filename if not provided
  if (is.null(output_filename)) {
    output_filename <- glue("skill_scores_comparison_{pollutant}.pdf")
  }
  
  output_file <- file.path(output_dir, output_filename)
  
  # Save plot as PDF
  ggsave(output_file, p, width = pdf_width, height = pdf_height, 
         units = "in", dpi = 300)
  
  flog.info(glue("Plot saved: {output_file}"))
  flog.info(glue("Metrics plotted: {paste(valid_metrics, collapse = ', ')}"))
  
  return(p)
}

# =============================================================================
# Create the plot
# =============================================================================
plot <- create_selection_plot(
  data = skill_data,
  pollutant = pollutant,
  output_dir = output_dir,
  nrows = nrows,
  metrics_to_plot = metrics_to_plot,
  reference_values = reference_values,
  facet_scales = facet_scales,
  plot_median = plot_median,
  pdf_width = pdf_width,
  pdf_height = pdf_height,
  output_filename = output_filename
)

# =============================================================================
# Print summary statistics
# =============================================================================
flog.info("Summary statistics for plotted metrics:")
for (metric in metrics_to_plot) {
  if (metric %in% names(skill_data)) {
    metric_values <- skill_data[[metric]]
    flog.info(glue("  {metric}:"))
    flog.info(glue("    Mean: {round(mean(metric_values, na.rm = TRUE), 3)}"))
    flog.info(glue("    Median: {round(median(metric_values, na.rm = TRUE), 3)}"))
    flog.info(glue("    Min: {round(min(metric_values, na.rm = TRUE), 3)}"))
    flog.info(glue("    Max: {round(max(metric_values, na.rm = TRUE), 3)}"))
  }
}

flog.info("Skill scores plotting completed successfully!")

# =============================================================================
# End of script
# =============================================================================