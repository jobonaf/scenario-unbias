# =============================================================================
# title          :upset_scores
# description    :Two UpSet plots for BCM selection analysis:
#                 1. Classic diagram with size bars (ordered by numerosity)
#                 2. Score-based diagram with boxplots (ordered by median RMSE)
# author         :Giovanni Bonafe'
# date           :20250209
# version        :1.4
# notes          :Requires ComplexUpset, ggplot2, dplyr
# R_version      :3.5.2
# =============================================================================

library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(readr)
library(glue)
library(scales)

# =============================================================================
# 1. Data preparation: read prepared BCM dataset
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript upset_scores.R <INQUINANTE>")
}

pollutant <- args[1]
#pollutant <- "O3"   # same convention as in preparation script

input_path <- glue(
  "data/criteria-assessment/bcm_upset_data_{pollutant}.csv"
)

bcm_data <- read_csv(input_path, show_col_types = FALSE)

stopifnot(nrow(bcm_data) > 0)

# Automatically detect criteria columns (boolean columns)
criteria <- names(bcm_data)[sapply(bcm_data, is.logical)]

# Explicit criterion for BCMs satisfying NO selection criteria
bcm_data <- bcm_data %>%
  mutate(
    `no criterion fulfilled` = !Reduce(`|`, across(all_of(criteria)))
  )

# Update criteria list
criteria <- c(criteria, "no criterion fulfilled")

cat("Detected criteria:", paste(criteria, collapse = ", "), "\n")
cat("Number of BCMs:", nrow(bcm_data), "\n\n")

# =============================================================================
# 1.5. Define consistent set order for both plots
# =============================================================================

# Calculate set sizes to order by frequency (most common first)
set_sizes <- sapply(criteria, function(crit) sum(bcm_data[[crit]], na.rm = TRUE))
criteria_ordered <- names(sort(set_sizes, decreasing = TRUE))

cat("Set order (by frequency):\n")
print(data.frame(criterion = criteria_ordered, size = set_sizes[criteria_ordered]))
cat("\n")

# =============================================================================
# 2. UpSet Plot 1: Classic diagram (ordered by intersection size)
# =============================================================================

cat("Creating classic UpSet plot...\n")

upset_classic <- upset(
  bcm_data,
  criteria_ordered,  # Use ordered criteria
  
  # --- Layout ---
  width_ratio = 0.15,
  min_size = 1,
  keep_empty_groups = TRUE,
  
  # --- Mode: inclusive intersection (at least these criteria) ---
  mode = "inclusive_intersection",
  
  # --- Ordering: by intersection size (default = descending) ---
  sort_intersections_by = 'cardinality',
  sort_intersections = 'descending',
  sort_sets = FALSE,  # Respect our ordering
  
  # --- Remove connection lines in intersection matrix ---
  matrix = intersection_matrix(
    geom = geom_point(size = 3),
    segment = geom_segment(linetype = 'blank')  # Remove connecting lines
  )
)

# Save classic plot
output_classic <- glue("data/criteria-assessment/upset_sizes_{pollutant}.pdf")
if (!dir.exists(dirname(output_classic))) {
  dir.create(dirname(output_classic), recursive = TRUE)
}
ggsave(output_classic, upset_classic, width = 10, height = 8, units = "in", dpi = 300)
cat("✓ Classic UpSet plot saved:", output_classic, "\n\n")

# =============================================================================
# 3. UpSet Plot 2: With score panels ordered by median RMSE
# =============================================================================

cat("Creating UpSet plot with score panels ordered by RMSE...\n")

# Compute intersection-level statistics
bcm_with_intersection <- bcm_data %>%
  mutate(
    intersection_id = apply(.[criteria_ordered], 1, function(row) {
      paste(as.integer(row), collapse = "|")
    })
  )

# Compute intersection-level statistics for inclusive mode
# For each intersection, include all BCMs that have AT LEAST those criteria
intersection_stats <- lapply(unique(bcm_with_intersection$intersection_id), function(id) {
  # Parse the intersection pattern
  values <- as.integer(strsplit(id, "\\|")[[1]])
  required_criteria <- criteria_ordered[values == 1]
  
  # Find all BCMs that have AT LEAST these criteria active
  if (length(required_criteria) == 0) {
    # Empty intersection - all BCMs
    matching_bcms <- bcm_data
  } else {
    # Filter BCMs that have all required criteria = TRUE
    matching_bcms <- bcm_data
    for (crit in required_criteria) {
      matching_bcms <- matching_bcms[matching_bcms[[crit]] == TRUE, ]
    }
  }
  
  data.frame(
    intersection_id = id,
    median_rmse = median(matching_bcms$RMSE, na.rm = TRUE),
    size = nrow(matching_bcms),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  arrange(median_rmse)

cat("Intersection ordering by median RMSE:\n")
print(head(intersection_stats, 10))
cat("\n")

# Create list of intersections in RMSE order
# Each element is a vector of CRITERIA names that define that intersection
intersections_ordered <- lapply(intersection_stats$intersection_id, function(id) {
  # Split "1|0|1|0" into c(1, 0, 1, 0)
  values <- as.integer(strsplit(id, "\\|")[[1]])  # [[1]] to extract from list, escape |
  # Get criteria where value is 1 (TRUE)
  active_criteria <- criteria_ordered[values == 1]
  
  return(active_criteria)
})

cat("\nTotal intersections to plot:", length(intersections_ordered), "\n")

cat("Example intersections (first 5):\n")
print(head(intersections_ordered, 5))
cat("\n")

upset_scores <- upset(
  bcm_data,
  intersect = criteria_ordered,           # Use same ordered criteria
  intersections = intersections_ordered,  # custom order
  sort_intersections = FALSE,             # respect the order we provided
  sort_sets = FALSE,                      # Respect our ordering
  
  # --- Layout ---
  width_ratio = 0.35,
  min_size = 1,
  keep_empty_groups = TRUE,
  set_sizes = FALSE,
  
  # --- Mode: inclusive intersection (at least these criteria) ---
  mode = "inclusive_intersection",
  
  # --- Remove connection lines in intersection matrix ---
  matrix = intersection_matrix(
    geom = geom_point(size = 2),
    segment = geom_segment(linetype = 'blank')  # Remove connecting lines
  ),
  
  # --- Score panels as base_annotations ---
  base_annotations = list(
    
    # Panel 1: RMSE (lower is better)
    "RMSE" = ggplot(mapping = aes(x = intersection, y = RMSE)) +
      geom_boxplot(outlier.shape = NA, fill = "lightcoral", alpha = 0.6) +
      geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.2, color = "darkred") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      labs(y = "RMSE") +
      scale_y_continuous(trans="pseudo_log") +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()),
    
    # Panel 2: ME (closer to 0 is better)
    "Bias" = ggplot(mapping = aes(x = intersection, y = ME)) +
      geom_boxplot(outlier.shape = NA, fill = "lightskyblue", alpha = 0.6) +
      geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.2, color = "steelblue") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      labs(y = "Bias") +
      scale_y_continuous(trans="pseudo_log") +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()),
    
    # Panel 3: Correlation (higher is better)
    "Correlation" = ggplot(mapping = aes(x = intersection, y = correlation)) +
      geom_boxplot(outlier.shape = NA, fill = "lightgreen", alpha = 0.6) +
      geom_jitter(width = 0.2, height = 0, alpha = 0.5, size = 1.2, color = "darkgreen") +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      labs(y = "Correlation") +
      scale_y_continuous(trans=logit_trans(),
                         breaks=c(0.2,0.4,0.6,0.8,0.9,0.95,
                                  0.99,0.995,0.999)) +
      theme_bw() +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  )
)

# Save scores plot
output_scores <- glue("data/criteria-assessment/upset_scores_{pollutant}.pdf")
ggsave(output_scores, upset_scores, width = 10, height = 8, units = "in", dpi = 300)
cat("✓ UpSet plot with scores saved:", output_scores, "\n\n")

# =============================================================================
# 4. Summary
# =============================================================================

cat("=============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=============================================================================\n")
cat("Two plots created:\n")
cat("  1. Classic UpSet (by size):", output_classic, "\n")
cat("  2. UpSet with scores (by RMSE):", output_scores, "\n")
cat("\nKey message:\n")
cat("  - Left side of scores plot = lowest RMSE (best performance)\n")
cat("  - BCMs selected by more criteria tend to perform better\n")
cat("=============================================================================\n")

# =============================================================================
# End of script
# =============================================================================
