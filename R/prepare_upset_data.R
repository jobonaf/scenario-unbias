# =============================================================================
# file           : prepare_upset_data.R
# description    : Prepare BCM dataset for UpSet analysis (criteria + skill scores)
# author         : Giovanni Bonafe'
# created        : 2025-02-09
# version        : 1.2
# dependencies   : dplyr, readr, stringr, glue
# =============================================================================

library(dplyr)
library(readr)
library(stringr)
library(glue)
library(tidyr)

# =============================================================================
# 1. Configuration
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript prepare_upset_data.R <INQUINANTE>")
}

pollutant <- args[1]
#pollutant <- "O3"   # e.g. PM25, NO2, O3, ...
if (pollutant == "PM25") {
  selected_clusters <- c(9, 10, 11)            # Jaccard PM25
} else if (pollutant == "NO2") {
  selected_clusters <- c(2, 8, 9, 10, 11)       # Jaccard NO2
} else if (pollutant == "O3") {
  selected_clusters <- c(2, 7)               # Jaccard O3
} 

# =============================================================================
# 2. Input paths (using glue)
# =============================================================================

path_scores <- glue(
  "data/models_verification/fairmode_exercise/skill_scores_{pollutant}.csv"
)

path_clustering <- glue(
  "data/clustering/clustering_{pollutant}.csv"
)

path_consistency <- "data/criteria-assessment/basecase_consistent.csv"
path_pca <- glue("data/clustering-with-pca/scores_pca_{pollutant}.csv")

# =============================================================================
# 3. Read skill scores
# =============================================================================

scores <- read_csv(path_scores, show_col_types = FALSE) %>%
  rename(BCM = model)

stopifnot(nrow(scores) > 0)

# =============================================================================
# 4.1 Read clustering results
# =============================================================================

clustering <- read_csv(path_clustering, show_col_types = FALSE) %>%
  rename(BCM = member)

stopifnot(nrow(clustering) > 0)

# =============================================================================
# 4.2 Read consistency results
# =============================================================================

consistent <- read_csv(path_consistency, show_col_types = FALSE) %>%
  filter(pollutant==pollutant) %>%
  mutate(BCM = paste(sequence,calibration,adjustment,spatialization,sep="."),
         basecase_consistency = TRUE)

stopifnot(nrow(consistent) > 0)

# =============================================================================
# 4.3 Read PCA results
# =============================================================================

pca <- read_csv(path_pca, show_col_types = FALSE) %>%
  separate(Scenario, sep = "_", into = c("pollutant","BCM")) %>%
  filter(pollutant==pollutant) %>%
  transmute(
	    BCM, 
	    score_PCA_weighted,
	    score_cluster,
	    score_consenso,
	    score_PCchamp
  )

stopifnot(nrow(consistent) > 0)

# =============================================================================
# 5. Merge scores + clustering + other
# =============================================================================
# NOTE:
# - BCM names are kept AS IS
# - suffixes (e.g. .ok, .ked, etc.) are meaningful and preserved

bcm_data <- scores %>%
  left_join(
    clustering %>%
      select(
        BCM,
        cluster_index
      ),
    by = "BCM"
  ) %>%
  separate(
    col = "BCM", sep = "\\.",
    into = c("sequence","calibration","adjustment","spatialization"), 
    fill = "right", remove = FALSE
  ) %>%
  separate(
    col = "spatialization", sep = "\\+",
    into = c(NA,"extra"),
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(cross_validation = !is.na(extra) & extra=="xv") %>%
  left_join(
    consistent %>%
      select(BCM, basecase_consistency)
  ) %>%
  left_join(pca)

# =============================================================================
# 6. Criterion 1: BCM naming / methodological family
# =============================================================================
# Definition:
#   CSA OR ((SCA OR CAS) AND ked)

bcm_data <- bcm_data %>%
  mutate(
    concentration_variability =
      !cross_validation &
      !(calibration %in% c("All","Grid")) &
      (sequence == "CSA" |
      (
        sequence %in% c("SCA", "CAS") &
          spatialization == "ked"
      ))
  )

# =============================================================================
# 7.1 Criterion 2: Membership in selected clusters
# =============================================================================
# (example definition, easily adjustable)


bcm_data <- bcm_data %>%
  mutate(
    Jaccard_clustering =
      !is.na(cluster_index) & cluster_index %in% selected_clusters
  )

# =============================================================================
# 7.2 Criterion 3: Basecase consistency
# =============================================================================

bcm_data <- bcm_data %>%
  mutate(basecase_consistency = replace_na(basecase_consistency, FALSE))

# =============================================================================
# 7.3 Criterion 4: PCA
# =============================================================================

bcm_data <- bcm_data %>%
  mutate(
    pca_informativeness = score_PCA_weighted>=quantile(score_PCA_weighted,2/3,na.rm=T),
    cluster_centrality = score_cluster>=quantile(score_cluster,2/3,na.rm=T),
    consensus_score = score_consenso>=quantile(score_consenso,2/3,na.rm=T),
    pc_champion = score_PCchamp>0
  ) %>%
  mutate(
    pca_informativeness = replace_na(pca_informativeness, FALSE),
    cluster_centrality = replace_na(cluster_centrality, FALSE),
    consensus_score = replace_na(consensus_score, FALSE),
    pc_champion = replace_na(pc_champion, FALSE)
  )

# =============================================================================
# 8. Final dataset for UpSet
# =============================================================================

bcm_upset <- bcm_data %>%
  transmute(
    BCM,
    
    # --- criteria ---
    concentration_variability,
    Jaccard_clustering,
    cross_validation,
    basecase_consistency,
    pca_informativeness,
    cluster_centrality,
    consensus_score,
    pc_champion,
    
    # --- skill scores ---
    ME,
    MAE,
    RMSE,
    R2,
    IOA,
    correlation,
    n_pixels
  ) %>%
  distinct()

bcm_upset_grouped <- bcm_upset%>%
  group_by(across(where(is.logical))) %>%
  summarise(
    across(where(is.numeric), ~ round(median(.x, na.rm = TRUE),3)),
    across(where(is.character), ~ paste(unique(.x), collapse = " ")),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(RMSE)

# =============================================================================
# 9. Output
# =============================================================================

output_path <- glue("data/criteria-assessment/bcm_upset_data_{pollutant}.csv")
output_grouped_path <- glue("data/criteria-assessment/grouped_upset_data_{pollutant}.csv")

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

write_csv(bcm_upset, output_path)
write_csv(bcm_upset_grouped, output_grouped_path)

cat("Prepared UpSet dataset written to:\n", output_path, "\nand\n", output_grouped_path,"\n")
