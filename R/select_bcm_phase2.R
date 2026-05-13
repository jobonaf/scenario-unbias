library(dplyr)

# ---- LOAD PHASE 1 SCORES ----
scores_ph1 <- bind_rows(
  read.csv("data/models_verification/fairmode_exercise/skill_scores_NO2.csv")  %>% mutate(pollutant = "NO2"),
  read.csv("data/models_verification/fairmode_exercise/skill_scores_O3.csv")   %>% mutate(pollutant = "O3"),
  read.csv("data/models_verification/fairmode_exercise/skill_scores_PM25.csv") %>% mutate(pollutant = "PM25")
) %>%
  filter(!grepl("xv$", model))

# ---- CRITERIO 1: top 5 per IOA in fase 1 ----
top5_ph1 <- scores_ph1 %>%
  group_by(pollutant) %>%
  slice_max(IOA, n = 5) %>%
  select(pollutant, model, IOA, RMSE, correlation) %>%
  arrange(pollutant, desc(IOA))
print(top5_ph1)

# output condensato
top5_ph1_compact <- top5_ph1 %>% ungroup() %>%
  summarise(BCM = paste(model, collapse = " "), .by = pollutant)
print(top5_ph1_compact)

# ---- LOAD UPSET DATA ----
upset_data <- bind_rows(
  read.csv("data/criteria-assessment/grouped_upset_data_NO2.csv") %>% mutate(pollutant = "NO2"),
  read.csv("data/criteria-assessment/grouped_upset_data_O3.csv")  %>% mutate(pollutant = "O3"),
  read.csv("data/criteria-assessment/grouped_upset_data_PM25.csv") %>% mutate(pollutant = "PM25")
) %>%
  filter(!grepl("xv", BCM))

# ---- CRITERIO 2: righe con IOA desc / RMSE asc fino ad avere almeno 4 BCM ----
top_upset <- upset_data %>%
  group_by(pollutant) %>%
  arrange(desc(IOA), RMSE, .by_group = TRUE) %>%
  mutate(nsum = cumsum(n)) %>%
  filter(replace(is.na(lag(nsum)), 1, TRUE) | lag(nsum) < 4) %>%
  summarise(BCM = paste(BCM, collapse = " "), n_bcm = sum(n)) %>%
  ungroup()

print(top_upset)


# ---- LOAD PHASE 2 SCORES ----
scores_ph2 <- bind_rows(
  read.csv("data/models_verification/fairmode_exercise_phase2/skill_scores_phase2_NO2.csv")  %>% mutate(pollutant = "NO2"),
  read.csv("data/models_verification/fairmode_exercise_phase2/skill_scores_phase2_O3.csv")   %>% mutate(pollutant = "O3"),
  read.csv("data/models_verification/fairmode_exercise_phase2/skill_scores_phase2_PM25.csv") %>% mutate(pollutant = "PM25")
) %>%
  filter(!grepl("xv$", model))

# ---- CRITERIO 0: top 5 per IOA in fase 2 ----
top5_ph2 <- scores_ph2 %>%
  group_by(pollutant) %>%
  slice_max(IOA, n = 5) %>%
  select(pollutant, model, IOA, RMSE, correlation) %>%
  arrange(pollutant, desc(IOA))
print(top5_ph2)

top5_ph2_compact <- top5_ph2 %>% ungroup() %>%
  summarise(BCM = paste(model, collapse = " "), .by = pollutant)
print(top5_ph2_compact)


# ---- SELEZIONE FINALE ----
pollutants <- c("NO2", "O3", "PM25")

selection <- lapply(setNames(pollutants, pollutants), function(pol) {
  ph1_apriori     <- strsplit(top_upset$BCM[top_upset$pollutant == pol],     " ")[[1]]
  ph1_aposteriori <- strsplit(top5_ph1_compact$BCM[top5_ph1_compact$pollutant == pol], " ")[[1]]
  ph2_aposteriori <- strsplit(top5_ph2_compact$BCM[top5_ph2_compact$pollutant == pol], " ")[[1]]
  list(
    ph1_apriori     = ph1_apriori,
    ph1_aposteriori = ph1_aposteriori,
    ph2_aposteriori = ph2_aposteriori,
    all             = unique(c(ph1_apriori, ph1_aposteriori, ph2_aposteriori))
  )
})

for (pol in pollutants) {
  cat("\n====", pol, "====\n")
  cat("ph1_apriori    :", selection[[pol]]$ph1_apriori,     "\n")
  cat("ph1_aposteriori:", selection[[pol]]$ph1_aposteriori, "\n")
  cat("ph2_aposteriori:", selection[[pol]]$ph2_aposteriori, "\n")
  cat("ALL (no dups)  :", selection[[pol]]$all,             "\n")
}


# ---- COPIA NETCDF SELEZIONATI ----
nc_dir  <- "data/fairmode-wg5-exercise-phase2-output"
out_dir <- "output/fairmode_phase2"
years   <- c(2015, 2022, 2023, 2024)

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (pol in pollutants) {
  for (bcm in selection[[pol]]$all) {
    for (yr in years) {
      src <- file.path(nc_dir, glue("Scen_{yr}_ITAWG_{pol}_{bcm}_CORR_YEARLY.nc"))
      if (!file.exists(src)) {
        cat("MISSING:", basename(src), "\n")
        next
      }
      file.copy(src, file.path(out_dir, basename(src)))
    }
  }
}
cat("\nCopia completata in", out_dir, "\n")