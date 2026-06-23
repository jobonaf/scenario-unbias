# =============================================================================
# title          :compare_clustering
# description    :Compare two different clustering approaches for air quality scenarios
#                 Calculates Adjusted Rand Index, confusion matrices, and creates
#                 synthetic summary tables
# author         :Giovanni Bonafe'
# date           :20250130
# version        :1.0
# notes          :Requires mclust, ggplot2, dplyr, tidyr, ggrepel
# R_version      :3.5.2
# =============================================================================

# Load required libraries
library(mclust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

# =============================================================================
# Helper Functions
# =============================================================================

# Funzione per caricare e preparare i dati
carica_clustering <- function(file1, file2, inquinante) {
  # Carica primo clustering
  clust1 <- read.csv(file1, stringsAsFactors = FALSE)
  
  # Carica secondo clustering  
  clust2 <- read.csv(file2, stringsAsFactors = FALSE)
  
  # Estrai identificatori degli scenari
  clust1 <- clust1 %>% 
    select(scenario = member, cluster1 = cluster_index)
  clust2 <- clust2 %>% 
    separate(col = "Scenario", sep = "_", into = c(NA, "scenario")) %>% 
    transmute(scenario, cluster2 = letters[Cluster])
  
  # Unisci i due clustering
  merged <- inner_join(clust1, clust2, by = "scenario")
  merged$inquinante <- inquinante
  
  return(merged)
}

# Funzione per calcolare ARI
calcola_ari <- function(data) {
  ari <- adjustedRandIndex(data$cluster1, data$cluster2)
  return(ari)
}

# Funzione per creare matrice di confusione
crea_confusion_matrix <- function(data) {
  table(Clustering1 = data$cluster1, Clustering2 = data$cluster2)
}

# Funzione per creare CSV sintetico
crea_csv_sintetico <- function(data, inquinante, output_dir = ".") {
  # Raggruppa per cluster e crea lista di scenari
  sintesi <- data %>%
    group_by(cluster1, cluster2) %>%
    summarise(scenari = paste(scenario, collapse = ","), .groups = 'drop') %>%
    arrange(cluster1, cluster2)
  
  # Salva CSV
  output_file <- file.path(output_dir, paste0("clustering_sintesi_", inquinante, ".csv"))
  write.csv(sintesi, output_file, row.names = FALSE)
  
  cat("CSV sintetico salvato:", output_file, "\n")
  
  return(sintesi)
}

# Funzione per trovare lo scenario più rappresentativo di un gruppo
trova_rappresentativo <- function(scenari) {
  # Splitta tutti i nomi in token (separati da . o _)
  tokens_list <- lapply(scenari, function(x) unlist(strsplit(x, "[._]")))
  
  # Per ogni scenario, conta quanti token condivide con gli altri
  overlap_scores <- sapply(seq_along(scenari), function(i) {
    miei_tokens <- tokens_list[[i]]
    # Conta overlap con tutti gli altri scenari
    sum(sapply(tokens_list[-i], function(altri_tokens) {
      length(intersect(miei_tokens, altri_tokens))
    }))
  })
  
  # Trova il massimo overlap
  max_overlap <- max(overlap_scores)
  
  # In caso di parità, pesca a caso
  candidati <- scenari[overlap_scores == max_overlap]
  set.seed(42)
  return(sample(candidati, size = 1))
}

# Funzione per visualizzare il confronto
plot_confronto <- function(data, inquinante, ari_value) {
  # Conta scenari per combinazione di cluster
  count_data <- data %>%
    group_by(cluster1, cluster2) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Prepara dati per etichette: solo gruppi con almeno 4 scenari
  # e scegli lo scenario più rappresentativo
  label_data <- data %>%
    group_by(cluster1, cluster2) %>%
    mutate(count = n()) %>%
    filter(count >= 2) %>%
    summarise(
      scenario = trova_rappresentativo(scenario),
      .groups = 'drop'
    )
  
  # Determina i limiti degli assi per le linee della griglia
  n_cluster1 <- length(unique(data$cluster1))
  n_cluster2 <- length(unique(data$cluster2))
  
  # Crea il plot con punti jitterati
  p <- ggplot(data, aes(x = factor(cluster1), y = factor(cluster2))) +
    # Tile di sfondo colorato per densità
    geom_tile(data = count_data, aes(fill = count), 
              color = "transparent", size = 1, alpha = 0.6) +
    # Punti per ogni scenario con jitter
    geom_point(position = position_jitter(width = 0.25, height = 0.25, seed = 42),
               size = 1.5, color = "grey40", alpha = 0.6) +
    # Linee della griglia che separano le celle
    geom_vline(xintercept = seq(0.5, n_cluster1 + 0.5, 1), 
               color = "gray70", size = 0.5) +
    geom_hline(yintercept = seq(0.5, n_cluster2 + 0.5, 1), 
               color = "gray70", size = 0.5) +
    # Etichette con ggrepel solo per gruppi grandi
    ggrepel::geom_text_repel(data = label_data, 
                             aes(label = scenario),
                             size = 2.8,
                             max.overlaps = 20,
                             box.padding = 0.6,
                             point.padding = 0.4,
                             segment.size = 0.2,
                             segment.alpha = 0.5,
                             min.segment.length = 0.2) +
    scale_fill_gradientn(colours=alpha(RColorBrewer::brewer.pal(9,"YlGnBu")[1:7],0.6),
                         breaks = scales::pretty_breaks(n = 5)) +
    labs(
      title = paste("Clustering comparison -", inquinante),
      subtitle = paste("ARI =", round(ari_value, 3)),
      x = "Jaccard clustering",
      y = "Clustering with PCA",
      fill = "No. of BCMs"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(ratio = 1)  +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(p)
}

# =============================================================================
# Main Analysis
# =============================================================================

cat("\n")
cat("=============================================================================\n")
cat("CLUSTERING COMPARISON ANALYSIS\n")
cat("=============================================================================\n\n")

# Create output directory if it doesn't exist
if (!dir.exists("data/clusters-comparison")) {
  dir.create("data/clusters-comparison", recursive = TRUE)
}

inquinanti <- c("O3", "PM2.5", "NO2")
risultati <- list()

for (inq in inquinanti) {
  cat("\n========================================\n")
  cat("Analisi per", inq, "\n")
  cat("========================================\n")
  
  # Definisci percorsi file
  file_standard <- paste0("data/clustering/clustering_", inq, ".csv")
  file_pca <- paste0("data/clustering-with-pca/", inq, "/cluster_scenari.csv")
  
  # Verifica esistenza file
  if (!file.exists(file_standard)) {
    cat("ATTENZIONE: File non trovato:", file_standard, "\n")
    next
  }
  if (!file.exists(file_pca)) {
    cat("ATTENZIONE: File non trovato:", file_pca, "\n")
    next
  }
  
  # Carica dati
  dati <- carica_clustering(file_standard, file_pca, inq)
  
  # Calcola ARI
  ari <- calcola_ari(dati)
  cat("\nAdjusted Rand Index:", round(ari, 4), "\n")
  
  # Crea matrice di confusione
  conf_matrix <- crea_confusion_matrix(dati)
  cat("\nMatrice di Confusione:\n")
  print(conf_matrix)
  
  # Statistiche aggiuntive
  cat("\nN. cluster standard:", length(unique(dati$cluster1)), "\n")
  cat("N. cluster PCA:", length(unique(dati$cluster2)), "\n")
  cat("N. scenari totali:", nrow(dati), "\n\n")
  
  # Crea CSV sintetico
  sintesi <- crea_csv_sintetico(dati, inq, output_dir = "data/clusters-comparison")
  
  # Salva risultati
  risultati[[inq]] <- list(
    dati = dati,
    ari = ari,
    confusion_matrix = conf_matrix,
    sintesi = sintesi,
    plot = plot_confronto(dati, inq, ari)
  )
}

# =============================================================================
# Visualization
# =============================================================================

cat("\n========================================\n")
cat("Creazione grafici...\n")
cat("========================================\n\n")

# Mostra tutti i plot
for (inq in names(risultati)) {
  print(risultati[[inq]]$plot)
}

# Crea plot comparativo degli ARI
ari_values <- data.frame(
  Inquinante = names(risultati),
  ARI = sapply(risultati, function(x) x$ari)
)

p_ari <- ggplot(ari_values, aes(x = Inquinante, y = ARI)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.6) +
  geom_text(aes(label = round(ARI, 3)), vjust = -0.5, size = 5) +
  ylim(0, 1) +
  labs(
    title = "Adjusted Rand Index Comparison",
    y = "Adjusted Rand Index",
    x = "Pollutant"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

print(p_ari)

# =============================================================================
# Export Results
# =============================================================================

cat("\n========================================\n")
cat("Esportazione risultati...\n")
cat("========================================\n\n")

# Salva ARI in un file
write.csv(ari_values, "data/clusters-comparison/risultati_ari.csv", row.names = FALSE)
cat("File salvato: data/clusters-comparison/risultati_ari.csv\n")

# Salva i plot
for (inq in names(risultati)) {
  filename <- paste0("data/clusters-comparison/confronto_clustering_", inq, ".pdf")
  ggsave(
    filename = filename,
    plot = risultati[[inq]]$plot,
    width = 6.5,
    height = 6
  )
  cat("Grafico salvato:", filename, "\n")
}

ggsave(
  filename = "data/clusters-comparison/confronto_ari_totale.pdf",
  plot = p_ari,
  width = 8,
  height = 6
)
cat("Grafico salvato: data/clusters-comparison/confronto_ari_totale.pdf\n")

cat("\n")
cat("=============================================================================\n")
cat("ANALISI COMPLETATA!\n")
cat("=============================================================================\n")
cat("File creati:\n")
cat("  - CSV sintetici: data/clusters-comparison/clustering_sintesi_*.csv\n")
cat("  - ARI values: data/clusters-comparison/risultati_ari.csv\n")
cat("  - Grafici: data/clusters-comparison/confronto_*.pdf\n")
cat("=============================================================================\n\n")

# =============================================================================
# End of script
# =============================================================================