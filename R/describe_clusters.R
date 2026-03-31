library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(glue)


plot_cluster <- function(poll) {
  clu <- read_csv(glue("data/clustering/clustering_{poll}.csv"))
  
  # Calcolo delle frequenze per le combinazioni richieste
  pdat1 <- clu %>%
    group_by(cluster_index, spatialization, sequence) %>%
    summarise(n = n(), .groups = "drop") %>%
    ungroup() %>%
    group_by(spatialization, sequence) %>%
    mutate(p = n==sum(n)) %>%
    ungroup() %>%
    filter(n > 0) %>%
    mutate(cluster_label = glue("{sprintf('%02i',cluster_index)}"))  # Aggiunge etichetta per i pannelli
  
  pdat2 <- clu %>%
    group_by(cluster_index, calibration, correction) %>%
    summarise(n = n(), .groups = "drop") %>%
    ungroup() %>%
    group_by(calibration, correction) %>%
    mutate(p = n==sum(n)) %>%
    ungroup() %>%
    filter(n > 0) %>%
    mutate(cluster_label = glue("{sprintf('%02i',cluster_index)}"))  # Aggiunge etichetta per i pannelli
  
  # Primo plot: spatialization vs sequence con etichette
  p1 <- ggplot(pdat1, aes(x = spatialization, y = sequence)) +
    geom_point(data = filter(pdat1, p),  # Filtra solo le righe con p == TRUE
               size = 6, 
               color = "firebrick", 
               shape = 21,
               show.legend = FALSE  # Disabilita la legenda per size
    ) +
    geom_text(aes(label = n), color = "black", size = 4) +
    facet_wrap(~ cluster_label, nrow = 2) +  # Usa le etichette personalizzate
    labs(x = "Spatialization Method",
         y = "Sequence") +
    theme_bw() +
    theme(panel.spacing = unit(1, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),  # Rimuove lo sfondo dei titoli
          strip.text = element_text(face = "bold"),  # Testo in grassetto
          plot.title = element_blank())  # Rimuove il titolo principale
  
  # Secondo plot: calibration vs correction con etichette
  p2 <- ggplot(pdat2, aes(x = calibration, y = correction)) +
    geom_point(data = filter(pdat2, p),  # Filtra solo le righe con p == TRUE
               size = 6, 
               color = "firebrick", 
               shape = 21,
               show.legend = FALSE  # Disabilita la legenda per size
    ) +
    geom_text(aes(label = n), color = "black", size = 4) +
    facet_wrap(~ cluster_label, nrow = 2) +  # Usa le etichette personalizzate
    labs(x = "Calibration Method",
         y = "Correction Method") +
    theme_bw() +
    theme(panel.spacing = unit(1, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),  # Rimuove lo sfondo dei titoli
          strip.text = element_text(face = "bold"),  # Testo in grassetto
          plot.title = element_blank())  # Rimuove il titolo principale
  
  # Salvataggio dei plot
  pdf(glue("clusters_description_{poll}.pdf"), width = 10, height = 3)
  print(p1)
  print(p2)
  dev.off()
}

plot_cluster("NO2")
plot_cluster("PM25")
plot_cluster("O3")
