# =============================================================================
# title          :identify_homogeneous_zones
# description    :Identify homogeneous and compact zones based on annual mean
#                 concentrations of PM2.5, PM10, NO2, O3 
# author         :Giovanni Bonafe'
# date           :20250904
# version        :1.0
# notes          :Requires NetCDF files with air pollution data
# R_version      :3.5.2
# =============================================================================

# Load required libraries
library(terra)         # For raster operations
library(sf)            # For vector operations
library(futile.logger) # For logging
library(glue)          # For string interpolation
library(ncdf4)         # For NetCDF handling
library(tidyr)
library(ggplot2)

# Initialize logging (stdout only)
flog.threshold(INFO)
flog.appender(appender.console())  # scrive su standard output
flog.info("Starting homogeneous zones identification analysis")

# Define file paths
base_path <- "data/esercizio-dominio-italiano/BaseCase_Reference_Gridded/"
italy_boundary_path <- "/atlas/arpa/bonafeg/data/geo/LimitiAmministrativi/Italy_WGS84_LatLong/ITA_adm0.shp"

file_pm25 <- paste0(base_path, "SCEN_BASE_c_PM25_YEARLY.nc")
file_pm10 <- paste0(base_path, "SCEN_BASE_c_PM10_YEARLY.nc")
file_no2  <- paste0(base_path, "SCEN_BASE_c_NO2_YEARLY.nc")
file_o3   <- paste0(base_path, "SCEN_BASE_c_O3_YEARLY.nc")

flog.info(glue("Base path: {base_path}"))
flog.info(glue("Italy boundary path: {italy_boundary_path}"))

# Function to read a NetCDF file as SpatRaster, adjusting extent if needed
source("R/read_netcdf_as_raster.R")

flog.info("Loading pollution data...")
rast_pm25 <- read_nc(file_pm25)
rast_pm10 <- read_nc(file_pm10)
rast_no2  <- read_nc(file_no2)
rast_o3   <- read_nc(file_o3)

flog.info("Loading Italy boundary shapefile...")
italy_boundary <- st_read(italy_boundary_path, quiet = TRUE) %>% vect()
flog.info(glue("Italy boundary CRS: {crs(italy_boundary)}"))

# Data preprocessing
names(rast_pm25) <- "PM25"
names(rast_pm10) <- "PM10"
names(rast_no2)  <- "NO2"
names(rast_o3)   <- "O3"
raster_stack <- c(rast_pm25, rast_pm10, rast_no2, rast_o3)

# Crop and mask to Italy boundary
flog.info("Cropping and masking to Italy boundary...")
raster_stack_italy <- crop(raster_stack, italy_boundary)
raster_stack_italy <- mask(raster_stack_italy, italy_boundary)
flog.info(glue("Final raster stack dimensions: {nrow(raster_stack_italy)} x {ncol(raster_stack_italy)}"))

# Data cleaning
flog.info("Removing NA values...")
values_df <- as.data.frame(raster_stack_italy, na.rm = TRUE)
flog.info(glue("Number of valid pixels: {nrow(values_df)}"))
if(any(!is.finite(as.matrix(values_df)))) {
  flog.warn("Non-finite values detected, removing...")
  values_df <- values_df[apply(values_df, 1, function(x) all(is.finite(x))), ]
  flog.info(glue("Remaining valid pixels after cleaning: {nrow(values_df)}"))
}

# Variable standardization
flog.info("Standardizing variables...")
values_scaled <- scale(values_df)

# Principal Component Analysis
flog.info("Performing Principal Component Analysis...")
pca_result <- prcomp(values_scaled, center = FALSE, scale. = FALSE)
pca_variance <- summary(pca_result)$importance[2, ]
flog.info("Explained variance by principal components:")
for(i in seq_along(pca_variance)) {
  flog.info(glue("PC{i}: {round(pca_variance[i] * 100, 1)}%"))
}
cumulative_var <- cumsum(pca_variance)
n_components <- which(cumulative_var >= 0.90)[1]
flog.info(glue("Number of selected principal components (90% variance): {n_components}"))
pca_scores <- pca_result$x[, 1:n_components]

# K-means clustering
n_clusters <- 5
flog.info(glue("Performing K-means clustering on {n_components} PCA components with {n_clusters} clusters..."))
set.seed(123)
kmeans_result <- kmeans(pca_scores, centers = n_clusters, iter.max = 100, nstart = 25)
flog.info(glue("K-means clustering completed with total within-cluster SS: {round(kmeans_result$tot.withinss, 2)}"))

# Assign clusters to raster
flog.info("Assigning clusters to raster...")
cluster_raster <- rast(raster_stack_italy[[1]])
full_values <- rep(NA, ncell(cluster_raster))
non_na_indices <- which(!is.na(values(raster_stack_italy[[1]])[]))
full_values[non_na_indices] <- kmeans_result$cluster
values(cluster_raster) <- full_values
cluster_raster_categ <- as.factor(cluster_raster)
levels(cluster_raster_categ) <- data.frame(value = 1:n_clusters, label = 1:n_clusters)

# Map clusters
flog.info("Saving cluster map to PDF...")
pdf("cluster_map.pdf", width = 6, height = 6)
plot(cluster_raster_categ,
     col = fields::tim.colors(n_clusters),
     legend = TRUE)
dev.off()
flog.info("Cluster map saved as cluster_map.pdf")

# Save results
flog.info("Saving results...")
zones_file <- "homogeneous_zones.tif"
writeRaster(cluster_raster, zones_file, overwrite = TRUE, datatype = "INT1U")
stack_file <- "original_data_italy.tif"
writeRaster(raster_stack_italy, stack_file, overwrite = TRUE)
pca_file <- "pca_results.rds"
saveRDS(pca_result, pca_file)
flog.info(glue("Output files created: {zones_file}, {stack_file}, {pca_file}"))

# Boxplot analysis of zones
flog.info("Creating boxplots for each pollutant by cluster...")
values_df$cluster <- factor(kmeans_result$cluster)
cluster_colors <- fields::tim.colors(n_clusters)
values_long <- values_df %>%
  pivot_longer(cols = c(PM25, PM10, NO2, O3),
               names_to = "pollutant",
               values_to = "value")
pdf("zone_boxplots.pdf", width = 8, height = 6)
ggplot(values_long, aes(x = cluster, y = value, fill = cluster)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ pollutant, scales = "free_y") +
  scale_fill_manual(values = cluster_colors) +
  #scale_y_continuous(trans="pseudo_log") +
  theme_bw(base_size = 12) +
  labs(x = "cluster",
       y = expression("concentration"~(mu*g/m^3)))
dev.off()
flog.info("Boxplots saved to zone_boxplots.pdf")


# PCA scree plot
flog.info("Creating PCA scree plot...")
pdf("pca_scree_plot.pdf", width = 8, height = 6)
par(mar = c(5, 5, 4, 2))
plot(cumulative_var * 100, type = "b", pch = 19, col = "blue",
     xlab = "Number of Principal Components", 
     ylab = "Cumulative Variance Explained (%)",
     main = "PCA Scree Plot - Cumulative Variance Explained",
     ylim = c(0, 100),
     cex.lab = 1.2, cex.axis = 1.1)
abline(h = 90, col = "red", lty = 2, lwd = 2)
text(x = n_components, y = 92, labels = glue("90% variance\n({n_components} components)"), 
     col = "red", pos = 3)
grid()
dev.off()
flog.info("PCA scree plot saved as pca_scree_plot.pdf")

# Final summary
flog.info(glue("Analysis completed successfully. Identified {n_clusters} homogeneous zones."))
flog.info(glue("Main output file: {zones_file}"))

