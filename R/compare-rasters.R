library(terra)
library(leaflet)
library(ggplot2)
library(ggmap)
library(ncdf4)
library(classInt)
library(RColorBrewer)

# Function to read a NetCDF file as SpatRaster
read_nc <- function(file) {
  nc_data <- nc_open(file)
  var_data <- ncvar_get(nc_data, nc_data$var[[1]]$name)
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat")
  raster_data <- rast(t(var_data)[ncol(var_data):1,], crs="EPSG:4326")
  ext(raster_data) <- c(min(lon), max(lon), min(lat), max(lat))
  nc_close(nc_data)
  return(raster_data)
}

# Function to compare two raster scenarios
display_raster_comparison <- function(file1, file2, title1 = basename(file1), title2 = basename(file2), ...) {
  # Read raster files
  raster1 <- if (grepl(".nc$", file1)) read_nc(file1) else rast(file1)
  raster2 <- if (grepl(".nc$", file2)) read_nc(file2) else rast(file2)
  
  # Ensure rasters have the same extent and resolution
  if (!compareGeom(raster1, raster2, stopOnError = FALSE)) {
    stop("The rasters have different extents or resolutions.")
  }
  
  # Compute the difference
  diff_raster <- raster1 - raster2
  
  # Extract values for visualization
  diff_values <- values(diff_raster, mat = FALSE, na.rm = TRUE)
  
  # Find min and max values and their coordinates
  min_val <- min(diff_values, na.rm = TRUE)
  max_val <- max(diff_values, na.rm = TRUE)
  min_coords <- xyFromCell(diff_raster, which.min(diff_values))
  max_coords <- xyFromCell(diff_raster, which.max(diff_values))
  
  # Create quantile-based color bins
  bins <- classIntervals(diff_values, n = 6, style = "quantile")$brks
  bins <- unique(bins)
  palette <- colorBin("RdBu", domain = diff_values, bins = bins, reverse = TRUE)
  
  # Major European cities with population > 1M and capitals > 200k
  cities <- data.frame(
    name = c("London", "Berlin", "Madrid", "Rome", "Paris", "Bucharest", "Vienna", 
             "Hamburg", "Budapest", "Warsaw", "Barcelona", "Munich", "Milan",
             "Stockholm", "Prague", "Sofia", "Copenhagen", "Brussels", "Amsterdam", 
             "Athens", "Oslo", "Helsinki", "Lisbon", "Dublin"),
    lat = c(51.5074, 52.5200, 40.4168, 41.9028, 48.8566, 44.4268, 48.2082, 53.5511, 
            47.4979, 52.2298, 41.3851, 48.1351, 45.4642,
            59.3293, 50.0755, 42.6975, 55.6761, 50.8503, 52.3676, 37.9838, 59.9139, 
            60.1695, 38.7169, 53.3498),
    lon = c(-0.1278, 13.4050, -3.7038, 12.4964, 2.3522, 26.1025, 16.3738, 9.9937, 
            19.0402, 21.0122, 2.1734, 11.5820, 9.1900,
            18.0686, 14.4378, 23.3242, 12.5683, 4.3517, 4.9041, 23.7275, 10.7522, 
            24.9354, -9.1399, -6.2603)
  )
  
  # Extract raster values at city locations
  city_values <- extract(diff_raster, cities[, c("lon", "lat")], ID=FALSE)
  names(city_values) <- "value"
  
  # Create Leaflet map
  map <- leaflet() %>%
    addTiles() %>%
    addRasterImage(diff_raster, colors = palette, opacity = 0.7) %>%
    addControl(html = paste0("<b>Comparison:</b><br>", title1, " - ", title2), position = "topright") %>%
    addLegend(pal = palette, values = diff_values, title = "Difference", labFormat = labelFormat(...)) %>%
    addCircleMarkers(lng = min_coords[1], lat = min_coords[2], radius = 7, color = "black", fillOpacity = 0, weight = 3, label = paste("Min Value:", round(min_val, 2))) %>%
    addCircleMarkers(lng = max_coords[1], lat = max_coords[2], radius = 7, color = "black", fillOpacity = 0, weight = 3, label = paste("Max Value:", round(max_val, 2))) %>%
    addCircleMarkers(data = cities, ~lon, ~lat, radius = 5, color = "black", fillOpacity = 0, weight = 1, label = ~paste(name, ":", round(city_values$value, 2)))
  
  # Create histogram of differences
  hist_plot <- ggplot(data.frame(Difference = diff_values), aes(x = Difference)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    theme_minimal() +
    labs(title = "Distribution of Differences", x = "Difference", y = "Frequency")
  
  # Create static map using ggplot2 with binned scale
  names(diff_raster) <- "value"
  df_raster <- as.data.frame(diff_raster, xy = TRUE, na.rm = TRUE)
  df_raster$bins <- cut(df_raster$value, breaks = bins, include.lowest = TRUE)
  static_map <- ggplot(df_raster, aes(x = x, y = y, fill = bins)) +
    geom_tile() +
    scale_fill_manual(values = rev(brewer.pal(length(bins) - 1, "RdBu"))) +
    theme_minimal() +
    labs(title = paste("Difference Map:", title1, "vs", title2), fill = "Difference")
  
  # Print results
  print(hist_plot)
  print(static_map)
  return(list(map=map, hist_plot=hist_plot, static_map=static_map))
}
