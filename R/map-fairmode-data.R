

# Function to create the leaflet map with a layer control widget and a binned color key (quantiles)
create_leaflet_map <- function(parameter) {
  library(terra)
  library(dplyr)
  library(leaflet)
  library(glue)
  
  parameter_description <- case_when(
    parameter == "NO2" ~ "nitrogen dioxide",
    parameter == "O3" ~ "ozone",
    parameter == "PM25" ~ "PM2.5"
  )

    # Read all the data
  source("/u/arpa/bonafeg/src/scenario-unbias/R/read-fairmode-data.R")
  data_list  <- read_data(parameter)

    # Combine values from all sources 
  all_values <- c(
    values(data_list$base_case),
    values(data_list$scenario),
    data_list$observed_data$value
  )
  all_values <- na.omit(all_values)
  
  # Define a color palette based on quantiles
  bins <- signif(quantile(data_list$observed_data$value, 
                 probs = seq(0, 1, length.out = 8), na.rm = TRUE), 2)
  bins[1] <- min(all_values)
  bins[length(bins)] <- max(all_values)
  color_palette <- colorBin(
    palette = "Spectral", reverse = T,
    domain = all_values, 
    bins = bins,
    na.color = "transparent"
  )
  
  # Start building the map
  map <- leaflet() %>%
    addProviderTiles("CartoDB.Positron") %>%
    
    # Add base case raster layer (as an image overlay)
    addRasterImage(
      data_list$base_case, 
      colors = color_palette, 
      opacity = 0.6, 
      project = TRUE, 
      group = "Base Case"
    ) %>%
    
    # Add scenario raster layer (as an image overlay)
    addRasterImage(
      data_list$scenario, 
      colors = color_palette, 
      opacity = 0.6, 
      project = TRUE, 
      group = "Scenario"
    ) %>%
    
    # Add observed data points as circles, colored by value
    addCircleMarkers(
      data = data_list$observed_data, 
      lat = ~y, 
      lng = ~x, 
      weight = 1,         # Border weight
      radius = 3,         # Point size
      opacity = 1,        # Border opacity
      color = "black",    # Border color
      fillColor = ~color_palette(value), 
      fillOpacity = 0.8, 
      group = "Observed Data", 
      popup = ~glue("{parameter_description}: {round(value, 1)}")
    ) %>%
    
    # Add Layers Control to toggle between base case, scenario, and observed data
    addLayersControl(
      overlayGroups = c("Observed Data"), 
      baseGroups = c("Base Case", "Scenario"),
      options = layersControlOptions(collapsed = FALSE)
    ) %>%
    
    # Add a color key (legend) to the map
    addLegend(
      position = "bottomright",
      pal = color_palette,
      values = all_values, 
      title = parameter_description,
      opacity = 1
    )
  
  # Return the map object
  return(map)
}

# Example
if(FALSE) {
  map <- create_leaflet_map("NO2")
  map
  map <- create_leaflet_map("PM25")
  map
  map <- create_leaflet_map("O3")
  map
}
