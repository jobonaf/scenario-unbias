# Load required package
library(terra)

# Function to extract raster metadata
get_tiff_metadata <- function(file_list) {
  # Check if file_list is empty
  if (length(file_list) == 0) {
    stop("The file list is empty.")
  }
  
  # Initialize an empty list to store results
  metadata_list <- lapply(file_list, function(file) {
    if (!file.exists(file)) {
      warning(paste("File not found:", file))
      return(NULL)
    }
    
    # Load the raster
    raster_data <- rast(file)
    
    # Extract metadata
    data.frame(
      file = file,
      ncol = ncol(raster_data),
      nrow = nrow(raster_data),
      xmin = ext(raster_data)$xmin,
      xmax = ext(raster_data)$xmax,
      ymin = ext(raster_data)$ymin,
      ymax = ext(raster_data)$ymax,
      res_x = res(raster_data)[1],
      res_y = res(raster_data)[2]
    )
  })
  
  # Remove NULL entries (files that were not found)
  metadata_list <- Filter(Negate(is.null), metadata_list)
  
  # Combine results into a single data.frame
  metadata_df <- do.call(rbind, metadata_list)
  
  return(metadata_df)
}

# Example usage:
file_list <- Sys.glob("data/processed/*.tif")
metadata_table <- get_tiff_metadata(file_list)
rownames(metadata_table) <- NULL
View(metadata_table)
