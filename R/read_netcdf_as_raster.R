library(ncdf4)

# Function to read a NetCDF file as SpatRaster, adjusting extent if needed
read_nc <- function(file, loncode="lon", latcode="lat") {
  nc_data <- nc_open(file)
  var_data <- ncvar_get(nc_data, nc_data$var[[1]]$name)
  lon <- ncvar_get(nc_data, loncode)
  lat <- ncvar_get(nc_data, latcode)
  res_x <- mean(diff(lon))
  res_y <- mean(diff(lat))
  xmin_new <- min(lon) - res_x / 2
  xmax_new <- max(lon) + res_x / 2
  ymin_new <- min(lat) - res_y / 2
  ymax_new <- max(lat) + res_y / 2
  raster_data <- rast(t(var_data)[ncol(var_data):1,], crs="EPSG:4326")
  ext(raster_data) <- c(xmin_new, xmax_new, ymin_new, ymax_new)
  nc_close(nc_data)
  return(raster_data)
}
