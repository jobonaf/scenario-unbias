source("R/spatialize-points-to-grid.R")

# Load necessary libraries
library(terra)
library(fields)
library(gstat)

# Example usage
# Load example data
r <- rast(system.file("ex/elev.tif", package = "terra"))
ra <- aggregate(r, 10)
xy <- as.data.frame(xyFromCell(ra, 1:ncell(ra)))
v <- values(ra)

# Filter valid points
valid_idx <- !is.na(v)
xy <- xy[valid_idx, ]
v <- v[valid_idx]

# Thin Plate Spline interpolation
result_tps <- spatialize(points_data = data.frame(x = xy[, 1], y = xy[, 2], values = v), 
                               grid_data = r, method = "tps")
plot(result_tps, main = "Thin Plate Spline")

# IDW interpolation
result_idw <- spatialize(points_data = data.frame(x = xy[, 1], y = xy[, 2], values = v), 
                               grid_data = r, method = "idw", 
                               params = list(idp = 2, nmax = 5))
plot(result_idw, main = "IDW")

# Ordinary Kriging interpolation
result_ok <- spatialize(points_data = data.frame(x = xy[, 1], y = xy[, 2], values = v), 
                              grid_data = r, method = "ok",
                              params = list(psill = NA, model = "Sph", range = 0.1, nugget = 10))
plot(result_ok, main = "Ordinary Kriging")

# Kriging with External Drift interpolation
result_ked <- spatialize(points_data = data.frame(x = xy[, 1], y = xy[, 2], values = v), 
                               grid_data = r, method = "ked")
plot(result_ked, main = "Kriging with External Drift")
