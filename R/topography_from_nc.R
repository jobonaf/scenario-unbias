# Function to load topography from NetCDF
topography_from_nc <- function(nc_file, varname = "topography", grid) {
  nc <- nc_open(nc_file)
  topo_vals <- ncvar_get(nc, varname)
  nc_close(nc)
  topo_rast <- rast(t(topo_vals))
  ext(topo_rast) <- ext(grid)
  return(topo_rast)
}
