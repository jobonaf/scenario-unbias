library(terra)
library(leaflet)

# read the data
r18 <- rast("data/pm25_avg_18.tif")
r22 <- rast("data/pm25_avg_22.tif")

# focus on Po Valley, delta 2022 vs 2018
e <- ext(4050000, 4650000, 2300000, 2600000)
rc18 <- crop(r18, e)
rc22 <- crop(r22, e)
rcd <- rc22-rc18

# palette
ceil <- values(rcd) %>% na.omit() %>% abs() %>% max() %>% ceiling() 
pal <- colorNumeric(palette = "RdBu", domain=c(-ceil, ceil), reverse = T)

# plot the map
plot(rcd, range=c(-ceil, ceil), col=pal(seq(-ceil,ceil,.1)))

