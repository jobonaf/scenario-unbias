# scenario-unbias
removing bias from air quality scenarios

## R code

|script|description|
|------|-----------|
|prepare-aq-obs-data.R|prepare air quality data for script testing|
|check-aq-raster.R|check the prepared AQ data|
|unbias-aq-scenario.R|function to process the unbiasing of a scenario|
|interpolate-points-to-grid.R|interpolate sparse data (or coefficients) to a grid|
|test_interpolation.R|test and plot results of the interpolation|
|calibrate-unbias-coefficient.R|calibrate the coefficients for correction, comparing base case vs observed data|
|apply-unbiasing.R|apply the correction coefficients to a scenario|

## data

|file name|description|source|
|---------|-----------|------|
|pm25_avg_18.tif|PM2.5 2018 annual mean, interpolated on 1km grid|downloaded from https://www.eea.europa.eu/en/datahub/datahubitem-view/938bea70-07fc-47e9-8559-8a09f7f92494|
|pm25_avg_22.tif|PM2.5 2022 annual mean, interpolated on 1km grid|downloaded from https://www.eea.europa.eu/en/datahub/datahubitem-view/938bea70-07fc-47e9-8559-8a09f7f92494|
|eea_obs_pm25_bkg_2022.csv|PM2.5 2022 annual mean, measured in European background stations|prepared with prepare-aq-obs-data.R|
|eea_obs_pm25_bkg_2023.csv|PM2.5 2023 annual mean, measured in European background stations|prepared with prepare-aq-obs-data.R|
