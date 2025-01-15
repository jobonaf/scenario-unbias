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

The dataset, located in `data/fairmode-wg5-exercise-202501/`, supports FAIRMODE WG5 exercise for AQ scenario unbiasing. It contains annual data for NO₂, O₃, and PM₂.₅ (humidity-adjusted at 50%) in gridded and point formats.

### `YEARLY`
#### 1. `BaseCase_Perturbed_Gridded`
Gridded NetCDF files for perturbed base case:

- `BaseCase_PERT_NO2_YEARLY.nc`: Annual NO₂ (µg/m³).
- `BaseCase_PERT_O3_YEARLY.nc`: Annual O₃ (ppb).
- `BaseCase_PERT_PM25_rh50_YEARLY.nc`: Annual PM₂.₅ (µg/m³).

#### 2. `BaseCase_Reference_Points`
Point-based reference data in CSV:

- `yearly_SURF_ppb_O3.csv`: O₃ (ppb).
- `yearly_SURF_ug_NO2.csv`: NO₂ (µg/m³).
- `yearly_SURF_ug_PM25_rh50.csv`: PM₂.₅ (µg/m³).

#### 3. `Scenario_Perturbed_Gridded`
Gridded NetCDF files for perturbed scenario:

- `SCEN_PERT_NO2_YEARLY.nc`: Annual NO₂ (µg/m³).
- `SCEN_PERT_O3_YEARLY.nc`: Annual O₃ (ppb).
- `SCEN_PERT_PM25_rh50_YEARLY.nc`: Annual PM₂.₅ (µg/m³).

**Formats:** NetCDF for gridded data, CSV for point data.
