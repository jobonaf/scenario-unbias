# scenario-unbias <img src='inst/unbias2.png' align="right" height="138" />

The repository contains code designed to test different approaches for reducing bias in air quality (AQ) scenarios, specifically within the context of the FAIRMODE WG5 exercise. While the primary focus is on this particular exercise, the methods and tools provided can also be useful for unbiasing AQ scenarios more generally, even outside of the FAIRMODE framework. 

## FAIRMODE WG5 Exercise

The FAIRMODE WG5 bias projection exercise aims to benchmark methodologies for removing bias from air quality model simulations, particularly in future policy scenarios. It focuses on deriving bias correction fields and projecting biases into future scenarios using synthetic datasets. Participants apply their preferred methodologies to post-process provided data and benchmark their results against known synthetic truths. The exercise includes annual data for PM2.5, NO2, and O3 in gridded (NetCDF) and point (CSV) formats for both reference and future projections.

## R code

|script|description|
|------|-----------|
|prepare-aq-obs-data.R|prepare air quality data for script testing|
|check-aq-raster.R|check the prepared AQ data|
|unbias-aq-scenario.R|function to process the unbiasing of a scenario|
|spatialize-points-to-grid.R|spatialize sparse data (or coefficients) to a grid|
|test_spatialization.R|test and plot results of the spatialization|
|calibrate-unbias-coefficient.R|calibrate the coefficients for correction, comparing base case vs observed data|
|apply-unbiasing.R|apply the correction coefficients to a scenario|
|read-fairmode-data.R|read the dataset for FAIRMODE WG5 exercise|
|map-fairmode-data.R|plot on a map FAIRMODE data|
|dashboard-fairmode-data.Rmd|show on a dashboard maps of FAIRMODE data|

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
