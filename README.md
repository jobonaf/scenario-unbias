# scenario-unbias <img src='inst/unbias2.png' align="right" height="138" />

The repository contains code designed to test different approaches for reducing bias in air quality (AQ) scenarios, specifically within the context of the FAIRMODE WG5 exercise. While the primary focus is on this particular exercise, the methods and tools provided can also be useful for unbiasing AQ scenarios more generally, even outside of the FAIRMODE framework. 

For an overview of the methodology, refer to the [Unbiasing Air Quality Scenarios Documentation](https://github.com/jobonaf/scenario-unbias/blob/main/docs/unbiasing-aq-scenarios.md).


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
|process-fairmode-data.R|read and process FAIRMODE data|

## Data

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

---

## Getting Started

### Clone the Repository

To clone the repository, use the following command:

```bash
git clone https://github.com/jobonaf/scenario-unbias.git
```

### Download as ZIP

If you prefer to download the repository as a ZIP file, follow these steps:  
1. Visit the repository page: [scenario-unbias](https://github.com/jobonaf/scenario-unbias)  
2. Click the green **Code** button.  
3. Select **Download ZIP** and extract the contents to your desired location.

## Processing FAIRMODE Data

To process all the FAIRMODE data and apply various unbiasing methods, use the `R/process-fairmode-data.R` script. This script allows you to specify pollutants, unbiasing sequences, calibration methods, correction algorithms, and spatialization techniques through command-line options.

### Example Usage

```bash
Rscript R/process-fairmode-data.R
```

By default, the script processes the following:
- **Pollutants:** NO₂, O₃, PM₂.₅
- **Unbias sequences:** SCA, CSA, CAS, CA
- **Calibration methods:** Point, Grid, Cell
- **Correction algorithms:** Add, Mult, Lin
- **Spatialization methods:** TPS, IDW, OK, KED
- **Output directory:** `data/processed`

### Customizing the Processing

You can customize the parameters by passing options to the script. Below is an example of how to specify different values:

```bash
Rscript R/process-fairmode-data.R \
  --pollutants=NO2,O3 \
  --output_dir=custom_output \
  --unbias_sequences=SCA,CSA \
  --calibration_methods=Point,Grid \
  --correction_algorithms=Add,Mult \
  --spatialization_methods=tps,idw
```

### Command-Line Options

To see all available options, run the help command:

```bash
Rscript R/process-fairmode-data.R -h
```

#### Available Options

```text
Usage: R/process-fairmode-data.R [options]

Options:
        -p POLLUTANTS, --pollutants=POLLUTANTS
                Comma-separated list of pollutants to process [default: NO2,O3,PM25]

        -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
                Output directory for processed data [default: data/processed]

        -u UNBIAS_SEQUENCES, --unbias_sequences=UNBIAS_SEQUENCES
                Comma-separated list of unbias sequences [default: SCA,CSA,CAS,CA]

        -c CALIBRATION_METHODS, --calibration_methods=CALIBRATION_METHODS
                Comma-separated list of calibration methods [default: Point,Grid,Cell]

        -a CORRECTION_ALGORITHMS, --correction_algorithms=CORRECTION_ALGORITHMS
                Comma-separated list of correction algorithms [default: Add,Mult,Lin]

        -s SPATIALIZATION_METHODS, --spatialization_methods=SPATIALIZATION_METHODS
                Comma-separated list of spatialization methods [default: tps,idw,ok,ked]

        -h, --help
                Show this help message and exit
```

### Output

Processed data will be saved in the specified output directory. By default, this is `data/processed`.
This flexible script allows for batch processing and comparison of multiple unbiasing methods, ensuring comprehensive analysis of the FAIRMODE data.
