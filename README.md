# scenario-unbias <img src='inst/unbias2.png' align="right" height="138" />

The repository contains code designed to test different approaches for reducing bias in air quality (AQ) scenarios, specifically within the context of the FAIRMODE WG5 exercise. While the primary focus is on this particular exercise, the methods and tools provided can also be useful for unbiasing AQ scenarios more generally, even outside of the FAIRMODE framework.

## Repository Structure

```
scenario-unbias/
├── R/                      # Core R scripts
├── data/                   # Data directory (input/output)
├── inst/                   # Installation-related files
├── docs/                   # Generated documentation
└── README.md              # This file
```

## FAIRMODE WG5 Exercise – Phase 1 (Synthetic Data)

The FAIRMODE WG5 bias projection exercise aims to benchmark methodologies for removing bias from air quality model simulations, particularly in future policy scenarios. It focuses on deriving bias correction fields and projecting biases into future scenarios using synthetic datasets. Participants apply their preferred methodologies to post-process provided data and benchmark their results against known synthetic truths. The exercise includes annual data for PM2.5, NO2, and O3 in gridded (NetCDF) and point (CSV) formats for both reference and future projections.

## Italian National Exercise

A complementary national exercise has been established to test and benchmark unbias methodologies specifically for Italian air quality planning. Mirroring the structure of the FAIRMODE WG5 Phase 1 exercise, it provides a synthetic dataset focused on the Italian domain, featuring annual data for PM10, PM2.5, NO₂, and O₃.

This initiative aims to evaluate the performance of different correction techniques within the complex orography and diverse pollution climates characteristic of Italy. The data structure and formats (NetCDF for gridded model data, CSV for point observations) are identical to the FAIRMODE exercise, ensuring methodological consistency and allowing for direct comparison of results between the European and national contexts. The provided tools and scripts in this repository are fully compatible for processing this additional dataset.

## FAIRMODE WG5 Exercise – Phase 2 (Real-World Data)

A follow-up intercomparison exercise based on real-world observations (EEA stations) and EMEP simulations for years 2015, 2022, 2023, and 2024. Participants define a bias correction using 2015 data and project it to future years, validating against observed concentrations. Daily and annual data are available for PM₂.₅, NO₂, and O₃.

## Installation

### Prerequisites

- **R version:** 3.5.2 or higher (tested with 3.5.2)

### Required R Packages

Based on the actual scripts in this repository:

```r
# Core spatial and data manipulation
install.packages(c(
  "raster",       # Spatial raster data handling
  "sp",           # Spatial data classes
  "ncdf4",        # NetCDF file I/O
  "rgdal",        # GDAL bindings for raster/vector I/O
  "dplyr",        # Data manipulation
  "tidyr",        # Data tidying
  "readr"         # CSV reading
))

# Spatial interpolation
install.packages(c(
  "fields",       # Thin plate spline (TPS)
  "gstat"         # IDW, OK, KED
))

# Visualization (optional, for dashboards and maps)
install.packages(c(
  "leaflet",      # Interactive maps
  "shiny",        # Dashboard framework
  "rmarkdown",    # Rendering RMarkdown dashboards
  "ggplot2",      # Static plots
  "RColorBrewer"  # Color palettes
))

# Clustering and analysis
install.packages(c(
  "cluster"       # Clustering algorithms
))
```

### Clone the Repository

```bash
git clone https://github.com/jobonaf/scenario-unbias.git
cd scenario-unbias
```

## Code Overview

### Data Ingestion
| Script | Description |
|--------|-------------|
| `read-fairmode-data.R` | Reads and preprocesses the dataset for the FAIRMODE WG5 exercise. |
| `read-italian-data.R` | Reads and preprocesses the dataset for the Italian case study. |
| `read_netcdf_as_raster.R` | Imports a NetCDF file and converts it to a Raster* object. |

### Core Unbiasing Methods (Based on the BCM Classification Framework)

| Script | Description |
|--------|-------------|
| `unbias-aq-scenario.R` | Core algorithm implementing bias correction methodologies (BCM) according to the four‑dimensional classification framework (sequence, adjustment, calibration strategy, spatialization). |
| `spatialize-points-to-grid.R` | Spatializes point data (correction coefficients) to a continuous grid using selected spatialization approaches. |
| `calibrate-unbias-coefficients.R` | Calibrates correction coefficients based on a chosen calibration strategy (per station, per grid cell, or global) by comparing a base scenario to observed data. |
| `apply-unbiasing.R` | Applies the pre‑calibrated correction coefficients to a target scenario using a specified type of adjustment. |
| `process-fairmode-data.R` | End‑to‑end workflow that applies the full BCM framework: reads data, calibrates, spatializes (when applicable), and applies adjustment according to the defined sequence of operations. |

### Visualization & Dashboards
| Script | Description |
|--------|-------------|
| `map-fairmode-data.R` | Generates static maps of the FAIRMODE data. |
| `dashboard-fairmode-data.Rmd` | Creates an interactive dashboard for visualizing FAIRMODE data on maps. |
| `dashboard-italian-data.Rmd` | Creates an interactive dashboard for visualizing the Italian case study data. |

### Evaluation & Analysis
| Script | Description |
|--------|-------------|
| `scenario_boxplot.R` | Summarizes and compares unbiasing output results using boxplots. |
| `clustering_output.R` | Performs cluster analysis on the model output data. |
| `distance_scenarios.R` | Calculates the Jaccard distance matrix between scenarios for clustering. |
| `describe_clusters.R` | Post-processing and statistical description of identified clusters. |
| `model_verification.R` | Verifies unbiasing performance by comparison against a gridded reference field. Also suitable for general AQ model evaluation. |

### Utilities
| Script | Description |
|--------|-------------|
| `test_spatialization.R` | Tests and visualizes the results of the spatialization procedure. |
| `tiff2netcdf.R` | Converts model output from GeoTIFF to NetCDF format. |
| `compare-rasters.R` | Compares two raster-based scenarios pixel-by-pixel. |
| `identify_homogeneous_zones.R` | Identifies homogeneous zones based on spatial patterns of annual mean concentrations. |

## FAIRMODE WG5 Exercise Data Structure

The dataset supports the FAIRMODE WG5 exercise for AQ scenario unbiasing. It contains **annual** data for NO₂, O₃, and PM₂.₅ (humidity-adjusted at 50%) in gridded and point formats.

### Gridded Data (NetCDF)
- **Base Case (Perturbed):** Simulated baseline scenario with bias.
  - `BaseCase_PERT_NO2.nc`: Annual NO₂ (µg/m³)
  - `BaseCase_PERT_O3.nc`: Annual O₃ (ppb)
  - `BaseCase_PERT_PM25_rh50.nc`: Annual PM₂.₅ (µg/m³)
- **Future Scenario (Perturbed):** Projected AQ scenario with bias.
  - `SCEN_PERT_NO2.nc`: Annual NO₂ (µg/m³)
  - `SCEN_PERT_O3.nc`: Annual O₃ (ppb)
  - `SCEN_PERT_PM25_rh50.nc`: Annual PM₂.₅ (µg/m³)

### Point-Based Data (CSV)
- **Reference Observations:** Surface measurements used for bias correction.
  - `yearly_SURF_ppb_O3.csv`: Annual O₃ (ppb)
  - `yearly_SURF_ug_NO2.csv`: Annual NO₂ (µg/m³)
  - `yearly_SURF_ug_PM25_rh50.csv`: Annual PM₂.₅ (µg/m³)

## BCM Framework Details

The implementation follows the classification proposed by Bonafè et al. (2025) with four key dimensions:

| Dimension | Options Implemented | Description |
|-----------|---------------------|-------------|
| **Sequence of operations** | `SCA`, `CSA`, `CAS`, `CA` | Order of Spatialization, Calibration, and Application steps. |
| **Type of adjustment** | `Add` (additive), `Mult` (multiplicative), `Lin` (linear regression) | Mathematical form of the correction applied to model outputs. |
| **Calibration strategy** | `Point` (per station), `Grid` (per grid cell), `Cell` (cell‑by‑cell) | Level at which correction coefficients are estimated. |
| **Spatialization approach** | `TPS` (thin‑plate spline), `OK` (ordinary kriging), `KED` (kriging with external drift), `IDW` (inverse distance weighting) | Method used to interpolate point‑based coefficients to the grid. |

These dimensions can be combined flexibly to define a wide range of bias correction methodologies, allowing systematic comparison and reproducibility.

## Getting Started

### Quick Start Example

Here's a minimal example to process data for a single pollutant:

```r
# Source the main processing script
source("R/process-fairmode-data.R")

# Process NO2 with default settings (SCA sequence, Point calibration,
# Additive adjustment, TPS spatialization)
process_pollutant(
  pollutant = "NO2",
  unbias_seq = "SCA",
  calib_method = "Point",
  corr_algorithm = "Add",
  spat_method = "TPS"
)
```

### Processing FAIRMODE Data

To process all the FAIRMODE data and apply different bias correction methodologies, use `R/process-fairmode-data.R`. The script accepts command‑line arguments to specify the four BCM dimensions.

**Command-line usage:**
```bash
# Process all pollutants with default settings (all combinations)
Rscript R/process-fairmode-data.R

# Process only NO₂ using specific dimensions
Rscript R/process-fairmode-data.R --pollutant NO2 --sequence SCA,CSA --adjustment Add --calibration Point --spatialization TPS,OK
```

**Default processed combinations:**
- **Pollutants:** NO₂, O₃, PM₂.₅
- **Sequences:** SCA, CSA, CAS, CA
- **Adjustment types:** Add, Mult, Lin
- **Calibration strategies:** Point, Grid, Cell
- **Spatialization approaches:** TPS, IDW, OK, KED

### Output Structure

Processed data is saved in the `data/processed/` directory with filenames encoding the four dimensions:

```
{POLLUTANT}_{SEQUENCE}_{CALIBRATION}_{ADJUSTMENT}_{SPATIALIZATION}.nc
```

**Example:** `NO2_SCA_Point_Add_TPS.nc`

### Running Dashboards

To launch the interactive dashboard for the Italian case study:

```r
rmarkdown::run("R/dashboard-italian-data.Rmd")
```

For the FAIRMODE dashboard:
```r
rmarkdown::run("R/dashboard-fairmode-data.Rmd")
```

## How to Contribute

We welcome contributions to this repository! To ensure a smooth collaboration, please follow these steps:

1. **Fork the Repository** (if you don't have write access).
2. **Clone Your Fork:**
   ```bash
   git clone https://github.com/jobonaf/scenario-unbias.git
   cd scenario-unbias
   ```
3. **Create a New Branch from `devel`:**
   ```bash
   git checkout -b devel-<feature_name>-YYYYMMDD origin/devel
   ```
4. **Develop Your Changes** and test them locally.
5. **Commit Your Changes:**
   ```bash
   git add .
   git commit -m "Add feature X or Fix bug Y"
   ```
6. **Push to Your Branch:**
   ```bash
   git push origin devel-<feature_name>-YYYYMMDD
   ```
7. **Open a Pull Request** from your branch to `devel` on GitHub.
8. **Code Review and Merge:** The maintainers will review your pull request and merge it if everything is in order.

## Submitting an Issue

If you encounter a bug, have a question, or want to propose a new feature, please submit an issue:

1. Navigate to the **Issues** tab on GitHub
2. Click **New Issue**
3. Select the appropriate issue template (if available)
4. Provide:
   - A clear title and description
   - Steps to reproduce bugs
   - Relevant code snippets or error messages
   - Your R version and package versions (`sessionInfo()` output)
5. Add appropriate labels
6. Submit the issue

## Citation

If you use this code in your research, please cite:

Bonafè, G., et al. (2025). *A Classification of Bias Correction Methodologies used in Air Quality Scenarios Modelling Applications*. Zenodo. https://doi.org/10.5281/zenodo.17183678

## Contact

For questions or feedback, please open an issue or contact the maintainers.
