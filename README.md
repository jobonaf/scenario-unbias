# scenario-unbias <img src='inst/unbias2.png' align="right" height="138" />

The repository contains code designed to test different approaches for reducing bias in air quality (AQ) scenarios, specifically within the context of the FAIRMODE WG5 exercise. While the primary focus is on this particular exercise, the methods and tools provided can also be useful for unbiasing AQ scenarios more generally, even outside of the FAIRMODE framework.

## FAIRMODE WG5 Exercise

The FAIRMODE WG5 bias projection exercise aims to benchmark methodologies for removing bias from air quality model simulations, particularly in future policy scenarios. It focuses on deriving bias correction fields and projecting biases into future scenarios using synthetic datasets. Participants apply their preferred methodologies to post-process provided data and benchmark their results against known synthetic truths. The exercise includes annual data for PM2.5, NO2, and O3 in gridded (NetCDF) and point (CSV) formats for both reference and future projections.

## Italian National Exercise

A complementary national exercise has been established to test and benchmark unbias methodologies specifically for Italian air quality planning. Mirroring the structure of the FAIRMODE WG5 exercise, it provides a synthetic dataset focused on the Italian domain, featuring annual data for PM10, PM2.5, NO₂, and O₃.

This initiative aims to evaluate the performance of different correction techniques within the complex orography and diverse pollution climates characteristic of Italy. The data structure and formats (NetCDF for gridded model data, CSV for point observations) are identical to the FAIRMODE exercise, ensuring methodological consistency and allowing for direct comparison of results between the European and national contexts. The provided tools and scripts in this repository are fully compatible for processing this additional dataset.

## Code Overview

| Group | Script | Description |
| :--- | :--- | :--- |
| **Data Ingestion** | `read-fairmode-data.R` | Reads and preprocesses the dataset for the FAIRMODE WG5 exercise. |
| | `read-italian-data.R` | Reads and preprocesses the dataset for the Italian case study. |
| | `read_netcdf_as_raster.R` | Imports a NetCDF file and converts it to a SpatRaster object. |
| **Core Unbiasing** | `unbias-aq-scenario.R` | Core algorithm for unbiasing air quality concentration scenarios. |
| | `spatialize-points-to-grid.R` | Spatializes sparse point data (or correction coefficients) to a continuous grid. |
| | `calibrate-unbias-coefficients.R` | Calibrates unbiasing coefficients by comparing a base scenario to observed data. |
| | `apply-unbiasing.R` | Applies the pre-calibrated correction coefficients to a target scenario. |
| | `process-fairmode-data.R` | End-to-end workflow: reads, processes, and applies unbiasing methods to the FAIRMODE data. |
| **Visualization** | `map-fairmode-data.R` | Generates static maps of the FAIRMODE data. |
| | `dashboard-fairmode-data.Rmd` | Creates an interactive dashboard for visualizing FAIRMODE data on maps. |
| | `dashboard-italian-data.Rmd` | Creates an interactive dashboard for visualizing the Italian case study data. |
| **Evaluation & Analysis** | `scenario_boxplot.R` | Summarizes and compares unbiasing output results using boxplots. |
| | `clustering_output.R` | Performs cluster analysis on the model output data. |
| | `distance_scenarios.R` | Calculates the Jaccard distance matrix between scenarios for clustering. |
| | `describe_clusters.R` | Post-processing and statistical description of identified clusters. |
| | `model_verification.R` | Verifies unbiasing performance by comparison against a gridded reference field. Also suitable for general AQ model evaluation. |
| **Utilities** | `test_spatialization.R` | Tests and visualizes the results of the spatialization procedure. |
| | `tiff2netcdf.R` | Converts model output from GeoTIFF to NetCDF format. |
| | `compare-rasters.R` | Compares two raster-based scenarios pixel-by-pixel. |
| | `identify_homogeneous_zones.R` | Identifies homogeneous zones based on spatial patterns of annual mean concentrations. |

## FAIRMODE WG5 Exercise Data Structure

The dataset supports the FAIRMODE WG5 exercise for AQ scenario unbiasing. It contains **annual** data for NO₂, O₃, and PM₂.₅ (humidity-adjusted at 50%) in gridded and point formats.

### Gridded Data (NetCDF)
- **Base Case (Perturbed):** Simulated baseline scenario with bias.
  - `BaseCase_PERT_NO2.nc`: Annual NO₂ (µg/m³).
  - `BaseCase_PERT_O3.nc`: Annual O₃ (ppb).
  - `BaseCase_PERT_PM25_rh50.nc`: Annual PM₂.₅ (µg/m³).
- **Future Scenario (Perturbed):** Projected AQ scenario with bias.
  - `SCEN_PERT_NO2.nc`: Annual NO₂ (µg/m³).
  - `SCEN_PERT_O3.nc`: Annual O₃ (ppb).
  - `SCEN_PERT_PM25_rh50.nc`: Annual PM₂.₅ (µg/m³).

### Point-Based Data (CSV)
- **Reference Observations:** Surface measurements used for bias correction.
  - `yearly_SURF_ppb_O3.csv`: Annual O₃ (ppb).
  - `yearly_SURF_ug_NO2.csv`: Annual NO₂ (µg/m³).
  - `yearly_SURF_ug_PM25_rh50.csv`: Annual PM₂.₅ (µg/m³).

### Formats
- **NetCDF:** Used for spatially continuous gridded data.
- **CSV:** Used for point-based reference observations.

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

## How to Contribute

We welcome contributions to this repository! To ensure a smooth collaboration, please follow these steps:

1. **Fork the Repository** (if you don't have write access).
2. **Clone Your Fork or the Main Repository:**
   ```bash
   git clone https://github.com/jobonaf/scenario-unbias.git
   cd scenario-unbias
   ```
3. **Create a New Branch from `devel`:**
   ```bash
   git checkout -b <new_branch_name> origin/devel
   ```
   *Tip*: Use the naming convention `devel-<feature_name>-YYYYMMDD` for your development branches.
4. **Develop Your Changes** and test them locally.
5. **Commit Your Changes:**
   ```bash
   git add .
   git commit -m "Add feature X or Fix bug Y"
   ```
6. **Push to Your Branch:**
   ```bash
   git push origin devel-<cosa_fa>-YYYYMMDD
   ```
7. **Open a Pull Request** from your branch to `devel` on GitHub.
8. **Code Review and Merge:** The maintainers will review your pull request and merge it if everything is in order.

## Submitting an Issue

If you encounter a bug, have a question, or want to propose a new feature, you can submit an issue to help improve the repository. Here's how:

1. **Navigate to the Issues Tab**  
   Go to the repository on GitHub and click on the **"Issues"** tab.

2. **Click "New Issue"**  
   Click the **"New Issue"** button to start creating your report.

3. **Choose an Issue Type**  
   Depending on the repository's setup, you might see templates for different types of issues (e.g., bug reports, feature requests). Select the most appropriate one.

4. **Provide a Clear Title and Description**  
   - Write a concise title summarizing the issue.
   - In the description, include:
     - A detailed explanation of the problem or suggestion.
     - Steps to reproduce the bug (if applicable).
     - Relevant files, code snippets, or screenshots to help illustrate the issue.

5. **Assign Labels (Optional)**  
   If you have permission, add labels to categorize the issue (e.g., `bug`, `enhancement`, `question`).

6. **Submit the Issue**  
   Once everything is filled out, click **"Submit New Issue"** to post it.

7. **Stay Engaged**  
   Be prepared to answer follow-up questions or clarify details as contributors review your issue.

By submitting a well-documented issue, you'll help the maintainers address problems or implement improvements more effectively!
