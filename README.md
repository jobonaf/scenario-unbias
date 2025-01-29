# scenario-unbias <img src='inst/unbias2.png' align="right" height="138" />

The repository contains code designed to test different approaches for reducing bias in air quality (AQ) scenarios, specifically within the context of the FAIRMODE WG5 exercise. While the primary focus is on this particular exercise, the methods and tools provided can also be useful for unbiasing AQ scenarios more generally, even outside of the FAIRMODE framework.

For an overview of the methodology, refer to the [Unbiasing Air Quality Scenarios Documentation](https://github.com/jobonaf/scenario-unbias/blob/main/docs/unbiasing-aq-scenarios.md).

## FAIRMODE WG5 Exercise

The FAIRMODE WG5 bias projection exercise aims to benchmark methodologies for removing bias from air quality model simulations, particularly in future policy scenarios. It focuses on deriving bias correction fields and projecting biases into future scenarios using synthetic datasets. Participants apply their preferred methodologies to post-process provided data and benchmark their results against known synthetic truths. The exercise includes annual data for PM2.5, NO2, and O3 in gridded (NetCDF) and point (CSV) formats for both reference and future projections.

## Code Overview

| Script | Description |
|--------|------------|
| `prepare-aq-obs-data.R` | Prepares air quality observation data for testing. |
| `check-aq-raster.R` | Checks the prepared AQ data and ensures consistency. |
| `unbias-aq-scenario.R` | Core function for unbiasing AQ scenarios. |
| `spatialize-points-to-grid.R` | Converts sparse point data (or coefficients) into a gridded format. |
| `test_spatialization.R` | Tests and visualizes the spatialization results. |
| `calibrate-unbias-coefficient.R` | Calibrates correction coefficients by comparing base case vs observed data. |
| `apply-unbiasing.R` | Applies the correction coefficients to an AQ scenario. |
| `read-fairmode-data.R` | Reads the dataset for the FAIRMODE WG5 exercise. |
| `map-fairmode-data.R` | Visualizes FAIRMODE data on maps. |
| `dashboard-fairmode-data.Rmd` | Displays interactive maps of FAIRMODE data in a dashboard format. |
| `process-fairmode-data.R` | Reads, processes, and applies unbiasing methods to FAIRMODE data. |

## Data Structure

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
   git checkout -b my-feature-branch origin/devel
   ```
4. **Develop Your Changes** and test them locally.
5. **Commit Your Changes:**
   ```bash
   git add .
   git commit -m "Add feature X or Fix bug Y"
   ```
6. **Push to Your Branch:**
   ```bash
   git push origin my-feature-branch
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
