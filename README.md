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

---

## How to Contribute

We welcome contributions to this repository! To ensure a smooth collaboration, please follow these steps based on your needs:

### For Personal Use
If you want to work independently on a stable copy of the repository:

1. Clone the repository and work directly on the `main` branch.
2. Pull updates periodically to stay aligned with the latest stable version:
   ```bash
   git pull origin main
   ```
3. Avoid making pull requests.

### For Contributors
If you want to contribute code to the repository:

1. **Clone the Repository**
   Clone the repository to your local machine:
   ```bash
   git clone https://github.com/jobonaf/scenario-unbias.git
   cd scenario-unbias
   ```

2. **Switch to the Development Branch**
   Always work on the `devel` branch to keep the `main` branch stable:
   ```bash
   git checkout -b devel origin/devel
   ```

3. **Make Your Changes**
   Develop your features or fix bugs locally. Make sure your changes are properly tested.

4. **Commit and Push**
   Commit your changes with a descriptive message:
   ```bash
   git add .
   git commit -m "Add feature X or Fix bug Y"
   git push origin devel
   ```

5. **Open a Pull Request**
   - Navigate to the repository on GitHub.
   - Open a pull request to align `devel` with `main`.
   - Add a clear description of your changes and why they are needed.

6. **Code Review and Merge**
   The codeowner will review your pull request and merge it into the `main` branch if everything is in order.

### Staying Updated
To keep your local copy of `devel` in sync with the latest changes, periodically run:
```bash
git fetch origin
git merge origin/devel
```

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

### How to Contribute

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

### Submitting an Issue

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
