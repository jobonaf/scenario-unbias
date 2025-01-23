# Unbiasing Methods for Air Quality Scenarios
## Proposal for a Classification Scheme  
**Author:** G. Bonafè  
**Date:** 2025-01-15  

---

## Introduction  
An unbiasing method can be described by defining:  
1. the **unbiasing sequence**,  
2. the **correction algorithm**,  
3. its **calibration method**, and  
4. the **spatialization algorithm**.  

These elements provide a systematic approach to deriving and applying corrections to air quality scenarios by comparing a simulated base case with observed data, identifying biases, and using this information to remove them in the scenario.

---

## Unbiasing Sequences  
Each sequence below defines a unique order for calibration, correction, and spatialization processes, determining the workflow.  

- **SCA (Spatialize - Calibrate - Apply):**  
  Spatializes the observed data to the entire grid, then calibrates and applies the correction algorithm.  

- **CSA (Calibrate - Spatialize - Apply):**  
  First calibrates the correction coefficients, then spatializes them across the grid, and finally applies corrections.  

- **CAS (Calibrate - Apply - Spatialize):**  
  The correction algorithm is calibrated and applied to monitoring sites, and then sparse data are spatialized across the grid.  

- **CA (Calibrate - Apply):**  
  A correction is applied after calibration, either globally over the entire grid or locally at monitoring sites.

---

## Correction Algorithms  
- **Additive (Add):** Adds a constant value to correct biases.  
- **Multiplicative (Mult):** Multiplies by a factor to adjust values.  
- **Rescaled Additive (Resc):** Adds corrections scaled by the ratio between base case and scenario.  
- **Quantile-based (Quant):** Corrects values based on quantile-specific factors.  
- **Linear (Lin):** Applies linear corrections across the range of values.  

These algorithms define how the adjustments are computed by analyzing biases between the simulated base case and observed data, ensuring that modeled scenarios more accurately reflect real-world conditions.

---

## Calibration Methods  
- **Point-based (Point):** Calibration of the coefficients of functions to remove bias, performed at monitoring sites.  
- **Grid-based (Grid):** Calibration of the coefficients applied over the entire grid to correct for bias.  
- **Cell vs Cell (Cell):** Calibration considers each cell individually, adjusting coefficients to remove bias at the cell level.  
- **Cell Neighborhood (Neigh):** Calibration uses information from surrounding cells to adjust coefficients for bias removal.

---

## Spatialization Algorithms  
- **Thin Plate Spline (TPS):** Produces smooth surfaces by minimizing bending energy.  
- **Inverse Distance Weighted (IDW):** Estimates values using a weighted average of nearby points.  
- **Ordinary Kriging (OK):** Geostatistical spatialization based on a variogram model.  
- **Kriging with External Drift (KED):** Extends Kriging by incorporating external variables.  
- **Successive Correction Method (SCM):** An iterative method following Bratseth's approach.  

Spatialization algorithms determine how information is propagated from points to the grid. Depending on the sequence of unbiasing, the propagated information can either be air quality indicators (for sequences SCA, CAS) or the coefficients of correction algorithms (for sequence CSA).

---

## Overview of Unbiasing Methods  
This table summarizes the combinations of sequences, correction algorithms, calibration methods, and spatialization techniques that are consistent and appropriate.  

| **Sequence** | **Correction Algorithm** | **Calibration** | **Spatialization Algorithm** |
|--------------|---------------------------|-----------------|------------------------------|
| **SCA**      | Lin, Quant               | Grid            | KED, SCM                     |
|              | Add, Mult, Resc          | Cell            | KED, SCM                     |
|              | Add, Mult, Lin           | Neigh           | KED, SCM                     |
| **CSA**      | Add, Mult, Resc          | Point           | TPS, IDW, OK, KED            |
| **CAS**      | Add, Mult, Lin, Resc     | Point           | KED, SCM                     |
| **CA**       | Add, Mult, Lin, Quant    | Point           | -                            |
