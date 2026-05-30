# inconel-processing-maps-dmm
Bivariate spline modeling and execution of the Dynamic Materials Model (DMM) to compute processing maps and Ziegler instability criteria for nickel-based superalloys.
# Inconel Superalloys Processing Maps via Dynamic Materials Model (DMM)

## Project Overview
This repository contains a Materials Informatics workflow developed to model and evaluate the hot workability of nickel-based superalloys (**IN718, IN625, IN100, and IN600**). Utilizing raw flow stress data across varying strain rates and temperatures, the project implements the **Dynamic Materials Model (DMM)** in Python to synthesize continuous power dissipation efficiency contours and predict metallurgical flow instabilities. 

The generated processing windows are cross-validated against standard aerospace and handbook reference data to optimize high-temperature industrial forging and manufacturing parameters.

## Core Metallurgy & DMM Theory
The Dynamic Materials Model treats a deforming metal as a power-dissipating system where total power ($P$) is divided into temperature/microstructural evolution paths:
* **Strain Rate Sensitivity ($m$):** Evaluated locally as:
  $$m = \frac{\partial(\ln \sigma)}{\partial(\ln \dot{\epsilon})}$$
* **Power Dissipation Efficiency ($\eta$):** Represents the percentage of power consumed by microstructural changes (e.g., Dynamic Recrystallization vs. simple dynamic recovery):
  $$\eta = \frac{2m}{m + 1} \times 100\%$$
* **Ziegler Instability Criterion ($\xi$):** Defines a negative physical regime ($\xi < 0$) where plastic flow becomes unstable, risking defects like adiabatic shear bands or localized cracking:
  $$\xi = \frac{\partial \left(\ln \left(\frac{m}{m+1}\right)\right)}{\partial(\ln \dot{\epsilon})} + m < 0$$

## Python Architecture & Data Analysis
To transition from coarse experimental discrete data points into continuous visual maps, the computational pipeline executes the following steps:
1. **Data Transformation:** Logs raw stress and strain rate arrays to achieve log-linear normalization.
2. **High-Resolution Mesh Construction:** Generates an intensive $300 \times 300$ resolution high-density grid domain.
3. **Bivariate Stress Surface Modeling:** Fits an advanced `RectBivariateSpline` (Cubic Spline interpolation) to calculate exact local derivatives smoothly.
4. **Criterion Mapping:** Computes matrices for $\eta$ and $\xi$, plotting efficiency curves atop shaded flow instability regions.

## Consolidated Alloy Findings

### 1. IN718 Superalloy (Starting Grain Size: 25 µm, Target Strain: 0.4)
* **Domain 1 (910–1010 °C):** Controlled strictly by the $\delta$-phase ($Ni_3Nb$). Precipitate pinning limits extensive grain growth during hot work. However, the model correctly flags a structural instability risk below a critical strain rate threshold ($\sim 10^{-1.7} \text{ s}^{-1}$), capturing cracking windows.
* **Domain 2 (1040–1100 °C):** Above the $\delta$-solvus temperature ($\sim 1025\text{ °C}$), total matrix dissolution occurs, enabling effortless dynamic recrystallization (DRX) and maximum work efficiencies exceeding 32%.

### 2. IN625 Superalloy (Target Strain: 0.5)
* **Optimal Window discovered:** $1150\text{ °C}$ at a strain rate of $1.0\text{ s}^{-1}$ yielding reliable Dynamic Recrystallization.
* **Danger Zones:** High strain rates at lower temperatures trigger intense flow localization (shear banding). Conversely, high-temperature, low-rate zones (Domain 3) risk severe wedge cracking despite registering false high mathematical efficiency values.

### 3. IN600 & IN100 Superalloys
* Continuous processing boundaries map microstructural shifts, isolating the precise effects of chemical additions (like Iron improving global diffusion properties) on lowering deformation activation energies.
* Documented discrepancies between the code output and literature benchmarks are systematically traced to interpolation edge constraints, local slope calculations, and starting experimental point limitations.

## Repository Contents
* INCONEL_ALLOYS_PROCESSING_MAPS.pdf: Detailed project presentation summarizing microstructural phases, DMM workflows, and comparative map sets.
* MM208-code.py: Complete Python script housing raw arrays, scipy.interpolate mechanics, calculations, and matplotlib contour visualization properties.
