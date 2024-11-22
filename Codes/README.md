# Code for Microclimate Modeling and Bias-Correction

This directory contains the code associated with the manuscript "From Big Data to Small Scales: Using Machine Learning to Correct Biases in Microclimate Predictions" by Itzkovitch et al. The repository is organized into subdirectories based on the different stages of the study, from data preparation to analysis and visualization.

## Directory Structure

### 1. `Figures/`
Scripts and resources for creating the figures in the manuscript and supplementary materials.

- Contains scripts in both Python and R for generating all figures, including main text and supplementary figures.
- Example: Scripts for generating the results visualizations in Figure S3.

### 2. `Machine_Learning_model/`
Code for training and testing the machine learning bias-correction model.

- Includes the Random Forest implementation (`RF_model.py`) for predicting and correcting microclimate model biases.

### 3. `Meteorological_data_preparation/`
Scripts for preparing meteorological input data required by the physical model.

- Example: `download_GLDAS_data.R` for downloading and organizing GLDAS data.

### 4. `Physical_framework/`
Scripts related to the physical microclimate modeling framework.

- Includes the main Julia script for running the heat balance model with additional optimizations and a loop for automation.

### 5. `Statistical_Analysis/`
Code for statistical analyses presented in the manuscript.

- Example: Linear regression models and Bayesian analysis tools to assess the machine learning model's impact on prediction accuracy.
- Includes both R and Python scripts for the statistical models.

### 6. `maps_preparation/`
Scripts for preparing input maps used in the microclimate and machine learning models.

- Includes processes for creating vegetation, height, solar radiation, shade, and skyview maps.
- Example: `README.md` in this directory provides detailed instructions on generating and cropping the maps.

## Usage Instructions

1. **Prepare Input Maps**: Start with the `maps_preparation/` scripts to generate input maps from raw drone imagery.
2. **Meteorological Data**: Use the scripts in `Meteorological_data_preparation/` to download and format meteorological data.
3. **Run Physical Model**: Navigate to `Physical_framework/` and run the Julia scripts to execute the physical microclimate model.
4. **Train Machine Learning Model**: Use the Python code in `Machine_Learning_model/` to train and test the bias-correction model.
5. **Statistical Analysis**: Analyze the results and evaluate the models using the scripts in `Statistical_Analysis/`.
6. **Create Figures**: Use the scripts in `Figures/` to reproduce the visualizations in the manuscript.

## Requirements

- Python 3.x
- Julia 1.x
- R version 4.x or higher
- Required Python, Julia, and R packages are listed in their respective scripts. Install them as needed.
