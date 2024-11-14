# Code for Microclimate Modeling and Bias-Correction

This folder contains the code used in our study on microclimate modeling and machine learning-based bias correction. The code is organized into three main stages, each with its own set of scripts and instructions.

## Overview

The code integrates physical modeling and machine learning to improve the accuracy of ground temperature predictions. This folder includes:
1. Scripts for setting up and running the physical microclimate model.
2. A machine learning model to correct biases in the physical model's predictions.
3. Scripts for data analysis and generating publication-quality figures.

## Repository Structure

- **`Physical framework`**: Code for setting up and running the physical model. 

- **`Statistical framework`**: Code for the machine learning bias-correction model. 

- **`Figures`**: Scripts for generating analysis and figures.

## Getting Started

### Prerequisites
To run the code, you need the following software environments:

- **R**: Version 4.2.x or lower to ensure compatibility with the `rgdal` and raster packages. 
- **Julia**: Required packages are automatically downloaded at the start of the script.
- **Python**: Ensure that all necessary packages (listed in `requirements.txt` if available) are installed in your Python environment.

### Running the Code

1. **Stage 1: Physical Framework**
   - Navigate to the `Physical framework` folder.
   - Run the R script for data preparation, then execute the Julia script to run the physical model.

2. **Stage 2: Statistical Framework**
   - Navigate to the `Statistical framework` folder.
   - Run the Python script to train the random forest model and correct the physical model’s biases.

3. **Stage 3: Analysis and Figures**
   - Navigate to the `Figures` folder.
   - Run the relevant scripts in R and Python to reproduce all figures for the paper.
