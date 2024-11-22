# Statistical Analysis

This folder contains the scripts used for conducting Bayesian statistical analyses. The analyses focus on comparing model errors (ME, MAE, and MSE) before and after machine learning correction, and on understanding the relationships between environmental parameters and microclimate model biases using PyMC for Bayesian modeling.

## Folder Contents

### 1. `Stat Model 1 - ME_differences_with_variance_structure.py`, `Stat Model 1 - MAE_differences.py`, and `Stat Model 1 - MSE_differences.py`
- **Purpose**: These scripts analyze the differences in prediction errors before and after machine learning correction, focusing on Mean Error (ME), Mean Absolute Error (MAE), and Mean Squared Error (MSE).
- **Output**:
  - Summarized posterior distributions for error differences.
  - Credible intervals and statistical summaries showing improvements across all metrics.

### 2. `Stat Model 2 - bias vs features model.py`
- **Purpose**: Investigates the relationships between microclimate model bias and environmental features (e.g., solar radiation, TGI, height, and skyview) before and after machine learning correction.
- **Output**:
  - Posterior distributions and credible intervals for the relationships between features and bias.
  - Comparisons of slopes before and after machine learning correction.

