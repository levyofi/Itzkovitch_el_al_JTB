import pandas as pd
from patsy import dmatrix
import numpy as np
import pymc as pm
#import pymc_bart as pmb
import pytensor
import pytensor.tensor as at

#import pymc.sampling.jax as pyjax #- need cuda 11.8 and above for that
#from statsmodels.stats.stattools import durbin_watson as dwtest
pytensor.config.exception_verbosity = 'high'
import arviz as az
az.rcParams["plot.matplotlib.show"] = True  # bokeh plots are automatically shown by default


data = pd.read_csv("/home/ofir/Dropbox/eclipse workspace/lab/Alon/regression_table.csv")
data['map'] = pd.Categorical(data['sys_error']).codes

# Generate the design matrix with main effects and all 2-way and 3-way interactions
formula = "TGI + height + real_solar + skyview + shade + shade:(TGI + height + real_solar + skyview)"
X = dmatrix(formula, data, return_type='dataframe')
map = data['map']
real_solar = data['real_solar']
height = data['height']
phyPE = data['phyPE']
mlPE = data['mlPE']

TGI = data['TGI']
# Number of observations
N = X.shape[0]

# Number of predictors in the design matrix (main effects + 2-way and 3-way interactions)
n_predictors = X.shape[1]

# Number of unique groups for the random effect (based on unique values in sys_error)
n_maps = len(np.unique(map))

with pm.Model() as combined_model:
    # Model for physPE
    intercept_physPE = pm.Normal("intercept_physPE", mu=0, sigma=5)
    betas_physPE = pm.Normal("betas_physPE", mu=0, sigma=100, shape=n_predictors)

    # Calculate slopes for shade in physPE
    betas_shade_physPE = pm.Deterministic("betas_shade_physPE", betas_physPE[1:5] + betas_physPE[6:10])

    # Random effects for 'sys_error' in physPE
    sigma_u_physPE = pm.HalfNormal("sigma_u_physPE", sigma=1)
    u_map_physPE = pm.Normal("u_map_physPE", mu=0, sigma=sigma_u_physPE, shape=n_maps)

    sigma_s_physPE = pm.HalfNormal("sigma_s_physPE", sigma=1)
    s_map_physPE = pm.Normal("s_map_physPE", mu=0, sigma=sigma_s_physPE, shape=n_maps)

    # Heteroscedasticity in physPE
    alpha_0_physPE = pm.Normal("alpha_0_physPE", mu=0, sigma=10)
    alpha_1_physPE = pm.Normal("alpha_1_physPE", mu=0, sigma=10)
    alpha_2_physPE = pm.Normal("alpha_2_physPE", mu=0, sigma=10)

    # Linear predictor for physPE
    mu_physPE = pm.math.dot(betas_physPE, X.T) + TGI.values * s_map_physPE[map] + u_map_physPE[map]
    phyPE_pred = pm.Deterministic("phyPE_pred", mu_physPE)

    sigma_physPE = pm.Deterministic("sigma_physPE",
                                    pm.math.exp(alpha_0_physPE + alpha_1_physPE * height + alpha_2_physPE * real_solar))

    # Likelihood for physPE
    phyPE_obs = pm.Normal("phyPE_obs", mu=mu_physPE, sigma=sigma_physPE, observed=phyPE)

    # Model for mlPE
    intercept_mlPE = pm.Normal("intercept_mlPE", mu=0, sigma=5)
    betas_mlPE = pm.Normal("betas_mlPE", mu=0, sigma=100, shape=n_predictors)

    # Calculate slopes for shade in mlPE
    betas_shade_mlPE = pm.Deterministic("betas_shade_mlPE", betas_mlPE[1:5] + betas_mlPE[6:10])

    # Random effects for 'sys_error' in mlPE
    sigma_u_mlPE = pm.HalfNormal("sigma_u_mlPE", sigma=1)
    u_map_mlPE = pm.Normal("u_map_mlPE", mu=0, sigma=sigma_u_mlPE, shape=n_maps)

    sigma_s_mlPE = pm.HalfNormal("sigma_s_mlPE", sigma=1)
    s_map_mlPE = pm.Normal("s_map_mlPE", mu=0, sigma=sigma_s_mlPE, shape=n_maps)

    # Heteroscedasticity in mlPE
    alpha_0_mlPE = pm.Normal("alpha_0_mlPE", mu=0, sigma=10)
    alpha_1_mlPE = pm.Normal("alpha_1_mlPE", mu=0, sigma=10)
    alpha_2_mlPE = pm.Normal("alpha_2_mlPE", mu=0, sigma=10)

    # Linear predictor for mlPE
    mu_mlPE = pm.math.dot(betas_mlPE, X.T) + real_solar.values * s_map_mlPE[map] + u_map_mlPE[map]
    mlPE_pred = pm.Deterministic("mlPE_pred", mu_mlPE)

    sigma_mlPE = pm.Deterministic("sigma_mlPE", pm.math.exp(alpha_0_mlPE + alpha_1_mlPE * height + alpha_2_mlPE * TGI))

    # Likelihood for mlPE
    mlPE_obs = pm.Normal("mlPE_obs", mu=mu_mlPE, sigma=sigma_mlPE, observed=mlPE)

    # Calculate the differences between physPE and mlPE betas
    betas_diff = pm.Deterministic("betas_diff", betas_physPE - betas_mlPE)
    betas_shade_diff = pm.Deterministic("betas_shade_diff", betas_shade_physPE - betas_shade_mlPE)

    # Calculate the differences between physPE and mlPE betas
    abs_betas_diff = pm.Deterministic("abs_betas_diff", pm.math.abs(betas_physPE) - pm.math.abs(betas_mlPE))
    abs_betas_shade_diff = pm.Deterministic("abs_betas_shade_diff", pm.math.abs(betas_shade_physPE) - pm.math.abs(betas_shade_mlPE))

    # Sample from the combined model
    trace = pm.sample(tune=2000, draws=2000, return_inferencedata=True, idata_kwargs={"log_likelihood": True})

# List of variable names for the summary, including the differences
var_names = [
    "betas_physPE", "betas_shade_physPE", "alpha_0_physPE", "alpha_1_physPE", "alpha_2_physPE",
    "betas_mlPE", "betas_shade_mlPE", "alpha_0_mlPE", "alpha_1_mlPE", "alpha_2_mlPE",
    "betas_diff", "betas_shade_diff", "abs_betas_diff", "abs_betas_shade_diff"
]

# Generate the summary with the new variable names
summary = az.summary(trace, var_names=var_names, hdi_prob=0.95)
print(summary)
summary.to_csv("Alon_trace_both models.csv")

# Plot the posterior trace
#az.plot_trace(trace)

trace.to_netcdf("alon_trace_both_models.nc")

print(az.waic(trace))
print(az.loo(trace))

az.plot_ppc(trace)
