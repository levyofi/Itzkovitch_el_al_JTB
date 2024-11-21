import pymc as pm
import numpy as np
import pandas as pd
import arviz as az
from numpy.core.defchararray import lower

import pymc.sampling.jax as pyjax #- need cuda 11.8 and above for that

az.rcParams["plot.matplotlib.show"] = True  # bokeh plots are automatically shown by default

data = pd.read_csv("/home/ofir/Dropbox/eclipse workspace/lab/Alon/M1_ME.csv")
# Prepare categorical indices for 'Map'
maps = pd.Categorical(data['Map']).codes
# Prepare stages, treating 'M1_ME' as the baseline
stage_effect = pd.Categorical(data['Stage'], categories=['M1_ME', 'M2_ME', 'M3_ME'], ordered=True).codes - 1
data["stage_effect"] = stage_effect

with pm.Model() as model:
    # Intercept (baseline mean for M1_ME)
    intercept = pm.Normal('intercept', mu=0, sigma=10)

    # Effects for M2_ME and M3_ME relative to M1_ME
    delta = pm.Normal('delta', mu=0, sigma=10, shape=2)  # Only two deltas since M1_ME is the reference

    # Standard deviation for the random effects
    sigma_u = pm.HalfNormal('sigma_u', sigma=1)

    # Random effects for each map
    u = pm.Normal('u', mu=0, sigma=sigma_u, shape=len(np.unique(maps)))

    # Stage-specific variances for the observation model
    #sigma_y = pm.Uniform('sigma_y', lower=0.05, upper=5, shape=3)
    sigma_y = pm.HalfNormal('sigma_y', sigma=5, shape=3)
    #tau = pm.Gamma('tau', alpha = 0.01, beta = 0.01, shape = 3)
    betas_diff = pm.Deterministic("betas_diff", delta[0]-delta[1])
    # Expected value calculation:
    # If stage_effect is -1 (M1_ME), just use mu_M1_ME. Otherwise, use mu_M1_ME + delta[stage_effect]
    y_hat = intercept + pm.math.switch(pm.math.eq(stage_effect, -1), 0, delta[stage_effect]) + u[maps]

    # Likelihood
    y = pm.Normal('y', mu=y_hat, sigma=sigma_y[stage_effect+1], observed=data['Value'])  # +1 to adjust index for sigma_y

    # Sampling from the model
    #trace = pm.sample(tune=10000, draws=10000, return_inferencedata=True, idata_kwargs={"log_likelihood": True})
    trace = pyjax.sample_numpyro_nuts(5000, tune=2000, chains=3,
                                      idata_kwargs={"log_likelihood": True})  # , chain_method="sequential")

# Summarize the posterior
summary = az.summary(trace, var_names=["delta", "intercept", "betas_diff"], hdi_prob=0.95)
print(summary)
summary.to_csv("Alon_ME_betas.csv")

# Plot the posterior trace
#az.plot_trace(trace)
trace.to_netcdf("Alon_ME_betas.nc")

az.plot_trace(trace, var_names=["delta", "intercept", "betas_diff"])

#create posterior predictions for model evaluation
with model:
    pp_samples = pm.sample_posterior_predictive(trace, var_names=["delta", "intercept", "betas_diff", "sigma_y"], extend_inferencedata=True)

summary = az.summary(trace, var_names=["delta", "intercept", "betas_diff", "sigma_y"], hdi_prob=0.95)
print(summary)

#save model traces
#trace.to_netcdf("alon_full.nc")
pp_samples.to_netcdf("Alon_ME_betas_pp.nc")

print(az.waic(trace))
print(az.loo(trace))
