import pymc as pm
import numpy as np
import pandas as pd
import arviz as az

az.rcParams["plot.matplotlib.show"] = True  # bokeh plots are automatically shown by default

data = pd.read_csv("/home/ofir/Dropbox/eclipse workspace/lab/Alon/M1_MSE.csv")
# Prepare categorical indices for 'Map'
maps = pd.Categorical(data['Map']).codes
# Prepare stages, treating 'M1_ME' as the baseline
stage_effect = pd.Categorical(data['Stage'], categories=['M1_MSE', 'M2_MSE', 'M3_MSE'], ordered=True).codes - 1

with pm.Model() as model:
    # Intercept (baseline mean for M1_ME)
    intercept = pm.Normal('intercept', mu=0, sigma=10)

    # Effects for M2_ME and M3_ME relative to M1_ME
    delta = pm.Normal('delta', mu=0, sigma=10, shape=2)  # Only two deltas since M1_ME is the reference
    exp_delta = pm.Deterministic("exp_delta", 1-pm.math.exp(delta))

    # Standard deviation for the random effects
    sigma_u = pm.HalfNormal('sigma_u', sigma=10)

    # Random effects for each map
    u = pm.Normal('u', mu=0, sigma=sigma_u, shape=len(np.unique(maps)))

    # Stage-specific variances for the observation model
    sigma_y = pm.HalfNormal('sigma_y', sigma=1, shape=1)

    betas_diff = pm.Deterministic("betas_diff", exp_delta[0]-exp_delta[1])
    # Expected value calculation:
    # If stage_effect is -1 (M1_ME), just use mu_M1_ME. Otherwise, use mu_M1_ME + delta[stage_effect]
    y_log = intercept + pm.math.switch(pm.math.eq(stage_effect, -1), 0, delta[stage_effect]) + u[maps]
    y_hat = pm.Deterministic("p", pm.math.exp(y_log))
    # Likelihood
    y = pm.Gamma('y', mu=y_hat, sigma=sigma_y, observed=data['Value'])  # +1 to adjust index for sigma_y

    # Sampling from the model
    trace = pm.sample(tune=5000, draws=5000, return_inferencedata=True, idata_kwargs={"log_likelihood": True})

# Summarize the posterior
summary = az.summary(trace, var_names=["delta", "intercept", "betas_diff", "exp_delta"], hdi_prob=0.95)
print(summary)
summary.to_csv("Alon_MSE_betas.csv")

# Plot the posterior trace
#az.plot_trace(trace)
trace.to_netcdf("Alon_MSE_betas.nc")

az.plot_trace(trace, var_names=["delta", "intercept", "betas_diff"])

#create posterior predictions for model evaluation
with model:
    pp_samples = pm.sample_posterior_predictive(trace, var_names=["delta", "intercept", "betas_diff", "sigma_y"], extend_inferencedata=True)

summary = az.summary(trace, var_names=["delta", "intercept", "betas_diff", "sigma_y"], hdi_prob=0.95)
print(summary)

#save model traces
#trace.to_netcdf("alon_full.nc")
pp_samples.to_netcdf("Alon_MSE_betas_pp.nc")

print(az.waic(trace))
print(az.loo(trace))
