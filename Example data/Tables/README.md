In this folder all relevant tables and dataframes that are used in the codes attached.
- `met_data.csv:` meteorological data file.
- `RF_CorrectionModel_summary.csv:` dataframe that summarize the results of all models (pysical and machine learning) for every cropped map.
- `SampledPixels.csv:` dataframe contains 1,000 sampled pixels for every test map, in order to visualize the results.
   - `SampledPixels_[1-2][2-3].csv` are once again tables that sample test map, to visualize differences between models.
- `online_met_dbs_comparison.csv:` a table that summarize the results of before_Ml and after_ML models, while the meteorological data derived from online databases (GLDAS and ERA5).
- `ML_NN_comparison.xlsx:` a table that compare the results of random forest and DNN corrections.
- in Figure S4, two more files are used. The two files are larger than the maximum volume allowed to upload in Github.
   - path_before_ml <- `'/home/ofir/Dropbox/pycharm_projects/pymc_models/Alon/alon_trace_phys.nc'`
   - path_after_ml  <- `'/home/ofir/Dropbox/pycharm_projects/pymc_models/Alon/alon_trace_ml.nc'`
