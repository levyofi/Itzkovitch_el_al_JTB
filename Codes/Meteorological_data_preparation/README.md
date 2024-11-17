# Meteorological Data Generation

This folder contains the code used to generate a table of meteorological data required for the physical model. The process consists of two stages:

1. **`download_GLDAS_data.R`**  
   This script downloads the relevant GLDAS data for each flight. The downloaded files are saved in the appropriate input folder.

2. **`collect_meteorological_data.R`**  
   This script compiles the meteorological data for each flight. It combines data from two sources:
   - The flight information table, located at `Example data/Input data/flights_info_final.csv`.
   - The GLDAS data downloaded in the previous step.

At the end of this process, a consolidated table with all the required meteorological data is saved. For an example, see the file `Example data/Input data/input_meteorological_data.csv`, which contains the meteorological data for the example flight provided in this repository.
