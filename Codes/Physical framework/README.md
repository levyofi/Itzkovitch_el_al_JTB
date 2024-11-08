### Here we can found scripts to run physical deterministic model (before ML)
To use the deterministic physical model, two steps are needed. 
The first is to run R scripts to create input data (with or without meteorological data, depend on the user), and the second is to run julia script to get physical predictions.

## Step 1 - convert drone images to maps and create input data (python, R and Pix4DMapper)
1. The first step is to create the model input data. In order to do so, one needs to convert the raw drone maps into input maps (tif format), as well as creating a meteorological dataframe.

We used Pix4DMapper software as the first image proccessing stage.  
2. `create_input_maps.R:` Next, we convert the maps from Pix4DMapper software (DSM, DTM, RGB) into input maps (TGI, height, real solar, skyview and shade) using the attached script.  

The code need three parameters to run properly: 
  - folder in which images are saved
  - date of flight
  - time of flight

3. `resample_images.py:` cropped the map into 5 sub-maps with  1,000*1,000 shape.

4. `create_meteorological_df_WithStation.R:` the script stands for taking meteorological data recorded in a mobile meteorological station (saved in dataframe, each row related to one flight), merge it with an online meteorological data (taken from GLDAS dataset), and save it in csv file. The file contains one row, and each column related to different meteorological parameter.
    
The code need three parameters to run properly: 
  - folder in which images are saved
  - time of flight (format DD.MM.YY_HHMM)
  - timezone

4.1 `create_meteorological_df_NoStation.R:` this script does the same as the above script ("withStation"), and used for cases where online data are accessed instead of mobile station.

##  Step 2 - Physical model (in julia)
To run the physical model script, julia environment need to be installed.
1. `physical_model.jl:` this script runs the physical model.
2. To run the model on your own data, one change need to be done - "ARGS" variable.
   
ARGS = [input, output, soil type (default = 0)]. 
ARGS variable need to be set as following:  
- fit soil type: we fitted to 8 (desert). Default is 0.
- input = input folder.
- output = output folder.

The script generates the ground temperature predictions, and saves it in the output folder.
