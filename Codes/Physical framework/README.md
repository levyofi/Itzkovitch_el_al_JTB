To use the deterministic physical model, two steps are needed. 
The first is to run R scripts to create input data (with or without meteorological data, depend on the user), and the second is to run julia script to get physical predictions.

input data
-
The first step is to create the model input data. In order to do so, one needs to convert the raw drone maps into input maps (tif format), as well as creating a meteorological dataframe.
We used Pix4DMapper software as the first image proccessing stage.  
Next, we convert the maps from Pix4DMapper software (DSM, DTM, RGB) into input maps (TGI, height, real solar, skyview and shade) using the attached script create_input_maps.R .  
The code need three parameters to run properly: 
  - folder in which images are saved
  - date of flight
  - time of flight
    
In parallell, the script named create_meteorological_df.R stands for taking meteorological data taken from a mobile meteorological station (saved in dataframe, each row related to one flight), merge it with an online meteorological data (taken from GLDAS dataset), and save it in csv file. The file contains one row, and each column related to different meteorological parameter.   
The code need three parameters to run properly: 
  - folder in which images are saved
  - time of flight (format DD.MM.YY_HHMM)
  - timezone

Physical model
-
To run the physical model script, julia environment need to be installed. 
To run the model on your own data, one change need to be done - on "ARGS" variable.   
ARGS variable need to be set as following:  
ARGS = [input, output, soil type (default = 0)].   
The script generates the ground temperature predictions, and saves it in the output folder.
