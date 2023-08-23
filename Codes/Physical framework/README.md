To use the deterministic physical model, two steps are needed.

  The first step is to create the model input data. In order to do so, one needs to convert the raw drone maps into input maps (tif format), as well as creating a meteorological dataframe.
We used Pix4DMapper software as the first image proccessing stage. 
Next, we convert the maps from Pix4DMapper software (DSM, DTM, RGB) into input maps (TGI, height, real solar, skyview and shade) using the attached script create_input_maps.R .  The code need three parameters to run properly: 
  - folder in which images are saved
  - date of flight
  - time of flight
In parallell, the script nameג create_meteorological_df.R XXXX

the second step
-
