# The physical mircoclimate model

To run the physical model script, julia environment need to be installed.
1. `physical_model.jl:` this script runs the physical model.
2. To run the model on your own data, one change need to be done - "ARGS" variable.
   
ARGS = [input, output, soil type (default = 0)]. 
ARGS variable need to be set as following:  
- fit soil type: we fitted to 8 (desert). Default is 0.
- input = input folder.
- output = output folder.

The script generates a raster with the ground temperature predictions, and saves it in the output folder.
