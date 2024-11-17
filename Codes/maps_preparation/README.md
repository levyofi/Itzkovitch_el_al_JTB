# Maps generation
The code in this folder is used to create the input maps of the physical and machine learning models. The preparation consists of 3 stages:
1) Generating 5 random maps from each flight.
2) Generating height maps.
3) Generating solar, TGI, shade, and skyview maps.
   
## Generating 5 random maps 
The code in `main.py` generates random 1100x1100 pixel maps of RGB, Digital Surface Model (DSM), and ground temperature (IR). 

The code reads the original maps from Pix4DMapper, and saves the output maps in a folder with the name of the flight. In the example of this repository, the maps will be saved at `Example data/Input data/Zeelim_31.05.21_1516`. The code uses the table in `Example data/Input data/flights_info_final.csv` to find the files in our server.

More information and explanation about using the code can be found in `Resampling Usage.md`. 

## Generating height maps



