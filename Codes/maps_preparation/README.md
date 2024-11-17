## Generating random square subimages from an image in `random_subimages.py`
This code is used to generate random square maps from a larger map. These subimages are created with minimal congruence to ensure variety. In this study, we used this code to generate 1100x1100 pixel maps of RGB, Digital Surface Model (DSM), and ground temperature (IR). 

The maps is this study can be generated using the `main.py` file. The code reads the original maps from Pix4DMapper, and saves the output maps in a folder with the name of the flight. In the example of this repository, the maps will be saved at `Example data/Input data/Zeelim_31.05.21_1516`

More information and explanation about using the code can be found in `Usage.md`.
