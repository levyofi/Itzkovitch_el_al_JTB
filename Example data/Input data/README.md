### Input Data for Physical Model, Bias Calculation, and Machine Learning Bias-Correction

This folder contains the data used to run the physical model, calculate bias, and train the machine learning bias-correction model. The maps below provide essential input for generating physical model predictions and are also used as features for the machine learning model to correct biases in those predictions.

**Files for Physical Model Predictions:**
1. `TGI.tif`: A vegetation index map based on the Triangular Greenness Index (TGI), indicating the presence and density of vegetation.
2. `height.tif`: A map of the height above ground for each pixel, which influences exposure to environmental factors.
3. `shade.tif`: A binary map where each pixel is labeled as either shaded or exposed to sunlight, used to model solar radiation effects.
4. `skyview.tif`: A skyview factor map, representing the fraction of visible sky from each pixel, which affects longwave radiation exposure.
5. `real_solar.tif`: A map of real solar radiation values, adjusted for local conditions, that helps model direct and diffuse solar radiation effects.

The above maps serve as input for the physical model in `Codes/Physical framework/physical_model.jl`, guiding predictions of ground temperature. They are also used as features in the machine learning (ML) bias-correction model in `Codes/Statistical framework/RF_model.py` to help identify and reduce biases in the physical model's predictions.

**Observed Temperature Map for Validation:**
`IR.tif`: A temperature map recorded with an infrared camera, providing real ground temperature values for model validation. The map is used to calculate the bias before and after the ML correction. 

**RGB Orthophoto map:**
3. `RGB.tif`: An RGB image for visualization, aiding in the interpretation of map data and comparison with real-world conditions.
