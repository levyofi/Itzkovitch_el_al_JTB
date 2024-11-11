### Input Data for Physical Model, Bias Calculation, and Machine Learning Bias-Correction

This folder contains the data used to run the physical model, calculate bias, and train the machine learning bias-correction model. The maps below serve as essential inputs for generating physical model predictions and are also used as features in the machine learning model to correct biases in those predictions.

**Files for Physical Model Predictions:**
1. `TGI.tif`: Vegetation index map based on the Triangular Greenness Index (TGI), indicating the presence and density of vegetation.
2. `height.tif`: Map of the height above ground for each pixel, influencing exposure to environmental factors.
3. `shade.tif`: Binary map where each pixel is labeled as either shaded or exposed to sunlight, used to model solar radiation effects.
4. `skyview.tif`: Skyview factor map, representing the fraction of visible sky from each pixel, affecting longwave radiation exposure.
5. `real_solar.tif`: Map of real solar radiation values, adjusted for local conditions, helping to model direct and diffuse solar radiation effects.

These maps serve as input for the physical model in `Codes/Physical framework/physical_model.jl`, guiding predictions of ground temperature. They are also used as features in the machine learning (ML) bias-correction model in `Codes/Statistical framework/RF_model.py` to help identify and reduce biases in the physical model's predictions.

**Observed Temperature Map for Validation:**

- `IR.tif`: Temperature map recorded with an infrared camera, providing real ground temperature values for model validation. This map is used in `Codes/Statistical framework/RF_model.py` to calculate the bias before and after the ML correction. 

**RGB Orthophoto Map:**

- `RGB.tif`: RGB image used to calculate `TGI.tif`. It also serves for visualization, aiding in the interpretation of map data and comparison with real-world conditions.
