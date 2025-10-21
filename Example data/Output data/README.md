### Model Output Data

* **`Zeelim_31.05.21_1516_1_Physical_prediction.tif`:**
  Contains the temperature predictions from the physical model (before any machine learning correction).
* **`Zeelim_31.05.21_1516_1_Correction_map.npy`:**
  A NumPy file containing the temperature predictions **after correction** using the random forest model.

**Note:**
These files represent the temperature predictions before and after ML correction. `Physical_prediction.tif` provides the initial predictions from the microclimate model, while `correction_map.npy` contains the final, corrected temperatures after applying the random forest model.
