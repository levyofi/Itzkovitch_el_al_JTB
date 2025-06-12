## Contents

This folder contains two main scripts for running and comparing predictive models:

1. **`RF_model.py`**
   Implements a random forest model to predict the errors of a physical model.
2. **`DNN_Correction_Model_Class.py`**
   Provides a neural network–based correction model that can be used in place of the random forest in `RF_model.py` to compare the performance of a simple DNN. More details are available in the Supplementary Information.

---

## Pipeline: How to Run the Models

To use these models, follow these steps:

1. **Specify the Data Folder:**
   Indicate the directory where your data is stored.
2. **Training and Testing Split:**
   The model is trained on 80% of the data. To prevent information leakage, all samples from the same day (i.e., all flights from a single day) are assigned to either the training or testing set, never both.
3. **Prediction:**
   The trained model predicts results on the remaining 20% of the data.
4. **Saving Results:**
   The output predictions are saved, generating a prediction map for each input file.

By default, the pipeline uses a random forest model, but you can easily switch to the DNN-based model if you wish to compare performance.
