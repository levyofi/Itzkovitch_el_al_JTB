### Here we can find
Two scripts to run the following models:
1. `RF_model.py:` The random forest model to predict prediction error of physical model.
2. `DNN_Correction_Model_Class.py:` A class that can replace random-forest modeling class from the `RF_model.py` script, in order to examine the performances of simple neural-network. Further elaborate in Supplementary.
## Pipeline (How to run the model)
The following steps are required and processed:
1. the user should indicate in which folder the data is saved.
2. the model trained on 80% of the data. In our case, there is no leakage of information within days, meaning all flights from  the same day are in one set (train or test).
3. the model predict results on the other 20%.
4. data is saved, meaning one prediction map for each file.

The base model (class) is random forest, which can be switched to DNN if wanted.
