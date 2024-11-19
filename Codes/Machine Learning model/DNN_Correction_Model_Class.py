import glob
import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import r2_score

class DNN_Correction_Model():

    def __init__(self, X_path, y_path, name, map_size = 1_024):
        '''
        Create model class.
        At initiation - create dataframe that contain all paths to input maps
        '''
        self.name = name
        self.N = map_size
        self.features = ['TGI', 'height', 'shade', 'real_solar', 'skyview']

        files_dict = {'IR':[], 'TGI':[], 'height':[], 'shade':[], 'real_solar':[], 'skyview':[]}
        temp_list = glob.glob(f'{y_path}/physical_model/*.tif')
        temp_list.sort()
        files_dict['M1'] = temp_list
        flights = [f.split('/')[-1][:-6] for f in temp_list]
        files_dict['Flight'] = flights
        for name in np.unique(flights):
            for feature in self.features:
              temp_list = glob.glob(f'{X_path}/{name}/{feature}_*.tif')
              temp_list.sort()
              files_dict[feature] += temp_list
            temp_list = glob.glob(f'{y_path}/IR_fixed/{name}*.npy')
            temp_list.sort()
            files_dict['IR'] += temp_list

        self.files_df = pd.DataFrame(files_dict)

    def split_data(self, pixels, train_set_maps):
        '''
        Create dataframe for random forest model.
        For each map within the train set maps, and for each crop at those maps,
        the script takes [pixels] random pixels.
        '''
        self.train_flights = pd.Series(train_set_maps)
        self.len_of_random_pix = pixels

        masks = {}
        dctNN = {'PredM1': [], 'PredErrorM1':[], 'IR': [],
               'TGI': [], 'Height': [], 'Shade': [], 'RealSolar': [], 'Skyview': [],
               }

        for flight in self.train_flights:
            subset = self.files_df[self.files_df['Flight'] == flight].reset_index()
            for crop in range(len(subset)):
                rand = np.random.randint(0,self.N**2, pixels) # list of random indexes from the map
                m = np.zeros((self.N,self.N)).flatten()
                m[rand] = 1
                masks[f'{flight}_{crop}'] = m.reshape((self.N, self.N))

                dctNN['TGI'] += list(gdal.Open(subset['TGI'][crop]).ReadAsArray().flatten()[rand])
                dctNN['Height'] += list(gdal.Open(subset['height'][crop]).ReadAsArray().flatten()[rand])
                dctNN['Shade'] += list(gdal.Open(subset['shade'][crop]).ReadAsArray().flatten()[rand])
                dctNN['RealSolar'] += list(gdal.Open(subset['real_solar'][crop]).ReadAsArray().flatten()[rand])
                dctNN['Skyview'] += list(gdal.Open(subset['skyview'][crop]).ReadAsArray().flatten()[rand])

                pred_m1 = gdal.Open(subset['M1'][crop]).ReadAsArray().flatten()[rand]
                dctNN['PredM1'] += list(pred_m1)

                ir = np.load(subset['IR'][crop]).flatten()[rand]
                dctNN['IR'] += list(ir)

                pred_error_m1 = gdal.Open(subset['M1'][crop]).ReadAsArray().flatten()[rand] - \
                    np.load(subset['IR'][crop]).flatten()[rand] - 273.16
                dctNN['PredErrorM1'] += list(pred_error_m1 - np.nanmean(pred_error_m1)) # centralized to calculate the residuals

        train_df_NN = pd.DataFrame(dctNN)
        train_df_NN = train_df_NN.dropna()
        train_df_NN = train_df_NN.reset_index()
        self.NN_train_df = train_df_NN

        self.mask = masks

    def trainNN(self, y_var, plot = False):
        '''
        run basic MLP pipeline on the training data.
        '''
        assert y_var in ['IR', 'PredErrorM1']

        X = self.NN_train_df.drop(['index', 'PredErrorM1' ,'IR'],1)
        if y_var == 'IR':
          X = self.NN_train_df.drop(['index', 'PredErrorM1' ,'IR', 'PredM1'],1)
        y = self.NN_train_df[y_var]

        nn_model = MLPRegressor(hidden_layer_sizes=(10, 10), max_iter=100, random_state=42)
        nn_model.fit(X, y)
        y_pred = nn_model.predict(X)
        r2_NN = r2_score(y, y_pred)
        self.NNModel = nn_model
        self.trainedR2_NN = r2_NN

        print(f'r^2:\tNN = {r2_NN:.3f}')

        if plot:
            rand = np.random.randint(0, len(y), 1000)
            sns.regplot(x = y[rand], y = y_pred[rand])
            plt.ylabel('Predicted')
            plt.xlabel('Real')
            plt.title('NN model')
            plt.show()

    def test(self, path):
        '''
        for maps in test set (not in train set), the function create dataframe from
        input maps and run the model.
        '''
        for flight in self.files_df['Flight'].unique():
            if flight in list(self.train_flights):
                continue
            subset = self.files_df[self.files_df['Flight'] == flight].reset_index()
            for crop, sub_flight in enumerate(subset['M1']):
                name = sub_flight.split('/')[-1][:-4]
                tgi = gdal.Open(subset['TGI'][crop]).ReadAsArray().flatten()
                height = gdal.Open(subset['height'][crop]).ReadAsArray().flatten()
                shade = gdal.Open(subset['shade'][crop]).ReadAsArray().flatten()
                real_solar = gdal.Open(subset['real_solar'][crop]).ReadAsArray().flatten()
                skyview = gdal.Open(subset['skyview'][crop]).ReadAsArray().flatten()
                pred_m1 = gdal.Open(subset['M1'][crop]).ReadAsArray().flatten()

                dataNN = pd.DataFrame({#'PredM1':pred_m1, #
                                     'TGI':tgi, 'Height':height,
                                     'Shade':shade, 'RealSolar':real_solar, 'Skyview':skyview})

                m2_map = self.NNModel.predict(dataNN).reshape((self.N, self.N))
                np.save(f'{path}/{name}_Model_NN.npy', m2_map)