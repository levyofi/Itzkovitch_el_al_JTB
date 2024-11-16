import rasterio as rs
from rasterio.windows import Window
import os
from random_subimages import RandomSubimages
from resample_images import ResampleImages

              
class RandomSubimagesAllMaps:
    # Class for creating random subimages for all the maps in the given directory
    
    def __init__(self, flights_info, output_path):
        # Initializes the input and output paths
        self.flights = flights_info
        self.output_path = output_path
        
    def save_subimages(self, paths, flight, coordinates, size):
        # Saves a list of subimages in a given flight directory in self.output_path
        if not os.path.exists(os.path.join(self.output_path, flight)):
            os.makedirs(os.path.join(self.output_path, flight), exist_ok=True)
        for path in paths:
            image_filename = os.path.basename(path)[:-4]
            with rs.open(path) as image:
                for index, coordinate in enumerate(coordinates):
                    window = Window(*coordinate[::-1], size, size)
                    transform = image.window_transform(window)
                    profile = image.profile
                    profile.update({'height': size, 'width': size, 'transform': transform})
                    subimage_path = os.path.join(self.output_path, flight, f'{image_filename}_{index + 1}.tif')
                    with rs.open(subimage_path, 'w', **profile) as subimage:
                        subimage.write(image.read(window=window))
    
    def map_path(self, f, map_type):
        # Returns the map path of a given flight
        cropped_path = map_type + '_cropped'
        if self.flights[cropped_path][f] == self.flights[cropped_path][f]:
            return self.flights[cropped_path][f]
        return self.flights[map_type][f]
    
    def create_subimages_all_dir(self, num, size, resample_type='downsample'):
        # Creates subimages with RandomSubimages for all the maps in self.flights_info
        resamples_dir = os.path.join(os.path.dirname(self.output_path), 'resampled_maps')
        if not os.path.exists(resamples_dir):
            os.mkdir(resamples_dir)
        for f in range(self.flights.shape[0]):
            rgb_path = self.map_path(f, 'rgb')
            dsm_path = self.map_path(f, 'dsm')
            ir_path = self.map_path(f, 'ir')
            print(rgb_path.split('/')[4])
            paths = ResampleImages([rgb_path, dsm_path, ir_path], resample_type).resample_images(resamples_dir)
            coordinates = RandomSubimages(paths[0]).create_subimages(num, size)
            self.save_subimages(paths, self.flights['zone'][f], coordinates, size)
        print('Done!')
