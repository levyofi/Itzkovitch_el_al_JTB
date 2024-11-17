from random_subimages_all_maps import RandomSubimagesAllMaps
import pandas as pd

# Load the flight information file and define the output path
flights_info = pd.read_csv('flights_info_final.csv')
output_path = 'Example data/Input maps/'

# Initialize the subimage generator
subimage_generator = RandomSubimagesAllMaps(flights_info, output_path)

# Generate and save subimages
subimage_generator.create_subimages_all_dir(num=6, size=1100, resample_type='downsample')
