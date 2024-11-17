# Instructions
## Generating random square subimages from an image in `random_subimages.py`
This script is designed to generate random square subimages from a larger image. These subimages are created with minimal congruence to ensure variety.
### Steps to Use:
#### 1. Initialize the RandomSubimages Class:
- Begin by creating an instance of the `RandomSubimages` class, providing the path to the image file:
  ```python
  from random_subimages import RandomSubimages
  subimage_generator = RandomSubimages('path/to/your/image.tif')
#### 2. Create Subimages:
- Use the `create_subimages` method to generate a specified number of random subimages. You need to define the size of the subimages and the allowable congruence percentage (`cong_percent`).
- Optionally, enable plotting to visualize the generated subimages on the original image:
  ```python
  coordinates = subimage_generator.create_subimages(num=10, size=100, cong_percent=0.2, plot=True)
#### 4. Plot the Subimages:
- If you did not plot the subimages in the previous step, you can use the `plot_subimages` method to display the subimages generated:
  ```python
  subimage_generator.plot_subimages(size=100)
### Key Methods:
- `__init__(self, path)`: Initializes the class with the image path and prepares the image for subimage extraction.
- `path_to_image(self)`: Reads the image from the provided path.
- `is_expandable(self, coordinate, size)`: Checks if a given coordinate can be expanded to a square of the specified size without exceeding the image boundaries.
- `is_congruential(self, coordinate, size, cong_percent)`: Ensures that the new subimage has a minimal congruence area with previously created subimages.
- `resize(self, size, image)`: Adjusts the subimage size to fit within the image dimensions if necessary.
- `create_subimages(self, num, size, cong_percent, plot)`: Generates the desired number of subimages, adhering to the congruence constraint.
- `plot_subimages(self, size)`: Visualizes the generated subimages on the original image.

## Resampling images for a given flight in `resample_images.py`
This script provides functionality for resampling a set of images associated with a flight. This can be useful for upscaling or downscaling images to match the resolution of other datasets or to meet specific analysis requirements.
### Steps to Use:
#### 1. Initialize the `ResampleImages` Class:
- Create an instance of the `ResampleImages` class by providing a list of paths to the images and specifying whether you want to `upsample` or `downsample`:
  ```python
  from resample_images import ResampleImages
  
  paths = [
      'path/to/RGB_image.tif',
      'path/to/DSM_image.tif',
      'path/to/thermal_image.tif'
  ]
  resample_type = 'downsample'  # or 'upsample'
  resampler = ResampleImages(paths, resample_type)
#### 2. Resample the Images:
- Use the `resample_images` method to resample the images based on the initialization settings. Provide a directory path where the resampled images should be saved.
- This method returns a list of paths to the updated (resampled) images:
  ```python
  resampled_paths = resampler.resample_images('path/to/save/resampled_images')
## Key Methods:
- `__init__(self, paths, resample_type)`: Initializes the class with the image paths and the desired resampling type (`'upsample'` or `'downsample'`).
- `resample_images_state(self)`: Determines which images need to be resampled based on the specified resampling type. This method returns a boolean dictionary where each path corresponds to a `True` or `False` value, indicating whether the image should be resampled.
- `get_scale_factor(self)`: Computes the scale factor for resampling based on the ratio of image dimensions. This factor is used to adjust the image resolution.
- `resample_images(self, resamples_dir)`: Performs the actual resampling of images. The method reads the images, resamples them according to the scale factor, and saves the resampled images to the specified directory. It returns a list of paths to the updated (resampled) images.
### Example Usage:
```python
# Initialize the ResampleImages class with paths and the resampling type
resampler = ResampleImages(
    paths=['path/to/RGB_image.tif', 'path/to/DSM_image.tif', 'path/to/thermal_image.tif'],
    resample_type='downsample'
)

# Resample the images and save them to the specified directory
resampled_paths = resampler.resample_images('path/to/save/resampled_images')

# The resampled_paths list now contains the paths to the resampled images
```
## Generate random subimages across multiple maps using `create_subimages_all_maps.py`
This script is designed to generate random subimages across various maps from a specified directory. This can be particularly useful for extracting smaller image patches for analysis or training machine learning models. The script integrates the functionalities of the `RandomSubimages` and `ResampleImages` classes to achieve this.
### Steps to Use:
#### 1. Initialize the `RandomSubimagesAllMaps` Class:
- Create an instance of the `RandomSubimagesAllMaps` class by providing information about the flight maps and the output directory where the subimages will be saved:
  ```python
  from random_subimages_all_maps import RandomSubimagesAllMaps
  
  flights_info = {
      # Add paths to your flight maps here
      'rgb': [...],
      'dsm': [...],
      'ir': [...],
      'rgb_cropped': [...],
      'dsm_cropped': [...],
      'ir_cropped': [...],
      'zone': [...],  # List of zones corresponding to the flight maps
  }
  output_path = 'path/to/save/subimages'
  
  subimage_generator = RandomSubimagesAllMaps(flights_info, output_path)
#### 2. Generate and Save Subimages:
- Use the `create_subimages_all_dir` method to generate random subimages for all maps in the provided directory. Specify the number of subimages, the size of each subimage, and the resampling type (`'downsample'` or `'upsample'`):
  ```python
  num_subimages = 100  # Number of subimages to generate per map
  subimage_size = 32   # Size of each subimage (e.g., 32x32 pixels)
  resample_type = 'downsample'  # or 'upsample'
  
  subimage_generator.create_subimages_all_dir(num_subimages, subimage_size, resample_type)
## Key Methods:
- `__init__(self, flights_info, output_path)`: Initializes the class with flight information (paths to maps) and the output directory where the subimages will be saved.
- `save_subimages(self, paths, flight, coordinates, size)`: Saves a list of subimages generated from the provided paths in the specified flight directory within the output path.
- `map_path(self, f, map_type)`: Returns the path of the specified map type (e.g., `'rgb'`, `'dsm'`, `'ir'`) for a given flight. If a cropped version of the map exists, it returns the path to the cropped map.
- `create_subimages_all_dir(self, num, size, resample_type='downsample')`: Creates random subimages for all maps in the directory. This method integrates the `RandomSubimages` class to generate subimages and the `ResampleImages` class to handle resampling if needed.
### Example Usage:
```python
from random_subimages_all_maps import RandomSubimagesAllMaps

# Define the flight information and output path
flights_info = {
    'rgb': ['path/to/RGB_map1.tif', 'path/to/RGB_map2.tif'],
    'dsm': ['path/to/DSM_map1.tif', 'path/to/DSM_map2.tif'],
    'ir': ['path/to/IR_map1.tif', 'path/to/IR_map2.tif'],
    'zone': ['zone1', 'zone2']
}
output_path = 'path/to/save/subimages'

# Initialize the subimage generator
subimage_generator = RandomSubimagesAllMaps(flights_info, output_path)

# Generate and save subimages
subimage_generator.create_subimages_all_dir(num=100, size=32, resample_type='downsample')
```

