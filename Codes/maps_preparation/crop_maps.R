# Load the raster package for working with spatial raster data
library(raster)

# Define working directories
wd = 'Example data/Input data'  # Directory containing the original raster files
cropped_wd = 'Example data/Input data/cropped_'  # Directory to save cropped raster files

# Function to calculate the cropped extent of a raster
get_cropped_extent = function(image, ext) {
  # Calculate the length of the raster in the x-direction
  image_len = xmax(image) - xmin(image)
  
  # Calculate the amount to crop from each side of the raster
  residue_length = (image_len - (image_len / ncol(image) * ext)) / 2
  
  # Create a new extent object for the cropped raster
  cropped_extent = extent(image)
  xmin(cropped_extent) = xmin(cropped_extent) + residue_length  # Adjust minimum x-coordinate
  xmax(cropped_extent) = xmax(cropped_extent) - residue_length  # Adjust maximum x-coordinate
  ymin(cropped_extent) = ymin(cropped_extent) + residue_length  # Adjust minimum y-coordinate
  ymax(cropped_extent) = ymax(cropped_extent) - residue_length  # Adjust maximum y-coordinate
  
  return(cropped_extent)  # Return the cropped extent
}

# Function to crop a raster and save the cropped version
crop_and_save = function(image, center_extent, save_path) {
  # Get the cropped extent for the raster
  cropped_extent = get_cropped_extent(image, center_extent)
  
  # Crop the raster to the calculated extent
  cropped_image = crop(image, cropped_extent)
  
  # Save the cropped raster to the specified path
  writeRaster(cropped_image, file=save_path, overwrite=TRUE)
}

# Get the list of directories (zones) inside the working directory
zones = list.dirs(wd, full.names = FALSE)
zones = zones[2:length(zones)]  # Remove the root directory from the list

# Loop through each zone to process the raster files
for (path in zones) {
  zone = path  # Current zone name
  
  # Define the save path for the cropped rasters in this zone
  save_path = paste0(cropped_wd, zone)
  
  # Create the save directory if it does not exist
  if (!file.exists(save_path)) {
    dir.create(save_path)
  }
  
  # Print the current zone being processed
  print(zone)
  
  # Loop through all .tif raster files in the current zone directory
  for (image in list.files(file.path(wd, path), pattern='.tif$')) {
    # Load the raster file as a stack (to handle multi-layer rasters)
    image_raster = stack(file.path(wd, path, image))
    
    # Crop the raster and save the result
    crop_and_save(image_raster, 1024, file.path(paste0(cropped_wd, zone), image))
  }
}
