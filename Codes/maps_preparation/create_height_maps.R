library(sp)
library (raster)
library(tools)

wd = 'Example data/Input data'

create_height = function(dsm_path, dtm_path){
  dsm = raster(dsm_path)
  dtm = raster(dtm_path)
  cropped_dtm = resample(crop(dtm, extent(dsm)), dsm, 'bilinear')
  return(dsm - cropped_dtm)
}

flights_df = read.csv('Example data/Input data/flights_info_final.csv')

# For the GitHub example:
# If working on a specific map (e.g., for an example), set `is_example = TRUE`.
# Otherwise, set `is_example = FALSE` to process all maps.
is_example = TRUE  # Change to FALSE when working on all maps
if (is_example) {
  flights_df = flights_df[flights_df$map=="Zeelim_31.05.21_1516", ]  # Select the 31st row as an example (adjust index as needed)
}

#The code assumes that there is a folder with the name of each map and the original dtm and the cropped dsm maps are in this folder
for (flight in 1:nrow(flights_df)){
  dsm = flights_df[flight, 'dsm']
  map = flights_df[flight, 'map']
  base_dsm = file_path_sans_ext(basename(dsm))
  dtm = file.path(wd, map, basename(flights_df[flight, 'dtm']))
  print(map)
  for (i in 1:5){
    dsm = paste(wd, map, paste0(base_dsm, '_', 'cropped_downsampled_', i, '.tif'), sep='/')
    height = create_height(dsm, dtm)
    writeRaster(height, file=paste(wd, zone, paste0('height_', i, '.tif'), sep='/'), overwrite=TRUE)
  }
}