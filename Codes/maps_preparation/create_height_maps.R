library(sp)
library (raster)

wd = 'Example data/Input data'

create_height = function(dsm_path, dtm_path){
  dsm = raster(dsm_path)
  dtm = raster(dtm_path)
  cropped_dtm = resample(crop(dtm, extent(dsm)), dsm, 'bilinear')
  return(dsm - cropped_dtm)
}

flights_df = read.csv('Example data/Input data/flights_info_final.csv')

#The code assumes that there is a folder with the name of each map and the original dtm and the cropped dsm maps are in this folder
for (flight in 1:nrow(flights_df)){
  dsm = flights_df[flight, 'dsm']
  map = flights_df[flight, 'map']
  split_dsm = strsplit(dsm, '/')[[1]]
  dtm = flights_df[flight, 'dtm']
  print(zone)
  for (i in 1:5){
    dsm = paste(wd, map, paste0(split_dsm[5+extra_idx], '_dsm_', is_cropped, 'downsampled_', i, '.tif'), sep='/')
    height = create_height(dsm, dtm)
    writeRaster(height, file=paste(wd, zone, paste0('height_', i, '.tif'), sep='/'), overwrite=TRUE)
  }
}