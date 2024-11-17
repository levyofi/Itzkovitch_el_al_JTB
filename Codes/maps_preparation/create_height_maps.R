library(sp)
library (raster)

wd = '/big_data/idan/subimages'

create_height = function(dsm_path, dtm_path){
  dsm = raster(dsm_path)
  dtm = raster(dtm_path)
  cropped_dtm = resample(crop(dtm, extent(dsm)), dsm, 'bilinear')
  return(dsm - cropped_dtm)
}

flights_df = read.csv('/data/idan/flights_info_final_with_dtm.csv')

for (flight in 1:nrow(flights_df)){
  dsm = flights_df[flight, 'dsm']
  zone = flights_df[flight, 'zone']
  split_dsm = strsplit(dsm, '/')[[1]]
  dtm = flights_df[flight, 'dtm']
  if (flights_df[flight, 'dsm_cropped'] == ''){
    is_cropped = ''
  } else {
    is_cropped = 'cropped_'
  }
  if ('good_projects' %in% split_dsm){
    extra_idx = 1
  }
  else {
    extra_idx = 0
  }
  print(zone)
  for (i in 1:5){
    dsm = paste(wd, zone, paste0(split_dsm[5+extra_idx], '_dsm_', is_cropped, 'downsampled_', i, '.tif'), sep='/')
    height = create_height(dsm, dtm)
    writeRaster(height, file=paste(wd, zone, paste0('height_', i, '.tif'), sep='/'), overwrite=TRUE)
  }
}