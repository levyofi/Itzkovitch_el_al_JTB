library(raster)

wd = 'Example data/Input data'
cropped_wd = 'Example data/Input data/cropped_'

get_cropped_extent = function(image, ext){
  image_len = xmax(image) - xmin(image)
  residue_length = (image_len - (image_len / ncol(image) * ext)) / 2
  cropped_extent = extent(image)
  xmin(cropped_extent) = xmin(cropped_extent) + residue_length
  xmax(cropped_extent) = xmax(cropped_extent) - residue_length
  ymin(cropped_extent) = ymin(cropped_extent) + residue_length
  ymax(cropped_extent) = ymax(cropped_extent) - residue_length
  return(cropped_extent)
}

crop_and_save = function(image, center_extent, save_path){
  cropped_extent = get_cropped_extent(image, center_extent)
  cropped_image = crop(image, cropped_extent)
  writeRaster(cropped_image, file=save_path, overwrite=TRUE)
}

zones = list.dirs(wd, full.names = FALSE)
zones = zones[2:length(zones)]
for (path in zones){
  zone = path
  save_path = paste0(cropped_wd, zone)
  if (!file.exists(save_path)){
    dir.create(save_path)
  }
  print(zone)
  for (image in list.files(file.path(wd, path), pattern='.tif$')){
    image_raster = stack(file.path(wd, path, image))
    crop_and_save(image_raster, 1024, file.path(paste0(cropped_wd, zone), image))
  }
}