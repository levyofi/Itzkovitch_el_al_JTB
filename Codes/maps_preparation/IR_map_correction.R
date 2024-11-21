library(raster)  # For handling raster data

pwd = dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

flight = 'Zeelim_31.05.21_1516/'
map_name = 'IX-12-11646_0206_index_thermal_ir_cropped_5.tif'

#Load the thermal image file
thermal_image <- raster(paste0(pwd, "/Example data/Input data/", flight, map_name))

#Load the CSV file and extract the correction value
station_data <- fread(paste0(pwd, "/Example data/Tables/station_correction_values.csv"))
correction_value <- station_data$`station - IR`[120]

#Apply correction
corrected_image <- thermal_image + correction_value

#Save the corrected image
output_filename <- paste0(pwd, "/Example data/Input data/", flight, 'corrected_', map_name)
writeRaster(corrected_image, output_filename, format = "GTiff", overwrite = TRUE)

cat("Corrected image saved to:", output_filename, "\n")
