library(raster)  # For handling raster data

pwd = dirname(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

flight = 'Zeelim_31.05.21_1516/'
map_name = 'IX-12-11646_0206_index_thermal_ir_cropped_5.tif'
# Step 1: Load the thermal image file
thermal_image <- raster(paste0(pwd, "/Example data/Input data/", flight, map_name))

# Step 2: Load the CSV file and extract the 'station-map' value
station_data <- fread(paste0(pwd, "/Example data/Tables/station_correction_values.csv"))
correction_value <- station_data$`station - IR`[120]

# Step 3: Apply correction (example: subtracting the station-map value)
corrected_image <- thermal_image + correction_value

# Step 4: Save the corrected image
output_filename <- paste0(pwd, "/Example data/Input data/", flight, 'corrected_', map_name)
writeRaster(corrected_image, output_filename, format = "GTiff", overwrite = TRUE)

cat("Corrected image saved to:", output_filename, "\n")
