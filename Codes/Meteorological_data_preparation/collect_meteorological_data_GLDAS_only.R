# Load necessary libraries
library(terra)
library(lubridate)
library(stringr)
library(ncdf4)
library(humidity)
library(bigleaf)

# Define helper functions
specific_humidity <- function(RH, T) { # RH - relative humidity [0-100 %]; T - air temperature [Kelvin]
  qair <- SH(vapor_pressure(RH, T))
  qair
}

vapor_pressure <- function(RH, T) { # RH - relative humidity [0-100 %]; T - air temperature [Kelvin]
  SVP <- SVP(T) # return hPa 
  SVP * RH / 100 * 100 # return Pa
}

raster_fix_NA <- function(rs) {
  rs_median <- global(rs, median, na.rm = TRUE)$median
  rs[is.na(rs[])] <- rs_median
  return(rs)
}

get_value_lat_lon <- function(coord, varname, nc_name) {
  nc_file <- nc_open(nc_name)
  varlayer <- ncvar_get(nc_file, varid = varname)
  lat <- ncvar_get(nc_file, varid = "lat")
  lon <- ncvar_get(nc_file, varid = "lon")
  lon[lon > 180] <- lon[lon > 180] - 360
  fillvalue <- ncatt_get(nc_file, varid = 0, attname = "missing_value")
  nc_close(nc_file)
  varlayer[varlayer == fillvalue$value] <- NA
  varraster <- rast(t(varlayer), extent = ext(min(lon), max(lon), min(lat), max(lat)), crs = "EPSG:4326")
  varraster <- flip(varraster, direction = 'vertical')
  extract(varraster, coord)[,2]
}

# Read input data
data_files <- read.table("Example data/Input data/flights_info_final.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Initialize an empty dataframe to store results
meteo <- data.frame()

#For the github example - choose the relevant flight
is_example=TRUE #change to FALSE when working on all maps
if (is_example){
  data_files=data_files[data_files$map=="Zeelim_31.05.21_1516",]
}

# Loop through each row of the data_files table
for (ifile in 1:nrow(data_files)) { 
  map <- basename(data_files[ifile, ]$map)
  print(paste("Processing map:", map))
  
  # Extract datetime 
  dt <- as_datetime(word(map, 2, 3, '_'), format = '%d.%m.%y_%H%M', tz = paste0("Etc/GMT-2"))
  # Change time to UTC
  dt <- with_tz(dt, tzone = "UTC")
  # Round to the nearest 3-hour interval
  dt_rounded <- round_date(dt, unit = "3 hours")
  
  # Extract location name
  location <- word(map, 1, 1, '_')
  
  # Define working directory if running on all folders
  dir <- paste0("Example data/Input data/", map)
  setwd(dir)
  
  # Get relevant rasters for the calculation of cloud cover and coordinates
  solar <- rast("real_solar_1.tif")
  slope <- rast("slope_1.tif")
  cleansky_solar <- rast("solar_1.tif")
  shade <- rast("shade_1.tif")
  
  # Get the extent of the rast object
  e <- ext(solar)
  # Calculate the middle coordinates
  middle_x <- (e$xmin + e$xmax) / 2  # Midpoint of longitude
  middle_y <- (e$ymin + e$ymax) / 2  # Midpoint of latitude
  middle_coords <- vect(cbind(lon = middle_x, lat = middle_y), crs = "EPSG:32636")
  middle_coords_latlon <- project(middle_coords, "EPSG:4326")
  
  # Get the GLDAS file name
  nc_name <- list.files(path = "../", 
                        pattern = paste0("^.*", format(dt_rounded, "%Y%m%d.%H00"), ".*.nc4"), 
                        full.names = TRUE)
  
  #### Meteorological variables ####
  # Extract meteorological variables from GLDAS
  WIND <- get_value_lat_lon(middle_coords_latlon, "Wind_f_inst", nc_name)
  pressure <- get_value_lat_lon(middle_coords_latlon, "Psurf_f_inst", nc_name)
  TAIR <- get_value_lat_lon(middle_coords_latlon, "Tair_f_inst", nc_name)
  QAIR <- get_value_lat_lon(middle_coords_latlon, "Qair_f_inst", nc_name)
  
  # Additional radiation and solar variables
  SWD <- get_value_lat_lon(middle_coords_latlon, "SWdown_f_tavg", nc_name)
  SWN <- get_value_lat_lon(middle_coords_latlon, "Swnet_tavg", nc_name)
  real_solar <- SWD
  clean_solar <- max(cleansky_solar[slope < 0.1 & shade == 0], na.rm = TRUE)
  cloud_cover <- max(0, 1 - real_solar / clean_solar)
  
  #### Extract GLDAS soil and albedo data ####
  ALBEDO <- get_value_lat_lon(middle_coords_latlon, "Albedo_inst", nc_name)
  soil_temp10 <- get_value_lat_lon(middle_coords_latlon, "SoilTMP0_10cm_inst", nc_name)
  soil_temp40 <- get_value_lat_lon(middle_coords_latlon, "SoilTMP10_40cm_inst", nc_name)
  soil_temp100 <- get_value_lat_lon(middle_coords_latlon, "SoilTMP40_100cm_inst", nc_name)
  soil_temp200 <- get_value_lat_lon(middle_coords_latlon, "SoilTMP100_200cm_inst", nc_name)
  soil_mois10 <- get_value_lat_lon(middle_coords_latlon, "SoilMoi0_10cm_inst", nc_name) / 100
  soil_mois40 <- get_value_lat_lon(middle_coords_latlon, "SoilMoi10_40cm_inst", nc_name) / 300
  soil_mois100 <- get_value_lat_lon(middle_coords_latlon, "SoilMoi40_100cm_inst", nc_name) / 600
  soil_mois200 <- get_value_lat_lon(middle_coords_latlon, "SoilMoi100_200cm_inst", nc_name) / 1000
  
  #### Compile row for meteo dataframe ####
  meteo_row <- data.frame(
    map = map, Date = dt, 
    TG = data_files[ifile, ]$ground_temp + 273.15,  # Convert to Kelvin
    Albedo = ALBEDO,
    ST10 = soil_temp10, ST40 = soil_temp40, ST100 = soil_temp100, ST200 = soil_temp200,
    SM10 = soil_mois10, SM40 = soil_mois40, SM100 = soil_mois100, SM200 = soil_mois200,
    Wind = WIND, Pressure = pressure, TAIR = TAIR, QAIR = QAIR, Cloud_cover = cloud_cover
  )
  
  # Append to meteo dataframe
  meteo <- rbind(meteo, meteo_row)
}

# Save the consolidated meteorological data to a CSV file
write.csv(meteo, "input_meteorological_data_GLDAS_only.csv", row.names = FALSE)

print("Processing complete. Consolidated data saved.")
