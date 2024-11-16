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

get_value_lat_lon <- function(sp, varname, nc_name) {
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
  extract(varraster, geom(sp))
}

# Read input data
data_files <- read.table("Example data/Input data/flights_info_final.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Initialize an empty dataframe to store results
meteo <- data.frame()

#For the guthub example - choose the relevant flight
is_example=TRUE #change to FALSE when working on all maps
if (is_example){
  data_files=data_files[31,]
}

# Loop through each row of the data_files table
for (ifile in 1:nrow(data_files)) { 
  map <- basename(data_files[ifile, ]$map)
  print(paste("Processing map:", map))
  
  # Extract datetime 
  dt <- as_datetime(word(map, 2, 3, '_'), format = '%d.%m.%y_%H%M', tz = paste0("Etc/GMT-2"))
  #change time to UTC
  dt <- with_tz(dt, tzone = "UTC")
  
  # Extract location name
  location <- word(map, 1, 1, '_')
  
  # Define working directory if running on all folders
  if (!is_example){
    dir <- paste0("Example data/Input data/", map)
  } else {
    dir <- paste0("Example data/Input data/")
  }
  setwd(dir)
  
  # Load and process rasters
  temp <- rast("IR.tif")
  solar <- rast("real_solar.tif")
  shade <- rast("shade.tif")
  slope <- rast("SLP.tif")
  cleansky_solar <- rast("solar.tiff")
  
  # Meteorological station data
  station <- data_files[ifile, ]
  WIND <- station$wind_speed # m/s
  TG <- station$soil + 273.15 # Kelvin
  RH <- station$rh # %
  pressure <- station$pressure # hPa
  TAIR <- station$temp + 273.15 # Kelvin
  QAIR <- specific_humidity(RH, TAIR)
  
  # Cloud cover calculation
  clean_solar <- max(values(cleansky_solar[slope < 0.1 & shade == 0]), na.rm = TRUE)
  real_solar <- station$radiation
  cloud_cover <- max(0, 1 - real_solar / clean_solar)
  
  # NetCDF data integration
  sp1 <- vect(cbind(lon = station$longitude, lat = station$latitude), crs = "EPSG:4326")
  nc_name <- list.files(path = "../../rasters/", pattern = paste0("^.*", format(dt, "%Y%m%d.%H00"), ".*.nc4"), full.names = TRUE)
  
  ALBEDO <- get_value_lat_lon(sp1, "Albedo_inst", nc_name)
  soil_temp10 <- get_value_lat_lon(sp1, "SoilTMP0_10cm_inst", nc_name)
  soil_temp40 <- get_value_lat_lon(sp1, "SoilTMP10_40cm_inst", nc_name)
  soil_temp100 <- get_value_lat_lon(sp1, "SoilTMP40_100cm_inst", nc_name)
  soil_temp200 <- get_value_lat_lon(sp1, "SoilTMP100_200cm_inst", nc_name)
  soil_mois10 <- get_value_lat_lon(sp1, "SoilMoi0_10cm_inst", nc_name) / 100
  soil_mois40 <- get_value_lat_lon(sp1, "SoilMoi10_40cm_inst", nc_name) / 300
  soil_mois100 <- get_value_lat_lon(sp1, "SoilMoi40_100cm_inst", nc_name) / 600
  soil_mois200 <- get_value_lat_lon(sp1, "SoilMoi100_200cm_inst", nc_name) / 1000
  
  # Compile row for meteo dataframe
  meteo_row <- data.frame(
    map = map, Date = dt, Latitude = station$latitude, Longitude = station$longitude,
    TG = TG, Albedo = ALBEDO,
    ST10 = soil_temp10, ST40 = soil_temp40, ST100 = soil_temp100, ST200 = soil_temp200,
    SM10 = soil_mois10, SM40 = soil_mois40, SM100 = soil_mois100, SM200 = soil_mois200,
    Wind = WIND, Pressure = pressure, TAIR = TAIR, QAIR = QAIR, Cloud_cover = cloud_cover
  )
  
  # Append to meteo dataframe
  meteo <- rbind(meteo, meteo_row)
}

# Save the consolidated meteorological data to a CSV file
write.csv(meteo, "consolidated_meteorological_data.csv", row.names = FALSE)

print("Processing complete. Consolidated data saved.")
