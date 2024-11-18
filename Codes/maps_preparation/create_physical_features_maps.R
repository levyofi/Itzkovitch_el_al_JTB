# Load necessary libraries for data processing and spatial analysis
library(leaflet)
library(stringr)
library(lubridate)
library(humidity)
library(bigleaf)
library(suncalc)
library(sf)
library(rgrass)
library(terra)

# Set the home directory and working directory
homedir <- "Example data/Input data/"
setwd(homedir)

# Read flight information data from CSV
data_files <- read.table("flights_info_final.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)

# Use a specific map for an example (for GitHub demonstration purposes)
is_example <- TRUE  # Set to FALSE to process all maps
if (is_example) {
  data_files <- data_files[data_files$map == "Zeelim_31.05.21_1516", ]  # Filter data for a specific map example
}

# Iterate over each flight (row in the data table)
for (ifile in 1:nrow(data_files)) {
  
  # Extract map name and parse datetime information
  map <- basename(data_files[ifile, ]$map)
  dt <- as_datetime(word(map, 2, 3, '_'), format = '%d.%m.%y_%H%M', tz = paste0("Etc/GMT-", 2))
  dt <- with_tz(dt, tzone = "UTC")
  
  # Check if directory exists for the map
  dir <- map
  if (dir.exists(dir)) {
    print("DIR EXISTS!")
  }
  
  setwd(dir)  # Change working directory to the map's folder
  files <- list.files(".", full.names = TRUE)  # List all files in the folder
  
  # Filter and identify DSM and RGB files
  dsm_names <- files[grep("dsm.*\\.tif$", files, ignore.case = TRUE)]
  rgb_files <- files[grep("mosaic.*\\.tif$", files, ignore.case = TRUE)]
  
  # Load the first DSM file
  dsm <- rast(dsm_names[1])
  
  # Retrieve projection information from the DSM
  p <- crs(dsm, proj = TRUE, describe = TRUE, parse = TRUE)
  
  # Initialize GRASS GIS and set the projection
  gisBase <- "/usr/local/grass85"  # Path to GRASS GIS installation
  gisDbase <- map
  location <- "microclimate"
  mapset <- "PERMANENT"
  
  # Start GRASS GIS session and set the projection
  loc <- initGRASS(gisBase = gisBase, home = getwd(), gisDbase = "GRASS_R", override = TRUE, mapset = "PERMANENT")
  execGRASS("g.proj", flags = "c", parameters = list(proj4 = as.character(p$proj)))
  
  # Process multiple sampled maps (1 to 5 for this example)
  for (i in 1:5) {
    # Set output file names for various results
    filename_incidence <- paste0("incidence_", i, ".tif")
    filename_solar <- paste0("solar_", i, ".tif")
    filename_real_solar <- paste0("real_solar_", i, ".tif")
    filename_skyview <- paste0("skyview_", i, ".tif")
    filename_slope <- paste0("slope_", i, ".tif")
    filename_shade <- paste0("shade_", i, ".tif")
    filename_tgi <- paste0("TGI_", i, ".tif")
    
    # Calculate the solar time
    sp1 <- SpatialPoints(xyFromCell(dsm, 1), proj4string = CRS(as.character(p$proj)))
    sp1degrees <- spTransform(sp1, CRS("+proj=longlat +datum=WGS84"))
    lat <- coordinates(sp1degrees)[1, 2]
    lon <- coordinates(sp1degrees)[1, 1]
    
    solar_noon <- getSunlightTimes(date = date(dt), lat = lat, lon = lon, keep = c("solarNoon"), tz = "UTC")
    solar_time <- 12 - hour(solar_noon$solarNoon) - minute(solar_noon$solarNoon) / 60 + (hour(dt) + minute(dt) / 60)
    solar_time <- solar_time - 0.5  # Adjust solar time
    
    # Import DSM to GRASS and calculate slope and aspect, and save the slope raster to a file
    execGRASS("r.in.gdal", parameters = list(input = dsm_names[i], output = "dsm"), flags = c("quiet", "overwrite"))
    execGRASS("r.slope.aspect", flags = c("overwrite"), parameters = list(elevation = "dsm", slope = "slope", aspect = "aspect"))
    execGRASS("r.out.gdal", flags = "overwrite", parameters = list(createopt = "PROFILE=GeoTIFF,TFW=YES", output = filename_slope, input = "slope"))
    
    # Calculate radiation incidence angle
    execGRASS("r.sun", flags = "overwrite", parameters = list(elevation = "dsm", day = yday(date(dt)), time = solar_time, aspect = "aspect", slope = "slope", beam_rad = "tmpbeam", incidout = "tmpincidout"))
    execGRASS("r.out.gdal", flags = "overwrite", parameters = list(createopt = "PROFILE=GeoTIFF,TFW=YES", output = filename_incidence, input = "tmpincidout"))
    incident <- rast(filename_incidence)
    
    # Generate global solar radiation map
    execGRASS("r.sun", flags = "overwrite", parameters = list(elevation = "dsm", day = yday(date(dt)), time = solar_time, aspect = "aspect", slope = "slope", glob_rad = "tmpglob", albedo_value = 0))
    execGRASS("r.out.gdal", flags = "overwrite", parameters = list(createopt = "PROFILE=GeoTIFF,TFW=YES", output = filename_solar, input = "tmpglob"))
    
    # Create shade map from incidence data
    shade <- incident  # Create a copy of the incident raster
    shade[!is.na(shade)] <- 0  # Set all non-NA values to 0 (not shaded)
    shade[is.na(shade)] <- 1  # Set NA values to 1 (shaded)
    writeRaster(shade, filename_shade, overwrite = TRUE)
    
    # Generate vegetation index (TGI) map from RGB data
    r <- stack(rgb_files[1])
    R <- r[[1]]; G <- r[[2]]; B <- r[[3]]
    r_tgi <- (G - 0.39 * R - 0.61 * B) / max(R, max(G, B))
    writeRaster(r_tgi, file = filename_tgi, overwrite = TRUE)
    
    # Convert solar radiation to real solar radiation
    solar <- rast(filename_solar)
    real_solar <- data_files[ifile, ]$radiation
    slope <- rast(filename_slope)
    clean_solar <- max(solar[slope < 0.1 & shade == 0], na.rm = TRUE)
    
    if (clean_solar < real_solar) print("Warning: clean_solar is lower than measured solar")
    r_real_solar <- solar * real_solar / clean_solar
    writeRaster(r_real_solar, filename_real_solar, overwrite = TRUE)
    
    # Generate skyview map (requires r.skyview extension in GRASS GIS). 
    # At first usage, please install r.skyview by running "g.extension extension=r.skyview" in grass terminal
    execGRASS("r.skyview", flags = "overwrite", parameters = list(input = "dsm", output = "tmpskyview"))
    execGRASS("r.out.gdal", flags = "overwrite", parameters = list(createopt = "PROFILE=GeoTIFF,TFW=YES", output = filename_skyview, input = "tmpskyview"))
  }
  
  # Return to the parent directory after processing the map
  setwd("..")
}
print("Finished processing all maps")

