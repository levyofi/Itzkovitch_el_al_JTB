library(httr)
library(lubridate)
library(stringr)

# Before using the script:
# Set up your ~/.netrc file for authentication as described here:
# https://disc.gsfc.nasa.gov/information/howto?title=How%20to%20Generate%20Earthdata%20Prerequisite%20Files
netrc_path <- ".netrc"  # Path to the .netrc file containing authentication information
cookie_path <- ".urs_cookies"  # Path to the cookies file for session handling
gldas_str_for_file = "GLDAS_NOAH025_3H.A"  # Base string for constructing GLDAS file names
gldas_str_for_path = "GLDAS_NOAH025_3H.2.1"  # Base path for GLDAS data on the server

# Read input data:
# The input data file contains flight information, including map names and datetime metadata.
data_files <- read.table("Example data/Input data/flights_info_final.csv", 
                         sep = ",", header = TRUE, stringsAsFactors = FALSE)

# For the GitHub example:
# If working on a specific map (e.g., for an example), set `is_example = TRUE`.
# Otherwise, set `is_example = FALSE` to process all maps.
is_example = TRUE  # Change to FALSE when working on all maps
if (is_example) {
  data_files = data_files[31, ]  # Select the 31st row as an example (adjust index as needed)
}

# Loop through each row of the data_files table:
# This processes each map specified in the input data file.
for (ifile in 1:nrow(data_files)) { 
  # Extract the map name from the current row
  map <- basename(data_files[ifile, ]$map)
  print(paste("Processing map:", map))  # Log the current map being processed
  
  # Extract datetime from the map name (assumes specific format with underscores)
  dt <- as_datetime(word(map, 2, 3, '_'), format = '%d.%m.%y_%H%M', tz = paste0("Etc/GMT-2"))
  
  # Convert datetime to UTC time zone
  dt <- with_tz(dt, tzone = "UTC")
  
  # Extract date and Julian day from the datetime
  date = date(dt)
  julian_day = yday(date)
  hour_3H = hour(dt)-hour(dt)%%3
  # Construct the GLDAS file name based on the extracted datetime
  str_file = paste0(gldas_str_for_file, year(date), 
                    str_pad(month(date), width = 2, pad = "0"), 
                    str_pad(day(dt), width = 2, pad = "0"), ".", 
                    str_pad(hour_3H, width = 2, pad = "0"), "00.021.nc4")
  
  # Construct the path string for the file on the server
  file_str = paste(year(date), str_pad(julian_day, width = 3, pad = "0"), str_file, sep = "/")
  
  # Set the local file path for downloading
  downloaded_file_path <- str_file
  
  # Configure the session for downloading:
  # Enables session handling using .netrc and cookies.
  set_config(config(followlocation = 1, netrc = 1, 
                    netrc_file = netrc_path, 
                    cookie = cookie_path, 
                    cookiefile = cookie_path, 
                    cookiejar = cookie_path))
  
  # Download the file from the server using the constructed URL and save it locally
  httr::GET(url = paste0("https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/", gldas_str_for_path, "/", file_str),
            write_disk(downloaded_file_path, overwrite = TRUE))
  
  # Log the download success
  print(paste("Downloaded", downloaded_file_path))
}

# Indicate completion of all downloads
print("Done.")
