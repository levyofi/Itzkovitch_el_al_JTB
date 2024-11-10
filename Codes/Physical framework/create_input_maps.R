library(httr)
library(lubridate)
library(stringr)
library(sp)
library(raster)
library(ncdf4)
library(lubridate)
library(rgdal)
library(raster)
library(sf)
library(rgrass7)
library(lubridate)
library(suncalc)

working_path <- '../' # Change it to wanted directory - NTBC

date = "18/09/2019" # args[1] # d/m/y format for the relevant data - NTBC
time = "13:00" # args[2] # hh:mm format for the relevant data - NTBC

date = dmy(date)
julian_day= yday(date)
time = hm(time)
gldas_str_for_file = ifelse(year(date)>=2022,"GLDAS_NOAH025_3H_EP.A", "GLDAS_NOAH025_3H.A")
str_file = paste0(gldas_str_for_file,year(date),str_pad(month(date), width = 2, pad = "0"),str_pad(day(date), width = 2, pad = "0"),".",str_pad(hour(date), width = 2, pad = "0"),"00.021.nc4" )
file_str=paste(year(date),str_pad(julian_day, width = 3, pad = "0"), str_file, sep="/")
netrc_path <- paste0(working_path, ".netrc")
cookie_path <- paste0(working_path, ".urs_cookies")
downloaded_file_path <- paste0(working_path, 'Data/', str_file)
# Before using the script
#Set up your ~/.netrc file as listed here: https://wiki.earthdata.nasa.gov/display/EL/How+To+Access+Data+With+cURL+And+Wget
set_config(config(followlocation=1,netrc=1,netrc_file=netrc_path,cookie=cookie_path,cookiefile=cookie_path,cookiejar=cookie_path))
gldas_str_for_path = ifelse(year(date)>=2022, "GLDAS_NOAH025_3H_EP.2.1","GLDAS_NOAH025_3H.2.1")
httr::GET(url = paste0("https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/", gldas_str_for_path,"/",file_str), ,
          write_disk(downloaded_file_path, overwrite = TRUE))

# Get coordinates from DSM image in degree ----
dsm = raster(paste0(working_path, 'Data/DSM.tif'))
xmin = extent(dsm)@xmin
ymin = extent(dsm)@ymin

sp1 <- SpatialPoints(matrix(c(xmin, ymin), ncol=2), proj4string = CRS(as.character(dsm@crs)) )
#transform to meters lat and lon
sp1Transformed <- spTransform(sp1, CRS("+proj=longlat +datum=WGS84"))

# get meteorological data from GLDAS
lat = extent(sp1Transformed)@ymin 
lon = extent(sp1Transformed)@xmin 
sp1deg = sp1Transformed
sp1 = sp1Transformed
coordinates(sp1deg)
nc_name = downloaded_file_path

get_value_lat_lon = function(sp, varname, nc_name){
  nc_file = nc_open(nc_name)
  varlayer = ncvar_get(nc_file, varid=varname)
  lat = ncvar_get(nc_file, varid ="lat")
  lon = ncvar_get(nc_file, varid="lon")
  lon[lon > 180] <- lon[lon > 180] - 360
  fillvalue <- ncatt_get(nc_file, varid = 0, attname =  "missing_value")
  nc_close(nc_file)
  varlayer[varlayer == fillvalue$value] <- NA
  varraster <- raster(t(varlayer), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
  varraster <- flip(varraster, direction='y')
  latlonvalue = extract (varraster, coordinates(sp))
  latlonvalue
}

ALBEDO = get_value_lat_lon(sp1deg, "Albedo_inst", nc_name)
SWR =  read.table("Example data/Input data/flight_met_station.csv", sep=",", header=T, stringsAsFactors = F)[66,]$radiation #get_value_lat_lon(sp1deg, "SWdown_f_tavg", nc_name) if no meteorological station data provided

loc <- initGRASS("/Users/user/Desktop/ofir_lab/MicroDrone/GRASS78", home=working_path, gisDbase="GRASS_R", override=TRUE ,mapset="PERMANENT")#
use_sp()

#### create slope and aspect maps ####
#set grass projection
dsm_file = paste0(working_path, "Data/DSM.tif")
r_dsm = raster(dsm_file)

execGRASS("g.proj", flags="c", parameters=list(proj4=as.character(r_dsm@crs)) )

#create slope and aspect maps
#load the dsm map to grass
execGRASS("r.in.gdal", flags=c("o", "overwrite"), parameters=list(input=dsm_file, output="dsm"))
execGRASS("g.region", parameters=list(raster="dsm") ) 

pars <- list(elevation="dsm", slope="slope", aspect="aspect")
execGRASS("r.slope.aspect", flags=c("overwrite"), parameters=pars)

#output the maps as tif files
filename_slope = paste0(working_path, "/Data/SLP.tif")
filename_aspect = paste0(working_path, "/Data/aspect.tif")
execGRASS("r.out.gdal", flags="overwrite", parameters=list(createopt = "PROFILE=GeoTIFF,TFW=YES", output=filename_slope, input="slope"))
execGRASS("r.out.gdal", flags="overwrite", parameters=list(createopt = "PROFILE=GeoTIFF,TFW=YES", output=filename_aspect, input="aspect"))

#### create vegetation index (TGI) and then vegetation (yes/no) map ####
rgb_file = paste0(working_path, "/Data/RGB.tif")
r <- stack(rgb_file)
R = r[[1]]; G = r[[2]]; B= r[[3]]
r_tgi = (G-0.39*R-0.61*B)/max(R, max(G, B))
writeRaster(r_tgi, file=paste0(working_path, "/Data/TGI.tif"), overwrite=T)
#plot(r_tgi)

r_veg = r_tgi
r_veg[r_tgi>0.04]=1
r_veg[r_tgi<0.04]=0
writeRaster(r_veg, file=paste0(working_path, "/Data/veg_yesno.tif"), overwrite=T)

#### create solar radiation and shade map ####
#get time
dt = ymd_hms("2019-11-11 13:00:00 UTC") #note that the time zone is UTC not local time
#get coordinate
sp1 <- SpatialPoints(xyFromCell(r_dsm, 1), proj4string = CRS(as.character(r_dsm@crs))) 
#convert to lat lon degrees
sp1degrees <- spTransform(sp1, CRS("+proj=longlat +datum=WGS84")) 
lat = coordinates(sp1degrees)[1,2]
lon = coordinates(sp1degrees)[1,1]

#calculate solar time needed for the r.sun grass command
solar_noon = getSunlightTimes(date = date(dt), lat = lat, lon = lon, keep = c("solarNoon"), tz = "UTC")
# solar_time = 12 - solar noon + current time
solar_time = 12 - hour(solar_noon$solarNoon) - minute(solar_noon$solarNoon)/60 + (hour(dt)+minute(dt)/60)

# Calculate radiation, assuming zero albedo (can be corrected later), also calculate incidence angle (the angle of the radiation hitting the ground)
albedo_value = ALBEDO/100 # 0.3 * length(r_veg[r_veg == 0])/ length(r_veg) + 0.2 * length(r_veg[r_veg == 1])/ length(r_veg) # 0.3 for soil. 0.2 for vegs
execGRASS("r.sun", flags="overwrite", parameters=list(elevation="dsm", day=yday(date(dt)), time=solar_time, aspect = "aspect", slope="slope",glob_rad="tmpglob", albedo_value=albedo_value, incidout="tmpincidout"))
execGRASS("r.out.gdal", flags="overwrite", parameters=list(createopt = "PROFILE=GeoTIFF,TFW=YES", output= paste0(working_path, "Data/solar.tiff"), input="tmpglob"))
execGRASS("r.out.gdal", flags="overwrite", parameters=list(createopt = "PROFILE=GeoTIFF,TFW=YES", output= paste0(working_path, "Data/incidence.tif"), input="tmpincidout"))

#create the shade map
shade = raster(paste0(working_path, "Data/incidence.tif"))
#plot(shade)
shade[!is.na(shade)] = 0
shade[is.na(shade)] = 1
#plot(shade)
writeRaster(shade, paste0(working_path, 'Data/shade.tif'), overwrite = T)

#create real_solar
solar = raster(paste0(working_path, 'Data/solar.tiff'))
#correct solar by real solar
real_solar = SWR
slope = raster(paste0(working_path, 'Data/SLP.tif'))
#get clean solar on sunlit (no shade) and flat surface (small slope)
clean_solar = max(solar[slope<0.1 & shade==0], na.rm = T)
if (clean_solar<real_solar) print ("warning: clean_solar is lower than measured solar")
r_real_solar = solar*real_solar/clean_solar
#plot(r_real_solar, add=T)
writeRaster(r_real_solar, paste0(working_path, 'Data/real_solar.tif'), overwrite = T)

# create heights map
dsm = raster(paste0(working_path, "Data/DSM.tif"))
dtm = raster(paste0(working_path, "Data/DTM.tif"))
writeRaster(dsm - dtm, paste0(working_path, 'Data/height.tif'), overwrite = T)


#create skyview map
alt = as(dsm, 'SpatialGridDataFrame')
write_RAST(alt, "dsm", flags=c("quiet", "overwrite"))
execGRASS("r.info", map="dsm") 
execGRASS("g.region", parameters=list(raster="dsm") ) 
execGRASS("r.skyview", flags="overwrite", parameters=list(input="dsm", output="tmpskyview"))
execGRASS("r.out.gdal", flags="overwrite", parameters=list(createopt = "PROFILE=GeoTIFF,TFW=YES", output=paste0(working_path, 'Data/skyview.tiff'), input="tmpskyview"))
