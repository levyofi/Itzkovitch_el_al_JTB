#create data file for physical microclimate model
#setwd("~/eclipse workspace/downscaling/for_paper")
#library(RgoogleMaps)
library(rgdal)
library(knitr)
library(leaflet)

library(lubridate)

library(plotKML)
library(stringr)

library(raster)
library(humidity)
library(bigleaf)

specific_humidity = function(RH, T){#RH - relative humidity [0-100 %]; T - air temperature [Kelvin]
  qair = SH(vapor_pressure(RH, T))
  qair
}
vapor_pressure = function(RH, T){#RH - relative humidity [0-100 %]; T - air temperature [Kelvin]
  SVP = SVP(T) #return hPa 
  SVP*RH/100*100 #return Pa
}

########### set of input data inludes:
### meteorological data from station 
# TG, radiation
### data from drone
# TGI, real_solar, shade, height

data_files = read.table("double_cuts.csv", sep=",", header=T, stringsAsFactors = F)
ifile=66
for (ifile in 1:nrow(data_files)){ 
  zone <- basename(data_files[ifile,]$zone)
  print(zone)
  dt <- as_datetime(word(zone, 2, 3, '_'), format = '%d.%m.%y_%H%M', tz=paste0("Etc/GMT-2"))
  location = word(zone, 1, 1, '_')
  dt = with_tz(dt, tzone="UTC")
  #create shade and solar files
  #get dsm raster
  dir = paste0("/home/ofir/eclipse workspace/downscaling/Double_Cuts_Complete/", zone)
  setwd(dir)
  
  raster_fix_NA = function(rs){
    # browser()
    rs_median = median(rs[], na.rm=T)
    rs[is.na(rs[])] = rs_median
    return(rs)
  }
  temp = raster_fix_NA(raster("IR.tif"))
  #r <- crop(resample(stack("RGB.tif"), temp), bbox)
  solar = raster_fix_NA(raster("real_solar.tif"))
  shade = raster_fix_NA(raster("shade.tif"))
  skyview = raster_fix_NA(raster("skyview.tiff"))
  r_tgi = raster_fix_NA(raster("TGI.tif"))
  heights = raster_fix_NA(raster("height.tif"))
  heights2 = heights
  heights2[heights<1.5 & r_tgi <0]=NA
  #plot(heights2, add=T)
  
  
  IR_emiss = 0.95 #crop(resample(raster(dir(pattern="^.+_transparent_reflectance_.+\\.tif")), solar), bbox)
  
  station_data = data_files
  station = station_data[ifile,]
  WIND = station$wind_speed #m/s
  TG = station$soil + 273.15 #this is actually the IR_Temp_Avg but Ronny changed the name to soil
  ir_emiss_ground = 0.95 #1-mean(IR_emiss[ shade==0 & r_tgi< 0], na.rm=T)
  ir_emiss_canopy = 0.95 #1-mean(IR_emiss[ shade==0 & r_tgi> 0], na.rm=T)
  pressure = data_files[ifile,]$pressure #hPa
  RH = station$rh #%
  TAIR = station$temp + 273.15 #Kelvin
  QAIR = specific_humidity(RH, TAIR)
  RHOAIR = air.density(TAIR - 273.15, pressure/10)
  canopy_temp = temp
  canopy_temp[r_tgi<0.1]=NA
  TV = cellStats(canopy_temp, stat = "mean") + 273
  TAH = (TV+TAIR)/2
  ir_emiss_sky = 1.72*( (vapor_pressure(RH, TAIR)/1000)/TAIR)^(1/7)
  
  #set cloud cover
  cleansky_solar = raster("solar.tiff")
  #plot(solar)
  #calculate clouds cover (%)
  slope = raster("SLP.tif")
  real_solar = station$radiation
  clean_solar = max(cleansky_solar[slope<0.1 & shade==0], na.rm = T)
  
  cloud_cover = max(0, 1-real_solar/clean_solar) 
  
  #set data from netcdf
  #set data from netcdf
  library(ncdf4)
  sp1 <- SpatialPoints(xyFromCell(slope, 1), proj4string = CRS(as.character(slope@crs))) 
  
  #convert to lat lon degrees
  sp1degrees <- spTransform(sp1, CRS("+proj=longlat +datum=WGS84")) 
  lat = coordinates(sp1degrees)[1,2]
  lon = coordinates(sp1degrees)[1,1]
  
  #for alona
  # lat = 32.57495 
  # lon = 34.98929 
  
  dt.num = sprintf("%04d%02d%02d%02d%02d", year(dt),month(dt), day(dt), hour(dt),minute(dt),second(dt)) #201911110810
  dt.str.for.GLDAS = sprintf("%04d%02d%02d.%02d00", year(dt),month(dt), day(dt), hour(dt)%/%3*3) #201911110810
  sp1 <- SpatialPoints(matrix(c(lon,lat), ncol=2), proj4string = CRS("+proj=longlat +datum=WGS84")) 
  coordinates(sp1)
  #change this line to look for the file
  nc_name = list.files(path="../../rasters/", pattern = paste0("^.*", dt.str.for.GLDAS,".*.nc4"),full.names = T)
  
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
    latlonvalue = extract (varraster, coordinates(sp1))
    latlonvalue
  }
  
  ALBEDO = get_value_lat_lon(sp1, "Albedo_inst", nc_name)
  soil_temp10 =  get_value_lat_lon(sp1, "SoilTMP0_10cm_inst", nc_name)
  soil_temp40 =  get_value_lat_lon(sp1, "SoilTMP10_40cm_inst", nc_name)
  soil_temp100 =  get_value_lat_lon(sp1, "SoilTMP40_100cm_inst", nc_name)
  soil_temp200 =  get_value_lat_lon(sp1, "SoilTMP100_200cm_inst", nc_name)
  soil_mois10 =  get_value_lat_lon(sp1, "SoilMoi0_10cm_inst", nc_name)
  soil_mois40 =  get_value_lat_lon(sp1, "SoilMoi10_40cm_inst", nc_name)
  soil_mois100 =  get_value_lat_lon(sp1, "SoilMoi40_100cm_inst", nc_name)
  soil_mois200 =  get_value_lat_lon(sp1, "SoilMoi100_200cm_inst", nc_name)
  #convert soil moisture from kg m-2 to m3 m-3 using the info in https://ldas.gsfc.nasa.gov/faq/nca-ldas
  #thickness of               1000 kg      1 m        layer in mm         kg
  #volumetric soil moisture  X  -------  X  -------  X  -----------   =   --
  #                               m^3       1000 mm          1            m^2
  #so basically we need to divide by the width of the layer in mm:
  soil_mois10 = soil_mois10/100
  soil_mois40 = soil_mois40/300
  soil_mois100 = soil_mois100/600
  soil_mois200 = soil_mois200/1000
  
  NSOIL=4

  meteo = data.frame(size = length(shade), NSOIL = NSOIL, width = ncol(shade), height = nrow(shade),
                   Latitude = lat, Longitude = lon, Date = date,
                   TG = TG, Albedo = ALBEDO, 
                   ST10 = soil_temp10, ST40 = soil_temp40, ST100 = soil_temp100, ST200 = soil_temp200,
                   SM10 = soil_mois10, SM40 = soil_mois40, SM100 = soil_mois100, SM200 = soil_mois200,
                   Wind = WIND, IR_em_G = ir_emiss_ground, IR_em_C = ir_emiss_canopy,
                   P = pressure, TAIR = TAIR, QAIR = QAIR, RHOAIR = RHOAIR,
                   TV = TV, TAH = TAH, IR_em_S = ir_emiss_sky, Cloud_cover = cloud_cover)

write.csv(meteo, paste0(wd, "Data/input_meteorological_data.csv"))
