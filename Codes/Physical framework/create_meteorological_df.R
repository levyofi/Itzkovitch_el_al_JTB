
library(RgoogleMaps)
library(rgdal)
library(httr)
library(knitr)
library(leaflet)
library(ncdf4)
library(lubridate)
library(stringr)
library(raster)
library(humidity)
library(bigleaf)

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

raster_fix_NA = function(rs){
  rs_median = median(rs[], na.rm=T)
  rs[is.na(rs[])] = rs_median
  return(rs)
}

vapor_pressure = function(RH, T){#RH - relative humidity [0-100 %]; T - air temperature [Kelvin]
  SVP = SVP(T) #return hPa 
  SVP*RH/100*100 #return Pa
}



wd = "/Users/user/Desktop/ofir_lab/MicroDrone/" # change to relevant directory

time = '11.11.19_1010' # change to relevant time (string)
timezone = 2 # change to relevant timezone
dt <- as_datetime(time, format = '%d.%m.%y_%H%M', tz=paste0("Etc/GMT-", timezone))

dt.num = sprintf("%04d%02d%02d%02d%02d%02d", year(dt),month(dt), day(dt), hour(dt),minute(dt),second(dt)) #201911110810
dt.str.for.GLDAS = sprintf("%04d%02d%02d.%02d00", year(dt),month(dt), day(dt), hour(dt)%/%3*3) #201911110810

date = dmy(paste0(day(dt), '/', month(dt), '/', year(dt)))
julian_day= yday(date)
time = hm(paste0(hour(dt), ':', minute(dt)))

gldas_str_for_file = ifelse(year(date)>=2022,"GLDAS_NOAH025_3H_EP.A", "GLDAS_NOAH025_3H.A")
str_file = paste0(gldas_str_for_file,year(date),str_pad(month(date), width = 2, pad = "0"),str_pad(day(date), width = 2, pad = "0"),".",str_pad(hour(date), width = 2, pad = "0"),"00.021.nc4" )
file_str=paste(year(date),str_pad(julian_day, width = 3, pad = "0"), str_file, sep="/")
netrc_path <- "Data/.netrc"
cookie_path <- "Data/.urs_cookies"
downloaded_file_path <- paste0(wd, 'Data/', str_file)

set_config(config(followlocation=1,netrc=1,netrc_file=netrc_path,cookie=cookie_path,cookiefile=cookie_path,cookiejar=cookie_path))
gldas_str_for_path = ifelse(year(date)>=2022, "GLDAS_NOAH025_3H_EP.2.1","GLDAS_NOAH025_3H.2.1")
httr::GET(url = paste0("https://hydro1.gesdisc.eosdis.nasa.gov/data/GLDAS/", gldas_str_for_path,"/",file_str), ,
          write_disk(downloaded_file_path, overwrite = TRUE))

nc_name = downloaded_file_path

shade = raster_fix_NA(raster(paste0(wd, "Data/shade.tif")))
r_tgi = raster_fix_NA(raster(paste0(wd, "Data/tgi.tif")))

sp1 <- SpatialPoints(xyFromCell(shade, 1), proj4string = CRS(as.character(shade@crs))) 
sp1degrees <- spTransform(sp1, CRS("+proj=longlat +datum=WGS84")) 
lat = coordinates(sp1degrees)[1,2]
lon = coordinates(sp1degrees)[1,1]
sp1 <- SpatialPoints(matrix(c(lon,lat), ncol=2), proj4string = CRS("+proj=longlat +datum=WGS84"))

TG = get_value_lat_lon(sp1, "AvgSurfT_inst", nc_name)
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

WIND = get_value_lat_lon(sp1, "Wind_f_inst", nc_name)
pressure = get_value_lat_lon(sp1, "Psurf_f_inst", nc_name)
TAIR = as.numeric(readline("Enter Air Temp: ")) #get_value_lat_lon(sp1, "Tair_f_inst", nc_name) # Need to be taken from GUI 
QAIR = get_value_lat_lon(sp1, "Qair_f_inst", nc_name)
RHOAIR = air.density(TAIR - 273.15, pressure/1000)

SWN = get_value_lat_lon(sp1, "Swnet_tavg", nc_name); SWD = get_value_lat_lon(sp1, "SWdown_f_tavg", nc_name)
ir_emiss_ground = 2*SWD / (2*SWD - SWN)
if (is.na(ir_emiss_ground) | is.infinite(ir_emiss_ground)){ir_emiss_ground = 0.95}

LWN = get_value_lat_lon(sp1, "Lwnet_tavg", nc_name); LWD = get_value_lat_lon(sp1, "LWdown_f_tavg", nc_name)
ir_emiss_sky = 2*LWD / (2*LWD - LWN)
if (is.na(ir_emiss_sky) | is.infinite(ir_emiss_sky)){ir_emiss_sky = 0.95}

ir_emiss_canopy = 0.95

TV = TAIR
TAH = (TV+TAIR)/2

if (SWD == 0) {SWD = -1*(LWN - 2*LWD)/ir_emiss_ground}
cleansky_solar = raster(paste0(wd, "Data/solar.tiff"))
slope = raster(paste0(wd, "Data/SLP.tif"))
real_solar = SWD
clean_solar = max(cleansky_solar[slope<0.1 & shade==0], na.rm = T)
cloud_cover = max(0, 1-real_solar/clean_solar)

meteo = data.frame(size = length(shade), NSOIL = 4, width = ncol(shade), height = nrow(shade),
                   Latitude = lat, Longitude = lon, Date = date,
                   TG = TG, Albedo = ALBEDO, 
                   ST10 = soil_temp10, ST40 = soil_temp40, ST100 = soil_temp100, ST200 = soil_temp200,
                   SM10 = soil_mois10, SM40 = soil_mois40, SM100 = soil_mois100, SM200 = soil_mois200,
                   Wind = WIND, IR_em_G = ir_emiss_ground, IR_em_C = ir_emiss_canopy,
                   P = pressure, TAIR = TAIR, QAIR = QAIR, RHOAIR = RHOAIR,
                   TV = TV, TAH = TAH, IR_em_S = ir_emiss_sky, Cloud_cover = cloud_cover)

write.csv(meteo, paste0(wd, "Data/input_meteorological_data.csv"))

