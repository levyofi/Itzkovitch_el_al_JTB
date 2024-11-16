###############################################################################
################################ PACKAGES #####################################
###############################################################################
try
    using CSV
    using DataFrames
    using ArchGDAL 
    using Statistics
    using Colors
    using Plots
    using DelimitedFiles

catch
    using Pkg
    Pkg.add("CSV")
    Pkg.add("DataFrames")
    Pkg.add("ArchGDAL")
    Pkg.add("Statistics")
    Pkg.add("Colors")
    Pkg.add("Plots")

    using CSV
    using DataFrames
    using ArchGDAL 
    using Statistics
    using Colors
    using Plots
    using DelimitedFiles
end

const AG = ArchGDAL;
###############################################################################
############################# DECLARATION #####################################
###############################################################################

# Variable          Description                                         Units
# ------------------------------------------------------------------------------
# Albedo            ground Albedo                                       dec. %
# EMG               ground emissivity                                   dec. %
# EMISS             surface emissivity                                  dec. %
# EMV               vegetation emissivity                               dec. %
# GLW               Near-IR (LW) downward flux                          W m^-2
# PSFC              surface pressure                                    Pa
# QAIR              specific humidity at surface (2m height)            Kg Kg^-1 
# RHOAIR            air density                                         Kg m^-3
# SMOIS             soil moisture                                       m^3 m^-3
# SWDOWN            Visible (SW) downward flux                          W m^-2
# T2                ground surface (2m) temprature                      K
# T2B               Bare ground surface (2m) temprature                 K
# T2V               2m temprature over canopy part                      K
# TAH               capony air temprature                               K
# TG (*)            Bulk ground temprature                              K
# TSLB (*)          soil temprature                                     K
# TV                vegetation leaf temprature                          K
# U10               u-component of wind speed (10m)                     m s^-1
# V10               v-component of wind speed (10m)                     m s^-1
# TGV (#)           ground temprature under canopy                      K
# TGB (#)           bare ground temprature                              K
# FVEG (#)          greeness vegetation fraction                        dec. %
# BGAP (#)          between canopy gap                                  dec. %
# WGAP (#)          within canopy gap                                   dec. %

### (*) - for initialization
### (#) - for output in ecological models


# ------------------------------------------------------------------------
# Physical Constansts:
# ------------------------------------------------------------------------
GRAV   = 9.80616   ; #acceleration due to gravity (m/s2)
SB     = 5.67E-08  ; # Stefan-Boltzmann constant (w/m2/k4)
VKC    = 0.40      ; # von Karman constant
TFRZ   = 273.16    ; # freezing/melting point (k)
HSUB   = 2.8440E06 ; # latent heat of sublimation (j/kg)
HVAP   = 2.5104E06 ; # latent heat of vaporization (j/kg)
HFUS   = 0.3336E06 ; # latent heat of fusion (j/kg)
CWAT   = 4.188E06  ; # specific heat capacity of water (j/m3/k)
CICE   = 2.094E06  ;  # specific heat capacity of ice (j/m3/k)
CPAIR  = 1004.64   ;   # heat capacity dry air at const pres (j/kg/k)
TKWAT  = 0.6       ;   # thermal conductivity of water (w/m/k)
TKICE  = 2.2       ;   # thermal conductivity of ice (w/m/k)
TKAIR  = 0.023     ;   # thermal conductivity of air (w/m/k)
RAIR   = 287.04    ;   # gas constant for dry air (j/kg/k)
RW     = 461.269   ;   # gas constant for  water vapor (j/kg/k)
DENH2O = 1000.0    ;   # density of water (kg/m3)
DENICE = 917.0     ;   # density of ice (kg/m3)
M      = 1.0       ;   # melting factor (-)
Z0SNO  = 0.002     ;   # snow surface roughness length (m) (0.002)
SSI    = 0.03      ;   # liquid water holding capacity for snowpack (m3/m3) (0.03)
SWEMX  = 1.00      ;  # new snow mass to fully cover old snow (mm)
CSOIL = 2.00E+6    ;   # vol. soil heat capacity (j/m3/K)

# ---------------------------------------------------------
# Noah-mp configuration variables:
# ---------------------------------------------------------
OPT_FRZ  = 1.0  ;
OPT_STC = 1.0   ;
OPT_SFC = 1.0   ;
OPT_TBOT = 1.0  ;
wrfNSNOW = 3.0  ; # maximum no. of snow layers [=3]
wrfNSOIL = 4.0  ; # No. of soil layers [=4]
IST = 1.0       ; # surface type: 1->soil; 2->lake [=1]
ZLVL = 2.0      ; # reference height for 2m temperature (m)
dt = 300.0      ; # land model time step (sec) [=3600 seconds]

# -----------------------------------------------------------
# location specific parameters: 
# -----------------------------------------------------------
## B parameter
EXPB   = [2.79,   4.26,   4.74,   5.33,   5.33,   5.25,   6.66,   8.72,   8.17,   10.73,   10.39,   11.55,   5.25,   0.0,   2.79,   4.26,   11.55,   2.79,   2.79] ; #soil specific b parameter
## porosity, saturated value of soil moisture (volumetric)
MAXSMC = [0.339,   0.421,   0.434,   0.476,   0.476,   0.439,   0.404,   0.464,   0.465,   0.406,   0.468,   0.468,   0.439,   1.0,   0.20,   0.421,   0.468,   0.200,   0.339] ; 
## saturated soil matric potential
SATPSI = [0.069,   0.036,   0.141,   0.759,   0.759,   0.355,   0.135,   0.617,   0.263,   0.098,   0.324,   0.468,   0.355,   0.0,   0.069,   0.036,   0.468,   0.069,   0.069] ;
## wilting point soil moisture (volumetric)
WLTSMC = [0.010,   0.028,   0.047,   0.084,   0.084,   0.066,   0.067,   0.120,   0.103,   0.100,   0.126,   0.138,   0.066,   0.0,   0.006,   0.028,   0.030,   0.006,   0.01] ;
## soil quartz content
QTZ = [0.92,   0.82,   0.60,   0.25,   0.10,   0.40,   0.60,   0.10,   0.35,   0.52,   0.10,   0.25,   0.05,   0.60,   0.07,   0.25,   0.60,   0.52,   0.92] ;
## momentum roughness length (m)
Z0MVT = [1.00,  0.06,  0.06,  0.06,  0.06,  0.15,  0.06,  0.06,  0.06,  0.86,  0.80,  0.85,  1.10,  1.09,  0.80,  0.0001,  0.06,  0.05, 0.001,  0.04,  0.06,  0.06,  0.03,  0.001,  0.01,  0.00,  0.00] ;

saved_layers = [1,2,3,4,5,6,7,8,9,10,16, 22, 28, 34, 40, 46, 52, 58, 66] ;

#########################################################
################## PreProccessing #######################
#########################################################

function get_raster_data(name)
    global working_directory
    raster = 1
    try
        raster = AG.read(working_directory * "input_data/" * name * ".tif"); 
    catch
        raster = AG.read(working_directory * "input_data/" * name * ".tiff");
    end
    band = AG.getband(raster, 1);
    map = AG.read(band);
    return map
end
function fix_raster_NA(data)
    d = replace(data, NaN => missing);
    # mean = sum(skipmissing(d)) / (length(data) - sum(isnan.(data)));
    median = Statistics.median(skipmissing(d))
    replace(data, NaN => median);
end
get_map_data(name) = fix_raster_NA(get_raster_data(name)) ;

# ARGS = [input, output, soiltype (default = 0)]
ARGS = ["../Example data/Input data", 
    "../Example data/Output data/", 
    8] ;

args = map(x -> string(x), ARGS) ; 
working_directory = String(args[1]) ;
meteoro = CSV.read(working_directory * "input_data/input_meteorological_data.csv", DataFrame) ;

timeValues_num = meteoro.Date[1] ;
TG = meteoro.TG[1] ;
ALBEDO = meteoro.Albedo[1] ;
TSLB = [meteoro.ST10[1], meteoro.ST40[1], meteoro.ST100[1], meteoro.ST200[1]] ; 
SMOIS = [meteoro.SM10[1], meteoro.SM40[1], meteoro.SM100[1], meteoro.SM200[1]] ; 
WIND = meteoro.Wind[1] ;
EMG = meteoro.IR_em_G[1] ;
EMV = meteoro.IR_em_C[1] ;
PSFC = meteoro.P[1] ;
TAIR = meteoro.TAIR[1] ;
QAIR = meteoro.QAIR[1] ;
RHOAIR = meteoro.RHOAIR[1] ;
TV = meteoro.TV[1] ;
TAH = meteoro.TAH[1] ;
SKYEMISS = meteoro.IR_em_S[1] ;
CLD  = meteoro.Cloud_cover[1] ;

HVEG = get_map_data("height") ;
SWDOWN = get_map_data("solar") ;
SHADE = get_map_data("shade") ;
SKYVIEW = get_map_data("skyview") ;
TGI = get_map_data("tgi") ;

soiltype = parse(Int, (args[3])) ; 
if soiltype == 0
    BEXP = Statistics.mean(EXPB) ;
    SMCMAX = Statistics.mean(MAXSMC) ;
    PSISAT = Statistics.mean(SATPSI) ;
    SMCWLT = Statistics.mean(WLTSMC) ;
    QUARTZ = Statistics.mean(QTZ) ;
else
    BEXP = EXPB[soiltype] ;
    SMCMAX = MAXSMC[soiltype] ;
    PSISAT = SATPSI[soiltype] ;
    SMCWLT = WLTSMC[soiltype] ;
    QUARTZ = QTZ[soiltype] ;
end

########################################################################
########################### MODELING FUNCTIONS #########################
########################################################################

function surface_temp(dz, dair, perm)
    # previously called surfacet. 
    # compute surface temperature (?)
    #dz -> soil layer thickness (m)
    #dair -> air layer thickness (m)  
    
    # global all - how to?
    n = meteoro.width[1] ; # number of coords - one row each permutation

    # create depth array
    NSOIL = UInt8(floor(2/dz)) ;
    NAIR = UInt8(floor(2/dair)) ;

    z = fill(0.0, NSOIL); # soil depth
    for i in 1:NSOIL
        z[i] *= i ;
    end

    # assume no snow
    # snow parameters
    NSNOW = 0 ; 
    ISNOW = 0 ;
    SNOWH = 0 ;
    SNEQV = 0 ;
    SNICE = 0 ;
    SNLIQ = 0 ;
    
    # create matrices in dimentions time step * depth 
    # when data will be allocated to this matrices, it will need to be NSNOW + i, in order to handle the snow layers.
    # in our model we assume no snow, so this is not affect anything.
    tsoil = fill(0.0, (n,NSNOW + NSOIL)) ; # matrix of soil temperature
    tsurf = fill(0.0, (n)) ; # calculated surface temperature
    TAIRS = fill(0.0, (n,NAIR)) ; # calculated air temperature at various heights
    UVH = fill(0.0, (n,NAIR)) ; # calculated wind speed at various heights
    xDF = fill(0.0, (n,NSNOW + NSOIL)) ; # matrix of soil thermal conductivity (w/m/kelvin)
    xHCPCT = fill(0.0, (n,NSNOW + NSOIL)) ; # matrix of soil heat capacity (j/m**3/kelvin)
    xSMOIS = fill(0.0, (n,NSOIL)) ; # matrix of total soil water (m**3/m**3)
    SH2O = fill(0.0, (1,NSOIL)) ; # vector of soil liquid water (m**3/m**3)
    xZSNSO = fill(0.0, (NSNOW + NSOIL)) ; # vector of snow/soil layer-bottom depth from snow/ground surface [m]
    DZSNSO = fill(0.0, (NSNOW + NSOIL)) ; # snow/soil layer thickness [m]
    SNICE = fill(0.0, (NSNOW + NSOIL)) ; 
    SNLIQ = fill(0.0, (NSNOW + NSOIL)) ;  
    xSAG = fill(0.0, (n)) ; 
    xSH = fill(0.0, (n)) ;  
    xEV = fill(0.0, (n)) ; 
    xIR = fill(0.0, (n)) ; 
    xGH = fill(0.0, (n)) ; 

    for i in 1:NSOIL
        depth = i*dz
        if depth < 0.1
            level = 1
        elseif depth < 0.4
            level = 2
        elseif depth < 1
            level = 3
        else
            level = 4
        end
        xSMOIS[:,i] .= SMOIS[level] ;
        tsoil[:,i+NSNOW] .= TSLB[level] ;
    end
    xSMOIS[:,1] = xSMOIS[:,1]/3.0 ;

    # snow/soil layer thickness
    for iz in ISNOW+1 : NSOIL
        if iz <= 0
            xZSNSO[iz] = ZSNSO[iz] # NTBA - ZSNOW dont defined!!
        elseif iz == ISNOW + 1
            xZSNSO[iz] = -dz
        else
            xZSNSO[iz] = xZSNSO[iz-1] - dz
        end
    end
    for iz in ISNOW+1 : NSOIL
        if iz == ISNOW+1
            DZSNSO[iz] = -xZSNSO[iz]
        else
            DZSNSO[iz] = xZSNSO[iz-1] - xZSNSO[iz]
        end
    end

    # start calculation for each coordinate
    for k in 1:n
        SH2O = xSMOIS[k,:] ;
        TAIRS[k,:] .= TAIR ;
        tsurf[k] = TG ;
        cur_shade = SHADE[k+n*(perm-1)] ;

        t = 0 ;
        while t < 600.0 # why twice?
            SNEQV, SH2O, xSMOIS[k,:], SNICE, SNLIQ, tsurf[k], tsoil[k,:], xHCPCT[k,:], xDF[k,:], 
            TAIRS[k,:], UVH[k,:], xSAG[k], xSH[k], xEV[k], xIR[k], xGH[k] = ENERGY(cur_shade, dair, SKYVIEW[k+n*(perm-1)], HVEG[k+n*(perm-1)],
            TGI[k+n*(perm-1)], NAIR, NSOIL, NSNOW, ISNOW, RHOAIR, PSFC, QAIR, TAIR, WIND, SWDOWN[k+n*(perm-1)], xZSNSO,
            CSOIL, DZSNSO, EMG, EMV, SKYEMISS, SNOWH, TAH,  TV, SNEQV, SH2O, xSMOIS[k,:], SNICE, SNLIQ,
            ALBEDO, tsurf[k], tsoil[k,:], xHCPCT[k,:], xDF[k,:]) ;
            t += dt;           
        end
    end

    tsurf_all = tsurf ; # need to be changed for parallelization

    return tsurf_all
end

function ESAT(T)
    # use polynomials to calculate saturation vapor pressure and derivative with
    # respect to temperature: over water when t > 0 c and over ice when t <= 0 c
    A0=6.107799961 ; A1=4.436518521E-01 ; A2=1.428945805E-02 ; A3=2.650648471E-04 ; 
    A4=3.031240396E-06 ; A5=2.034080948E-08 ; A6=6.136820929E-11 ;

    B0=6.109177956 ; B1=5.034698970E-01 ; B2=1.886013408E-02 ; B3=4.176223716E-04 ;
    B4=5.824720280E-06 ; B5=4.838803174E-08 ; B6=1.838826904E-10 ;

    C0= 4.438099984E-01 ; C1=2.857002636E-02 ; C2= 7.938054040E-04 ; C3=1.215215065E-05 ;
    C4= 1.036561403E-07 ; C5=3.532421810e-10 ; C6=-7.090244804E-13 ;

    D0=5.030305237E-01 ; D1=3.773255020E-02 ; D2=1.267995369E-03 ; D3=2.477563108E-05 ; 
    D4=3.005693132E-07 ; D5=2.158542548E-09 ; D6=7.131097725E-12 ;
    
    ESW  = 100.0*(A0+T*(A1+T*(A2+T*(A3+T*(A4+T*(A5+T*A6)))))) ;
    ESI  = 100.0*(B0+T*(B1+T*(B2+T*(B3+T*(B4+T*(B5+T*B6)))))) ;
    DESW = 100.0*(C0+T*(C1+T*(C2+T*(C3+T*(C4+T*(C5+T*C6)))))) ;
    DESI = 100.0*(D0+T*(D1+T*(D2+T*(D3+T*(D4+T*(D5+T*D6)))))) ;
    return ESW, ESI, DESW, DESI
end

function TDFCND(SMC, SH2O)
    # Calculate thermal diffusivity and conductivity of the soil.
    # Peters-Lidard approach (Peters-Lidard et al., 1998)
    SATRATIO = SMC / SMCMAX  ; # SATURATION RATIO
    THKW = 0.57 ; # water thermal conductivity
    THKO = 2.0 ; # thermal conductivity for other soil components
    THKQTZ = 7.7 ; # thermal conductivity for quartz
    THKS = (THKQTZ ^ QUARTZ)* (THKO ^ (1.0 - QUARTZ)) ; # thermal conductivity for the solids
    XUNFROZ = SH2O / SMC ; # UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
    XU = XUNFROZ * SMCMAX ; # SATURATED THERMAL CONDUCTIVITY
    THKSAT = THKS ^ (1.0 - SMCMAX) * TKICE ^ (SMCMAX - XU) * THKW ^ (XU) ; # DRY DENSITY IN KG/M3
    GAMMD = (1.0 - SMCMAX)*2700.0 ; 
    THKDRY = (0.135* GAMMD+ 64.7)/ (2700.0 - 0.947* GAMMD); # DRY THERMAL CONDUCTIVITY IN W.M-1.K-1

    if (SH2O + 0.0005) < SMC
        AKE = SATRATIO # RANGE OF VALIDITY FOR THE KERSTEN NUMBER (AKE)
    else
        if SATRATIO > 0.1
            AKE = log10(SATRATIO) + 1.0
        else
            AKE = 0.0
        end
    end
    DF = AKE * (THKSAT - THKDRY) + THKDRY # thermal conductivity [w/m/k]
    return DF
end

function  Ground_Flux(SHADE, TGI, HVEG, SKYVIEW, ISNOW, NAIR, DZSNSO, SWDOWN, ALBEDO, LWDN, UR, 
    QAIR, RHOAIR, SNOWH, Z0M, dair, EMG, EMV, EMISS, STC, DF, RSURF,
    LATHEA, GAMMA, RHSUR, T2M , TAH, TV, TG)

    # use newton-raphson iteration to solve ground temperature (TG) that balances the surface energy budgets.
    
    function TDC(T) # Kalvin to Celsius with limit -50 to +50
        T  = T - TFRZ ;
        if T < -50 || T > 50
            return T # NTBA should be empty to limit the result
        end
        return T
    end
    MPE = 1E-6 ;
    DTG = 0.0 ;
    MOZSGN = 0.0 ;
    MOZOLD = 0.0 ;
    H      = 0.0 ;
    FV     = 0.1 ;
    T2 = T2M ;

    tairs = fill(0.0, NAIR) ;
    UVH = fill(0.0, NAIR) ;

    SAG=(1.0-ALBEDO/100.0)*SWDOWN ; # Solar Radiation Sn

    if TGI<0.01 || HVEG<1.0
        L_g_sky = EMG*SKYVIEW*LWDN ;
        L_g_canopy = (1-SKYVIEW)*(EMG*(1.0-EMV)*LWDN + EMG*SB*TV^4.0) ;
    else
        L_g_sky  = 0 ;
        L_g_canopy = (EMG*(1.0-EMV)*LWDN + EMG*SB*TV^4.0) ; #assuming sky view of 20% under the green vegetation
        SAG = 0 ;
    end

    L_c =  -L_g_sky- L_g_canopy ;
    CIR = EMG*SB ;

    # coefficient for ground heat flux (Eq. 10 in metadata)
    if ISNOW < 0
        CGH = 2.0*DF[ISNOW+1]/DZSNSO[ISNOW+1] ;
    else # no snow layer present  (combine the first two soil layers for stability)
        CGH = 2.0*DF[1]/0.1 ;
    end
    
    MOZ = 0 ;
    MOZSGN = 0 ;
    FV = 0 ;
    FM = 0 ;
    FH = 0 ;
    Z0H = Z0M ;
    for ITER in 1:10 # previously called loop3
        Z0H = Z0M ;

        if OPT_SFC == 1
            MOZ, MOZSGN, FM, FH, CM, CH, FV = SFCDIF1(ITER, T2, RHOAIR, H, QAIR, ZLVL, Z0M, Z0H, UR, 
            MPE, MOZ, MOZSGN, FV, FM, FH) ;
        end

        RAMB = max(1.0,1.0/(CM*UR)) ;
        RAHB = max(1.0,1.0/(CH*UR)) ;
        RAWB = RAHB ;
        EHB = 1.0/RAHB ;

        T = TDC(TG);
        ESATW, ESATI, DSATW, DSATI = ESAT(T) ;
        if T > 0
            ESTG  = ESATW ;
            DESTG = DSATW ;
        else 
            ESTG  = ESATI ;
            DESTG = DSATI ;
        end

        CSH = RHOAIR*CPAIR/RAHB ; # sensible heat coefficient
        CEV = RHOAIR*CPAIR/GAMMA/(RSURF+RAWB) ; # evaporation heat flux coefficient
        IR   = CIR * TG^4 + L_c ; # LW radiation flux

        if TGI < 0
            SH = CSH * (TG - T2M) ; # SW heat flux
        else
            SH = CSH * (TG - TAH) ;
        end

        EAIR = QAIR*PSFC / (0.622+0.378*QAIR) ;
        EAH = EAIR ; # no significant evaporation from leafs
        EV = CEV * (ESTG*RHSUR - EAIR*(1.0-SHADE) - EAH*SHADE) ; # latent radiation flux
        GH = CGH * (TG - STC[ISNOW+1]) ; # ground heat flux
        
        # newton-raphson parameters
        B = SAG-IR-SH-EV-GH ;
        A = 4.0*CIR*TG^3 + CSH + CEV*DESTG + CGH ;
        DTG   = B/A ;

        # update ground temperature and fluxes
        TG = TG + DTG ; 
        IR = IR + 4.0*CIR*TG^3*DTG ;
        SH = SH + CSH*DTG ;
        EV = EV + CEV*DESTG*DTG ;
        GH = GH + CGH*DTG ;

        T = TDC(TG);
        ESATW, ESATI, DSATW, DSATI = ESAT(T) ;
        if T > 0
            ESTG  = ESATW ;
        else 
            ESTG  = ESATI ;
        end

        # stability calculation
        if OPT_SFC == 1 || OPT_SFC ==2
            EHB2  = FV*VKC/log((2.0+Z0H)/Z0H) ;
            if EHB2 < 1e-5 
                T2  = TG ;
            else
                T2  = TG - SH/(RHOAIR*CPAIR) * 1.0/EHB2 ;
            end
        end

        CH = EHB ;
        # for M-O length
        H = CSH * (TG - T2)
    end # close loop3

    # if snow on ground and TG > TFRZ: reset TG = TFRZ. reevaluate ground fluxes.
    if OPT_STC == 1
        if SNOWH > 0.05 && TG > TFRZ
            TG = TFRZ ;
            IR = CIR * TG^4 - EMG*LWDN ;
            if TGI<0.00
                SH = CSH * (TG - T2M) ;
            else
                SH = CSH * (TG - TAH) ;
            end
            EV = CEV * (ESTG*RHSUR - EAIR) ;
            GH = SAG - (IR+SH+EV) ;
        end
    end

    # calculate air temperatures and wind speeds at various heights for output
    height = dair
    iheight=1
    while height<2.0
        tairs[iheight] = (log(height/Z0H+1)/log(2.0/Z0H+1))*(T2M-TG) + TG ; # NTBA tairs not defined?
        UVH[iheight] = max(FV/VKC*log(height/Z0H+1.),0.1) ;
        height = height + dair ;
        iheight = iheight + 1 ;
    end 

    if ISNOW<0
        GH = CGH * (TG - STC[ISNOW+1]) ;
    else
        GH = CGH * (TG - (STC[ISNOW+1]+STC[ISNOW+2])/2.0) ;
    end 

    xSAG = SAG ;
    xSH = SH ;
    xIR = IR ;
    xEV = EV ;
    xGH = GH ;

    return TG, GH, tairs, UVH, xSAG, xSH, xEV, xIR, xGH
end
    
function ENERGY(SHADE, dair, SKYVIEW, HVEG, TGI, NAIR, NSOIL, NSNOW, ISNOW, RHOAIR, SFCPRS, QAIR, TAIR, WIND, SWDOWN, ZSNSO, 
    CSOIL, DZSNSO, EMG, EMV, EMISS, SNOWH, TAH, TV, SNEQV, SH2O, SMC, SNICE, SNLIQ, ALBEDO, TG, STC, HCPCT, DF)

    UR = max(0.01, (1-SHADE)*WIND) ; # wind speed at reference height
    T2 = TAIR ;
    Z0 = 0.01 ; #Bare-soil roughness length (m) (i.e., under the canopy)

    FSNO = 0.0 # ground snow cover fraction
    if SNOWH > 0
        BDSNO = SNEQV / SNOWH ;
        FMELT = (BDSNO/100.0)^M ;
        FSNO = tanh( SNOWH /(2.5* Z0 * FMELT)) ;
    end

    # ground roughness length
    if IST == 2
        if TG <= TFRZ
            Z0MG = 0.01 * (1.0-FSNO) + FSNO * Z0SNO ;
        else
            Z0MG = 0.1 ;
        end
    else
        Z0MG = Z0 * (1.0-FSNO) + FSNO * Z0SNO ;
    end

    #roughness length and displacement height
    if TGI>=0.01
        Z0M  = max(Z0MG, 0.02*HVEG) ;
        UR = 0.1*UR ;
    else
        Z0M  = Z0MG ;
    end

    DF, HCPCT, FACT = Thermo_Prop(NSOIL, NSNOW, ISNOW,  DZSNSO, SNOWH, SNICE, SNLIQ, CSOIL, SMC, SH2O, TG, STC, UR, Z0M) ;

    # soil surface resistance for ground evap.
    # RSURF based on Sakaguchi and Zeng, 2009
    # taking the "residual water content" to be the wilting point,           ! and correcting the exponent on the D term (typo in SZ09 ?)
    L_RSURF = (0.1) * ( exp((1.0 - min(1.0,SH2O[1]/SMCMAX)) ^ 5 ) - 1.0 ) / ( 2.71828 - 1.0 ) ;
    D_RSURF = 2.2E-5 * SMCMAX * SMCMAX * ( 1.0 - SMCWLT / SMCMAX ) ^ (2.0+3.0/BEXP) ;
    RSURF = L_RSURF / D_RSURF ;
    if SH2O[1] < 0.01 && SNOWH == 0.0
        RSURF = 1e6 ;
    end

    # calculate relative humidity in soil/snow air space
    PSI   = -PSISAT*(max(0.01,SH2O[1])/SMCMAX)^(-BEXP) ;
    RHSUR = FSNO + (1.0-FSNO) * exp(PSI*GRAV/(RW*TG)) ;

    #calculate the psychrometric constant (Eq. 7 in metadata)
    if T2 > TFRZ
        LATHEA = HVAP ;
    else
        LATHEA = HSUB ;
    end
    GAMMA = CPAIR * SFCPRS / (0.622*LATHEA) ;
    LWDN = SKYEMISS * SB * (1-CLD)*(TAIR)^4 + SB*CLD*(TAIR-2)^4 ;

    # Surface temperatures of the ground
    TG, GH, TAIRS, UVH, xSAG, xSH, xEV, xIR, xGH = Ground_Flux(SHADE, TGI, HVEG, SKYVIEW, ISNOW, NAIR, 
    DZSNSO, SWDOWN, ALBEDO, LWDN, UR, QAIR, RHOAIR, SNOWH, Z0M,
    dair, EMG, EMV, EMISS, STC, DF, RSURF, LATHEA, GAMMA, RHSUR, T2, TAH, TV, TG) ; 

    # 3L snow & 4L soil temperatures
    STC = TSNOSOI(NSOIL, NSNOW, ISNOW, ZSNSO, DZSNSO, GH, DF, HCPCT, STC) ;

    # Energy released or consumed by snow & frozen soil
    IMELT, QMELT, PONDING, STC, SH2O = Phase_Change(NSOIL, NSNOW, ISNOW, FACT, DZSNSO, HCPCT, 
    SNEQV, SNOWH, STC, SH2O, SMC, SNICE, SNLIQ ) ;
    return SNEQV, SH2O, SMC, SNICE, SNLIQ, TG, STC, HCPCT, DF, TAIRS, UVH, xSAG, xSH, xEV, xIR, xGH
end

function TSNOSOI(NSOIL, NSNOW, ISNOW, ZSNSO, DZSNSO, SSOIL, DF, HCPCT, STC)
    # Compute snow (up to 3L) and soil (4L) temperature. Note that snow temperatures
    # during melting season may exceed melting point (TFRZ) but later in Phase_Change
    # subroutine the snow temperatures are reset to TFRZ for melting snow.

    # snow/soil heat storage for energy balance check
    TBEG = STC ;

    # compute soil temperatures
    AI, BI, CI, RHSTS, EFLXB = HRT(NSOIL, ISNOW, ZSNSO, STC, DF, HCPCT, SSOIL) ;
    AI, BI, CI, RHSTS, STC = HSTEP(NSOIL, NSNOW, ISNOW, AI, BI, CI, RHSTS, STC) ;
    if OPT_TBOT == 1
        EFLXB2 = 0.0 ;
    end

    # energy balance check
    ERR_EST = 0.0 ;
    for iz in ISNOW+1:NSOIL
        ERR_EST += (STC[iz]-TBEG[iz]) * DZSNSO[iz] * HCPCT[iz] / dt ;
    end

    if OPT_STC == 1 # semi-implicit
        ERR_EST -= - (SSOIL +EFLXB) ;
    end

    if abs(ERR_EST) > 1.0
        #println("TSNOSOI is losing(-)/gaining(+) false energy ", ERR_EST, " W/m2")
    end

    return STC
end

function HRT(NSOIL, ISNOW, ZSNSO, STC, DF, HCPCT, SSOIL)

    DENOM = fill(0.0, (ISNOW+NSOIL)) ;
    DDZ = fill(0.0, (ISNOW+NSOIL)) ;
    DTSDZ = fill(0.0, (ISNOW+NSOIL)) ;
    EFLUX = fill(0.0, (ISNOW+NSOIL)) ;

    # left-hand coefficient
    AI = fill(0.0, (ISNOW+NSOIL)) ;
    BI = fill(0.0, (ISNOW+NSOIL)) ;
    CI = fill(0.0, (ISNOW+NSOIL)) ;

    # right-hand matrix
    RHSTS = fill(0.0, (ISNOW+NSOIL)) ;
    BOTFLX = 0.0
    for k in ISNOW+1:NSOIL
        if k == ISNOW+1
            DENOM[k]  = - ZSNSO[k] * HCPCT[k] ;
            TEMP1     = - ZSNSO[k+1] ;
            DDZ[k]    = 2.0 / TEMP1 ;
            DTSDZ[k]  = 2.0 * (STC[k] - STC[k+1]) / TEMP1 ;
            EFLUX[k]  = DF[k] * DTSDZ[k] - SSOIL ;
        elseif k < NSOIL
            DENOM[k] = (ZSNSO[k-1] - ZSNSO[k]) * HCPCT[k] ;
            TEMP1 = ZSNSO[k-1] - ZSNSO[k+1] ;
            DDZ[k] = 2.0 / TEMP1 ;
            DTSDZ[k] = 2.0 * (STC[k] - STC[k+1]) / TEMP1 ;
            EFLUX[k] = (DF[k]*DTSDZ[k] - DF[k-1]*DTSDZ[k-1]) ;
        elseif k == NSOIL
            DENOM[k] = (ZSNSO[k-1] - ZSNSO[k]) * HCPCT[k] ;
            TEMP1 = ZSNSO[k-1] - ZSNSO[k] ;
            BOTFLX = 0.0 ;
            EFLUX[k] = (-BOTFLX - DF[k-1]*DTSDZ[k-1]) ;
        end
    end
    
    for k in ISNOW+1:NSOIL
        if k == ISNOW+1
            AI[k] = 0 ;
            CI[k] = -DF[k] * DDZ[k] / DENOM[k] ;
            if OPT_STC == 1
                BI[k] = -CI[k]
            elseif OPT_STC == 2
                BI[k] = -CI[k] * DF[k] / (0.5*ZSNSO[k]*ZSNSO[k]*HCPCT[k]) ; 
            end
        elseif k < NSOIL
            AI[k] = -DF[k-1] * DDZ[k-1] / DENOM[k] ;
            CI[k] = -DF[k] * DDZ[k] / DENOM[k] ;
            BI[k] = - (AI[k] + CI[k]) ;
        elseif k == NSOIL
            AI[k] = -DF[k-1] * DDZ[k-1] / DENOM[k] ;
            CI[k] = 0.0 ;
            BI[k] = - (AI[k] + CI[k]) ;
        end

        RHSTS[k] = EFLUX[k]/ -DENOM[k] ;
    end

    return AI, BI, CI, RHSTS, BOTFLX
end

function HSTEP(NSOIL, NSNOW, ISNOW, AI, BI, CI, RHSTS, STC)
    # calculate/update the soil temperature field

    RHSTS *= dt ;
    AI *= dt ;
    BI *= dt ;
    CI *= dt ;

    RHSTSIN = RHSTS ;
    CIIN = CI ;

    # solve the tri-diagonal matrix equation
    CIIN, CI, RHSTS = ROSR12(CI, AI, BI, CIIN, RHSTSIN, RHSTS, floor(ISNOW+1), NSOIL, NSNOW) ;

    # update snow and soil temperature
    STC += CI ;  

    return AI, BI, CI, RHSTS, STC
end

function ROSR12(P,A,B,C,D,DELTA,NTOP,NSOIL,NSNOW)
    # ----------------------------------------------------------------------
    # INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
    # ###                                            ### ###  ###   ###  ###
    # #B(1), C(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
    # #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
    # # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
    # # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
    # # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
    # # .                                          .   # #  .   # = #   .  #
    # # .                                          .   # #  .   #   #   .  #
    # # .                                          .   # #  .   #   #   .  #
    # # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
    # # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
    # # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
    # ###                                            ### ###  ###   ###  ###
    # ----------------------------------------------------------------------
    
    C = fill(0.0, NSNOW + NSOIL) ;
    P = fill(0.0, NSNOW + NSOIL) ;
    DELTA = fill(0.0, NSNOW + NSOIL) ;

    C[NSOIL] = 0 ;
    P[NTOP] = -C[NTOP] / B[NTOP] ;

    # solve coefficient for first layer
    DELTA[NTOP] = D[NTOP] / B[NTOP] ;

    # solve coefficients for layer 2 to NSOIL
    for k in NTOP+1:NSOIL
        P[k] = -C[k] * (1.0 / (B[k] + A[k] * P[k-1])) ;
        DELTA[k] = (D[k] - A[k]* DELTA[k-1]) * (1.0 / (B[k] + A[k] * P[k-1])) ;
    end

    # set P to DELTA for lowest soil layer
    P[NSOIL] = DELTA[NSOIL] ;

    # adjust P for layer 2 to NSOIL
    for k in NTOP+1:NSOIL
        kk = NSOIL - k + (NTOP - 1) + 1 ;
        P[kk] = P[kk] * P[kk+1] + DELTA[kk]
    end

    return C, P, DELTA
end

function Thermo_Prop(NSOIL, NSNOW, ISNOW, DZSNSO, SNOWH, SNICE, SNLIQ, CSOIL,
    SMC, SH2O, TG, STC, UR, Z0M)

    TKSNO, CVSNO, SNICEV, SNLIQV, EPORE = CSNOW(NSOIL, NSNOW, ISNOW, DZSNSO, SNICE , SNLIQ) ;

    SICE = fill(0.0, NSOIL) ;
    DF = fill(0.0, (NSNOW+NSOIL)) ;
    HCPCT = fill(0.0, (NSNOW+NSOIL)) ;
    FACT = fill(0.0, (NSNOW+NSOIL)) ;

    for iz in 1:ISNOW
        DF[iz] = TKSNO[iz] ;
        HCPCT[iz] = CVSNO[iz] ;
    end 

    for iz in 1:NSOIL
        SICE[iz] = SMC[iz] - SH2O[iz] ;
        HCPCT[iz] = SH2O[iz] * CWAT + (1.0-SMCMAX)*CSOIL + (SMCMAX - SMC[iz]) * CPAIR + SICE[iz] * CICE ;
        DF[iz] = TDFCND(SMC[iz], SH2O[iz]) ;
    end

    for iz in ISNOW+1:NSOIL
        FACT[iz] = dt / (HCPCT[iz] * DZSNSO[iz]) ;
    end

    if ISNOW == 0
        DF[1] = (DF[1]*DZSNSO[1] + 0.35 * SNOWH) / (SNOWH + DZSNSO[1]) ;
    else
        DF[1] = (DF[1]*DZSNSO[1] + DF[0]*DZSNSO[0]) / (DZSNSO[0] + DZSNSO[1]) ;
    end

    return DF, HCPCT, FACT
end

function CSNOW(NSOIL, NSNOW, ISNOW, DZSNSO, SNICE , SNLIQ)

    CVSNO = fill(0.0, NSNOW+1) ; #volumetric specific heat (j/m3/k)
    TKSNO = fill(0.0, NSNOW+1) ; #thermal conductivity (w/m/k)
    SNICEV = fill(0.0, NSNOW+1) ; #partial volume of ice [m3/m3]
    SNLIQV = fill(0.0, NSNOW+1) ; #partial volume of liquid water [m3/m3]
    EPORE = fill(0.0, NSNOW+1) ; #effective porosity [m3/m3]
    BDSNOI = fill(0.0, NSNOW+1) ; #bulk density of snow(kg/m3)

    #thermal capacity of snow
    for iz in ISNOW+1:NSNOW+1
        SNICEV[iz] = min(1., SNICE[iz]/(DZSNSO[iz]*DENICE) ) ;
        EPORE[iz] = 1.0 - SNICEV[iz] ;
        SNLIQV[iz] = min(EPORE[iz],SNLIQ[iz]/(DZSNSO[iz]*DENH2O)) ;

        BDSNOI[iz] = (SNICE[iz]+SNLIQ[iz])/DZSNSO[iz] ;
        CVSNO[iz] = CICE*SNICEV[iz]+CWAT*SNLIQV[iz] ;
    end

    # thermal conductivity of snow
    for iz in ISNOW+1:NSNOW+1
        TKSNO[iz] = 3.2217e-6*BDSNOI[iz]^2.0  ;
    end
    
    return TKSNO, CVSNO, SNICEV, SNLIQV, EPORE
end

function Phase_Change(NSOIL, NSNOW, ISNOW, FACT, DZSNSO, HCPCT, SNEQV, SNOWH, STC, SH2O, SMC, SNICE, SNLIQ )
    # melting/freezing of snow water and soil water

    QMELT   = 0.0 ; # snowmelt rate [mm/s]
    PONDING = 0.0 ; # snowmelt when snow has no layer [mm]
    XMF     = 0.0 ; # total latent heat of phase change
    SUPERCOOL = fill(0.0, NSOIL+NSNOW) ; # supercooled water in soil (kg/m2)
    MICE = fill(0.0, NSOIL+NSNOW) ; # soil/snow ice mass [mm]
    MLIQ = fill(0.0, NSOIL+NSNOW) ; # soil/snow liquid water mass [mm]

    for j in 1: ISNOW+NSOIL
        if j < ISNOW
            MICE[j] = SNICE[j] ;
            MLIQ[j] = SNLIQ[j] ;
        else
            SH2O[j] = min(SH2O[j], SMC[j]) ;
            MLIQ[j] = SH2O[j] * DZSNSO[j] * 1000.0 ;
            MICE[j] = (SMC[j] - SH2O[j])  * DZSNSO[j] * 1000.0 ;
        end
    end

    IMELT = fill(0.0, (NSNOW+NSOIL)) ; # phase change index
    HM = fill(0.0, (NSNOW+NSOIL)) ; # energy residual [w/m2]
    XM = fill(0.0, (NSNOW+NSOIL)) ; # melting or freezing water [kg/m2]
    WICE0 = MICE ;
    WLIQ0 =MLIQ ;
    WMASS0 = MICE + MLIQ ;

    if IST == 1
        for j in NSNOW+1:NSOIL
            if OPT_FRZ == 1
                if STC[j] < TFRZ
                    SMP = HFUS*(TFRZ-STC[j])/(GRAV*STC[j]) ; # [m]
                    SUPERCOOL[j] = SMCMAX*(SMP/PSISAT)^(-1.0/BEXP) ;
                    SUPERCOOL[j] = SUPERCOOL[j]*DZSNSO[j]*1000.0 ; #[mm]
                end
            end
        end
    end

    for j in ISNOW+1:NSOIL
        if MICE[j] > 0.0 && STC[j] >= TFRZ
            IMELT[j] = 1 ; # melting
        end
        if MLIQ[j] > SUPERCOOL[j] && STC[j] < TFRZ
            IMELT[j] = 2 ;
        end

        if ISNOW == 0 && SNEQV > 0.0 && j == 1
            if STC[j] >= TFRZ
                IMELT[j] = 1 ; # If snow exists, but its thickness is not enough to create a layer
            end
        end
    end

    # Calculate the energy surplus and loss for melting and freezing
    for j in ISNOW+1:NSOIL
        if IMELT[j] > 0
            HM[j] = (STC[j]-TFRZ) / FACT[j] ;
            STC[j] = TFRZ ; 
        end
        if IMELT[j] == 1 && HM[j] < 0.0
            HM[j] = 0.0 ;
            IMELT[j] = 0.0 ;
        end
        if IMELT[j] == 2 && HM[j] > 0.0
            HM[j] = 0.0 ;
            IMELT[j] = 0.0 ;
        end

        XM[j] = HM[j] * dt / HFUS
    end

    # The rate of melting and freezing for snow without a layer
    if ISNOW == 0 && SNEQV > 0.0 && XM[1] > 0
        TEMP1 = SNEQV ;
        SNEQV = max(0.0, TEMP1-XM[1]) ;
        PROPOR = SNEQV/TEMP1 ;
        SNOWH = max(0, PROPOR*SNOWH) ;
        HEATR = HM[1] - HFUS*(TEMP1-SNEQV)/dt ;
        if HEATR > 0.0
            XM[1] = HEATR*dt/HFUS ;
            HM[1] = HEATR ;
        else
            XM[1] = 0.0 ;
            HM[1] = 0.0 ;
        end

        QMELT = max(0.0,(TEMP1-SNEQV))/dt ;
        XMF = HFUS*QMELT ;
        PONDING = TEMP1-SNEQV ;
    end

    # The rate of melting and freezing for snow and soil
    for j in ISNOW+1:NSOIL
        if IMELT[j] > 0.0 && abs(HM[j]) > 0.0
            HEATR = 0.0 ;
            if XM[j] > 0.0
                MICE[j] = max(0.0, WICE0[j]-XM[j]) ;
                HEATR = HM[j] - HFUS* (WICE0[j]-MICE[j]) / dt ;
            elseif XM[j] < 0.0
                if j <= NSNOW + 1
                    MICE[j] = min(WMASS0[j], WICE0[0]-XM[j]) ;
                else
                    if WMASS0[j] < SUPERCOOL[j]
                        MICE[j] = 0.0 ;
                    else
                        MICE[j] = max(0.0, min(WMASS0[j]-SUPERCOOL[j], WICE0[j]-XM[j])) ;
                    end
                end
                HEATR = HM[j] - HFUS* (WICE0[j]-MICE[j]) / dt ;
            end
            
            MLIQ[j] = max(0.0, WMASS0[j]-MICE[j]) ; 
            if abs(HEATR) > 0.0
                STC[j] += FACT[j]*HEATR ;
                if j <= NSNOW + 1 # snow layers
                    if MLIQ[j]*MICE[j] > 0.0
                        STC[j] = TFRZ ;
                    end
                end
            end

            XMF += HFUS * (WICE0[j]-MICE[j])/dt ;

            if j < NSNOW+1
                QMELT += max(0, WICE0[j]-MICE[j]) / dt ;
            end
        end
    end

    for j in NSNOW+1:NSOIL # soil layers
        SH2O[j] = MLIQ[j] / (1000.0 * DZSNSO[j])
    end

    return IMELT, QMELT, PONDING, STC, SH2O
end

function SFCDIF1(ITER, SFCTMP, RHOAIR, H, QAIR, ZLVL, Z0M, 
    Z0H, UR, MPE, MOZ, MOZSGN, FV, FM, FH)
    # computing surface drag coefficient CM for momentum and CH for heat

    TMPCM = log((2.0 + Z0M) / Z0M) ; # temporary CM
    TMPCH = log((2.0 + Z0H) / Z0H) ; # temporary CH
    MOZOLD = MOZ ; # Monin-Obukhov stability parameter moz for next iteration

    if ITER == 1
            FV = 0.0 ;
            MOZ = 0.0 ;
            MOL = 0.0 ;
            FM = 0.0 ;
    else
        TVIR = (1.0 + 0.61*QAIR) * SFCTMP ;
        TMP1 = VKC * (GRAV/TVIR) * H/(RHOAIR*CPAIR) ;
        if abs(TMP1) <= MPE
            TMP1 = MPE ;
        end
        MOL = -1.0 * FV^3 / TMP1 ;
        MOZ = min((2.0+Z0H)/MOL, 1.0) ;
    end

    # accumulate number of times moz changes sign.
    if MOZOLD*MOZ < 0.0
        MOZSGN = MOZSGN+1 ;
    end
    if MOZSGN >= 2
        MOZ = 0.0 ;
        FM = 0.0 ;
        FH = 0.0 ;
    end

    # evaluate stability-dependent variables using moz from prior iteration
    if MOZ < 0.0
        TMP1 = (1.0 - 16.0*MOZ)^0.25 ;
        TMP2 = log((1.0+TMP1*TMP1)/2.0) ;
        TMP3 = log((1.0+TMP1)/2.0) ;
        FMNEW = 2.0*TMP3 + TMP2 - 2.0*atan(TMP1) + 1.5707963 ;
        FHNEW = 2*TMP2 ;
    else
        FMNEW = -5.0*MOZ ;
        FHNEW = FMNEW ;
    end

    # except for first iteration, weight stability factors for previous iteration to help avoid flip-flops from one iteration to the next
    if ITER == 1
        FM = FMNEW ;
        FH = FHNEW ;
    else
        FM = 0.5 * (FM+FMNEW) ;
        FH = 0.5 * (FH+FHNEW) ;
    end

    FH = min(FH,0.9*TMPCH) ;
    FM = min(FM,0.9*TMPCM) ;

    CMFM = TMPCM-FM ;
    CHFH = TMPCH-FH ;
    if abs(CMFM) <= MPE
        CMFM = MPE ;
    end
    if abs(CHFH) <= MPE
        CHFH = MPE ;
    end
    CM  = VKC*VKC/(CMFM*CMFM) ;
    CH  = VKC*VKC/(CMFM*CHFH) ;

    FV = UR * sqrt(CM) ; # friction velocity 

    return MOZ, MOZSGN, FM, FH, CM, CH, FV
end

########################################################################
############################## MODELING ################################
########################################################################

output_folder = args[2] ;

surf_temp_all = zeros(Float16, (meteoro.height[1], meteoro.width[1])) ;

println("Start Surface Temprature Calculation")
for h in 1: meteoro.height[1]
    global meteoro; 
    if h%(round(meteoro.height[1]/100)) == 0
        print("#")
    end
    surf_temp = surface_temp(0.03, 0.03, h) ; # h refer to the permutation number
    surf_temp_all[:,h] = surf_temp ;
end
println("")
println("End Surface Temprature Calculation")

# write matrix
writedlm(output_folder*"matrix.csv", surf_temp_all)

# writing the raster
dataset = AG.read("/"*working_directory*"/skyview.tiff") ;
ref = AG.getproj(dataset) ;
geotransform = AG.getgeotransform(dataset) ;

outfile = "/"*output_folder*"/predicted_thermo.tif" ; # output_folder * "prediction.tif" ; 

AG.create(
    outfile,
    driver = AG.getdriver("GTiff"),
    width = ArchGDAL.width(dataset),
    height = ArchGDAL.height(dataset),
    nbands = 1,
    dtype = Float64
) do raster
    AG.setgeotransform!(raster, geotransform) ;
    AG.setproj!(raster, ref) ;
    AG.write!(
        raster,
        convert(Array{Float64}, surf_temp_all), # image to "burn" into the raster
        1, # update band 1
    )
end 

# dispaly the raster  
display(heatmap(surf_temp_all))
