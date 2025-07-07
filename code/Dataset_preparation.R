# This script was used to prepare dataframes for ML models building and prediction.
# It took ~ 4 min to execute this script on our server.
# We have saved and uploaded the resulting dataframes, so you can
# skip this script execution and scroll down to the bottom.

library(tidyverse)
library(ncdf4)

# Functions `grid_lat` and `grid_lon` are used to 
# regrid longitude latitude to 1x1
grid_lat <- function(lat, bnds = 1, lat_min = -89.5, lat_max = 89.5){
  bnds_center <- seq(from = lat_min, to = lat_max, by = bnds)
  bnds_center_index <- which.min(abs(lat - bnds_center))
  bnds_center[bnds_center_index]
}

grid_lon <- function(lon, bnds = 1, lon_min = -179.5, lon_max = 179.5){
  bnds_center <- seq(from = lon_min, to = lon_max, by = bnds)
  bnds_center_index <- which.min(abs(lon - bnds_center))
  bnds_center[bnds_center_index]
}


# `Shao_Global_Tricho_nifH_existence.csv` is derived and modified from 
#  Shao et al.(2023, ESSD) to get dataframe `Tricho_nifH_existence_ds_regrided`

# Shao et al.(2023) have comprehensively collect Trichodesmium
# nifH abundance data. We transformed the abundance data into 
# presence-absence data. Because of the detection limit of qPCR,
# no detection of Tricho nifH genes should not be simply interpreted 
# as "absence". Instead, it indicates that Trichodesmium is absence or 
#of low abundance at the site.   

Tricho_nifH_existence_ds_regrided <- # regrid to 1X1 
  read_csv("data/Shao_Global_Tricho_nifH_existence.csv",
           col_types = cols(DATE = col_date(format = "%m/%d/%Y")),
           skip = 1) %>%
  filter(Data_Flag == 0) %>% # remove bad data points
  mutate(Sampling_month = month(DATE)) %>%
  rowwise() %>%
  mutate(lat_regrid = grid_lat(LATITUDE), # 
         lon_regrid = grid_lon(LONGITUDE) #
  ) %>%
  ungroup() %>%
  dplyr::select(-c(Station, Notes_nifH)) %>%
  group_by(lat_regrid, lon_regrid, Sampling_month) %>%
  summarise(Trichodesmium_nifH_existence_regrid = 
              if_else((mean(Tricho_exist, na.rm = T) * length(Tricho_exist)) > 0, T, F)) %>%
  ungroup()


# `Shao_Global_Tricho_nifH_integral.csv` is derived and modified from 
#  Shao et al.(2023, ESSD) to get dataframe `Tricho_nifH_ds_regrided`

Tricho_nifH_ds_regrided <-
  read_csv("data/Shao_Global_Tricho_nifH_integral.csv",
           col_types = cols(DATE = col_date(format = "%m/%d/%Y")),
           skip = 1) %>%
  filter(Data_Flag == 0) %>% # remove bad data points
  mutate(Sampling_month = month(DATE)) %>%
  rowwise() %>%
  mutate(lat_regrid = grid_lat(LATITUDE),
         lon_regrid = grid_lon(LONGITUDE)
  ) %>%
  ungroup() %>%
  dplyr::select(-c(Station, Notes_nifH)) %>%
  group_by(lat_regrid, lon_regrid, Sampling_month) %>%
  summarise(Trichodesmium_nifH_integral_regrid = mean(Trichodesmium_nifH_integral, na.rm = T)) %>% # 10^6 copies/m2
  ungroup()

# SST
## SST monthly climatology is downloaded from
## https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/
WOA2023_SST_monthly <- array(dim = c(360, 180, 12))
for (i in 1:12){
  if (i<10) {
    WOA_file_name <- paste0("data/environmental_variables/SST/WOA2023/woa23_decav_t0", i, "_01.nc")
  } else {
    WOA_file_name <- paste0("data/environmental_variables/SST/WOA2023/woa23_decav_t", i, "_01.nc")
  }
  Temp_WOA2023_nc <- 
    nc_open(WOA_file_name)
  
  Temp_WOA2023_nc_t_an <- 
    ncvar_get(Temp_WOA2023_nc, "t_an")  #[lon,lat,depth]
  
  Temp_WOA2023_nc_t_an_SST <- Temp_WOA2023_nc_t_an[ , , 1]
  WOA2023_SST_monthly[ , , i] <- Temp_WOA2023_nc_t_an_SST
}

## Function Match_SST is used to match SST and Trichodesmium nifH data
Match_SST <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index <- lat-(-89.5)+1
  lon_index <- lon-(-179.5)+1
  SST <- WOA2023_SST_monthly[lon_index, lat_index, mon]
  if(!is.na(SST)){
    SST
  } else{
    mean(c(WOA2023_SST_monthly[lon_index-1, lat_index, mon],
           WOA2023_SST_monthly[lon_index+1, lat_index, mon],
           WOA2023_SST_monthly[lon_index, lat_index-1, mon],
           WOA2023_SST_monthly[lon_index, lat_index+1, mon],
           WOA2023_SST_monthly[lon_index-1, lat_index-1, mon],
           WOA2023_SST_monthly[lon_index+1, lat_index+1, mon],
           WOA2023_SST_monthly[lon_index-1, lat_index+1, mon],
           WOA2023_SST_monthly[lon_index+1, lat_index-1, mon]),
         na.rm = T)
  }
}

# `SST_for_Tricho_nifH_regrided` contains SST data corresponding to
# Tricho nifH data in dataframes `Tricho_nifH_existence_ds_regrided` and
# `Tricho_nifH_ds_regrided`
SST_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_SST)


# DIP/Pi 
## DIP/Pi monthly climatology is downloaded from
## https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/
WOA2023_Pi_monthly <- array(dim = c(360, 180, 12))

for (i in 1:12){
  if (i<10) {
    WOA_file_name <- paste0("data/environmental_variables/DIP/WOA2023/woa23_all_p0", i, "_01.nc")
  } else {
    WOA_file_name <- paste0("data/environmental_variables/DIP/WOA2023/woa23_all_p", i, "_01.nc")
  }
  Pi_WOA2023_nc <- 
    nc_open(WOA_file_name)
  
  Pi_WOA2023_nc_p_an <- 
    ncvar_get(Pi_WOA2023_nc, "p_an")  #[lon,lat,depth]
  
  Pi_WOA2023_nc_p_an_0m <- Pi_WOA2023_nc_p_an[ , , 1]
  WOA2023_Pi_monthly[ , , i] <- Pi_WOA2023_nc_p_an_0m
}

Match_Pi <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index <- lat-(-89.5)+1
  lon_index <- lon-(-179.5)+1
  Pi <- WOA2023_Pi_monthly[lon_index, lat_index, mon]
  if(!is.na(Pi)){
    Pi
  } else{
    mean(c(WOA2023_Pi_monthly[lon_index-1, lat_index, mon],
           WOA2023_Pi_monthly[lon_index+1, lat_index, mon],
           WOA2023_Pi_monthly[lon_index, lat_index-1, mon],
           WOA2023_Pi_monthly[lon_index, lat_index+1, mon],
           WOA2023_Pi_monthly[lon_index-1, lat_index-1, mon],
           WOA2023_Pi_monthly[lon_index+1, lat_index+1, mon],
           WOA2023_Pi_monthly[lon_index-1, lat_index+1, mon],
           WOA2023_Pi_monthly[lon_index+1, lat_index-1, mon]),
         na.rm = T)
  }
}

Pi_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_Pi)

# DO (dissolved oxygen) minimum within 500 m
## DO monthly climatology is downloaded from 
## https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/
WOA2023_DO_500min_monthly <- array(dim = c(360, 180, 12)) # mininum DO within 500m
for (i in 1:12){
  if (i<10) {
    WOA_file_name <- paste0("data/environmental_variables/DO/WOA2023/woa23_all_o0", i, "_01.nc")
  } else {
    WOA_file_name <- paste0("data/environmental_variables/DO/WOA2023/woa23_all_o", i, "_01.nc")
  }
  DO_WOA2023_nc <- 
    nc_open(WOA_file_name)
  
  DO_WOA2023_nc_o_an <- 
    ncvar_get(DO_WOA2023_nc, "o_an")  #[lon,lat,depth]
  
  DO_WOA2023_nc_o_an_within_500m <- DO_WOA2023_nc_o_an[ , , 1:37]
  
  WOA2023_DO_500min_monthly[ , , i] <- apply(DO_WOA2023_nc_o_an_within_500m, 1:2, min, na.rm = T)
  WOA2023_DO_500min_monthly[WOA2023_DO_500min_monthly==Inf] <- NA
}

Match_DO_500min <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index <- lat-(-89.5)+1
  lon_index <- lon-(-179.5)+1
  DO_500min <- WOA2023_DO_500min_monthly[lon_index, lat_index, mon]
  if(!is.na(DO_500min)){
    DO_500min
  } else{
    mean(c(WOA2023_DO_500min_monthly[lon_index-1, lat_index, mon],
           WOA2023_DO_500min_monthly[lon_index+1, lat_index, mon],
           WOA2023_DO_500min_monthly[lon_index, lat_index-1, mon],
           WOA2023_DO_500min_monthly[lon_index, lat_index+1, mon],
           WOA2023_DO_500min_monthly[lon_index-1, lat_index-1, mon],
           WOA2023_DO_500min_monthly[lon_index+1, lat_index+1, mon],
           WOA2023_DO_500min_monthly[lon_index-1, lat_index+1, mon],
           WOA2023_DO_500min_monthly[lon_index+1, lat_index-1, mon]),
         na.rm = T)
  }
}

DO_500min_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_DO_500min)

# Nitrate
## Nitrate monthly climatology is downloaded from 
## https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/
WOA2023_Nitrate_monthly <- array(dim = c(360, 180, 12))
for (i in 1:12){
  if (i<10) {
    WOA_file_name <- paste0("data/environmental_variables/DIN/WOA2023/woa23_all_n0", i, "_01.nc")
  } else {
    WOA_file_name <- paste0("data/environmental_variables/DIN/WOA2023/woa23_all_n", i, "_01.nc")
  }
  Nitrate_WOA2023_nc <- 
    nc_open(WOA_file_name)
  
  Nitrate_WOA2023_nc_p_an <- 
    ncvar_get(Nitrate_WOA2023_nc, "n_an")  #[lon,lat,depth]
  
  Nitrate_WOA2023_nc_n_an_0m <- Nitrate_WOA2023_nc_p_an[ , , 1]
  WOA2023_Nitrate_monthly[ , , i] <- Nitrate_WOA2023_nc_n_an_0m
}

Match_Nitrate <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index <- lat-(-89.5)+1
  lon_index <- lon-(-179.5)+1
  Nitrate <- WOA2023_Nitrate_monthly[lon_index, lat_index, mon]
  if(!is.na(Nitrate)){
    Nitrate
  } else{
    mean(c(WOA2023_Nitrate_monthly[lon_index-1, lat_index, mon],
           WOA2023_Nitrate_monthly[lon_index+1, lat_index, mon],
           WOA2023_Nitrate_monthly[lon_index, lat_index-1, mon],
           WOA2023_Nitrate_monthly[lon_index, lat_index+1, mon],
           WOA2023_Nitrate_monthly[lon_index-1, lat_index-1, mon],
           WOA2023_Nitrate_monthly[lon_index+1, lat_index+1, mon],
           WOA2023_Nitrate_monthly[lon_index-1, lat_index+1, mon],
           WOA2023_Nitrate_monthly[lon_index+1, lat_index-1, mon]),
         na.rm = T)
  }
}

Nitrate_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_Nitrate)

# PAR
## PAR monthly climatology is downloaded from 
## https://oceandata.sci.gsfc.nasa.gov/l3/
MODIS_PAR_monthly <- array(dim = c(4320, 2160, 12))
for (i in 1:12){
  if (i<10) {
    file_name <- paste0("data/environmental_variables/PAR/AQUA_MODIS.0",i,".L3m.MC.PAR.par.9km.nc")
  } else {
    file_name <- paste0("data/environmental_variables/PAR/AQUA_MODIS.",i,".L3m.MC.PAR.par.9km.nc")
  }
  PAR_MODIS_nc <- 
    nc_open(file_name)
  
  PAR_MODIS_nc_par <- 
    ncvar_get(PAR_MODIS_nc, "par")  #[lon,lat]
  
  MODIS_PAR_monthly[ , , i] <- PAR_MODIS_nc_par
}

MODIS_lon <- ncvar_get(PAR_MODIS_nc, "lon")
MODIS_lat <- ncvar_get(PAR_MODIS_nc, "lat")

Match_PAR <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index_range <- (MODIS_lat <= lat+0.5)&(MODIS_lat >= lat-0.5)
  lon_index_range <- (MODIS_lon <= lon+0.5)&(MODIS_lon >= lon-0.5)
  PAR <- 
    MODIS_PAR_monthly[lon_index_range, lat_index_range, mon] %>%
    mean(na.rm = T)
  # lat_index <- which.min(abs(lat-MODIS_lat))
  # lon_index <- which.min(abs(lon-MODIS_lon))
  # PAR <- MODIS_PAR_monthly[lon_index, lat_index, mon]
  # if(!is.na(PAR)){
  #   PAR
  # } else{
  #   mean(c(MODIS_PAR_monthly[lon_index-1, lat_index, mon],
  #          MODIS_PAR_monthly[lon_index+1, lat_index, mon],
  #          MODIS_PAR_monthly[lon_index, lat_index-1, mon],
  #          MODIS_PAR_monthly[lon_index, lat_index+1, mon],
  #          MODIS_PAR_monthly[lon_index-1, lat_index-1, mon],
  #          MODIS_PAR_monthly[lon_index+1, lat_index+1, mon],
  #          MODIS_PAR_monthly[lon_index-1, lat_index+1, mon],
  #          MODIS_PAR_monthly[lon_index+1, lat_index-1, mon]),
  #        na.rm = T)
  # }
}
PAR_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_PAR)

# Kd490
## Kd490 monthly climatology is downloaded from 
## https://oceandata.sci.gsfc.nasa.gov/l3/
MODIS_Kd490_monthly <- array(dim = c(4320, 2160, 12))
for (i in 1:12){
  if (i<10) {
    file_name <- paste0("data/environmental_variables/Kd490/AQUA_MODIS.0", i, ".L3m.MC.KD.Kd_490.9km.nc")
  } else {
    file_name <- paste0("data/environmental_variables/Kd490/AQUA_MODIS.", i, ".L3m.MC.KD.Kd_490.9km.nc")
  }
  KD490_MODIS_nc <- 
    nc_open(file_name)
  
  Kd490_MODIS_nc_par <- 
    ncvar_get(KD490_MODIS_nc, "Kd_490")  #[lon,lat]
  
  MODIS_Kd490_monthly[ , , i] <- Kd490_MODIS_nc_par
}

MODIS_lon <- ncvar_get(KD490_MODIS_nc, "lon")
MODIS_lat <- ncvar_get(KD490_MODIS_nc, "lat")

Match_Kd490 <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index_range <- (MODIS_lat <= lat+0.5)&(MODIS_lat >= lat-0.5)
  lon_index_range <- (MODIS_lon <= lon+0.5)&(MODIS_lon >= lon-0.5)
  Kd490 <-
    MODIS_Kd490_monthly[lon_index_range, lat_index_range, mon] %>%
    mean(na.rm = T)
  # lat_index <- which.min(abs(lat-MODIS_lat))
  # lon_index <- which.min(abs(lon-MODIS_lon))
  # Kd490 <- MODIS_Kd490_monthly[lon_index, lat_index, mon]
  # if(!is.na(Kd490)){
  #   Kd490
  # } else{
  #   mean(c(MODIS_Kd490_monthly[lon_index-1, lat_index, mon],
  #          MODIS_Kd490_monthly[lon_index+1, lat_index, mon],
  #          MODIS_Kd490_monthly[lon_index, lat_index-1, mon],
  #          MODIS_Kd490_monthly[lon_index, lat_index+1, mon],
  #          MODIS_Kd490_monthly[lon_index-1, lat_index-1, mon],
  #          MODIS_Kd490_monthly[lon_index+1, lat_index+1, mon],
  #          MODIS_Kd490_monthly[lon_index-1, lat_index+1, mon],
  #          MODIS_Kd490_monthly[lon_index+1, lat_index-1, mon]),
  #        na.rm = T)
  # }
}
Kd490_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_Kd490)

# Chl a
## Chl a monthly climatology is downloaded from 
## https://oceandata.sci.gsfc.nasa.gov/l3/

MODIS_Chla_monthly <- array(dim = c(4320, 2160, 12))
for (i in 1:12){
  if (i<10) {
    file_name <- paste0("data/environmental_variables/Chla/AQUA_MODIS.0", i, ".L3m.MC.CHL.chlor_a.9km.nc")
  } else {
    file_name <- paste0("data/environmental_variables/Chla/AQUA_MODIS.", i, ".L3m.MC.CHL.chlor_a.9km.nc")
  }
  Chla_MODIS_nc <- 
    nc_open(file_name)

  Chla_MODIS_nc_chlor_a <- 
    ncvar_get(Chla_MODIS_nc, "chlor_a")  #[lon,lat]
  
  MODIS_Chla_monthly[ , , i] <- Chla_MODIS_nc_chlor_a
}

Match_Chla <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index_range <- (MODIS_lat <= lat+0.5)&(MODIS_lat >= lat-0.5)
  lon_index_range <- (MODIS_lon <= lon+0.5)&(MODIS_lon >= lon-0.5)
  Chla <-
    MODIS_Chla_monthly[lon_index_range, lat_index_range, mon] %>%
    mean(na.rm = T)
  # lat_index <- which.min(abs(lat-MODIS_lat))
  # lon_index <- which.min(abs(lon-MODIS_lon))
  # MODIS_Chla_monthly[lon_index, lat_index, mon]
}

Chla_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_Chla)

# SSS
## SSS monthly climatology is downloaded from 
## https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/

WOA2023_SSS_monthly <- array(dim = c(360, 180, 12))

for (i in 1:12){
  if (i<10) {
    WOA_file_name <- paste0("data/environmental_variables/SSS/WOA2023/woa23_decav_s0", i, "_01.nc")
  } else {
    WOA_file_name <- paste0("data/environmental_variables/SSS/WOA2023/woa23_decav_s", i, "_01.nc")
  }
  Sal_WOA2023_nc <- #Salinity
    nc_open(WOA_file_name)
  
  Sal_WOA2023_nc_s_an <- 
    ncvar_get(Sal_WOA2023_nc, "s_an")  #[lon,lat,depth]
  
  Sal_WOA2023_nc_s_an_SSS <- Sal_WOA2023_nc_s_an[ , , 1]
  WOA2023_SSS_monthly[ , , i] <- Sal_WOA2023_nc_s_an_SSS
}

Match_SSS <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index <- lat-(-89.5)+1
  lon_index <- lon-(-179.5)+1
  SSS <- WOA2023_SSS_monthly[lon_index, lat_index, mon]
  if(!is.na(SSS)){
    SSS
  } else{
    mean(c(WOA2023_SSS_monthly[lon_index-1, lat_index, mon],
           WOA2023_SSS_monthly[lon_index+1, lat_index, mon],
           WOA2023_SSS_monthly[lon_index, lat_index-1, mon],
           WOA2023_SSS_monthly[lon_index, lat_index+1, mon],
           WOA2023_SSS_monthly[lon_index-1, lat_index-1, mon],
           WOA2023_SSS_monthly[lon_index+1, lat_index+1, mon],
           WOA2023_SSS_monthly[lon_index-1, lat_index+1, mon],
           WOA2023_SSS_monthly[lon_index+1, lat_index-1, mon]),
         na.rm = T)
  }
}

SSS_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_SSS)

# MLD
## MLD monthly climatology is downloaded from 
## https://www.ncei.noaa.gov/access/world-ocean-atlas-2023/
WOA2018_MLD_monthly <- array(dim = c(360, 180, 12))

for (i in 1:12){
  if (i<10) {
    WOA_file_name <- paste0("data/environmental_variables/MLD/WOA2018/woa18_decav81B0_M020", i, "_01.nc")
  } else {
    WOA_file_name <- paste0("data/environmental_variables/MLD/WOA2018/woa18_decav81B0_M02", i, "_01.nc")
  }
  MLD_WOA2018_nc <- 
    nc_open(WOA_file_name)
  
  MLD_WOA2018_nc_M_an <- 
    ncvar_get(MLD_WOA2018_nc, "M_an")  #[lon,lat,depth]
  
  WOA2018_MLD_monthly[ , , i] <- MLD_WOA2018_nc_M_an
}

Match_MLD <- function(df){
  lat <- df$lat_regrid
  lon <- df$lon_regrid
  mon <- df$Sampling_month
  lat_index <- lat-(-89.5)+1
  lon_index <- lon-(-179.5)+1
  MLD <- WOA2018_MLD_monthly[lon_index, lat_index, mon]
  if(!is.na(MLD)){
    MLD
  } else{
    mean(c(WOA2018_MLD_monthly[lon_index-1, lat_index, mon],
           WOA2018_MLD_monthly[lon_index+1, lat_index, mon],
           WOA2018_MLD_monthly[lon_index, lat_index-1, mon],
           WOA2018_MLD_monthly[lon_index, lat_index+1, mon],
           WOA2018_MLD_monthly[lon_index-1, lat_index-1, mon],
           WOA2018_MLD_monthly[lon_index+1, lat_index+1, mon],
           WOA2018_MLD_monthly[lon_index-1, lat_index+1, mon],
           WOA2018_MLD_monthly[lon_index+1, lat_index-1, mon]),
         na.rm = T)
  }
}

MLD_for_Tricho_nifH_regrided <-
  Tricho_nifH_ds_regrided %>%
  rowwise() %>%
  group_split() %>%
  map_dbl(Match_MLD)

# Dataframes `DS_for_ML_building_Tricho_nifH` and 
# `DS_for_ML_building_Tricho_nifH_existence` were used to build ML models

DS_for_ML_building_Tricho_nifH <-
  Tricho_nifH_ds_regrided %>%
  mutate(SST = SST_for_Tricho_nifH_regrided,
         SSS = SSS_for_Tricho_nifH_regrided,
         DIP = Pi_for_Tricho_nifH_regrided,
         DIN = Nitrate_for_Tricho_nifH_regrided,
         DO_500min = DO_500min_for_Tricho_nifH_regrided,
         Chla = Chla_for_Tricho_nifH_regrided,
         PAR = PAR_for_Tricho_nifH_regrided,
         kd_490 = Kd490_for_Tricho_nifH_regrided,
         MLD = MLD_for_Tricho_nifH_regrided) %>% #rename the log transformed variable
  mutate(loc_A = sin(lat_regrid/180*pi),
         loc_B = sin(lon_regrid/180*pi)*cos(lat_regrid/180*pi),
         loc_C = -cos(lon_regrid/180*pi)*cos(lat_regrid/180*pi),
         month_cos = cos(Sampling_month/12*2*pi),
         month_sin = sin(Sampling_month/12*2*pi))

DS_for_ML_building_Tricho_nifH_existence <-
  Tricho_nifH_existence_ds_regrided %>%
  mutate(SST = SST_for_Tricho_nifH_regrided,
         SSS = SSS_for_Tricho_nifH_regrided,
         DIP = Pi_for_Tricho_nifH_regrided,
         DIN = Nitrate_for_Tricho_nifH_regrided,
         DO_500min = DO_500min_for_Tricho_nifH_regrided,
         Chla = Chla_for_Tricho_nifH_regrided,
         PAR = PAR_for_Tricho_nifH_regrided,
         kd_490 = Kd490_for_Tricho_nifH_regrided,
         MLD = MLD_for_Tricho_nifH_regrided) %>%
  mutate(loc_A = sin(lat_regrid/180*pi),
         loc_B = sin(lon_regrid/180*pi)*cos(lat_regrid/180*pi),
         loc_C = -cos(lon_regrid/180*pi)*cos(lat_regrid/180*pi),
         month_cos = cos(Sampling_month/12*2*pi),
         month_sin = sin(Sampling_month/12*2*pi))


# Dataframe `Predictors_for_model_prediction` were used to do the final 
# global predictions.

## MODIS data needed to be regrided to 1X1 degree

lat_grid <- seq(from = -89.5, to = 89.5, by = 1)
lon_grid <- seq(from = -179.5, to = 179.5, by = 1)
MODIS_lat <- ncvar_get(nc_open("data/environmental_variables/Chla/AQUA_MODIS.01.L3m.MC.CHL.chlor_a.9km.nc"), "lat")
MODIS_lon <- ncvar_get(nc_open("data/environmental_variables/Chla/AQUA_MODIS.01.L3m.MC.CHL.chlor_a.9km.nc"), "lon")

lat_lon_crossed <- crossing(lat_grid, lon_grid)

## Chla
Regrided_Chla_monthly <- array(dim = c(360, 180, 12))
Regrid_MODIS_Chla <- function(lat_lon_crossed) {
  lon_grid <- lat_lon_crossed$lon_grid
  lat_grid <- lat_lon_crossed$lat_grid
  lon_index_vec <- (MODIS_lon <= lon_grid+0.5) & (MODIS_lon >= lon_grid-0.5)
  lat_index_vec <- (MODIS_lat <= lat_grid+0.5) & (MODIS_lat >= lat_grid-0.5)
  MODIS_Chla_monthly[lon_index_vec, lat_index_vec, month] %>%
    mean(na.rm = T)
} 
for (month in 1:12) {
  Regrided_Chla_monthly[ , , month] <-
    lat_lon_crossed %>%
    rowwise() %>%
    group_split() %>%
    map_vec(Regrid_MODIS_Chla) %>%
    matrix(nrow = 360, ncol = 180)
}



## Kd490
Regrided_Kd490_monthly <- array(dim = c(360, 180, 12))

Regrid_MODIS_Kd490 <- function(lat_lon_crossed) {
  lon_grid <- lat_lon_crossed$lon_grid
  lat_grid <- lat_lon_crossed$lat_grid
  lon_index_vec <- (MODIS_lon <= lon_grid+0.5) & (MODIS_lon >= lon_grid-0.5)
  lat_index_vec <- (MODIS_lat <= lat_grid+0.5) & (MODIS_lat >= lat_grid-0.5)
  MODIS_Kd490_monthly[lon_index_vec, lat_index_vec, month] %>%
    mean(na.rm = T)
}
for (month in 1:12) {
  Regrided_Kd490_monthly[ , , month] <-
    lat_lon_crossed %>%
    rowwise() %>%
    group_split() %>%
    map_vec(Regrid_MODIS_Kd490) %>%
    matrix(nrow = 360, ncol = 180)
}

## PAR
Regrided_PAR_monthly <- array(dim = c(360, 180, 12))

Regrid_MODIS_PAR <- function(lat_lon_crossed) {
  lon_grid <- lat_lon_crossed$lon_grid
  lat_grid <- lat_lon_crossed$lat_grid
  lon_index_vec <- (MODIS_lon <= lon_grid+0.5) & (MODIS_lon >= lon_grid-0.5)
  lat_index_vec <- (MODIS_lat <= lat_grid+0.5) & (MODIS_lat >= lat_grid-0.5)
  MODIS_PAR_monthly[lon_index_vec, lat_index_vec, month] %>%
    mean(na.rm = T)
}

for (month in 1:12) {
  Regrided_PAR_monthly[ , , month] <-
    lat_lon_crossed %>%
    rowwise() %>%
    group_split() %>%
    map_vec(Regrid_MODIS_PAR) %>%
    matrix(nrow = 360, ncol = 180)
}


Sampling_month <- 1:12
lat_regrid <- seq(from = -89.5, to = 89.5, by = 1)
lon_regrid <- seq(from = -179.5, to = 179.5, by = 1)

Predictors_for_model_prediction <-
  crossing(Sampling_month, lat_regrid, lon_regrid) %>%
  mutate(SST = c(WOA2023_SST_monthly),
         SSS = c(WOA2023_SSS_monthly),
         DIP = c(WOA2023_Pi_monthly),
         DIN = c(WOA2023_Nitrate_monthly),
         DO_500min = c(WOA2023_DO_500min_monthly),
         MLD = c(WOA2018_MLD_monthly),
         Chla = c(Regrided_Chla_monthly),
         kd_490 = c(Regrided_Kd490_monthly),
         PAR = c(Regrided_PAR_monthly)) %>%
  #mutate(across(.cols = all_of(vars_to_log), ~log10(.x + .1))) %>%  #log transform those variables, shift a bit from 0 as well
  #rename_with( ~ str_c("Log_", vars_to_log), .cols = all_of(vars_to_log)) %>% #rename the log transformed variable
  mutate(loc_A = sin(lat_regrid/180*pi),
         loc_B = sin(lon_regrid/180*pi)*cos(lat_regrid/180*pi),
         loc_C = -cos(lon_regrid/180*pi)*cos(lat_regrid/180*pi),
         month_cos = cos(Sampling_month/12*2*pi),
         month_sin = sin(Sampling_month/12*2*pi)) %>%
  arrange(Sampling_month, lon_regrid, lat_regrid)

# Only keep the dataframes that were used for model building and prediction
# Remove the intermediate variables
objects_to_keep <- c("DS_for_ML_building_Tricho_nifH",
                     "DS_for_ML_building_Tricho_nifH_existence",
                     "Predictors_for_model_prediction")
rm(list = setdiff(ls(), objects_to_keep))

# The three dataframes `DS_for_ML_building_Tricho_nifH`,
# `DS_for_ML_building_Tricho_nifH_existence`,
# and `Predictors_for_model_prediction` were saved to the file 
# `DS_prepared_for_model_building.RData`

# save.image("data/data_R_exported/DS_prepared_for_model_building.RData")
