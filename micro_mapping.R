###
#' Desc: Map the microclimate with a linear mixed effects model
#' Created on 29.11.21 15:42
#' @author: sulmonie
###


library(sf)
library(terra)

setwd('ForestManagement')

## get predictor variables
vhm_sce_raw <- rast('aoi/adapt_chm/results/chm_random_0.25_fm1_buffer0.0m_10m.tif')
vhm_sce <- subst(vhm_sce_raw,NaN,0.02)
vhm <- aggregate(rast('aoi/adapt_chm/CHM.tif'),10)
temp <- rast('aoi/temp_aoi.tif')
rain <- rast('aoi/rain_aoi.tif')
twi <- rast('aoi/twi_aoi.tif')
tpi <- rast('aoi/tpi_aoi.tif')
slope <- rast('aoi/slope_aoi.tif')
aspect_n <- rast('aoi/aspect_n_aoi.tif')
skyview <- rast('aoi/skyview_aoi.tif')
potrad <- rast('aoi/potrad_aoi.tif')
trans_sce <- rast('aoi/transm_6_25_cut.tif')
transpot_sce <- potrad*trans_sce
trans <- rast('aoi/transm_6_orig.tif')
transpot <- potrad*trans


## stack the rasters
pred_sce <- c(temp,transpot_sce,vhm_sce,rain,twi,tpi,slope,aspect_n,skyview)
names(pred_sce) <- c('Tmax_meteo','trans_pot','vegh','Rain','topo_wetness','topo_index','slope','aspect_n','skyview')
pred <- c(temp,transpot,vhm,rain,twi,tpi,slope,aspect_n,skyview)
names(pred) <- c('Tmax_meteo','trans_pot','vegh','Rain','topo_wetness','topo_index','slope','aspect_n','skyview')

## get Microclimate model and dataset
load('lme.RData')

## Predict microclimate
micro_sce <- predict(pred_sce,lmes$low$lme_Tmax,re.form=NA)/100
writeRaster(micro_sce,'micro_maps/micro_map_sce.tif',overwrite=TRUE)
micro <- predict(pred,lmes$low$lme_Tmax,re.form=NA)/100
writeRaster(micro,'micro_maps/micro_map.tif',overwrite=TRUE)
diff <- (micro_sce-micro)
writeRaster(diff,'micro_maps/micro_map_diff.tif',overwrite=TRUE)
