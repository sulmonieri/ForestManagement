###
#' Desc: Crop rasters for area of interes (aoi) to change chm and compute microclimate
#' Created on 29.11.21 15:42
#' @author: sulmonie
###


library(sf)
library(terra)

setwd('ForestManagement')

## read nationwide rasters
vhm <- rast('nationwide/VHM.tif')
dtm <- rast('nationwide/DTM_5m.tif')
fmask <- rast('nationwide/ForestMask_10m.tif')
temp <- rast('nationwide/temperature.tif')
rain <- rast('nationwide/rain.tif')
twi <- rast('nationwide/TWI.tif')
tpi <- rast('nationwide/TPI.tif')
slope <- rast('nationwide/slope.tif')
aspect_n <- rast('nationwide/aspect_n.tif')
skyview <- rast('nationwide/skyview.tif')

## read shapfile
shp <- st_read('shp/BDM1.shp')

## crop rasters and write do disk
temp_crop <- crop(temp,shp,overwrite=TRUE,filename='aoi/temp_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'))
rain_crop <- crop(rain,shp,overwrite=TRUE,filename='aoi/rain_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'))
twi_crop <- crop(twi,shp,overwrite=TRUE,filename='aoi/twi_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'))
tpi_crop <- crop(tpi,shp,overwrite=TRUE,filename='aoi/tpi_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'))
slope_crop <- crop(slope,shp,overwrite=TRUE,filename='aoi/slope_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'))
aspect_n_crop <- crop(aspect_n,shp,overwrite=TRUE,filename='aoi/aspect_n_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'))
skyview_crop <- crop(skyview,shp,overwrite=TRUE,filename='aoi/skyview_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'))
fmask_crop <- crop(fmask,shp,overwrite=TRUE,filename='aoi/adapt_chm/fmask_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'))


## DTM needs to be disaggregated to 1m (VHM resolution)
dtm_crop <- disagg(crop(dtm,shp),5,method='bilinear',filename='aoi/adapt_chm/DTM_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'),overwrite=TRUE)
vhm_crop <- crop(vhm,shp,filename='aoi/adapt_chm/CHM_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'), overwrite=TRUE)

## generate digital surface model
dsm <- vhm_crop+dtm_crop
writeRaster(dsm,filename='aoi/adapt_chm/DSM_aoi.tif',gdal=c('COMPRESS=DEFLATE','PREDICTOR=3'),overwrite=TRUE)
