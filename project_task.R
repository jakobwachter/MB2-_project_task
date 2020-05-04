### analyzing scPDSI trends for different sites, bioclimatic divisions and Land Cover Classes in bavaria
### project work for MB1 - Programming and Geostatistics
### Jakob Wachter
### 18th February 2020

#########################################################################################
# in this script, scPDSI (self-calibrated Palmer Drought Severity Index) over bavaria   #
# should be analyzed. The scPDSI stack is applied to two vector layers, the landklif    #
# project quadrants and the bioclimatic regions of Bavaria. By this, regional patterns  #
# in the scPDSI trend should be recognized and can furthermore be interpreted.          #
# The input data is under https://github.com/jakobwachter/MB2-_project_task.            #
# The results of this script will be two updated versions of the input shapefiles       #
# for further analysis and graphs showing the scPDSI trends for the regions             #
# of the shapefiles that have been generated here.                                      #
#########################################################################################

## loading all required packages and setting working directory

if(!require(gdalUtils)){
  install.packages("gdalUtils")
  library(gdalUtils)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

if(!require(sp)){
  install.packages("sp")
  library(sp)
}

if(!require(rgdal)){
  install.packages("rgdal")
  library(rgdal)
}

if(!require(raster)){
  install.packages("raster")
  library(raster)
}

if(!require(snow)){
  install.packages("snow")
  library(snow)
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

setwd("D:\\JWachter\\project_task")

## loading input data

beginCluster(2)

# previously calculated scPDSI values in 100 fishnet polygons over germany - provided by Marius Philipp, added missing months of 2019 personally

scPDSI_stack<-brick("scPDSI_raster_JW.tif")

# only use monthly layers from 2000 to 2019 due to data amounts

scPDSI_stack2000<- brick(scPDSI_stack[[109:348]])

# Landklif project quadrants shapefile

landklif_quadrants<-readOGR("Final60Quadrants_epsg25832.shp")

# bioclimatic regions of bavaria shapefile

bioclim_regions<-readOGR("NatGlied_BAY_Meynen-Schmithuesen_UTM.shp")

## check the data for consistency
# raster data

plot(scPDSI_stack2000)
head(scPDSI_stack2000)
tail(scPDSI_stack2000)
mean(scPDSI_stack2000)

# vector data

## set projections 
# define wanted coordinate system
# in this case DHDN 3 Gauss-Krueger 4 

defaultproj<- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +datum=potsdam +units=m +no_defs+proj=tmerc +lat_0=0 +lon_0=12 +k=1 +x_0=4500000 +y_0=0 +ellps=bessel +datum=potsdam +units=m +no_defs"

# assign correct coordinate system to all data

scPDSI_stack2000@crs<-sp::CRS(defaultproj)
landklif_quadrants<-spTransform(landklif_quadrants,CRS(defaultproj))
bioclim_regions<-spTransform(bioclim_regions,CRS(defaultproj))

## after all data are in the same coordinate system, combine raster and vector data
# first, subset scPDSI data to the extent to bavaria --> not necessary for Landklif quadrants and bioclimatic regions

# get bavaria boundaries from the raster package

bnd <- raster::getData("GADM", country='DEU', level=1)
bnd <- spTransform(bnd, CRS(proj4string(scPDSI_stack2000)))
bnd_by <- bnd[bnd$NAME_1=="Bayern",]

plot(bnd_by)

# subset scPDSI stack

scPDSI_stack_by<-crop(scPDSI_stack2000,bnd_by)
scPDSI_stack_by<-mask(scPDSI_stack2000,bnd_by)

# extract scPDSI values at extent of different shapefiles

scPDSI_landklif_extracted<-as.matrix(extract(scPDSI_stack_by,landklif_quadrants))
scPDSI_bioclim_extracted<- as.matrix(extract(scPDSI_stack_by,bioclim_regions))

# calculate examplewise statistics --> scPDSI_landklif consists of 60 entities that represent the Landklif quadrants; scPDSI_bioclim consists of 96 elements representing the bioclimatic regions 

mean(scPDSI_landklif_extracted[[1]])
mean(scPDSI_landklif_extracted[[60]])

# apply statistical functions to both spatialPolygon dataframes using the lapply function 

scPDSI_landklif_extracted_mean<- lapply(scPDSI_landklif_extracted, mean)
scPDSI_landklif_extracted_min<- lapply(scPDSI_landklif_extracted, min)
scPDSI_landklif_extracted_max<- lapply(scPDSI_landklif_extracted, max)
scPDSI_landklif_extracted_sd<- lapply(scPDSI_landklif_extracted, sd)

scPDSI_bioclim_extracted_mean<- lapply(scPDSI_bioclim_extracted, mean)
scPDSI_bioclim_extracted_min<- lapply(scPDSI_bioclim_extracted, min)
scPDSI_bioclim_extracted_max<- lapply(scPDSI_bioclim_extracted, max)
scPDSI_bioclim_extracted_sd<- lapply(scPDSI_bioclim_extracted, sd)

## write the results back to the shapefiles extra columns and change their class types for writing back to the files

library(foreign)

landklif_quadrants$scPDSI_mean<-scPDSI_landklif_extracted_mean
landklif_quadrants$scPDSI_min<-scPDSI_landklif_extracted_min
landklif_quadrants$scPDSI_max<-scPDSI_landklif_extracted_max
landklif_quadrants$scPDSI_sd<-scPDSI_landklif_extracted_sd

class(landklif_quadrants$scPDSI_mean) = "integer"
class(landklif_quadrants$scPDSI_min) = "integer"
class(landklif_quadrants$scPDSI_max) = "integer"
class(landklif_quadrants$scPDSI_sd) = "integer"

bioclim_regions$scPDSI_mean<-scPDSI_bioclim_extracted_mean
bioclim_regions$scPDSI_min<-scPDSI_bioclim_extracted_min
bioclim_regions$scPDSI_max<-scPDSI_bioclim_extracted_max
bioclim_regions$scPDSI_sd<-scPDSI_bioclim_extracted_sd

class(bioclim_regions$scPDSI_mean) = "integer"
class(bioclim_regions$scPDSI_min) = "integer"
class(bioclim_regions$scPDSI_max) = "integer"
class(bioclim_regions$scPDSI_sd) = "integer"

# write the results back to the files

require(maptools)
writeSpatialShape(landklif_quadrants,"Final60Quadrants_epsg25832.shp")
writeSpatialShape(bioclim_regions,"NatGlied_BAY_Meynen-Schmithuesen_UTM.shp")

## plot results for nice visualization

timesteps<-seq(as.Date("2000/1/1"), by = "month",length.out = 240) 

scPDSI_landklif_extracted_mean_mat<-as.matrix(scPDSI_landklif_extracted_mean)
scPDSI_landklif_extracted_min_mat<-as.matrix(scPDSI_landklif_extracted_min)
scPDSI_landklif_extracted_max_mat<-as.matrix(scPDSI_landklif_extracted_max)
scPDSI_landklif_extracted_sd_mat<-as.matrix(scPDSI_landklif_extracted_sd)

DF_Landklif_mean<-cbind(scPDSI_landklif_extracted_mean,as.Date(timesteps))
DF_Landklif_mean<-as.data.frame(DF_Landklif_mean)
DF_Landklif_mean$timesteps<- as.Date(DF_Landklif_mean$timesteps)