### analyzing scPDSI trends for different sites, representing Landklif project quadrants in bavaria
### project work for MB1 - Programming and Geostatistics
### Jakob Wachter
### 11th May 2020

#########################################################################################
# in this script, scPDSI (self-calibrated Palmer Drought Severity Index) over bavaria   #
# should be analyzed. The scPDSI stack is applied to a vector layer, representing       #
# the landklif project quadrants. By this, regional patterns                            #
# in the scPDSI trend should be recognized and can furthermore be interpreted.          #
# The input data is under https://github.com/jakobwachter/MB2-_project_task.            #
# It is also available on the University M-Drive in the Landklif-Folder, or             #
# Carina Kübert as a responsible person for the project can provide it.                 #
# The results of this script will be an updated version of the input shapefile of the   #
# Landklif Quadrants for further analysis, including statistical scPDSI data for each   #
# quadrant, and graphs showing the scPDSI trends as a map as well as statistical        #
# illustrations.                                                                        #                                                                              
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

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

if(!require(rgeos)){
  install.packages("rgeos")
  library(rgeos)
}

if(!require(tidyr)){
  install.packages("tidyr")
  library(tidyr)
}

if(!require(sf)){
  install.packages("sf")
  library(sf)
}

if(!require(stringr)){
  install.packages("stringr")
  library(stringr)
}

working.dir<- choose.dir(caption = "Select working directory adapted to your computer")

### Change Directory to previously defined working directory
setwd(working.dir)

## loading input data
# create a cluster to optimize the available computing power for the script, this can be 
# adapted dependent on the number of cores of your computer.

beginCluster(2)

# previously calculated scPDSI values in 100 fishnet polygons over germany - 
# provided by Marius Philipp, added missing months of 2019 personally

scPDSI_stack<-brick("scPDSI_raster_JW.tif")

# only use monthly layers from 2000 to 2019 due to data amounts

scPDSI_stack2000<- brick(scPDSI_stack[[109:348]])

# Landklif project quadrants shapefile

landklif_quadrants<-readOGR("Final60Quadrants_epsg25832.shp")

## check the data for consistency
# raster data

head(scPDSI_stack2000)
tail(scPDSI_stack2000)
plot(scPDSI_stack2000)

# vector data

head(landklif_quadrants)
tail(landklif_quadrants)
plot(landklif_quadrants)

## set projections 
# define wanted coordinate system
# in this case DHDN 3 Gauss-Krueger 4 

defaultproj<- "+proj=tmerc +lat_0=0 +lon_0=9 +k=1 +x_0=3500000 +y_0=0 +ellps=bessel +datum=potsdam +units=m +no_defs+proj=tmerc +lat_0=0 +lon_0=12 +k=1 +x_0=4500000 +y_0=0 +ellps=bessel +datum=potsdam +units=m +no_defs"

# assign correct coordinate system to all data

scPDSI_stack2000@crs<-sp::CRS(defaultproj)
landklif_quadrants<-spTransform(landklif_quadrants,CRS(defaultproj))

## after all data are in the same coordinate system, combine raster and vector data
# first, subset scPDSI data to the extent to bavaria --> not necessary for Landklif quadrants

# get bavaria boundaries from the raster package

bnd <- raster::getData("GADM", country='DEU', level=1)
bnd <- spTransform(bnd, CRS(proj4string(scPDSI_stack2000)))
bnd_by <- bnd[bnd$NAME_1=="Bayern",]

plot(bnd_by)

# subset scPDSI stack

scPDSI_stack_by<-crop(scPDSI_stack2000,bnd_by)
scPDSI_stack_by<-mask(scPDSI_stack2000,bnd_by)

# extract scPDSI values at extent of different shapefiles
# note that data for quadrant @37 cannot be extracted here, because it is too close 
# to the border and partially reaching into the czech republic, where the scPDSI data is not available.
# It will be seen in the graphics as missing values later, constructing a buffer around the 
# bavarian border from which the data should be extracted did not solve the problem.

scPDSI_landklif_extracted<- raster::extract(scPDSI_stack_by, landklif_quadrants, df = T)

# group according to quadrant ID after mean

landklif_agg<-aggregate(scPDSI_landklif_extracted[,2:ncol(scPDSI_landklif_extracted)],list(scPDSI_landklif_extracted$ID),mean)
View(landklif_agg)

# invert in order to write all timesteps as one column

landklif_agg_inv<-landklif_agg %>% gather(timesteps,mean,scPDSI_raster_JW.109:scPDSI_raster_JW.348)
landklif_agg_inv$timesteps<-as.factor(landklif_agg_inv$timesteps)
View(landklif_agg_inv)

# save as csv

write.csv(landklif_agg,"landklif_agg.csv")
write.csv(landklif_agg_inv,"landklif_agg_inv.csv")

# loop to calculate statistics for every Quadrant over the timesteps

for(i in 1:nrow(landklif_agg)){
  landklif_quadrants@data[i,"timemean"]<- (sum(landklif_agg[i,2:241])/241)
  landklif_quadrants@data[i,"timemin"]<- min(landklif_agg[i,2:241])
  landklif_quadrants@data[i,"timemax"]<- max(landklif_agg[i,2:241])
  landklif_quadrants@data[i,"timesd"]<- sd(landklif_agg[i,2:241])
  if(i == 60){
    rm(i)
  }
}

# save back to the input shapefile for further analysis

writeOGR(landklif_quadrants, ".", layer = "landklif_quadrants_with_scPDSI_statistics", driver = "ESRI Shapefile", overwrite_layer = T)

## plot results for nice visualization

# first, we want to create an overview map with all plots over bavaria
# containing their scPDSI mean values over all timesteps

# reload shapefile as sf object

landklif_quadrants_sf<-st_read("landklif_quadrants_with_scPDSI_statistics.shp")

# save bavaria boundaries shapefile and reload as sf object for plotting

writeOGR(bnd_by, ".", layer = "bavariaboundaries", driver = "ESRI Shapefile")
bnd_by_sf<-st_read("bavariaboundaries.shp")


landklifmap<-ggplot()+
  geom_sf(data=bnd_by_sf,fill="white") +
  geom_sf(data=landklif_quadrants_sf,aes(fill=timemean,colour=landklif_quadrants_sf$Zones))+
  scale_fill_gradient(low="tomato", high="skyblue")+
  scale_color_discrete(l=40,c=350,name="climate zone(1-5) and LC class")+
  ggtitle("mean scPDSI values for Landklif plots 2000-2019")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_dark()

landklifmap


ggsave("landklifmap.png",landklifmap,width=210,height=160,unit="mm",device="png")

# secondly, we visualize scPDSI statistics of every single plot over time as a graph
# create sequence for all quadrants

QN<-seq(1,nrow(landklif_agg))

# for-loop to generate a dataframe for each Quadrant with the timesteps 
# to plot each Quadrant separately and saving them directly to the working directory

for(i in 1:length(QN)){               
  Quadrant<-landklif_agg_inv$mean[seq(from = QN[i], to = dim(landklif_agg_inv)[1], by = nrow(landklif_agg))]
  timesteps<-seq(as.Date("2000/1/1"), by = "month", length.out = 240)
  stat_df<-data.frame(Quadrant,timesteps)
  statplottime<-ggplot(stat_df)+
    geom_line(aes(x=timesteps,y=Quadrant,color="red"))+
    geom_smooth(aes(x=timesteps,y=Quadrant,color="blue"))+
    scale_color_discrete(name="scPDSI trend")+
    ggtitle("Development of scPDSI values over time",subtitle = paste0("Landklif-Quadrant: #", str_pad(QN[i],2,pad="0")))+
    ylab("mean scPDSI Value")+
    xlab("Time")+
    theme_bw()+
    theme(legend.position = "none")
  
  statplottime
  outname<-paste0(getwd(),"/scPDSI_timeseries_Quadrant",str_pad(QN[i],2,pad="0"),".png")
  ggsave(filename=outname,plot = statplottime,width=210,height=140,units="mm",device="png")
  
}

# third visualization: create a row of boxplots 
# representing the scPDSI value range of each Quadrant 

Quadrantboxplots<-ggplot(landklif_agg_inv)+
  geom_boxplot(aes(x=Group.1,y=mean,group=as.factor(Group.1),fill="red"))+
  ylab("scPDSI value range")+
  xlab("Quadrant")+
  ggtitle("Overview of scPDSI value ranges for all Landklif Quadrants")+
  theme(legend.position = "none")

ggsave("Quadrant_Boxplots.png",Quadrantboxplots,width=210,height=160,unit="mm",device="png")


