
# Handling Raster Data

# Script 3: SCCS-NY 2018
# Dt: October 21st 2018

# Using WorldClim data for this workshop

library(raster)
# library(sf)
library(sp)
library(rgdal)
library(mapview)
library(spatialEco)

setwd("C:\\Users\\vr235\\Desktop\\Spatial Analysis in R - SCCS 2018\\Data\\")

# Multiple options of loading data
# 1. Downloading layers for the entire globe at a coarse resolution

r <- getData('worldclim',var='bio', res= 2.5) 
?getData # Explore getData in detail to see what datasets one can download

# 2. Downloading a tile based on the location of your interest
# Simply add in a lat/lon for a location you want to look at (you generally get a MUCH bigger tile encompassing your location)

env<-getData(name="worldclim", var="bio", res=0.5, download=T, lat=10, lon=75)
plot(env) # Visualizing the 19 bioclimatic layers 
env

# Below section is optional depending on the need to stitch multiple tiles (if your locations don't fall within a single tile)
# Merge function from the raster package can be used to stitch multiple tiles if needed
# Try below example if you would like to see how merge works
# Downloaded multiple tiles for a different region of the world 
env1<-getData(name="worldclim", var="bio", res=0.5, download=T, lat=36, lon=-78)
env2<-getData(name="worldclim", var="bio", res=0.5, download=T, lat=36, lon=-100)
env3<-getData(name="worldclim", var="bio", res=0.5, download=T, lat=20, lon=-78)
env4<-getData(name="worldclim", var="bio", res=0.5, download=T, lat=20, lon=-100)

# Merge the four tiles into one rasterstack called envall
# Or simply use the mosaic function in the raster package

env12<-merge(env1, env2)
enva123<-merge(env12, env3)
envall<-merge(env123, env4)

# Define the extent of your study region / region of interest
# Using readOGR from rgdal 

WG <- readOGR("WG.shp")

# Note order as xmin, xmax, ymin, ymax
# Suppose I want to define an extent for the area I want to crop

e <- extent(c(72.87866, 78.07781,8.114791,21.25695))

WG.crop <- crop(env, extent(WG)) #Replace WG with e if you want to define a particular extent
plot(WG.crop)

WG.mask <- mask(WG.crop, WG)
plot(WG.mask)

# For interactive viewing
mapview(WG.mask[[1]])

# Manipulating Rasters
WG.mask # This tells you information about the raster, class, extent, coordinate system, res, dimensions and values

WG.mask <- stack(WG.mask) # Manipulating stacked rasters seem easier than manipulating bricks

# Dividing certain raster layers by a value of 10 as temp is multiplied by 10 in the original data
WG1 <-  stack()
a <- c(1,2,5,6,7,8,9,10,11)
for (i in unique(a)){
  b <- WG.mask@layers[[i]]/10
  WG1 <- stack(WG1,b)
}

# Dividing certain raster layers by a value of 100
WG2 <- stack()
d <- c(3,4)
for(i in unique(d)){
  e <- WG.mask@layers[[i]]/100
  WG2 <- stack(WG2,e)
}

WG_new <- stack(WG1, WG2, WG.mask@layers[[12]], WG.mask@layers[[13]],
                WG.mask@layers[[14]],WG.mask@layers[[15]], WG.mask@layers[[16]],
                WG.mask@layers[[17]], WG.mask@layers[[18]], WG.mask@layers[[19]])

WG_new # The new 19 layers are ready

# Setting a tranverse mercator projection for further calculation in the future
# Setting a resolution of 1km 

UTMCrs <- CRS("+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs ")
WG_new <- projectRaster(WG_new, crs = UTMCrs, res= 1000)
WG_new

# Extracting values from a raster for a given set of points
# Loading data from Ramesh et al., (2017)

barfly_pr <- read.csv("barfly_pr.csv")
head(barfly_pr)
coordinates(barfly_pr)<- ~LONGITUDE+LATITUDE
WGCrs<- CRS("+init=epsg:4326")
barfly_pr@proj4string <-WGCrs
barfly_pr

# Transform the coordinate system to match that of the rasters
# Here, we use spTransform for transforming the coordinate system of vectors
barfly_pr <- spTransform(barfly_pr, "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs")

# Using the raster::extract to obtain values for rasters at given locations of bird occurrence
Pres <- extract(WG_new,barfly_pr)
Pres <- data.frame(coordinates(barfly_pr), Pres)
head(Pres) # I have a dataframe ready for further analysis

# If you need help with preparing land cover rasters, feel free to reach out. 
?writeRaster
