# Handling Raster Data

# Script 4: SCCSNY-19
# Dt: 5th October 2019

# Using WorldClim data for this workshop

library(raster)
library(sf)
library(sp)
library(rgdal)
library(mapview)
library(spatialEco)

# setwd("C:\\Users\\vr235\\Desktop\\Spatial Analysis in R - SCCS 2019\\Data\\")

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

# Use data loaded previously
barfly

# Using the raster::extract to obtain values for rasters at given locations of bird occurrence
Pres <- extract(WG_new,barfly)
Pres <- data.frame(coordinates(barfly_pr), Pres)
head(Pres) # I have a dataframe ready for further analysis

# Working with Land-Cover Datasets
# How do you prepare the data? Aggregate classes? Resample ?

# Load Land Cover Raster
lulc_2015 <- raster("C:\\Users\\vr235\\Desktop\\Spatial Analysis in R - SCCS 2019\\Data\\2015_lulc")
lulc_2015
lulc_2015@data@attributes # Gives you a list of all land cover classes and their respective count

# Load matrix to be used for reclassification
rec <- read.csv("C:\\Users\\vr235\\Desktop\\Spatial Analysis in R - SCCS 2019\\Data\\LandCover_ReclassifyMatrix_2015.csv")

#Reclassifying the Raster to convert 70 LC to 15 LC (Since I have too many classes)
rc <- as.matrix(data.frame(from=rec$V2, to=rec$To))

New_LULC_2015 <- reclassify(lulc_2015,rc) # Reclassify it based on the matrix

New_LULC_2015 <- ratify(New_LULC_2015) # Ensure R knows that it is a categorical variable
rat <- data.frame(
  ID= 1:15,
  LandCover = c("Evergreen","Deciduous","Mixed Forest",
                "Sholas","Degraded Forest","Scrubland",
                "Dry Grassland","Wet Grassland","Plantation",
                "Settlement","High Elevation Plantation","Coastal",
                "Barren Land","Cropland","Wetlands")
)

levels(New_LULC_2015)<- rat

# Below lines of CODE need to be carefully carried out depending on your system (URGES CAUTION!)

# 1. Identify what is the difference between the resolution you want to sample to and your current resolution
# ie. 30m to 1000m, for example (by what factor)
fact <- round(dim(New_LULC_2015)[1:2]/dim(newelev)[1:2]) # Here newelev is a layer I want to resample my land cover rasters to

# 2. Aggregate your Land Cover layer to coarser resolution
# Remember that finer to coarser (aggregate) and coarser to finer (disaggregate)
# Below I aggregate by a factor and use the fun=modal, taking the mode of cells within a factor window and assigning that value
# to the land cover layer
agg_modal <- aggregate(New_LULC_2015,fact, modal)

# 3. Resample function to resample it to the coarser scale and using a nearest neighbour method 
Resam_LC <- resample(agg_modal, newelev, method='ngb')