
# Basics of Spatial Analyses - Data types, Mapping, Slots

# Script 1 - SCCS-NY 2018
# Dt: 21st October 2018

# A useful package to start learning Spatial Visualization  and Data Types in R is GISTools
# For further spatial analysis - packages such as sp, sf, raster, rgdal, rgeos and maptools are useful
# Code adapted from Brunsdon and Comber 2015

install.packages('GISTools', dependencies = T)
library(GISTools)

data(newhaven) # Pre-existing dataset from GISTools library on New Haven Crime Stats along with demographic info, railway lines, place names
# http://www.newhavencrimelog.org
ls()

# 3 common data objects - Spatial Points, Lines and Polygons (lack attributes)
# If attributes are associated with above, they are dataframes

getClass('SpatialPoints')
getClass('SpatialPointsDataFrame')
getClass('SpatialPolygons')
getClass('SpatialPolygonsDataFrame')

head(blocks) # See what's in the blocks data
blocks@data$OCCUPIED # Visualizing data within slots

plot(blocks) # Plot blocks
plot(roads,add=T, col='red') # One can visualize the roads over the blocks

# Similar to the newhaven dataset, we can load the georgia dataset to play with different spatial objects
# Data comprises info on county outlines along with 1990 census info
data(georgia)
plot(georgia) # Outlines of county boundaries

georgia.outline <- gUnaryUnion(georgia, id=NULL) # Suppose I want to merge different spatial polygons
plot(georgia.outline, lwd=3)
plot(georgia, col='red', border='blue', add=T)
axis(1)
axis(2)
box()
grid() # Optional when visualizing particular types of spatial datasets

# Working with attributes and spatial objects
# Going back to the NewHaven dataset for a bit

hist(blocks@data$P_VACANT) # One can plot histograms to visualize patterns

# Exploratory analyses such as Kernel Density Estimates
plot(tracts)
plot(breach, add=T)
breach.dens <- kde.points(breach, lims = tracts) # Computes Density
level.plot(breach.dens) # Creates a level plot
masker <- poly.outer(breach.dens, tracts, extend = 100) # Analogous to Extract by Mask in ArcMap
add.masking(masker) # Draws the mask over the image currently drawn
plot(tracts, add=T)

summary(breach.dens) # Notice the dataframe has been converted to a SpatialPixels DataFrame

# Choropleth Map
# Thematic Map in which areas are shaded in proportion to the attributes present

?choropleth # You need a SpatialPolygons or a SpatialPolygonsDataFrame Object

choropleth(blocks, blocks@data$P_VACANT)
shades <- auto.shading(blocks@data$P_VACANT) # Adding n = 5/6/7 will decide intervals for the legend
# Below function is interactive and you can click on locations where you would like to plot points
locator() # Use locator to decide where you want your legend
choro.legend(525175, 161853, shades)
axis(1)
axis(2)
box()
map.scale(575410, 155719, miles2ft(2), "Miles",4, 0.5)
north.arrow(577427, 178907, miles2ft(0.2), col='lightblue', cex.lab = 0.8)
title("Add meaningful title here")

# ColorBrewer package has a large number of color gradients one can choose from (for instance, Green)
shades <- auto.shading(blocks@data$P_VACANT, cols = brewer.pal(5, "Greens"))
choropleth(blocks, blocks@data$P_VACANT, shading = shades)
locator()
choro.legend(525175, 161853, shades)
axis(1)
axis(2)
box()

# Reading in Vector data
# Using the Example data from Ramesh et al 2017 [Preferred format are shapefiles, csv's]

install.packages("rgdal")
library(rgdal)

# Set working directory to where the data is located
setwd("C:\\Users\\vr235\\Desktop\\Spatial Analysis in R - SCCS 2018\\Data\\")

# An alternative to the below command is using readShapePoly from maptools package
WG <- readOGR("WG.shp")
WG # One can see if there are any attributes associated with this vector data
plot(WG, col='lightpink')
axis(1); axis(2) ; box(); grid()

# Converting dataframe to spatial object
barfly_pr <- read.csv("barfly_pr.csv")
head(barfly_pr)
coordinates(barfly_pr)<- ~LONGITUDE+LATITUDE # Using the Coordinates command converts it to a Spatial Object
summary(barfly_pr)

# However, we need to set the coordinate system
WGCrs<- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
barfly_pr@proj4string <-WGCrs
barfly_pr

# Visualizing data for bird occurrences on the boundary drawn
plot(WG)
plot(barfly_pr, add=T, col='red', pch=1) 
axis(1); axis(2); box()

# Showing interactive plotting with vector data
install.packages("mapview")
library(mapview)
m1 <- mapview(WG)
m2 <- mapview(barfly_pr)
m1+m2

# Transforming coordinate system
# For the below functions, one has to ensure that the coordinate systems are the same across polygons and shapefiles

WG <- spTransform(WG, "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs")
proj4string(WG)
barfly_pr <- spTransform(barfly_pr, "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs")
proj4string(barfly_pr)

# Drawing a buffer
WG_buff <- gBuffer(WG, byid = T, width=10000)
plot(WG)
plot(WG_buff, add=T, border='red')

# Counting the number of spatial points within a polygon
barfly_count <- poly.counts(barfly_pr, WG)
barfly_count 

# Calculating areas of polygon features
WG_area <- poly.areas(WG)
WG_area <- WG_area/(1000*1000) # In Square Kilometres
WG_area

# Using the Intersect and gIntersection function to check if all points are in-fact occurring within the Western Ghats polygon

# Trying gIntersection

# Testing if every feature in x, is intersecting with every feature in y and the intersecting feature is returned
# See ?over if you want to only check LOGICAL T/F regarding if geometries intersect

bar_inter_1 <- gIntersection(barfly_pr, WG, byid = T) # Error appears that there is self-intersection of geometries
WG_inter <- gBuffer(WG, byid = T, width = 0) # Adding a zero-width buffer to solve this issue

bar_inter_1 <- gIntersection(barfly_pr, WG_inter, byid = T)
bar_inter_1 # A Spatial Points Object is returned

# Trying Intersect (from the Raster package) # RECOMMENDED OPTION
bar_inter_2 <- intersect(barfly_pr, WG)
head(bar_inter_2)
bar_inter_2 # A Spatial Points DataFrame is returned with attributes of y merged to this dataset

# The visualization was helpful since we know that we have clustering of points (before further analysis is carried out)
# (Think about using spThin / hand-thinning approaches to reduce extent of spatial clustering)

# How close are these points to one another ?

# Distance Calculation
library(raster)

?pointDistance # Neat function in the raster package to calculate distance between every set of points

bar_dist <- pointDistance(barfly_pr@coords, lonlat = F)
bar_dist <- (bar_dist/1000) # Converting it to kilometres
head(bar_dist)

# One can use functions from other packages such as spThin (Spatial Thinning to remove records within a particular distance)
library(spThin)

# thinned_bar <- thin(barfly_pr,lat.col = "LATITUDE", long.col = "LONGITUDE",spec.col= "Barfly", thin.par=2, reps=10, out.dir="")
# Note above command requires a dataset with species column, latitude column and longitude column

# Other functions
# Do all points fall within the Western Ghats boundary? 
gWithin(barfly_pr, WG_inter) # Logical T/F

# Conversions - Vector to Raster (Rasterize) and Raster to Polygon (Polygonize)

# Rasterize
WG_rast <- raster()
extent(WG_rast) <- extent(WG) + 1 # Adding 1 degree of extent (0.5 to each end)
r  <- rasterize(WG, WG_rast)

# Polygonize
# One can use rasterToPolygons (takes a lot of TIME, BEWARE!)
# Instead, I recommend using the efficient helper code written by John Baumgartner

source('C:\\Users\\vr235\\Desktop\\Spatial Analysis in R - SCCS 2018\\Polygonize.R')

# NOTE - you need to have GDAL or OSGeo4W installed on your system to run the above function
# Follow instructions from Albert Kochaphum - https://sandbox.idre.ucla.edu/sandbox/tutorials/installing-gdal-for-windows
# EDIT : Please go to ProgramFiles/GDAL/gdalplugins and remove ogrMSSQLSpatial.dll (not needed)
# Above link tells you how to install GDAL and the entire process for Windows

system.time(p <- polygonizer(r)) # Takes 7.53 seconds for the entire Western Ghats
p

# Using John's example for the entire planet to see how long it takes

library(rasterVis)
download.file('https://www.dropbox.com/s/tk3kg2oce4h2snd/NEcountries.asc.zip?dl=1',
              destfile={f <- tempfile()}, quiet=TRUE, cacheOK=FALSE,
              mode='wb')
unzip(f, exdir = d <- tempdir() )
t <- raster(file.path(d, 'NEcountries.asc'), crs='+proj=longlat')
levelplot(t, margin=FALSE, col.regions=rainbow)

system.time(s <- polygonizer(t)) # Takes < 20 seconds!
s

# Fancy plots using functions from the rasterVis library
spplot(s, col.regions=rainbow(200))



