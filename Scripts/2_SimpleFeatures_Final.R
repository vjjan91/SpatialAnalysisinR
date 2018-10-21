
# Basics of Simple Features Package in R

# Script 2: SCCS-NY 2018
# Dt: October 21st 2018

# Trial using the sf package in R
library(spData)
library(spDataLarge)
library(sf)

# Using data provided in the sf vignette

nc <- st_read(system.file("shape/nc.shp", package="sf")) # Reading the NorthCarolina feature from the sf package
class(nc)
attr(nc,"sf_column")
print(nc[9:15],n=8) # Printing data from 9th to 15th column and first 8 features

methods(class="sf") # One can view all the different methods and operations you can carry out

# Viewing the geometry (either call the name of the column or use a function)
(nc_geom <- st_geometry(nc))

# Plotting
plot(nc[1])
plot(nc[1,1],col="grey",add=T) # Selecting a single feature of the entire set and greying it and adding it back to map

# Multipolygon structure in sf package is a list of lists of lists of matrices
nc_geom[[4]][[2]][[1]][1:3,] #4th feature, second exterior ring, first ring is always exterior and first 3 coordinate pairs

# Geometry columns have their own class
class(nc_geom)
methods(class="sfc")

# Geometry and Geometry collection are important to distinguish
# Writing a shapefile:
st_write(nc, "nc.shp", delete_layer = T) # Delete_layer ensures that layers are overwritten 

# CRUD - Create, Read, Update and Delete options

# Coordinate Reference Systems and Transformations
# sfc objects haev 2 attributes to store a CRS: epsg and proj4string

nc.web.mercator <- st_transform(nc, 3857)

# Calculate distance between geometries
dist <- st_distance(nc.web.mercator) # Provides a Matrix of distances between combinations of polygons

# Drawing Buffers around selected geomtery of choice
sel <- c(1,5,14)
geom <- st_geometry(nc.web.mercator[sel,])
buf <- st_buffer(geom, dist = 30000)
plot(buf, border = 'red')
plot(geom, add = TRUE)
plot(st_buffer(geom, -5000), add = TRUE, border = 'blue')

# Calculating areas
st_area(nc)

# A field can take 3 possibilities : constant, aggregate annd identity 
# Constant - variable has a constant value over a spatial extent - example - land-use, soil type, climate
# Aggregate - variable holds values that are summary values over the geometry - pop density
# Identity - variable identifies the geometry

data(meuse, package="sp")
meuse_sf <- st_as_sf(meuse, coords = c("x","y"), crs=28992, agr = "constant")
st_agr(meuse_sf)

plot(meuse_sf[2])

# Lines, multilines, multipolygons can all be included in a geometry collection

# Example for Coordinate Systems and Transformations
st_crs(nc)

# Transforming the Coordinate System
nc.web.mercator <- st_transform(nc, 3857)
st_crs(nc.web.mercator)

# Manipulating simple features
nc.transform <- st_transform(nc, 2264)
nc.transform[1,] # Selects the first record and prints the rest of the features

# Suppose I want to select a single column and the first two records, I can use the pipe operator
library(dplyr)
nc.transform %>% select(NWBIR74) %>% head(2)

# If I need to coerce the above to a dataframe and exclude the geometry and just need the information associated
nc.transform %>% as.data.frame %>% select(NWBIR74) %>% head(2)

# Suppose I want to know what counties intersect with other counties - Really useful tool
# Use square brackets to subset data, just the way you do it in base R

Ashe <- nc[nc$NAME=="Ashe",] # Gives all records that have the name=Ashe associated with it
nc[Ashe,] # This will give you all records that INTERSECT with Ashe 
nc[Ashe, op=st_touches] # Removing Self-intersection (Ashe with Ashe)

# The above can be done using dplyr as well
nc %>% filter(lengths(st_touches(.,Ashe))>0)


