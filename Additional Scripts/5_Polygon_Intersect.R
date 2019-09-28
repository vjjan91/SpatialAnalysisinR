# Script that identifies the number of polygons that intersect with a particular polygon
# Author: Vijay Ramesh
# Date : June 2017

#### Please use st_intersects / st_intersection - Outdated script ####

# Load libraries
library(rgdal)
library(GISTools)
library(sp)
library(raster)

# List all shapefiles for every country
shapes <- list.files(path = "C:\\Users\\rameshv\\Desktop\\2_AllCountries\\", pattern=".shp$", full.names=TRUE)
head(shapes)

# Read in the big spatial polygons dataframe - this contains ALL MAMMALS
Mammals <-readOGR("C:\\Users\\rameshv\\Desktop\\Terrestrial_mammals.shp")

for (i in 1:length(shapes)){
  a <- readOGR(shapes[i])
  o <- over(a, Mammals,returnList = T)
  s <- o[[1]]
  count <- length(s$binomial)
  a$NAME <- as.character(a$NAME)
  results <- data.frame(a$NAME,count)
  write.table(results,file="C:\\Users\\rameshv\\Desktop\\Overlap_Countries.csv", row.names=F,
              append = TRUE, col.names = F, sep = ",")
}


