
# Basics of Simple Features Package in R

# Script 1: SCCS-NY 2019
# Dt: October 5th 2019

# Trial using the sf package in R
library(sf)

# Using data provided in the data folder:
# Let's read in a .csv of points (Potential problems at this point?)

barfly <- read.csv("C:\\Users\\vr235\\Desktop\\SCCS2019\\SpatialAnalysisinR\\Data\\barfly_pr.csv",header=TRUE)
head(barfly)

# Convert it to a shapefile 
barfly <- st_as_sf(barfly, coords=c("LONGITUDE","LATITUDE"),crs=4326)
barfly # Let's take a look at it
class(barfly)

# A basic visualization
plot(barfly$geometry)
axis(1)
axis(2)
box()
grid()

# Let's try reading in a polygon
WG <- st_read("C:\\Users\\vr235\\Desktop\\SCCS2019\\SpatialAnalysisinR\\Data\\WG.shp")
class(WG)
st_crs(WG) <- 4326

# A quick reminder on the different types of functions one can use in sf
methods(class = "sf")

# Plotting
plot(WG[,1], reset = FALSE) # reset = FALSE: we want to add to a plot with a legend
plot(WG[2,1], col = 'grey', add = TRUE)

# Transforming coordinate systems
st_crs(WG)
WG_proj <- st_transform(WG, 32643)
WG_proj

# Let's read in an example from the sf package
nc <- st_read(system.file("shape/nc.shp", package="sf")) # Reading the North Carolina feature from the sf package
class(nc)
attr(nc,"sf_column")
print(nc[9:15],n=8) # Printing data from 9th to 15th column and first 8 features

# Plotting
plot(nc[1], reset=F)
plot(nc[1,1],col="grey",add=T) # Selecting a single feature of the entire set and greying it and adding it back to map

# Geometry and Geometry collection are important to distinguish
# Writing a shapefile:
st_write(nc, "nc.shp", delete_layer = T) # Delete_layer ensures that layers are overwritten 

# CRUD - Create, Read, Update and Delete options

# Coordinate Reference Systems and Transformations
# sfc objects have 2 attributes to store a CRS: epsg and proj4string

nc.web.mercator <- st_transform(nc, 3857)

# Calculate distance between geometries (If one wants distance between counties?)
dist <- st_distance(nc.web.mercator) # Provides a Matrix of distances between combinations of polygons

# Drawing Buffers around selected geomtery of choice
sel <- c(1,5,14)
geom <- st_geometry(nc.web.mercator[sel,])
buf <- st_buffer(geom, dist = 30000)
plot(buf, border = 'red')
plot(geom, add = TRUE)
plot(st_buffer(geom, -5000), add = TRUE, border = 'blue')

# Calculating areas
st_area(nc.web.mercator)

# A field can take 3 possibilities : constant, aggregate annd identity 
# Constant - Variable has a constant value over a spatial extent - Example - Land-use, Soil type, Climate
# Aggregate - Variable holds values that are summary values over the geometry - Example - Pop density
# Identity - Variable identifies the geometry

data(meuse, package="sp")
meuse_sf <- st_as_sf(meuse, coords = c("x","y"), crs=28992, agr = "constant")
st_agr(meuse_sf)

plot(meuse_sf[2])

### Tapping into tidyverse's capabilities for visualization
# Let's use ggplot2() on the nc shapefile 
# Code adapted from Matt-Strimas Mackey

library(ggplot2)

nc <- st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE)
ggplot(nc) +
  geom_sf(aes(fill = AREA)) +
  scale_fill_viridis_c("Area") +
  ggtitle("Area of counties in North Carolina") +
  theme_bw()

# Using the coord_sf() function to essentially plot the same using a different projection
ggplot(nc) +
  geom_sf(aes(fill = AREA)) +
  scale_fill_viridis_c("Area") +
  coord_sf(crs = st_crs(102003)) +
  ggtitle("Area of counties in North Carolina (Albers)") +
  theme_bw()

# Let's look at dplyr - the gold standard for data manipulation 
# Couple of functions that are very useful for manipulation of simple features:
# The following verbs operate only on the attribute data and leave the geometries untouched:
  
# select() keeps the specified variables, possibly renaming them
# rename() renames a variable and leaves all others unchanged
# filter() returns the rows that match the given conditions
# mutate() adds new variables based on existing variables
# transmute() creates new variables and drops existing variables
# arrange() sorts by the given variables
# slice() selects rows based on row number
# sample_n() samples n features randomly

# Suppose I want to know what counties intersect with other counties - Really useful tool
# Use square brackets to subset data, just the way you do it in base R

Ashe <- nc[nc$NAME=="Ashe",] # Gives all records that have the name=Ashe associated with it
nc[Ashe,] # This will give you all records that INTERSECT with Ashe 
nc[Ashe, op=st_touches] # Removing Self-intersection (Ashe with Ashe)

## Examples of using tidy ##
## Example 1
nc %>% 
  # salculate Area in km^2
  mutate(area_km2 = AREA * 10000) %>% 
  # select desired columns, note geometry column not explicitly selected
  select(name = NAME, area_km2) %>% 
  # filter to counties over 1,000 km^2
  filter(area_km2 > 2000) %>% 
  # arrange in descending order of area
  arrange(desc(area_km2)) %>% 
  # select first three rows
  slice(1:3)

## Example 2
## One could also use functions from sf within dplyr to calculate area (eg. st_area())
library(lwgeom)
nc %>% 
  mutate(area_m2 = st_area(geometry)) %>% 
  select(name = NAME, area_m2, area = AREA) %>% 
  head() %>% 
  as_tibble()

## Example 3
## Suppose I want to group data based on a given variable and then visualize the same
## For instance, I want to group information based on area

# Let's create a new arbitrary grouping variable and add it as a column to the nc dataset
nc_groups <- nc %>% 
  mutate(group = sample(LETTERS[1:3], nrow(.), replace = TRUE))
# Lets' average area by the groups created (So we can use this as a plotting parameter)
nc_mean_area <- nc_groups %>% 
  group_by(group) %>% 
  summarise(area_mean = mean(AREA))
# Plotting
ggplot(nc_mean_area) +
  geom_sf(aes(fill = area_mean)) +
  scale_fill_distiller("Area", palette = "Greens") +
  ggtitle("Mean area by group") +
  theme_bw()

## For more such examples, I recommend visiting: http://strimas.com/r/tidy-sf/

## Static Mapping Examples (Code borrowed from useR 2017 conference - Bhaskar K)

library(maps)
library(rnaturalearth) # Useful package for visualizations at large spatial scales

data(world.cities)

world.cities <- world.cities[world.cities$pop>1000000,]
world.cities <- st_as_sf(world.cities, coords = c("long","lat"), crs=4326)

world <- countries110
world <- world[world$name != 'Antarctica',]
grid.lines.mj <- gridlines(world,easts = seq(-180,180,by=30), norths = seq(-90,90,by=30))
grid.lines.mi <- gridlines(world,easts = seq(-165,195,by=15), norths = seq(-90,90,by=15))

# Using Base R visualization tools
par(mar = c(8, 0.1, 0.1, 0.1))
plot(methods::as(world, 'Spatial'), expandBB=c(0,0,0.05,0.05)) # Learned about this recently
plot(world, add=TRUE, border=grey(0.2))
plot(grid.lines.mi, col=grey(0.95), add=T)
plot(grid.lines.mj, col=grey(0.9), add=T)
text(labels(grid.lines.mj, side=1:2, col = grey(.6), offset=0.3))
plot(world.cities, add=TRUE, col='#FF5A0088', pch=20,
          cex=world.cities$pop/2000000)
v = c(1,4,8,12)
legend("topright", legend = v, pch = 20, pt.cex = v/2,
            text.col =grey(.7), box.col=grey(0.9), col = '#FF5A0088',
            title='Pop. (Millions)', horiz =T)

# It's not bad, but it's not perfect and can get a little clunky
## Question usually raised at this point: 
## How do I know what projection to use for my data. Fear not as proejctionWizard is here to help - http://projectionwizard.org/
## Also highly recommend using: https://spatialreference.org/

## Let's try the above using ggplot2()
library(rnaturalearthdata) # Wonderful dataset with boundaries, polygons 

world <- ne_countries(scale='medium',returnclass = 'sf')
class(world)

# Goal is to produce several maps - either side by side / or on a grid
# Below code was sourced from r-spatial.org

# Basics:
ggplot(data = world) +
  geom_sf() + theme_bw()

# Basics 2:
ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$name)), " countries)")) + theme_bw()

# Creating a map of the planet, with a rectangle around the gulf of mexico
(gworld <- ggplot(data = world) +
    geom_sf(aes(fill = region_wb)) +
    geom_rect(xmin = -102.15, xmax = -74.12, ymin = 7.65, ymax = 33.97,  # You need the coordinates for this (Unfortunately)
              fill = NA, colour = "black", size = 1.5) +
    scale_fill_viridis_d(option = "plasma") +
    theme(panel.background = element_rect(fill = "azure"),
          panel.border = element_rect(fill = NA)))

# Let's centre a map on the Gulf of Mexico
(ggulf <- ggplot(data = world) +
    geom_sf(aes(fill = region_wb)) +
    annotate(geom = "text", x = -90, y = 26, label = "Gulf of Mexico", 
             fontface = "italic", color = "grey22", size = 6) + 
    coord_sf(xlim = c(-102.15, -74.12), ylim = c(7.65, 33.97), expand = FALSE) + # Note: coord_sf() uses coordinates from the rectangle used earlier
    scale_fill_viridis_d(option = "plasma") +
    theme(legend.position = "none", axis.title.x = element_blank(), 
          axis.title.y = element_blank(), panel.background = element_rect(fill = "azure"), 
          panel.border = element_rect(fill = NA)))

# Two ways to arrange them in the same grid:
# Using ggplot()
ggplot() +
  coord_equal(xlim = c(0, 3.3), ylim = c(0, 1), expand = FALSE) +
  annotation_custom(ggplotGrob(gworld), xmin = 0, xmax = 1.5, ymin = 0, 
                    ymax = 1) +
  annotation_custom(ggplotGrob(ggulf), xmin = 1.5, xmax = 3, ymin = 0, 
                    ymax = 1) +
  theme_void()

# Much easier approach - using library(cowplot)
library(cowplot)
plot_grid(gworld, ggulf, nrow = 1, rel_widths = c(2.3, 1))

# To save the plots: 
# ggsave("grid.pdf", width = 15, height =  5)

## Some examples of interactive mapping
library(mapview)
m1 <- mapview(WG)
m1

m2 <- mapview(barfly)
m1+m2

# Edit map interactively:
library(mapedit)
a <- mapview(WG) %>% editMap('WG')
mapview(a$drawn)

### Finishing up with some basics of raster package
library(raster)

# Conversions - Vector to Raster (Rasterize) and Raster to Polygon (Polygonize)

# Rasterize
WG_rast <- raster()
extent(WG_rast) <- extent(WG) + 1 # Adding 1 degree of extent (0.5 to each end)
WG_m <- st_zm(WG)
r  <- rasterize(WG_m, WG_rast)


