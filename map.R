



## ---- Map showing structure similarity between islands
library(tmaptools)
library(rgdal)
library(sf)
library(raster)
library(OpenStreetMap)
library(maptiles)
library(sp)
library(glue)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(colorspace)
library(sf)

##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")






# Define the bounding box for the Canary Islands
canary_islands_bbox <- c(left = -18.542076, bottom = 26, right = -12.58351, top = 30.323990)

# Fetch OpenStreetMap data for the Canary Islands
canary_osm <- tm_basemap(canary_islands_bbox, server = "terrain")

tmap_mode("view")
sf_data <- as_tmap_sf(canary_osm)
sp_points <- SpatialPointsDataFrame(canary_osm, data = canary_osm)
str(canary_osm)
canary_osm = sf::st_as_sf(canary_osm) #make it Spatial format

tm_shape(canary_osm) +
  tm_rgb() +
  tm_borders() +
  tm_layout(frame = FALSE)


library(maps)
library(mapdata)
coord <- maps:::map.poly("worldHires", "Australia", exact=FALSE,
                         xlim=c(110,160),ylim=c(-45,-5),
                         boundary=FALSE,
                         interior=TRUE, fill=TRUE, as.polygon=TRUE)
coord.sp <- map2SpatialPolygons(coord)
par(mar=c(0,0,0,0), xaxs="i",yaxs="i")
plot(coord.sp, col="grey")

# Create a static map with OpenStreetMap as the basemap
tmap_mode("plot")

tm_basemap()

dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990)) 
data("dryad_location")


tm_shape(World) +
  tm_polygons("name")

tm_shape(dryad_location) +
  tm_raster("elevation", palette = terrain.colors(10)) +
  tm_shape(World) +
  tm_borders("white", lwd = .5) +
  tm_text("iso_a3", size = "AREA") +
  tm_shape(metro) +
  tm_symbols(col = "red", size = "pop2020", scale = .5) +
  tm_legend(show = FALSE)

# Add points to the map
my_map <- my_map +
  tm_symbols(size = 1, col = "red", legend.show = FALSE, data = YOUR_DATA_FRAME)

# Display the map
plot(my_map)


