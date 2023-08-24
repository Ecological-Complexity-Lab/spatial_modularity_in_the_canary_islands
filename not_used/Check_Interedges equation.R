library(infomapecology)
library(igraph)
library(bipartite)
library(tidyverse)
library(magrittr)
library(betalink)
library(readxl)
library(ggalluvial)
library(scatterpie)
library(reshape2)
library(ggforce)
library(ggmap)
library(ggraph)
library(ggpubr)
library(reshape2)



##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")

distances <- read.csv("./csvs/distances_file.csv")

old_names <- c("WesternSahara1", "WesternSahara2",
               "Fuerteventura1", "Fuerteventura2",
               "GranCanaria1", "GranCanaria2",
               "TenerifeSouth1", "TenerifeSouth2",
               "TenerifeTeno1", "TenerifeTeno2",
               "Gomera1", "Gomera2",
               "Hierro1", "Hierro2")

new_names <- c("WesternSahara", "WesternSahara",
               "Fuerteventura", "Fuerteventura",
               "GranCanaria", "GranCanaria",
               "TenerifeSouth", "TenerifeSouth",
               "TenerifeTeno", "TenerifeTeno",
               "Gomera", "Gomera",
               "Hierro", "Hierro")

distances$layer_from[distances$layer_from %in% old_names] <- 
  new_names[match(distances$layer_from, old_names)] #change to reflect layer = island

distances$layer_to[distances$layer_to %in% old_names] <- 
  new_names[match(distances$layer_to, old_names)] #change to reflect layer = island

distances_normalized <- distances %>% filter(layer_from != layer_to) %>% #delete distances between sites in the same island
  group_by(layer_to, layer_from) %>% #group will contain 4 sites- site 1 and 1 of layer from and site 1 and 2 or layer to
  summarise(mean_distance = mean(distance_in_meters)) %>%unique() #use an average distance of the 4 sites in 2 different islands to determine the distance between the islands

distances_normalized <- distances_normalized[c("layer_from", "layer_to", "mean_distance")]


all_distances <-distances_normalized$mean_distance

min(all_distances)

logdist <- function(d_ab, d_min= 52751.49){
  x=(d_ab)
  y=(d_min)
  return(y/x)
}

for (i in all_distances){
  print(logdist(i))
}

