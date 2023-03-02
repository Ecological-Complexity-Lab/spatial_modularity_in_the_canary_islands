##---load_libraries-------------------------------------------------------------------------------------------------
options(rgl.useNULL = TRUE)
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
#library(taxize)
#library(taxizedb)
#library(ggtree)
library(ggpubr)

#this portion of code turns the data into a multilayer network with 14 layers and does
#modularity analysis


##----get_data--------------------------------------------------------------------------------------------------------
#setwd("/Users/maya/Desktop/plant_pollinator_data/dryad_network")
#setwd("/Users/mayagoldstein/Desktop/project")
#getwd()

dryad_intralayer <- read.csv("./csvs/intralayer_file.csv")
#print(dryad_intralayer)
dryad_interlayer <- read.csv("./csvs/interlayer_file.csv") #already has inverted within
#print(dryad_interlayer)

##---- layers as islands and not sites -------------------------------------------------------------------------------
dryad_intralayer_islands <- dryad_intralayer

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

dryad_intralayer_islands$layer_from[dryad_intralayer_islands$layer_from %in% old_names] <- 
  new_names[match(dryad_intralayer_islands$layer_from, old_names)] #change to reflect layer = island

dryad_intralayer_islands$layer_to[dryad_intralayer_islands$layer_to %in% old_names] <- 
  new_names[match(dryad_intralayer_islands$layer_to, old_names)] #change to reflect layer = island

#if node_from, node_to, layer_from, layer_to are all the same need to sum the weight
dryad_intralayer_islands_grouped <- dryad_intralayer_islands %>% 
  group_by(layer_from, node_from, layer_to, node_to) %>% 
  summarise(sum_weight = sum(weight))

## ----dryad intralayer interlayer both ways-------------------------------------------------------------------------------------
intralayer_inverted <- tibble(values= dryad_intralayer_islands$layer_to, dryad_intralayer_islands$node_to, dryad_intralayer_islands$layer_from, 
                              dryad_intralayer_islands$node_from, dryad_intralayer_islands$weight) #create an inverted copy for directed intralayers
colnames(intralayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight")

## ---- normalize intralayer weights--------------------------------------------------------------------------
#plants in from
tot_plant <- dryad_intralayer_islands %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted <- dryad_intralayer_islands %>% left_join(tot_plant) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)


#pols in from
tot_pol <- intralayer_inverted %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted <- intralayer_inverted %>% left_join(tot_pol) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight) 

## ----interlayers as islands and not sites-----------------------------------------------------------------------
interlayers_new <- read.csv("./csvs/interlayers_new.csv")

interlayers_new_islands <- interlayers_new

interlayers_new_islands$layer_from[interlayers_new_islands$layer_from %in% old_names] <- 
  new_names[match(interlayers_new_islands$layer_from, old_names)] #change to reflect layer = island

interlayers_new_islands$layer_to[interlayers_new_islands$layer_to %in% old_names] <- 
  new_names[match(interlayers_new_islands$layer_to, old_names)] #change to reflect layer = island

interlayers_new_islands <- interlayers_new_islands %>% unique()

##---- create weight for interlayer edges --------------------------------------------------------
distances <- read.csv("./csvs/distances_file.csv")

distances$layer_from[distances$layer_from %in% old_names] <- 
  new_names[match(distances$layer_from, old_names)] #change to reflect layer = island

distances$layer_to[distances$layer_to %in% old_names] <- 
  new_names[match(distances$layer_to, old_names)] #change to reflect layer = island

distances <- distances %>% filter(layer_from != layer_to) %>% 
  group_by(layer_to, layer_from) %>% #add ave or group %>% unique()

shortest_distance <- min(distances$distance_in_meters)

interlayer_weight <- function(d){
  #recieves distance and normalizes it
  weight <- (1/log(d))/(1/log(shortest_distance))
  return(weight)
}

distances_with_weights <- distances %>% 
  mutate(weight = interlayer_weight(distance_in_meters)) %>% #add weight using the function
  subset(select = -distance_in_meters)

#write.csv(distances_with_weights, "./csvs/distances_with_weights.csv", row.names = FALSE)

interlayers_with_weights_islands <- interlayers_new_islands %>% inner_join(distances_with_weights, 
                                                           by = c("layer_from", "layer_to")) %>% unique()



## ----multilayer_extended_final--------------------------------------------------------------------------------------

edgelist_intralayers_both <- bind_rows(intralayer_weighted, intralayer_weighted_inverted) #combine weighted version of intra with inter

dryad_edgelist_complete <- bind_rows(edgelist_intralayers_both, dryad_interlayer) #combine inverted and non inverted verions

