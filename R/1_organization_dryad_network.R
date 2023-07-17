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
library(taxize)
library(taxizedb)

#this portion of the code is responsible for adding the interlayer weights 

##----get data and organize--------------------------------------------------------------------------------------------------------
#intralayer edges
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
interactions_csv <- read.csv("./csvs/interactions.csv")

intercations <- NULL

intercations <- interactions_csv %>% group_by(interaction_id) %>% #group 1 interaction at a time
  mutate(layer_from = value[2], node_from = value[1], layer_to = value[4], 
         node_to = value[3], weight = value[5]) %>% #create new columns based on data
  subset(select = -c(interaction_id, attribute, value)) %>% unique() #delete previous columns

#write.csv(intercations, "./csvs/intralayer_file.csv", row.names = FALSE)

#distances
distances_csv <- read.csv("./csvs/Distance_between_sites_Dryad.csv")

#change distances to not include aggregation
distances <- distances_csv %>% dplyr::rename(layer_from = From, 
                                             layer_to = To, 
                                             distance_in_meters = DISTANCE..m.)


#write.csv(distances, "./csvs/distances_file.csv", row.names = FALSE)

##---- interlayer edges-------------------------------------------------------
interlayer_edges_from <- intercations %>% group_by(node_from) %>% 
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3],
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in


interlayer_edges_to <- intercations %>% group_by(node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges <- rbind(interlayer_edges_from, interlayer_edges_to)

#combine species and distances for interlayer
interlayers <- NULL

for (i in 1:nrow(interlayer_edges)){ #for every species
  current_run <- interlayer_edges[i,]
  for (j in 3:ncol(interlayer_edges)){ #for all locations
    current_species <- current_run$node_from
    current_location <- current_run$layer_from
    current_location_to <- current_run[j]
    interlayers <- rbind(interlayers, tibble(layer_from = current_location,
                                             node_from = current_species, 
                                             layer_to = as.character(current_location_to), 
                                             node_to = current_species))
    
  }
}

interlayers_new <- interlayers %>% subset(layer_to != "NA")
interlayers_new <- interlayers_new %>% subset(layer_from != layer_to)



##---- create weight for interlayer edges --------------------------------------------------------
shortest_distance <- min(distances$distance_in_meters)

interlayer_weight <- function(d){
  #receieves distance and normalizes it
  weight <- (1/log(d))/(1/log(shortest_distance))
  return(weight)
}

distances_with_weights <- distances %>% 
  mutate(weight = interlayer_weight(distance_in_meters)) %>% #add weight using the function
  subset(select = -distance_in_meters)

#write.csv(distances_with_weights, "./csvs/distances_with_weights.csv", row.names = FALSE)

interlayers_with_weights <- interlayers_new %>% inner_join(distances_with_weights, 
                                                           by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights <- interlayers_with_weights[!duplicated(interlayers_with_weights[c(2,4,5)]),]

#write.csv(interlayers_with_weights, "./csvs/interlayer_file.csv", row.names = FALSE)



