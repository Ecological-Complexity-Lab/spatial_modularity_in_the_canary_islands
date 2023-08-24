# ---- NULL MODEL M1: SHUFFLING SPECIES BETWEEN ISLANDS ------------------------------------------------------------------------

# this portion of the code shuffles the networks between layers in one of 3 versions:
# 1. shuffling plants among themselves
# 2. shuffling pollinators among themselves
# 3. shuffling plants among themselves and then pollinators among themselves
# the shuffled networks are then compared the empirical network to determine whether 
# certain network properties are random or not

# ---- NULL MODEL M4: UNIFORM VALUE INTERLAYER LINKS BETWEEN ISLANDS ------------------------------------------------------------------------

# this portion of the code fixes the weight of all interlayer links between islands. The weight of all interlayer links to a uniform value 
#equal to the median of all the interlayer weights in the network. This null model maintains intralayer structure and the presence of interlayer links.


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
library(ggpubr)

##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")

dryad_intralayer <- read.csv("./csvs/intralayer_file.csv")
dryad_interlayer <- read.csv("./csvs/interlayer_file.csv") #already has inverted within

##---- layers as islands -------------------------------------------------------------------------------
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
  summarise(sum_weight = sum(weight)) #turn sums of sites to sum of island

#---- shuffle pollinators ---------------------------------------------------------------------------------------------
physical_nodes <- read.csv("./csvs/Islands/Jac/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs/Islands/Jac/layer_metadata_islands.csv")

intralayers_with_ids <- 
  dryad_intralayer_islands_grouped %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, sum_weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, sum_weight)

intralayers_with_ids$from <- paste(intralayers_with_ids$layer_from, "_", intralayers_with_ids$node_from) #combine layer and node from
intralayers_with_ids$to <- paste(intralayers_with_ids$layer_to, "_", intralayers_with_ids$node_to) #combine layer and node to 

intralayers_with_ids <- intralayers_with_ids %>% select(from, to, sum_weight) #select only new columns and weight

intralayer_matrix <- intralayers_with_ids %>%  #turn edge list to matrix
  select(from, to, sum_weight) %>%
  dcast(from ~ to, value.var = "sum_weight", fill = 0) %>%
  column_to_rownames(var="from")

intralayer_matrix <- t(intralayer_matrix) #change plants and pol- shuffling pols 

shuf_trial_matrix <- NULL
shuf_null_edge_list <- NULL

for (i in 1:1000){ #change from 500 to 1000 to shorten error bars
  trial <- intralayer_matrix
  print(i)
  random_num_list <- 1:nrow(trial) #list of all row numbers
  for (j in 1:3000){ 
    random_num_1 <- sample(random_num_list, 1, replace = FALSE) #pick 1st random number from list
    random_num_2 <- sample(random_num_list, 1, replace = FALSE) #pick 2nd random number from list
    random_species_1 <- str_split(rownames(trial)[random_num_1], " _ ")[[1]][2] #split the from to layer and node and only take node
    random_species_2 <- str_split(rownames(trial)[random_num_2], " _ ")[[1]][2] #split the to to layer and node and only take node
    random_layer_1 <- str_split(rownames(trial)[random_num_1], " _ ")[[1]][1] #split the from to layer and node and only take layer
    random_layer_2 <- str_split(rownames(trial)[random_num_2], " _ ")[[1]][1]#split the to to layer and node and only take layer
    
    wanted_value_layer_1 <- paste(random_layer_1 , "_" , random_species_2) #combination of layer 1 and node 2
    wanted_value_layer_2 <- paste(random_layer_2 , "_" , random_species_1) #combination of layer 2 and node 1
    
    if (random_species_1 == random_species_2) next #make sure we're not just switching the same species with itself
    if (wanted_value_layer_1 %in% rownames(trial)) next #make sure we don't already have the wanted combination in the layer
    if (wanted_value_layer_2 %in% rownames(trial)) next
    
    
    else{
      #if we dont have the combination:
      rownames(trial)[random_num_1] <- wanted_value_layer_1 #change the species to its new layer 
      rownames(trial)[random_num_2] <- wanted_value_layer_2 
    }
  }
  trial_with_shuf <- cbind(trial, i) #add trial number 
  shuf_trial_matrix <- rbind(shuf_trial_matrix, trial_with_shuf) #create big matrix of all matrices
  edge_list_version <- melt(as.matrix(trial)) %>% filter(value > 0) %>%
    select(from=Var1, to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
  edge_list_version$trial_number <- i #add trial number to the edge list
  shuf_null_edge_list <- rbind(shuf_null_edge_list, edge_list_version) #create mega edge list with all repetitions
}

shuf_null_edge_list$layer_from <- substr(shuf_null_edge_list$from, 1,2) #choose the layer portion of the mixed layer_node name
shuf_null_edge_list$node_from <- substr(shuf_null_edge_list$from, 5,8) #choose the node portion of the mixed layer_node name
shuf_null_edge_list$layer_to <- substr(shuf_null_edge_list$to, 1,2) #choose the layer portion of the mixed layer_node name
shuf_null_edge_list$node_to <- substr(shuf_null_edge_list$to, 5,8) #choose the node portion of the mixed layer_node name

shuf_null_edge_list <- shuf_null_edge_list %>% 
  select(layer_from, node_from, layer_to, node_to, weight, trial_number)

#write.csv(shuf_null_edge_list, "./csvs/Islands/Jac/shuf_null_edge_list_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix, "./csvs/Islands/Jac/shuf_trial_matrix_islands_as_layers.csv", row.names = TRUE)
#shuf_null_edge_list <- read.csv("./csvs/Islands/Jac/shuf_null_edge_list_islands_as_layers.csv")
#shuf_trial_matrix <- read.csv("./csvs/Islands/Jac/shuf_trial_matrix_islands_as_layers.csv")
