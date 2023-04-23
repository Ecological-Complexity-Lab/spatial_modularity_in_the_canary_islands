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
library(ggpubr)


# this portion of the code shuffles the networks between layers in one of 3 versions:
# 1. shuffling plants among themselves
# 2. shuffling pollinators among themselves
# 3. shuffling plants among themselves and then pollinators among themselves
# the shuffled networks are then compared the empirical network to determine whether 
# certain network properties are random or not

##---- null model rearrange matrices ----------------------------------------------------------------------------------
##---- shuffling pollinators ------------------------------------------------------------------------------------------
intralayers_with_ids <- 
  dryad_intralayer %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)

intralayers_with_ids$from <- paste(intralayers_with_ids$layer_from, "_", intralayers_with_ids$node_from) #combine layer and node from
intralayers_with_ids$to <- paste(intralayers_with_ids$layer_to, "_", intralayers_with_ids$node_to) #combine layer and node to 

intralayers_with_ids <- intralayers_with_ids %>% select(from, to, weight) #select only new columns and weight

#write.csv(intralayers_with_ids, "./csvs/intralayers_with_ids.csv")


intralayer_matrix <- intralayers_with_ids %>%  #turn edge list to matrix
  select(from, to, weight) %>%
  dcast(from ~ to, value.var = "weight", fill = 0) %>%
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

write.csv(shuf_null_edge_list, "./csvs/shuf_null_edge_list.csv", row.names = FALSE)
write.csv(shuf_trial_matrix, "./csvs/shuf_trial_matrix.csv", row.names = TRUE)
#shuf_null_edge_list <- read.csv("./csvs/shuf_null_edge_list.csv")
#shuf_trial_matrix <- read.csv("./csvs/shuf_trial_matrix.csv")

#here!!

# 1-39 to
# 40-288 from

#interlayer edges
interlayer_edges_from_shuf <- shuf_null_edge_list %>% group_by(trial_number, node_from) %>%
  select(trial_number, layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3], 
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in


interlayer_edges_to_shuf <- shuf_null_edge_list %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_shuf <- rbind(interlayer_edges_from_shuf, interlayer_edges_to_shuf) 


#write.csv(interlayer_edges_shuf, "./csvs/interlayer_edges_shuf.csv",row.names = FALSE)
#interlayer_edges_shuf <- read.csv("./csvs/interlayer_edges_shuf.csv")

##---- shuffling plants ----------------------------------------------------------------------------------------
intralayer_matrix_plants <- intralayer_matrix #change plants and pol- shuffling pols

shuf_trial_matrix_plants <- NULL
shuf_null_edge_list_plants <- NULL

for (i in 1:1000){ #change from 500 to 1000 to shorten error bars
  trial <- intralayer_matrix_plants
  print(i)
  random_num_list <- 1:ncol(trial) #list of all row numbers
  for (j in 1:3000){ 
    random_num_1 <- sample(random_num_list, 1, replace = FALSE) #pick 1st random number from list
    random_num_2 <- sample(random_num_list, 1, replace = FALSE) #pick 2nd random number from list
    random_species_1 <- str_split(colnames(trial)[random_num_1], " _ ")[[1]][2] #split the from to layer and node and only take node
    random_species_2 <- str_split(colnames(trial)[random_num_2], " _ ")[[1]][2] #split the to to layer and node and only take node
    random_layer_1 <- str_split(colnames(trial)[random_num_1], " _ ")[[1]][1] #split the from to layer and node and only take layer
    random_layer_2 <- str_split(colnames(trial)[random_num_2], " _ ")[[1]][1]#split the to to layer and node and only take layer
    
    wanted_value_layer_1 <- paste(random_layer_1 , "_" , random_species_2) #combination of layer 1 and node 2
    wanted_value_layer_2 <- paste(random_layer_2 , "_" , random_species_1) #combination of layer 2 and node 1
    
    if (random_species_1 == random_species_2) next #make sure we're not just switching the same species with itself
    if (wanted_value_layer_1 %in% colnames(trial)) next #make sure we dont already have the wanted combination in the layer
    if (wanted_value_layer_2 %in% colnames(trial)) next
    
    
    else{
      #if we dont have the combination:
      colnames(trial)[random_num_1] <- wanted_value_layer_1 #change the species to its new layer 
      colnames(trial)[random_num_2] <- wanted_value_layer_2 
    }
  }
  trial_with_shuf_plants <- cbind(trial, i) #add trial number 
  shuf_trial_matrix_plants <- rbind(shuf_trial_matrix_plants, trial_with_shuf_plants) #create big matrix of all matrices
  edge_list_version <- melt(as.matrix(trial)) %>% filter(value > 0) %>%
    select(from=Var1, to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
  edge_list_version$trial_number <- i #add trial number to the edge list
  shuf_null_edge_list_plants <- rbind(shuf_null_edge_list_plants, edge_list_version) #create mega edge list with all repetitions
}

shuf_null_edge_list_plants$layer_from <- substr(shuf_null_edge_list_plants$from, 1,2)
shuf_null_edge_list_plants$node_from <- substr(shuf_null_edge_list_plants$from, 5,8)
shuf_null_edge_list_plants$layer_to <- substr(shuf_null_edge_list_plants$to, 1,2)
shuf_null_edge_list_plants$node_to <- substr(shuf_null_edge_list_plants$to, 5,8)

shuf_null_edge_list_plants <- shuf_null_edge_list_plants %>% 
  select(layer_from, node_from, layer_to, node_to, weight, trial_number)

shuf_null_edge_list_plants$layer_from <- as.numeric(shuf_null_edge_list_plants$layer_from)
shuf_null_edge_list_plants$layer_to <- as.numeric(shuf_null_edge_list_plants$layer_to)
shuf_null_edge_list_plants$node_from <- as.numeric(shuf_null_edge_list_plants$node_from)
shuf_null_edge_list_plants$node_to <- as.numeric(shuf_null_edge_list_plants$node_to)

write.csv(shuf_null_edge_list_plants, "./csvs/shuf_null_edge_list_plants.csv", row.names = FALSE)
write.csv(shuf_trial_matrix_plants, "./csvs/shuf_trial_matrix_plants.csv", row.names = TRUE)
#shuf_null_edge_list_plants <- read.csv("./csvs/shuf_null_edge_list_plants.csv")
#shuf_trial_matrix_plants <- read.csv("./csvs/shuf_trial_matrix_plants.csv")

#interlayer edges
interlayer_edges_from_shuf_plants <- shuf_null_edge_list_plants %>% group_by(trial_number, node_from) %>% 
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3],
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in

interlayer_edges_to_shuf_plants <- shuf_null_edge_list_plants %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_shuf_plants <- rbind(interlayer_edges_from_shuf_plants, interlayer_edges_to_shuf_plants) 


#write.csv(interlayer_edges_shuf_plants, "./csvs/interlayer_edges_shuf_plants.csv",row.names = FALSE)
#interlayer_edges_shuf_plants <- read.csv("./csvs/interlayer_edges_shuf_plants.csv")

##---- shuffling both plants and pollinators -------------------------------------------------------------------
intralayer_matrix_both <- shuf_trial_matrix #copying and now we'll change columns and not rows


#test <- intralayer_matrix_both[intralayer_matrix_both[,"i"] == 500,]

shuf_trial_matrix_both <- NULL
shuf_null_edge_list_both <- NULL

for (i in 1:1000){ #change from 500 to 1000 to shorten error bars
  trial <- intralayer_matrix_both[intralayer_matrix_both[,"i"] == i,] #filter by trial number 1-1000
  print(i)
  random_num_list <- 1:(ncol(trial)-1) #list of all col numbers except for the "i" which is trial number
  for (j in 1:3000){ 
    random_num_1 <- sample(random_num_list, 1, replace = FALSE) #pick 1st random number from list
    random_num_2 <- sample(random_num_list, 1, replace = FALSE) #pick 2nd random number from list
    random_species_1 <- str_split(colnames(trial)[random_num_1], " _ ")[[1]][2] #split the from to layer and node and only take node
    random_species_2 <- str_split(colnames(trial)[random_num_2], " _ ")[[1]][2] #split the to to layer and node and only take node
    random_layer_1 <- str_split(colnames(trial)[random_num_1], " _ ")[[1]][1] #split the from to layer and node and only take layer
    random_layer_2 <- str_split(colnames(trial)[random_num_2], " _ ")[[1]][1]#split the to to layer and node and only take layer
    
    wanted_value_layer_1 <- paste(random_layer_1 , "_" , random_species_2) #combination of layer 1 and node 2
    wanted_value_layer_2 <- paste(random_layer_2 , "_" , random_species_1) #combination of layer 2 and node 1
    
    if (random_species_1 == random_species_2) next #make sure we're not just switching the same species with itself
    if (wanted_value_layer_1 %in% colnames(trial)) next #make sure we dont already have the wanted combination in the layer
    if (wanted_value_layer_2 %in% colnames(trial)) next
    
    
    else{
      #if we dont have the combination:
      colnames(trial)[random_num_1] <- wanted_value_layer_1 #change the species to its new layer 
      colnames(trial)[random_num_2] <- wanted_value_layer_2 
    }
  }
  trial_with_shuf_both <- cbind(trial, i) #add trial number 
  shuf_trial_matrix_both <- rbind(shuf_trial_matrix_both, trial_with_shuf_both) #create big matrix of all matrices
  edge_list_version <- melt(as.matrix(trial)) %>% filter(value > 0) %>%
    select(from=Var1, to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
  edge_list_version$trial_number <- i #add trial number to the edge list
  shuf_null_edge_list_both <- rbind(shuf_null_edge_list_both, edge_list_version) #create mega edge list with all repetitions
}


shuf_null_edge_list_both$layer_from <- substr(shuf_null_edge_list_both$from, 1,2)
shuf_null_edge_list_both$node_from <- substr(shuf_null_edge_list_both$from, 5,8)
shuf_null_edge_list_both$layer_to <- substr(shuf_null_edge_list_both$to, 1,2)
shuf_null_edge_list_both$node_to <- substr(shuf_null_edge_list_both$to, 5,8)

shuf_null_edge_list_both <- shuf_null_edge_list_both %>% 
  select(layer_from, node_from, layer_to, node_to, weight, trial_number)

shuf_null_edge_list_both <- subset(shuf_null_edge_list_both, layer_to != "i") #delete every row where "i" is found in layer_to
shuf_null_edge_list_both <- subset(shuf_null_edge_list_both, layer_from != "i") #delete every row where "i" is found in layer_from 
#write.csv(shuf_null_edge_list_both, "./csvs/shuf_null_edge_list_both.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_both, "./csvs/shuf_trial_matrix_both.csv", row.names = TRUE)


#interlayer edges
interlayer_edges_from_shuf_both <- shuf_null_edge_list_both %>% group_by(trial_number, node_from) %>% 
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3],
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in


interlayer_edges_to_shuf_both <- shuf_null_edge_list_both %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_shuf_both <- rbind(interlayer_edges_from_shuf_both, interlayer_edges_to_shuf_both) 

#write.csv(interlayer_edges_shuf_both, "./csvs/interlayer_edges_shuf_both.csv",row.names = FALSE)
#interlayer_edges_shuf_both <- read.csv("./csvs/interlayer_edges_shuf_both.csv")

##---- run oh HPC and get results for analysis -----------------------------------------------------------------
write.csv(interlayer_edges_shuf, "./HPC/shuf_between_layers/interlayer_edges_shuf.csv", row.names = FALSE) #create to run on HPC
write.csv(interlayer_edges_shuf_plants, "./HPC/shuf_between_layers/interlayer_edges_shuf_plants.csv", row.names = FALSE) #create to run on HPC
write.csv(interlayer_edges_shuf_both, "./HPC/shuf_between_layers/interlayer_edges_shuf_both.csv", row.names = FALSE) #create to run on HPC

#run on HPC and then come back with results
#each version has 3 code portions for the HPC:
# 1. HPC_network_x_shuffle.R (x being pollinators, plants or both)
# 2. 1_1000_x.sh (x being pollinators, plants or both)
# 3. i_x.sh (x being pollinators, plants or both)
# running 1_1000_x.sh manually in the cmd will make the other two run.
# all 1000 result csvs can be found in HPC/csvs_x (x being pollinators, plants or both) files

#both
files_both <- list.files("./HPC/shuf_between_layers/csvs_both/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_both <- read_csv(files_both) %>% bind_rows() #create a long edge list with all the csvs

#pollinators
files_pollinators <- list.files("./HPC/shuf_between_layers/csvs_pollinators/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_pol <- read_csv(files_pollinators) %>% bind_rows() #create a long edge list with all the csvs

#plants
files_plants <- list.files("./HPC/shuf_between_layers/csvs_plants/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_plants <- read_csv(files_plants) %>% bind_rows() #create a long edge list with all the csvs

#write.csv(my_merged_interlayer_shuf_both, "./csvs/my_merged_interlayer_shuf_both.csv", row.names = FALSE)
#write.csv(my_merged_interlayer_shuf_pol, "./csvs/my_merged_interlayer_shuf_pol.csv", row.names = FALSE)
#write.csv(my_merged_interlayer_shuf_plants, "./csvs/my_merged_interlayer_shuf_plants.csv", row.names = FALSE)


#---- interlayers with weights shuf version ------------------------------------------
distances_with_weights_ids <- distances_with_weights %>%
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, layer_to=layer_id.y, weight)

#write.csv(distances_with_weights_ids, "./csvs/distances_with_weights_ids.csv", row.names = FALSE)
#distances_with_weights_ids <- read.csv("./csvs/distances_with_weights_ids.csv")

#pols
interlayers_with_weights_shuf_pols <- my_merged_interlayer_shuf_pol %>% inner_join(distances_with_weights_ids, 
                                                                                   by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights_shuf_pols <- interlayers_with_weights_shuf_pols[!duplicated(interlayers_with_weights_shuf_pols[c(1,3,5,6)]),]

#plants
interlayers_with_weights_shuf_plants <- my_merged_interlayer_shuf_plants %>% inner_join(distances_with_weights_ids, 
                                                                                        by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights_shuf_plants <- interlayers_with_weights_shuf_plants[!duplicated(interlayers_with_weights_shuf_plants[c(1,3,5,6)]),]

#both
interlayers_with_weights_shuf_both <- my_merged_interlayer_shuf_both %>% inner_join(distances_with_weights_ids, 
                                                                                    by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights_shuf_both <- interlayers_with_weights_shuf_both[!duplicated(interlayers_with_weights_shuf_both[c(1,3,5,6)]),]

#write.csv(interlayers_with_weights_shuf_pols, "./csvs/interlayer_shuf_file_pols.csv", row.names = FALSE)
#write.csv(interlayers_with_weights_shuf_plants, "./csvs/interlayer_shuf_file_plants.csv", row.names = FALSE)
#write.csv(interlayers_with_weights_shuf_both, "./csvs/interlayer_shuf_file_both.csv", row.names = FALSE)


#create inter and intra for the 1000 shuf trials
dryad_interlayer_shuf_pols <- read.csv("./csvs/interlayer_shuf_file_pols.csv") #already has inverted
dryad_interlayer_shuf_plants <- read.csv("./csvs/interlayer_shuf_file_plants.csv") #already has inverted
dryad_interlayer_shuf_both <- read.csv("./csvs/interlayer_shuf_file_both.csv") #already has inverted

#----interlayer edges distribution------------------------------------------

dryad_interlayer_shuf_pols_dist <- dryad_interlayer_shuf_pols
dryad_interlayer_shuf_plants_dist <- dryad_interlayer_shuf_plants
dryad_interlayer_shuf_both_dist <- dryad_interlayer_shuf_both
inter_extended_dist <- inter_extended


#add type
dryad_interlayer_shuf_pols_dist$type <- "null_pols"
dryad_interlayer_shuf_plants_dist$type <- "null_plants"
dryad_interlayer_shuf_both_dist$type <- "null_both"
inter_extended_dist$type <- "empirical"
inter_extended_dist$trial_number <- NA
inter_extended_dist <- inter_extended_dist %>% select(trial_number, layer_from, node_from, layer_to, node_to, weight, type)

all_interlayers <- rbind(inter_extended_dist, dryad_interlayer_shuf_both_dist, 
                        dryad_interlayer_shuf_plants_dist, dryad_interlayer_shuf_pols_dist)

pdf('./graphs/shuffle_between_layers/interlayers_distribution_empirical_vs_shuf.pdf', 10, 6)
all_interlayers %>% ggplot(aes(x = weight, fill = type))+ geom_density(alpha = 0.4)+ theme_classic()+
  labs(x="Weight", y="Density")+ geom_vline(xintercept = 0.357602, linetype = "dashed", color = "tomato", size = 1)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()

#pols
dryad_intralayer_shuf_pols <- read.csv("csvs/shuf_null_edge_list.csv") 
dryad_intralayer_shuf_pols <- dryad_intralayer_shuf_pols[, c(6,1,2,3,4,5)]

#plants
dryad_intralayer_shuf_plants <- read.csv("csvs/shuf_null_edge_list_plants.csv") 
dryad_intralayer_shuf_plants <- dryad_intralayer_shuf_plants[, c(6,1,2,3,4,5)]

#both
dryad_intralayer_shuf_both <- read.csv("csvs/shuf_null_edge_list_both.csv") 
dryad_intralayer_shuf_both <- dryad_intralayer_shuf_both[, c(6,1,2,3,4,5)]

#view(dryad_interlayer_shuf)

#----create inverted versions--------------------------------------------
#pols
intralayer_inverted_shuf_pols <- tibble(values= dryad_intralayer_shuf_pols$layer_to, dryad_intralayer_shuf_pols$node_to, 
                                        dryad_intralayer_shuf_pols$layer_from, dryad_intralayer_shuf_pols$node_from, 
                                        dryad_intralayer_shuf_pols$weight, dryad_intralayer_shuf_pols$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_pols) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")


#plants
intralayer_inverted_shuf_plants <- tibble(values= dryad_intralayer_shuf_plants$layer_to, dryad_intralayer_shuf_plants$node_to, 
                                          dryad_intralayer_shuf_plants$layer_from, dryad_intralayer_shuf_plants$node_from, 
                                          dryad_intralayer_shuf_plants$weight, dryad_intralayer_shuf_plants$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_plants) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")

#both
intralayer_inverted_shuf_both <- tibble(values= dryad_intralayer_shuf_both$layer_to, dryad_intralayer_shuf_both$node_to, 
                                        dryad_intralayer_shuf_both$layer_from, dryad_intralayer_shuf_both$node_from, 
                                        dryad_intralayer_shuf_both$weight, dryad_intralayer_shuf_both$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_both) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")


## ----weighted intralayer edges--------------------------------------------------------------------------
##pols
#plants in from
tot_plant_shuf_pols <- intralayer_inverted_shuf_pols %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_pols <- intralayer_inverted_shuf_pols %>% left_join(tot_plant_shuf_pols) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_pols <- tot_plant_shuf_pols[, c(3,1,2,4)]

#pols in from
tot_pol_shuf_pols <- dryad_intralayer_shuf_pols %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_pols <- dryad_intralayer_shuf_pols %>% left_join(tot_pol_shuf_pols) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_pols <- tot_pol_shuf_pols[, c(3,1,2,4)]

##plants
#plants in from
tot_plant_shuf_plants <- intralayer_inverted_shuf_plants %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_plants <- intralayer_inverted_shuf_plants %>% left_join(tot_plant_shuf_plants) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_plants <- tot_plant_shuf_plants[, c(3,1,2,4)]

#pols in from
tot_pol_shuf_plants <- dryad_intralayer_shuf_plants %>%
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_plants <- dryad_intralayer_shuf_plants %>% left_join(tot_pol_shuf_plants) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_plants <- tot_pol_shuf_plants[, c(3,1,2,4)]

##both
#plants in from
tot_plant_shuf_both <- intralayer_inverted_shuf_both %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_both <- intralayer_inverted_shuf_both %>% left_join(tot_plant_shuf_both) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_both <- tot_plant_shuf_both[, c(3,1,2,4)]


#pols in from
tot_pol_shuf_both <- dryad_intralayer_shuf_both %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_both <- dryad_intralayer_shuf_both %>% left_join(tot_pol_shuf_both) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_both <- tot_pol_shuf_both[, c(3,1,2,4)]

## ----multilayer_extended_final--------------------------------------------------------------------------------------

edgelist_intralayer_shuf_pols <- bind_rows(intralayer_weighted_shuf_pols, intralayer_weighted_inverted_shuf_pols)
edgelist_intralayer_shuf_plants <- bind_rows(intralayer_weighted_shuf_plants, intralayer_weighted_inverted_shuf_plants)
edgelist_intralayer_shuf_both <- bind_rows(intralayer_weighted_shuf_both, intralayer_weighted_inverted_shuf_both)

dryad_edgelist_complete_shuf_pols <- bind_rows(edgelist_intralayer_shuf_pols, dryad_interlayer_shuf_pols) #combine inter and intra
dryad_edgelist_complete_shuf_plants <- bind_rows(edgelist_intralayer_shuf_plants, dryad_interlayer_shuf_plants) #combine inter and intra
dryad_edgelist_complete_shuf_both <- bind_rows(edgelist_intralayer_shuf_both, dryad_interlayer_shuf_both) #combine inter and intra


## ----distance decay in species shuf networks --------------------------------------------------------------------
#pols
all_species_all_layers_all_trials_pols <- rbind(tot_plant_shuf_pols, tot_pol_shuf_pols) %>% subset(select = -c(tot)) %>% 
  rename(node_id = node_from, layer_id = layer_from)

#plants
all_species_all_layers_all_trials_plants <- rbind(tot_plant_shuf_plants, tot_pol_shuf_plants) %>% subset(select = -c(tot)) %>% 
  rename(node_id = node_from, layer_id = layer_from)

#both
all_species_all_layers_all_trials_both <- rbind(tot_plant_shuf_both, tot_pol_shuf_both) %>% subset(select = -c(tot)) %>% 
  rename(node_id = node_from, layer_id = layer_from)

classic_layers_turnover_shuf_pols <- NULL
classic_layers_turnover_shuf_plants <- NULL
classic_layers_turnover_shuf_both <- NULL


#pols
classic_layers_turnover_shuf_output_pols <- distnace_decay_shuf(all_species_all_layers_all_trials_pols, 
                                                                classic_layers_turnover_shuf_pols) 

#write.csv(classic_layers_turnover_shuf_output_pols, "./csvs/classic_layers_turnover_shuf_output_pols.csv", row.names = FALSE)

#here!

#plants
classic_layers_turnover_shuf_output_plants <- distnace_decay_shuf(all_species_all_layers_all_trials_plants, 
                                                                  classic_layers_turnover_shuf_plants)

#write.csv(classic_layers_turnover_shuf_output_plants, "./csvs/classic_layers_turnover_shuf_output_plants.csv", row.names = FALSE)


#both
classic_layers_turnover_shuf_output_both <- distnace_decay_shuf(all_species_all_layers_all_trials_both, 
                                                                classic_layers_turnover_shuf_both)

#write.csv(classic_layers_turnover_shuf_output_both, "./csvs/classic_layers_turnover_shuf_output_both.csv", row.names = FALSE)


#pols
classic_layers_turnover_with_distances_shuf_pols <- right_join(classic_layers_turnover_shuf_output_pols, 
                                                               distances_with_ids, by= c("layer_from", "layer_to"))
classic_layers_turnover_with_distances_shuf_pols <- na.omit(classic_layers_turnover_with_distances_shuf_pols) #remove NA and delete layer name

#plants
classic_layers_turnover_with_distances_shuf_plants <- right_join(classic_layers_turnover_shuf_output_plants, 
                                                                 distances_with_ids, by= c("layer_from", "layer_to"))
classic_layers_turnover_with_distances_shuf_plants <- na.omit(classic_layers_turnover_with_distances_shuf_plants) #remove NA and delete layer name

#both
classic_layers_turnover_with_distances_shuf_both <- right_join(classic_layers_turnover_shuf_output_both, 
                                                               distances_with_ids, by= c("layer_from", "layer_to"))
classic_layers_turnover_with_distances_shuf_both <- na.omit(classic_layers_turnover_with_distances_shuf_both) #remove NA and delete layer name


#write.csv(classic_layers_turnover_with_distances_shuf_pols, "./csvs/classic_layers_turnover_with_distances_shuf_pols.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_with_distances_shuf_plants, "./csvs/classic_layers_turnover_with_distances_shuf_plants.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_with_distances_shuf_both, "./csvs/classic_layers_turnover_with_distances_shuf_both.csv", row.names = FALSE)
classic_layers_turnover_with_distances_shuf_pols <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_pols.csv")
classic_layers_turnover_with_distances_shuf_plants <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_plants.csv")
classic_layers_turnover_with_distances_shuf_both <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_both.csv")

#create an average for shuf with sd
#pols
ave_turnover_for_shuf_pols <- classic_layers_turnover_with_distances_shuf_pols %>% 
  group_by(layer_from, layer_to, distance_in_meters) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_pollinators") #create mean and sd for each point

#plants
ave_turnover_for_shuf_plants <- classic_layers_turnover_with_distances_shuf_plants %>% 
  group_by(layer_from, layer_to, distance_in_meters) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_plants") #create mean and sd for each point

#both
ave_turnover_for_shuf_both <- classic_layers_turnover_with_distances_shuf_both %>% 
  group_by(layer_from, layer_to, distance_in_meters) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_both") #create mean and sd for each point

#add the emprical classical turnover
empirical_turnover_for_shuf <- classic_layers_turnover_with_distances %>% 
  group_by(layer_from, layer_to, distance_in_meters) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#---- graphs for distance decay in species-----------------------------------------------------
turnover_shuf_and_empirical <- rbind(empirical_turnover_for_shuf, ave_turnover_for_shuf_pols, ave_turnover_for_shuf_plants,
                                     ave_turnover_for_shuf_both)

turnover_shuf_and_empirical <- turnover_shuf_and_empirical %>% mutate(distance_in_km = distance_in_meters/1000)

turnover_shuf_and_empirical %>% ggplot(aes(x= distance_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13),legend.text = element_text(size = 13))+
  labs(x="distance in km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400)+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.36, 0.34, 0.32, 0.30))

turnover_shuf_and_empirical %>% ggplot(aes(x= distance_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13),legend.text = element_text(size = 13))+
  labs(x="distance in km", y="Jaccard Similarity")

#-------making sure its significant---------------------------------------------

lm1 = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical,
                                     turnover_shuf_and_empirical$type=="empirical")) #in empirical
lm2 = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical,
                                     turnover_shuf_and_empirical$type=="null_pollinators")) #in null pols
lm3 = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical,
                                     turnover_shuf_and_empirical$type=="null_plants")) #in null plants
lm4 = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical,
                                     turnover_shuf_and_empirical$type=="null_both")) #in null pols


b1 <- summary(lm1)$coefficients[2,1] #coef for empirical
se1 <- summary(lm1)$coefficients[2,2] #sd for empirical

b2 <- summary(lm2)$coefficients[2,1] #coef for pols
se2 <- summary(lm2)$coefficients[2,2] #sd for pols

b3 <- summary(lm3)$coefficients[2,1] #coef for plants
se3 <- summary(lm3)$coefficients[2,2] #sd for plants

b4 <- summary(lm4)$coefficients[2,1] #coef for both
se4 <- summary(lm4)$coefficients[2,2] #sd for both

p_value_pol = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
p_value_pol

p_value_plants = 2*pnorm(-abs(compare.coeff(b1,se1,b3,se3)))
p_value_plants

p_value_both = 2*pnorm(-abs(compare.coeff(b1,se1,b4,se4)))
p_value_both

## ----multilayer_class-----------------------------------------------------------------------------------------------
# Input: An extended edge list.
dryad_edgelist_complete_shuf_pols <- dryad_edgelist_complete_shuf_pols[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_plants <- dryad_edgelist_complete_shuf_plants[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_both <- dryad_edgelist_complete_shuf_both[, c(1,2,3,4,6,5)] #make sure weight is #5

#write.csv(dryad_edgelist_complete_shuf_pols, "./csvs/dryad_edgelist_complete_shuf_pols.csv", row.names = FALSE)
#write.csv(dryad_edgelist_complete_shuf_plants, "./csvs/dryad_edgelist_complete_shuf_plants.csv", row.names = FALSE)
#write.csv(dryad_edgelist_complete_shuf_both, "./csvs/dryad_edgelist_complete_shuf_both.csv", row.names = FALSE)


dryad_multilayer_shuf_1000_pols <- NULL
dryad_multilayer_shuf_1000_plants <- NULL
dryad_multilayer_shuf_1000_both <- NULL


#pols
dryad_multilayer_shuf_1000_pols_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_pols, 
                                                              dryad_multilayer_shuf_1000_pols)

#plants
dryad_multilayer_shuf_1000_plants_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_plants, 
                                                                dryad_multilayer_shuf_1000_plants)


#both
dryad_multilayer_shuf_1000_both_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_both, 
                                                              dryad_multilayer_shuf_1000_both)

#write.csv(dryad_multilayer_shuf_1000_pols_output, "./csvs/dryad_multilayer_shuf_1000_pols_output.csv", row.names = FALSE)
#write.csv(dryad_multilayer_shuf_1000_plants_output, "./csvs/dryad_multilayer_shuf_1000_plants_output.csv", row.names = FALSE)
#write.csv(dryad_multilayer_shuf_1000_both_output, "./csvs/dryad_multilayer_shuf_1000_both_output.csv", row.names = FALSE)

dryad_multilayer_shuf_1000_pols_output <- read.csv("./csvs/dryad_multilayer_shuf_1000_pols_output.csv")
dryad_multilayer_shuf_1000_plants_output <- read.csv("./csvs/dryad_multilayer_shuf_1000_plants_output.csv")
dryad_multilayer_shuf_1000_both_output <- read.csv("./csvs/dryad_multilayer_shuf_1000_both_output.csv")

##---- distance decay of modules in islands shuf vs empirical----------------------------------------------------
#all_edge_list_island_combine_no_module_shuf_pols <- NULL
#all_edge_list_island_combine_no_module_shuf_plants <- NULL
#all_edge_list_island_combine_no_module_shuf_both <- NULL
#all_edge_list_layers_combine_no_module_shuf <- NULL


#islands_turnover_with_distnace_pols <- NULL
#islands_turnover_with_distnace_plants <- NULL
#islands_turnover_with_distnace_both <- NULL

#islands_turnover_with_distnaces_pols <- NULL
#islands_turnover_with_distnaces_plants <- NULL
#islands_turnover_with_distnaces_both <- NULL


#module_island_turnover_shuf <- NULL


##---- distance decay of modules in layers and ---------------------------------------------------------------
layers_turnover_with_distnace_pols <- NULL
layers_turnover_with_distnace_plants <- NULL
layers_turnover_with_distnace_both <- NULL

module_layer_turnover_shuf <- NULL


#---- outputs islands----------------------------------------------------------------------------
#pols
#all_edge_list_island_combine_no_module_shuf_pols_output <- module_distance_decay_func(dryad_multilayer_shuf_1000_pols_output,
#                                                                                      islands_turnover_with_distnace_pols,
#                                                                                      islands_turnover_with_distnaces_pols)


#write.csv(all_edge_list_island_combine_no_module_shuf_pols_output, 
#          "./csvs/all_edge_list_island_combine_no_module_shuf_pols_output.csv", 
#          row.names = FALSE)

#all_edge_list_island_combine_no_module_shuf_pols_output <- read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_pols_output.csv")

#plants
#all_edge_list_island_combine_no_module_shuf_plants_output <- module_distance_decay_func(dryad_multilayer_shuf_1000_plants_output,
#                                                                                        islands_turnover_with_distnace_plants,
#                                                                                        islands_turnover_with_distnaces_plants)

#write.csv(all_edge_list_island_combine_no_module_shuf_plants_output, 
#          "./csvs/all_edge_list_island_combine_no_module_shuf_plants_output.csv", 
#          row.names = FALSE)

#all_edge_list_island_combine_no_module_shuf_plants_output <- read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_plants_output.csv")

#both
#all_edge_list_island_combine_no_module_shuf_both_output <- module_distance_decay_func(dryad_multilayer_shuf_1000_both_output,
#                                                                                      islands_turnover_with_distnace_both,
#                                                                                      islands_turnover_with_distnaces_both)

#write.csv(all_edge_list_island_combine_no_module_shuf_both_output, 
#          "./csvs/all_edge_list_island_combine_no_module_shuf_both_output.csv", 
#          row.names = FALSE)

#all_edge_list_island_combine_no_module_shuf_both_output <- read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_both_output.csv")

#---- outputs layers--------------------------------------------------------------------------------------

all_edge_list_layer_combine_no_module_shuf_pols_output <- module_distance_decay_layer_func(dryad_multilayer_shuf_1000_pols_output,
                                                                                           layers_turnover_with_distnace_pols)


write.csv(all_edge_list_layer_combine_no_module_shuf_pols_output, 
          "./csvs/all_edge_list_layer_combine_no_module_shuf_pols_output.csv", 
          row.names = FALSE)

all_edge_list_layer_combine_no_module_shuf_pols_output <- read.csv("./csvs/all_edge_list_layer_combine_no_module_shuf_pols_output.csv")

#plants
all_edge_list_layer_combine_no_module_shuf_plants_output <- module_distance_decay_layer_func(dryad_multilayer_shuf_1000_plants_output,
                                                                                             layers_turnover_with_distnace_plants)

write.csv(all_edge_list_layer_combine_no_module_shuf_plants_output, 
          "./csvs/all_edge_list_layer_combine_no_module_shuf_plants_output.csv", 
          row.names = FALSE)

all_edge_list_layer_combine_no_module_shuf_plants_output <- read.csv("./csvs/all_edge_list_layer_combine_no_module_shuf_plants_output.csv")

#both
all_edge_list_layer_combine_no_module_shuf_both_output <- module_distance_decay_layer_func(dryad_multilayer_shuf_1000_both_output,
                                                                                           layers_turnover_with_distnace_both)

write.csv(all_edge_list_layer_combine_no_module_shuf_both_output, 
          "./csvs/all_edge_list_layer_combine_no_module_shuf_both_output.csv", 
          row.names = FALSE)

all_edge_list_layer_combine_no_module_shuf_both_output <- read.csv("./csvs/all_edge_list_layer_combine_no_module_shuf_both_output.csv")

##---------------------------------------------------------------------------------------------------------------------------------------
## number of modules
#pols
#pols_number_of_modules_output <- number_of_shared_modules_func(dryad_multilayer_shuf_1000_pols_output,
#                                                               all_edge_list_island_combine_no_module_shuf_pols)

#write.csv(pols_number_of_modules_output, "./csvs/pols_number_of_modules_output.csv", row.names = FALSE)
#pols_number_of_modules_output <- read.csv("./csvs/pols_number_of_modules_output.csv")

#plants
#plants_number_of_modules_output <- number_of_shared_modules_func(dryad_multilayer_shuf_1000_plants_output,
#                                                               all_edge_list_island_combine_no_module_shuf_plants)

#write.csv(plants_number_of_modules_output, "./csvs/plants_number_of_modules_output.csv", row.names = FALSE)
#plants_number_of_modules_output <- read.csv("./csvs/plants_number_of_modules_output.csv")


#both
#both_number_of_modules_output <- number_of_shared_modules_func(dryad_multilayer_shuf_1000_both_output,
#                                                               all_edge_list_island_combine_no_module_shuf_both)

#write.csv(both_number_of_modules_output, "./csvs/both_number_of_modules_output.csv", row.names = FALSE)
#both_number_of_modules_output <- read.csv("./csvs/both_number_of_modules_output.csv")


## ---------------------------------------------------------------------------------------------------------------------------------

#create an average for shuf with sd
#pols
#ave_turnover_for_module_shuf_pols <- pols_number_of_modules_output %>% 
#  group_by(layer_from, layer_to) %>%
#  summarise(ave=mean(number_of_modules), sd=sd(number_of_modules), ave_dist=mean(ave_distance)) %>% mutate(type="null_pollinators") #create mean and sd for each point


#plants
#ave_turnover_for_module_shuf_plants <- plants_number_of_modules_output %>% 
#  group_by(layer_from, layer_to) %>%
#  summarise(ave=mean(number_of_modules), sd=sd(number_of_modules), ave_dist=mean(ave_distance)) %>% mutate(type="null_plants") #create mean and sd for each point

#both
#ave_turnover_for_module_shuf_both <- both_number_of_modules_output %>% 
#  group_by(layer_from, layer_to) %>%
#  summarise(ave=mean(number_of_modules), sd=sd(number_of_modules), ave_dist=mean(ave_distance)) %>% mutate(type="null_both") #create mean and sd for each point


#add the emprical classical turnover
#empirical_turnover_for_module_shuf <- edge_list_island_combine_no_module %>% group_by(layer_from, layer_to) %>%
#  summarise(ave=mean(number_of_modules), sd=sd(number_of_modules), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#module_turnover_empirical_and_null <- rbind(empirical_turnover_for_module_shuf, ave_turnover_for_module_shuf_pols,
#                                            ave_turnover_for_module_shuf_plants, ave_turnover_for_module_shuf_both)

#module_turnover_empirical_and_null <- module_turnover_empirical_and_null %>% mutate(ave_dist_in_km = ave_dist/1000)

#---- create ave for jaccard islands -------------------------------------------------------------------------
#pols
#ave_module_island_turnover_shuf_pols <- all_edge_list_island_combine_no_module_shuf_pols_output %>% 
#  group_by(layer_from, layer_to) %>%
#  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_pollinators") #create mean and sd for each point

#plants
#ave_module_island_turnover_shuf_plants <- all_edge_list_island_combine_no_module_shuf_plants_output %>% 
#  group_by(layer_from, layer_to) %>%
#  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_plants") #create mean and sd for each point

#both
#ave_module_island_turnover_shuf_both <- all_edge_list_island_combine_no_module_shuf_both_output %>% 
#  group_by(layer_from, layer_to) %>%
#  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_both") #create mean and sd for each point

#add the empirical empirical
#empirical_turnover_for_module_island_shuf <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
#  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#empirical_turnover_for_module_island_shuf_no_self_loop <- empirical_turnover_for_module_island_shuf %>% subset(layer_from != layer_to) #for empirical only distance decay graph

#empirical_turnover_for_module_island_shuf_no_self_loop_km <- empirical_turnover_for_module_island_shuf_no_self_loop %>%
#  mutate(ave_dist_in_km = ave_dist/1000)

#write.csv(empirical_turnover_for_module_island_shuf_no_self_loop_km,
#          "./csvs/empirical_turnover_for_module_island_shuf_no_self_loop_km.csv", row.names = FALSE)

#---- create ave for jaccard layers ----------------------------------------------------------------------------------------
#pols
ave_module_layer_turnover_shuf_pols <- all_edge_list_layer_combine_no_module_shuf_pols_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_pollinators") #create mean and sd for each point

#plants
ave_module_layer_turnover_shuf_plants <- all_edge_list_layer_combine_no_module_shuf_plants_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_plants") #create mean and sd for each point

#both
ave_module_layer_turnover_shuf_both <- all_edge_list_layer_combine_no_module_shuf_both_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_both") #create mean and sd for each point

#add the empirical empirical
empirical_turnover_for_module_layer_shuf <- layers_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null


empirical_turnover_for_module_layer_shuf_no_self_loop <- empirical_turnover_for_module_layer_shuf %>% subset(layer_from != layer_to) #for empirical only distance decay graph

empirical_turnover_for_module_layer_shuf_no_self_loop_km <- empirical_turnover_for_module_layer_shuf_no_self_loop %>%
  mutate(ave_dist_in_km = ave_dist/1000)

#---- combine for jaccard analysis------------------------------------------------------------------------------------------------------
#combine all islands
#jaccard_similarity_empirical_and_null <- rbind(empirical_turnover_for_module_island_shuf, ave_module_island_turnover_shuf_pols,
#                                               ave_module_island_turnover_shuf_plants, ave_module_island_turnover_shuf_both)

#jaccard_similarity_empirical_and_null_no_self_loop <- jaccard_similarity_empirical_and_null %>% subset(layer_from != layer_to)

#jaccard_similarity_empirical_and_null_no_self_loop_km <- jaccard_similarity_empirical_and_null_no_self_loop %>% 
#  mutate(ave_dist_in_km = ave_dist/1000)

#combine all layers
jaccard_similarity_layer_empirical_and_null <- rbind(empirical_turnover_for_module_layer_shuf, ave_module_layer_turnover_shuf_pols,
                                                     ave_module_layer_turnover_shuf_plants, ave_module_layer_turnover_shuf_both)

jaccard_similarity_layer_empirical_and_null_no_self_loop <- jaccard_similarity_layer_empirical_and_null %>% subset(layer_from != layer_to)

jaccard_similarity_layer_empirical_and_null_no_self_loop_km <- jaccard_similarity_layer_empirical_and_null_no_self_loop %>% 
  mutate(ave_dist_in_km = ave_dist/1000)


#---- graphs for distance decay in modules shuf vs empirical--------------------------------
#island
#just emprical
empirical_turnover_for_module_island_shuf_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point(color = "#F47069")+ theme_classic()+ geom_smooth(color = "#F47069", method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400, color = "#F47069")+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = 0.61, color = "#F47069")

#empirical and null
jaccard_similarity_empirical_and_null_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400)+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.63, 0.60, 0.57, 0.54))


#version without trendline
jaccard_similarity_empirical_and_null_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ 
  geom_smooth(data = empirical_turnover_for_module_island_shuf_no_self_loop_km, method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")


#layer
#just emprical
empirical_turnover_for_module_layer_shuf_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="distance in km", y="Jaccard Similarity")

#empirical and null
jaccard_similarity_layer_empirical_and_null_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="distance in km", y="Jaccard Similarity")

# number of modules in common as func of distance between islands
module_turnover_empirical_and_null %>% ggplot(aes(x= ave_dist, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="similarity")


#------check if its significant for islands----------------------------------------------------------------
lm1_module = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop,
                                            jaccard_similarity_empirical_and_null_no_self_loop$type=="empirical")) #in empirical
lm2_module = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop,
                                            jaccard_similarity_empirical_and_null_no_self_loop$type=="null_pollinators")) #in null pols
lm3_module = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop,
                                            jaccard_similarity_empirical_and_null_no_self_loop$type=="null_plants")) #in null plants
lm4_module = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop,
                                            jaccard_similarity_empirical_and_null_no_self_loop$type=="null_both")) #in null both

#get equations
lm1_module_equation <- paste("y=", coef(lm1_module)[[1]], "+", coef(lm1_module)[[2]], "*x")
lm2_module_equation <- paste("y=", coef(lm2_module)[[1]], "+", coef(lm2_module)[[2]], "*x")
lm3_module_equation <- paste("y=", coef(lm3_module)[[1]], "+", coef(lm3_module)[[2]], "*x")
lm4_module_equation <- paste("y=", coef(lm4_module)[[1]], "+", coef(lm4_module)[[2]], "*x")

b1_module <- summary(lm1_module)$coefficients[2,1]
se1_module <- summary(lm1_module)$coefficients[2,2]
b2_module <- summary(lm2_module)$coefficients[2,1]
se2_module <- summary(lm2_module)$coefficients[2,2]
b3_module <- summary(lm3_module)$coefficients[2,1]
se3_module <- summary(lm3_module)$coefficients[2,2]
b4_module <- summary(lm4_module)$coefficients[2,1]
se4_module <- summary(lm4_module)$coefficients[2,2]

p_value_module_pols = 2*pnorm(-abs(compare.coeff(b1_module,se1_module,b2_module,se2_module)))
p_value_module_pols
p_value_module_plants = 2*pnorm(-abs(compare.coeff(b1_module,se1_module,b3_module,se3_module)))
p_value_module_plants
p_value_module_both = 2*pnorm(-abs(compare.coeff(b1_module,se1_module,b4_module,se4_module)))
p_value_module_both

##---- number of modules against shuffled networks -------------------------------------------------------
num_of_modules_pols <- NULL
num_of_modules_plants <- NULL
num_of_modules_both <- NULL
Num_of_modules <- NULL

num_of_modules_empirical <- modules_dryad_multilayer$modules %>% slice(which.max(module)) %>% select(module) #number of modules in the empirical network

Num_of_modules <- rbind(Num_of_modules, tibble(trial_num = NA ,module = num_of_modules_empirical, type = "empirical")) #add the number of modules for empirical
threshold <- unlist(Num_of_modules$module) #number of modules in empirical network

#pols
num_of_modules_pols <- dryad_multilayer_shuf_1000_pols_output %>% group_by(trial_num) %>% slice(which.max(module)) %>% select(module) #num of modules in each iteration of pols
num_of_modules_pols$type <- "null_pollinators"

#plants
num_of_modules_plants <- dryad_multilayer_shuf_1000_plants_output %>% group_by(trial_num) %>% slice(which.max(module)) %>% select(module) #num of modules in each iteration of plants
num_of_modules_plants$type <- "null_plants"

#both
num_of_modules_both <- dryad_multilayer_shuf_1000_both_output %>% group_by(trial_num) %>% slice(which.max(module)) %>% select(module) #num of modules in each iteration of both
num_of_modules_both$type <- "null_both"

Num_of_modules <- rbind(Num_of_modules, num_of_modules_pols, num_of_modules_plants, num_of_modules_both)
Num_of_modules$module <- unlist(Num_of_modules$module) #make sure that module column is not a list

#create graph for number of modules
Num_of_modules %>%
  ggplot(aes(x = module, fill = type))+
  geom_histogram(color = "black", alpha = 0.6, position = "identity")+
  geom_vline(aes(xintercept=threshold), color="red")+
  theme_classic()+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  geom_text(y=400, x=56, label="empirical network", color= "red", size=5, family= "Tahoma")+
  labs(x="number of modules", y="number of iterations")

#---- NMI against shuffled netwotworks----------------------------------------------------------------------------------------------------------------

I_es_pols <- NULL #I between empirical and pols
I_es_plants <- NULL #I between empirical and plants
I_es_both <- NULL #I between empirical and both

for(i in 1:1000){
  print(i)
  NMI_empirical <- modules_dryad_multilayer$modules %>% select(module, node_id, layer_id)  #NMI empirical
  NMI_pols <- dryad_multilayer_shuf_1000_pols_output %>% filter(trial_num == i) %>% select(module, node_id, layer_id) #pols
  NMI_plants <- dryad_multilayer_shuf_1000_plants_output %>% filter(trial_num == i) %>% select(module, node_id, layer_id) #plants
  NMI_both <- dryad_multilayer_shuf_1000_both_output %>% filter(trial_num == i) %>% select(module, node_id, layer_id) #both
  
  #pols    
  J_pols <- inner_join(NMI_empirical,NMI_pols,by=c('node_id','layer_id')) %>% 
    group_by(module.y) %>% 
    select(module.x) %>% table()
  I_es_pols <- rbind(I_es_pols, tibble(i=i, mutual_information=NMI(J_pols)))
  
  #plants
  J_plants <- inner_join(NMI_empirical,NMI_plants,by=c('node_id','layer_id')) %>% 
    group_by(module.y) %>% 
    select(module.x) %>% table()
  I_es_plants <- rbind(I_es_plants, tibble(i=i, mutual_information=NMI(J_plants)))
  
  #both
  J_both <- inner_join(NMI_empirical,NMI_both,by=c('node_id','layer_id')) %>% 
    group_by(module.y) %>% 
    select(module.x) %>% table()
  I_es_both <- rbind(I_es_both, tibble(i=i, mutual_information=NMI(J_both)))
}

I_es_pols$type <- "null_pollinators"
I_es_plants$type <- "null_plants"
I_es_both$type <- "null_both"
I_es_all <- rbind(I_es_pols, I_es_plants, I_es_both)

#create graph for species composition 
ggplot(I_es_all, aes(x = mutual_information, fill = type))+geom_histogram(color = "black", alpha = 0.6, position = "identity")+
  xlim(0,1)+ theme_classic()+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="mutual information", y="number of iterations")

##----partner comparison against shuffled networks ------------------------------------------------------------------------------------------------------
partner_comparison_pols <- NULL
partner_comparison_plants <- NULL
partner_comparison_both <- NULL


#pols
partner_comparison_pols <- partner_comparison_func(dryad_multilayer_shuf_1000_pols_output,
                                                   partner_comparison_pols)

partner_comparison_pols$node_id <- partner_comparison_pols$similarity$node %>% unlist() #unlist the stuff under similarity
partner_comparison_pols$jaccard_similarity <- partner_comparison_pols$similarity$similarity %>% unlist()
partner_comparison_pols$type <- partner_comparison_pols$similarity$type %>% unlist()
partner_comparison_pols <- partner_comparison_pols %>% select(trial_number, node_id, jaccard_similarity, type) #only save new columns
partner_comparison_pols$shuf_type <- "null_pollinators"

write.csv(partner_comparison_pols, "./csvs/partner_comparison_pols.csv", row.names = FALSE)

#plants
partner_comparison_plants <- partner_comparison_func(dryad_multilayer_shuf_1000_plants_output,
                                                     partner_comparison_plants)

partner_comparison_plants$node_id <- partner_comparison_plants$similarity$node %>% unlist() #unlist the stuff under similarity
partner_comparison_plants$jaccard_similarity <- partner_comparison_plants$similarity$similarity %>% unlist()
partner_comparison_plants$type <- partner_comparison_plants$similarity$type %>% unlist()
partner_comparison_plants <- partner_comparison_plants %>% select(trial_number, node_id, jaccard_similarity, type) #only save new columns
partner_comparison_plants$shuf_type <- "null_plants"

write.csv(partner_comparison_plants, "./csvs/partner_comparison_plants.csv", row.names = FALSE)

#both
partner_comparison_both <- partner_comparison_func(dryad_multilayer_shuf_1000_both_output,
                                                   partner_comparison_both)

partner_comparison_both$node_id <- partner_comparison_both$similarity$node %>% unlist() #unlist the stuff under similarity
partner_comparison_both$jaccard_similarity <- partner_comparison_both$similarity$similarity %>% unlist()
partner_comparison_both$type <- partner_comparison_both$similarity$type %>% unlist()
partner_comparison_both <- partner_comparison_both %>% select(trial_number, node_id, jaccard_similarity, type) #only save new columns
partner_comparison_both$shuf_type <- "null_both"

write.csv(partner_comparison_both, "./csvs/partner_comparison_both.csv", row.names = FALSE)

partner_comparison_shuf <- rbind(partner_comparison_pols, partner_comparison_plants, partner_comparison_both)


#create graph for partner similarity
partner_comparison_shuf %>%
  ggplot(aes(x=jaccard_similarity, fill = shuf_type))+geom_histogram(color = "black", alpha = 0.6, position = "identity")+
  theme_classic()+xlim(0,1)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Jaccard Similarity", y="number of iterations", fill = "type") 

##correlation and r sqaured between jaccard and distance for each run ----------------------------------------------------------------------------
#pols
iteration_correlation_pols <- NULL
iteration_correlation_data_pols <- all_edge_list_island_combine_no_module_shuf_pols_output %>% subset(layer_from != layer_to) 

for (i in 1:1000){
  trial_pols = iteration_correlation_data_pols %>% filter(trial == i)
  iteration_correlation_new_pols <- cor.test(trial_pols$turnover, trial_pols$ave_distance, method = "pearson")
  lm_val_pols <- lm(turnover ~ ave_distance, data = trial_pols)
  iteration_correlation_pols <- rbind(iteration_correlation_pols, tibble(estimate = iteration_correlation_new_pols$estimate, 
                                                               p_val = iteration_correlation_new_pols$p.value, 
                                                               statistic = iteration_correlation_new_pols$statistic, 
                                                               confidence_int_low = iteration_correlation_new_pols$conf.int[1],
                                                               confidence_int_high = iteration_correlation_new_pols$conf.int[2],
                                                               slope = lm_val_pols$coefficients[2],
                                                               intercept = lm_val_pols$coefficients[1],
                                                               rsquared = summary(lm_val_pols)$adj.r.squared,
                                                               trial_num = i))
}


#write.csv(iteration_correlation_pols, "./csvs/iteration_correlation_pols.csv", row.names = FALSE)

#correlation empirical
islands_turnover_with_distnace_empirical_no_loop <- islands_turnover_with_distnace_empirical %>% subset(layer_from != layer_to) 

correlation_empirical_data_pols <- cor.test(islands_turnover_with_distnace_empirical_no_loop$turnover, #ave is just the value of the turnover
                                            islands_turnover_with_distnace_empirical_no_loop$ave_distance, method = "pearson")

lm_val_empirical_pols <- lm(turnover ~ ave_distance, data = islands_turnover_with_distnace_empirical_no_loop)

correlation_empirical_pols <- tibble(estimate = correlation_empirical_data_pols$estimate, 
                                p_val = correlation_empirical_data_pols$p.value, 
                                statistic = correlation_empirical_data_pols$statistic, 
                                confidence_int_low = correlation_empirical_data_pols$conf.int[1],
                                confidence_int_high = correlation_empirical_data_pols$conf.int[2],
                                slope = lm_val_empirical_pols$coefficients[2],
                                intercept = lm_val_empirical_pols$coefficients[1], 
                                rsquared = summary(lm_val_empirical_pols)$adj.r.squared)

#write.csv(correlation_empirical_pols, "./csvs/correlation_empirical_pols.csv", row.names = FALSE) #so it can be used for classical shuffling
#correlation_empirical_pols <- read.csv("./csvs/correlation_empirical_pols.csv")

##distribution of estimate (r) and add empirical
iteration_correlation_pols %>% ggplot(aes(x = estimate))+ geom_density(fill = "#BE75FA", color = "#BE75FA", alpha = 0.4)+ 
  theme_classic()+ labs(x = "r")+
  geom_vline(xintercept = correlation_empirical_pols$estimate, linetype = "dashed", color = "#F47069") 

p_r_pols <- sum(iteration_correlation_pols$estimate < correlation_empirical_pols$estimate)/1000
p_r_pols

##distribution of p-val and add empirical
iteration_correlation_pols %>% ggplot(aes(x = p_val))+ geom_density(fill = "#BE75FA", color = "#BE75FA", alpha = 0.4)+ theme_classic()+ labs(x = "p-value")+
  geom_vline(xintercept = correlation_empirical_pols$p_val, linetype = "dashed", color = "#F47069") 

p_p_val_pols <- sum(iteration_correlation_pols$p_val < correlation_empirical_pols$p_val)/1000
p_p_val_pols

##distribution of slope and add empirical
iteration_correlation_pols %>% ggplot(aes(x = slope))+ geom_density(fill = "#BE75FA", color = "#BE75FA", alpha = 0.4)+ theme_classic()+ labs(x = "slope")+
  geom_vline(xintercept = correlation_empirical_pols$slope, linetype = "dashed", color = "#F47069") 

p_slope_pols <- sum(iteration_correlation_pols$slope > correlation_empirical_pols$slope)/1000
p_slope_pols

##distribution of intercept and add empirical
iteration_correlation_pols %>% ggplot(aes(x = intercept))+ geom_density(fill = "#BE75FA", color = "#BE75FA", alpha = 0.4)+ theme_classic()+ labs(x = "intercept")+
  geom_vline(xintercept = correlation_empirical_pols$intercept, linetype = "dashed", color = "#F47069") 

p_intercept_pols <- sum(iteration_correlation_pols$intercept > correlation_empirical_pols$intercept)/1000
p_intercept_pols

##distribution of rsquared and add empirical
iteration_correlation_pols %>% ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#BE75FA", color = "#BE75FA", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_pols$rsquared, linetype = "dashed", color = "#F47069")+
  theme(axis.title=element_text(size=22))

p_rsquared_pols <- sum(iteration_correlation_pols$rsquared > correlation_empirical_pols$rsquared)/1000
p_rsquared_pols

#---- distribution for each dot on graph and its corresponding partner -----------------------------------------------

correlation_empirical_data_partial_pols <- islands_turnover_with_distnace_empirical_no_loop %>% 
  select(layer_from, layer_to, turnover) %>% rename(empirical_turnover = turnover)

pair_distribution_pols <- merge(iteration_correlation_data_pols, correlation_empirical_data_partial_pols, 
                           by = c("layer_from", "layer_to")) #create new column where empirical turnover is found

pair_distribution_pols %>% group_by(layer_from, layer_to) %>%
  ggplot(aes(x = turnover))+ geom_density(fill = "#BE75FA", color = "#BE75FA", alpha = 0.4)+ theme_classic()+ labs(x = "turnover")+
  facet_wrap(~layer_from + layer_to)+ geom_vline(data = pair_distribution_pols %>% group_by(layer_from, layer_to), #group by same thing as rest of panel
                                                 aes(xintercept = empirical_turnover), #pair the correct empirical with the rest of the distribution
                                                 linetype = "dashed", color = "#F47069") 

##one tailed test p-value to determine if significantly different
p_values_for_pairs_pols <- NULL

for (i in 1:21){
  layer_from_emp <- correlation_empirical_data_partial_pols[i,]$layer_from
  layer_to_emp <- correlation_empirical_data_partial_pols[i,]$layer_to
  if (layer_from_emp == 1 && layer_to_emp == 6) next
  empirical_current_pair <- correlation_empirical_data_partial_pols[i,]$empirical_turnover #take empirical value
  current_pair <- iteration_correlation_data_pols %>% filter(layer_from == layer_from_emp, layer_to == layer_to_emp) #choose compatible p air to empirical
  current_pair_mean <- mean(current_pair$turnover)
  if (current_pair_mean > empirical_current_pair){
    p <- sum(current_pair$turnover < empirical_current_pair)/1000
  }
  else{
    p <- sum(current_pair$turnover > empirical_current_pair)/1000
  }
  p_values_for_pairs_pols <- rbind(p_values_for_pairs_pols, tibble(layer_from = layer_from_emp, layer_to = layer_to_emp, 
                                                         p_value = p))
}

#amount its smaller that 0.01
amount_pols <- sum(p_values_for_pairs_pols$p_value <= 0.01)/21
amount_pols
amount_0.001_pols <- sum(p_values_for_pairs_pols$p_value <= 0.001)/21
amount_0.001_pols



#----plants
iteration_correlation_plants <- NULL
iteration_correlation_data_plants <- all_edge_list_island_combine_no_module_shuf_plants_output %>% subset(layer_from != layer_to) 

for (i in 1:1000){
  trial_plants = iteration_correlation_data_plants %>% filter(trial == i)
  iteration_correlation_new_plants <- cor.test(trial_plants$turnover, trial_plants$ave_distance, method = "pearson")
  lm_val_plants <- lm(turnover ~ ave_distance, data = trial_plants)
  iteration_correlation_plants <- rbind(iteration_correlation_plants, tibble(estimate = iteration_correlation_new_plants$estimate, 
                                                                         p_val = iteration_correlation_new_plants$p.value, 
                                                                         statistic = iteration_correlation_new_plants$statistic, 
                                                                         confidence_int_low = iteration_correlation_new_plants$conf.int[1],
                                                                         confidence_int_high = iteration_correlation_new_plants$conf.int[2],
                                                                         slope = lm_val_plants$coefficients[2],
                                                                         intercept = lm_val_plants$coefficients[1],
                                                                         rsquared = summary(lm_val_plants)$adj.r.squared,
                                                                         trial_num = i))
}

#write.csv(iteration_correlation_plants, "./csvs/iteration_correlation_plants.csv", row.names = FALSE)


#correlation empirical is the same as pols
##distribution of estimate (r) and add empirical
iteration_correlation_plants %>% ggplot(aes(x = estimate))+ geom_density(fill = "#15B7BC", color = "#15B7BC", alpha = 0.4)+ 
  theme_classic()+ labs(x = "r")+
  geom_vline(xintercept = correlation_empirical_pols$estimate, linetype = "dashed", color = "#F47069") 

p_r_plants <- sum(iteration_correlation_plants$estimate > correlation_empirical_pols$estimate)/1000
p_r_plants

##distribution of p-val and add empirical
iteration_correlation_plants %>% ggplot(aes(x = p_val))+ geom_density(fill = "#15B7BC", color = "#15B7BC", alpha = 0.4)+ theme_classic()+ labs(x = "p-value")+
  geom_vline(xintercept = correlation_empirical_pols$p_val, linetype = "dashed", color = "#F47069") 

p_p_val_plants <- sum(iteration_correlation_plants$p_val > correlation_empirical_pols$p_val)/1000
p_p_val_plants

##distribution of slope and add empirical
iteration_correlation_plants %>% ggplot(aes(x = slope))+ geom_density(fill = "#15B7BC", color = "#15B7BC", alpha = 0.4)+ theme_classic()+ labs(x = "slope")+
  geom_vline(xintercept = correlation_empirical_pols$slope, linetype = "dashed", color = "#F47069") 

p_slope_plants <- sum(iteration_correlation_plants$slope < correlation_empirical_pols$slope)/1000
p_slope_plants

##distribution of intercept and add empirical
iteration_correlation_plants %>% ggplot(aes(x = intercept))+ geom_density(fill = "#15B7BC", color = "#15B7BC", alpha = 0.4)+ theme_classic()+ labs(x = "intercept")+
  geom_vline(xintercept = correlation_empirical_pols$intercept, linetype = "dashed", color = "#F47069") 

p_intercept_plants <- sum(iteration_correlation_plants$intercept < correlation_empirical_pols$intercept)/1000
p_intercept_plants

##distribution of rsquared and add empirical
iteration_correlation_plants %>% ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#15B7BC", color = "#15B7BC", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_pols$rsquared, linetype = "dashed", color = "#F47069")+
  theme(axis.title=element_text(size=22))

p_rsquared_plants <- sum(iteration_correlation_plants$rsquared < correlation_empirical_pols$rsquared)/1000
p_rsquared_plants

#---- distribution for each dot on graph and its corresponding partner -----------------------------------------------

pair_distribution_plants <- merge(iteration_correlation_data_plants, correlation_empirical_data_partial_pols, #empirical is the same as pols
                                by = c("layer_from", "layer_to")) #create new column where empirical turnover is found

pair_distribution_plants %>% group_by(layer_from, layer_to) %>%
  ggplot(aes(x = turnover))+ geom_density(fill = "#15B7BC", color = "#15B7BC", alpha = 0.4)+ theme_classic()+ labs(x = "turnover")+
  facet_wrap(~layer_from + layer_to)+ geom_vline(data = pair_distribution_plants %>% group_by(layer_from, layer_to), #group by same thing as rest of panel
                                                 aes(xintercept = empirical_turnover), #pair the correct empirical with the rest of the distribution
                                                 linetype = "dashed", color = "#F47069") 

##one tailed test p-value to determine if significantly different
p_values_for_pairs_plants <- NULL

for (i in 1:21){
  layer_from_emp <- correlation_empirical_data_partial_pols[i,]$layer_from
  layer_to_emp <- correlation_empirical_data_partial_pols[i,]$layer_to
  empirical_current_pair <- correlation_empirical_data_partial_pols[i,]$empirical_turnover #take empirical value
  current_pair <- iteration_correlation_data_plants %>% filter(layer_from == layer_from_emp, layer_to == layer_to_emp) #choose compatible p air to empirical
  current_pair_mean <- mean(current_pair$turnover)
  if (current_pair_mean > empirical_current_pair){
    p <- sum(current_pair$turnover < empirical_current_pair)/1000
  }
  else{
    p <- sum(current_pair$turnover > empirical_current_pair)/1000
  }
  p_values_for_pairs_plants <- rbind(p_values_for_pairs_plants, tibble(layer_from = layer_from_emp, layer_to = layer_to_emp, 
                                                                   p_value = p))
}

#amount its smaller that 0.01
amount_plants <- sum(p_values_for_pairs_plants$p_value <= 0.01)/21
amount_plants
amount_0.001_plants <- sum(p_values_for_pairs_plants$p_value <= 0.001)/21
amount_0.001_plants

#---- both
iteration_correlation_both <- NULL
iteration_correlation_data_both <- all_edge_list_island_combine_no_module_shuf_both_output %>% subset(layer_from != layer_to) 

for (i in 1:1000){
  trial_both = iteration_correlation_data_both %>% filter(trial == i)
  iteration_correlation_new_both <- cor.test(trial_both$turnover, trial_both$ave_distance, method = "pearson")
  lm_val_both <- lm(turnover ~ ave_distance, data = trial_both)
  iteration_correlation_both <- rbind(iteration_correlation_both, tibble(estimate = iteration_correlation_new_both$estimate, 
                                                                             p_val = iteration_correlation_new_both$p.value, 
                                                                             statistic = iteration_correlation_new_both$statistic, 
                                                                             confidence_int_low = iteration_correlation_new_both$conf.int[1],
                                                                             confidence_int_high = iteration_correlation_new_both$conf.int[2],
                                                                             slope = lm_val_both$coefficients[2],
                                                                             intercept = lm_val_both$coefficients[1],
                                                                             rsquared = summary(lm_val_both)$adj.r.squared,
                                                                             trial_num = i))
}

#write.csv(iteration_correlation_both, "./csvs/iteration_correlation_both.csv", row.names = FALSE)

#correlation empirical is the same as pols
##distribution of estimate (r) and add empirical
iteration_correlation_both %>% ggplot(aes(x = estimate))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ 
  theme_classic()+ labs(x = "r")+
  geom_vline(xintercept = correlation_empirical_pols$estimate, linetype = "dashed", color = "#F47069") 

p_r_both <- sum(iteration_correlation_both$estimate < correlation_empirical_pols$estimate)/1000
p_r_both

##distribution of p-val and add empirical
iteration_correlation_both %>% ggplot(aes(x = p_val))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "p-value")+
  geom_vline(xintercept = correlation_empirical_pols$p_val, linetype = "dashed", color = "#F47069") 

p_p_val_both <- sum(iteration_correlation_both$p_val < correlation_empirical_pols$p_val)/1000
p_p_val_both

##distribution of slope and add empirical
iteration_correlation_both %>% ggplot(aes(x = slope))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "slope")+
  geom_vline(xintercept = correlation_empirical_pols$slope, linetype = "dashed", color = "#F47069") 

p_slope_both <- sum(iteration_correlation_both$slope > correlation_empirical_pols$slope)/1000
p_slope_both

##distribution of intercept and add empirical
iteration_correlation_both %>% ggplot(aes(x = intercept))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "intercept")+
  geom_vline(xintercept = correlation_empirical_pols$intercept, linetype = "dashed", color = "#F47069") 

p_intercept_both <- sum(iteration_correlation_both$intercept > correlation_empirical_pols$intercept)/1000
p_intercept_both

##distribution of rsquared and add empirical
iteration_correlation_both %>% ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_pols$rsquared, linetype = "dashed", color = "#F47069")+
  theme(axis.title=element_text(size=22))

p_rsquared_both <- sum(iteration_correlation_both$rsquared > correlation_empirical_pols$rsquared)/1000
p_rsquared_both

#---- distribution for each dot on graph and its corresponding partner -----------------------------------------------

pair_distribution_both <- merge(iteration_correlation_data_both, correlation_empirical_data_partial_pols, #empirical is the same as pols
                                  by = c("layer_from", "layer_to")) #create new column where empirical turnover is found

pair_distribution_both %>% group_by(layer_from, layer_to) %>%
  ggplot(aes(x = turnover))+ geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ theme_classic()+ labs(x = "turnover")+
  facet_wrap(~layer_from + layer_to)+ geom_vline(data = pair_distribution_both %>% group_by(layer_from, layer_to), #group by same thing as rest of panel
                                                 aes(xintercept = empirical_turnover), #pair the correct empirical with the rest of the distribution
                                                 linetype = "dashed", color = "#F47069") 

##one tailed test p-value to determine if significantly different
p_values_for_pairs_both <- NULL

for (i in 1:21){
  layer_from_emp <- correlation_empirical_data_partial_pols[i,]$layer_from
  layer_to_emp <- correlation_empirical_data_partial_pols[i,]$layer_to
  empirical_current_pair <- correlation_empirical_data_partial_pols[i,]$empirical_turnover #take empirical value
  current_pair <- iteration_correlation_data_both %>% filter(layer_from == layer_from_emp, layer_to == layer_to_emp) #choose compatible p air to empirical
  current_pair_mean <- mean(current_pair$turnover)
  if (current_pair_mean > empirical_current_pair){
    p <- sum(current_pair$turnover < empirical_current_pair)/1000
  }
  else{
    p <- sum(current_pair$turnover > empirical_current_pair)/1000
  }
  p_values_for_pairs_both <- rbind(p_values_for_pairs_both, tibble(layer_from = layer_from_emp, layer_to = layer_to_emp, 
                                                                       p_value = p))
}

#amount its smaller that 0.01
amount_both <- sum(p_values_for_pairs_both$p_value <= 0.01)/21
amount_both
amount_0.001_both <- sum(p_values_for_pairs_both$p_value <= 0.001)/21
amount_0.001_both


