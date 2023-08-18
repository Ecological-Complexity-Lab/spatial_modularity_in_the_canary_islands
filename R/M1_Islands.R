# ---- NULL MODEL M1: SHUFFLING SPECIES BETWEEN ISLANDS ------------------------------------------------------------------------

# this portion of the code shuffles the networks between layers in one of 3 versions:
# 1. shuffling plants among themselves
# 2. shuffling pollinators among themselves
# 3. shuffling plants among themselves and then pollinators among themselves
# the shuffled networks are then compared the empirical network to determine whether 
# certain network properties are random or not

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

##---- create interedges  --------------------------------------------------------
# keep only species that occur in 2 or more layers
co_occurrence_from_shuf<- shuf_null_edge_list %>% 
  group_by(trial_number,node_from) %>%
  mutate(num_layers_from=n_distinct(layer_from)) %>% 
  filter(num_layers_from>="2")

# a for loop that calculates all the interlayer edges based on jaccard index
interlayers_with_weights_islands <- NULL

for (trial in 1:1000){
  print(trial) #to keep tab on how far along we are
  
  co_occurrence_from_shuf2 <- co_occurrence_from_shuf %>% filter(trial_number == trial) #take only 1 trial
  
for (i in unique(co_occurrence_from_shuf2$node_from)) {
  print(i)
  partners_sp <- 
    co_occurrence_from_shuf2 %>% 
    filter(node_from== i) %>%
    group_by(node_to) %>%
    select(c(node_to,layer_from)) %>%
    distinct() %>% 
    mutate(present=1) %>%
    spread(node_to, present, fill = 0) %>%
    column_to_rownames("layer_from")
  
  beta_layers_sp <- 1-as.matrix(vegdist(partners_sp, "jaccard"))
  beta_layers_sp_m <- melt(as.matrix(extRC::tril(beta_layers_sp)))
  inter_fid <- beta_layers_sp_m  %>%
    tibble() %>%
    filter(value!=0) %>%
    subset(Var1 != Var2) %>%
    mutate(node_from=i, node_to =i, trial_num=trial) %>%
    select(c(node_from,layer_from=Var1, layer_to=Var2,node_to, weight=value,trial_num))
  interlayers_with_weights_islands <- rbind(interlayers_with_weights_islands,inter_fid)
}
}

interlayer_edges_shuf_pols<-interlayers_with_weights_islands 
#write.csv(interlayer_edges_shuf_pols, "./csvs/Islands/Jac/interlayer_edges_shuf_pols_islands_as_layers.csv",row.names = FALSE)


##---- create intraedges inverted versions and put weight ---------------------------------
dryad_intralayer_shuf_pols <- read.csv("csvs/Islands/Jac/shuf_null_edge_list_islands_as_layers.csv") 
dryad_intralayer_shuf_pols <- dryad_intralayer_shuf_pols[, c(6,1,2,3,4,5)]

#----inverted versions
#pols
intralayer_inverted_shuf_pols <- tibble(values= dryad_intralayer_shuf_pols$layer_to, dryad_intralayer_shuf_pols$node_to, 
                                        dryad_intralayer_shuf_pols$layer_from, dryad_intralayer_shuf_pols$node_from, 
                                        dryad_intralayer_shuf_pols$weight, dryad_intralayer_shuf_pols$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_pols) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")

#--- weighted
## ----weighted intralayer edges--------------------------------------------------------------------------
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

## ----multilayer_extended_final--------------------------------------------------------------------------------------
edgelist_intralayer_shuf_pols <- bind_rows(intralayer_weighted_shuf_pols, intralayer_weighted_inverted_shuf_pols)#intraedges
edgelist_intralayer_shuf_pols<-edgelist_intralayer_shuf_pols[,c(1,2,3,4,6,5)]#change order columns
colnames(edgelist_intralayer_shuf_pols)[6]<-"trial_num"

interlayer_edges_shuf_pols<-read.csv("csvs/Islands/Jac/interlayer_edges_shuf_pols_islands_as_layers.csv")
dryad_interlayer_shuf_pols<-interlayer_edges_shuf_pols#interedges
dryad_interlayer_shuf_pols<-interlayer_edges_shuf_pols[,c(2,1,3,4,5,6)]#change order columns
dryad_interlayer_shuf_pols$node_from<-as.integer(dryad_interlayer_shuf_pols$node_from)
dryad_interlayer_shuf_pols$node_to<-as.integer(dryad_interlayer_shuf_pols$node_to)


dryad_edgelist_complete_shuf_pols <- bind_rows(edgelist_intralayer_shuf_pols, dryad_interlayer_shuf_pols) #combine inter and intra
#write.csv(dryad_edgelist_complete_shuf_pols, "./csvs/Islands/Jac/dryad_edgelist_complete_shuf_pols_islands_as_layers.csv", row.names = FALSE)



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

#write.csv(shuf_null_edge_list_plants, "./csvs/Islands/Jac/shuf_null_edge_list_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_plants, "./csvs/Islands/Jac/shuf_trial_matrix_plants_islands_as_layers.csv", row.names = TRUE)
#shuf_null_edge_list_plants <- read.csv("./csvs/Islands/Jac/shuf_null_edge_list_plants_islands_as_layers.csv")
#shuf_trial_matrix_plants <- read.csv("./csvs/Islands/Jac/shuf_trial_matrix_plants_islands_as_layers.csv")

##---- create interedges  --------------------------------------------------------
# keep only species that occur in 2 or more layers
co_occurrence_from_shuf<- shuf_null_edge_list_plants %>% 
  group_by(trial_number,node_from) %>%
  mutate(num_layers_from=n_distinct(layer_from)) %>% 
  filter(num_layers_from>="2")

# a for loop that calculates all the interlayer edges based on jaccard index
interlayers_with_weights_islands <- NULL

for (trial in 1:1000){
  print(trial) #to keep tab on how far along we are
  
  co_occurrence_from_shuf2 <- co_occurrence_from_shuf %>% filter(trial_number == trial) #take only 1 trial
  
  for (i in unique(co_occurrence_from_shuf2$node_from)) {
    print(i)
    partners_sp <- 
      co_occurrence_from_shuf2 %>% 
      filter(node_from== i) %>%
      group_by(node_to) %>%
      select(c(node_to,layer_from)) %>%
      distinct() %>% 
      mutate(present=1) %>%
      spread(node_to, present, fill = 0) %>%
      column_to_rownames("layer_from")
    
    beta_layers_sp <- 1-as.matrix(vegdist(partners_sp, "jaccard"))
    beta_layers_sp_m <- melt(as.matrix(extRC::tril(beta_layers_sp)))
    inter_fid <- beta_layers_sp_m  %>%
      tibble() %>%
      filter(value!=0) %>%
      subset(Var1 != Var2) %>%
      mutate(node_from=i, node_to =i, trial_num=trial) %>%
      select(c(node_from,layer_from=Var1, layer_to=Var2,node_to, weight=value,trial_num))
    interlayers_with_weights_islands <- rbind(interlayers_with_weights_islands,inter_fid)
  }
}

interlayer_edges_shuf_plants<-interlayers_with_weights_islands 
#write.csv(interlayer_edges_shuf_plants, "./csvs/Islands/Jac/interlayer_edges_shuf_plants_islands_as_layers.csv",row.names = FALSE)

##---- create intraedges inverted versions and put weight ---------------------------------
dryad_intralayer_shuf_plants <- read.csv("csvs/Islands/Jac/shuf_null_edge_list_plants_islands_as_layers.csv") 
dryad_intralayer_shuf_plants <- dryad_intralayer_shuf_plants[, c(6,1,2,3,4,5)]

#----inverted versions
#pols
intralayer_inverted_shuf_plants <- tibble(values= dryad_intralayer_shuf_plants$layer_to, dryad_intralayer_shuf_plants$node_to, 
                                          dryad_intralayer_shuf_plants$layer_from, dryad_intralayer_shuf_plants$node_from, 
                                          dryad_intralayer_shuf_plants$weight, dryad_intralayer_shuf_plants$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_plants) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")

#--- weighted
## ----weighted intralayer edges--------------------------------------------------------------------------
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

## ----multilayer_extended_final--------------------------------------------------------------------------------------
edgelist_intralayer_shuf_plants <- bind_rows(intralayer_weighted_shuf_plants, intralayer_weighted_inverted_shuf_plants)#intraedges
edgelist_intralayer_shuf_plants<-edgelist_intralayer_shuf_plants[,c(1,2,3,4,6,5)]#change order columns
colnames(edgelist_intralayer_shuf_plants)[6]<-"trial_num"

dryad_interlayer_shuf_plants<-interlayer_edges_shuf_plants#interedges
dryad_interlayer_shuf_plants<-interlayer_edges_shuf_plants[,c(2,1,3,4,5,6)]#change order columns
dryad_interlayer_shuf_plants$node_from<-as.integer(dryad_interlayer_shuf_plants$node_from)
dryad_interlayer_shuf_plants$node_to<-as.integer(dryad_interlayer_shuf_plants$node_to)

dryad_edgelist_complete_shuf_plants <- bind_rows(edgelist_intralayer_shuf_plants, dryad_interlayer_shuf_plants) #combine inter and intra
#write.csv(dryad_edgelist_complete_shuf_plants, "./csvs/Islands/Jac/dryad_edgelist_complete_shuf_plants_islands_as_layers.csv", row.names = FALSE)



##---- shuffling both plants and pollinators -------------------------------------------------------------------
intralayer_matrix_both <- shuf_trial_matrix #now we'll change columns and not rows

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
#write.csv(shuf_null_edge_list_both, "./csvs/Islands/Jac/shuf_null_edge_list_both_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_both, "./csvs/Islands/Jac/shuf_trial_matrix_both_islands_as_layers.csv", row.names = TRUE)

##---- create interedges  --------------------------------------------------------
# keep only species that occur in 2 or more layers
co_occurrence_from_shuf<- shuf_null_edge_list_both %>% 
  group_by(trial_number,node_from) %>%
  mutate(num_layers_from=n_distinct(layer_from)) %>% 
  filter(num_layers_from>="2")

# a for loop that calculates all the interlayer edges based on jaccard index
interlayers_with_weights_islands <- NULL

for (trial in 1:1000){
  print(trial) #to keep tab on how far along we are
  
  co_occurrence_from_shuf2 <- co_occurrence_from_shuf %>% filter(trial_number == trial) #take only 1 trial
  
  for (i in unique(co_occurrence_from_shuf2$node_from)) {
    print(i)
    partners_sp <- 
      co_occurrence_from_shuf2 %>% 
      filter(node_from== i) %>%
      group_by(node_to) %>%
      select(c(node_to,layer_from)) %>%
      distinct() %>% 
      mutate(present=1) %>%
      spread(node_to, present, fill = 0) %>%
      column_to_rownames("layer_from")
    
    beta_layers_sp <- 1-as.matrix(vegdist(partners_sp, "jaccard"))
    beta_layers_sp_m <- melt(as.matrix(extRC::tril(beta_layers_sp)))
    inter_fid <- beta_layers_sp_m  %>%
      tibble() %>%
      filter(value!=0) %>%
      subset(Var1 != Var2) %>%
      mutate(node_from=i, node_to =i, trial_num=trial) %>%
      select(c(node_from,layer_from=Var1, layer_to=Var2,node_to, weight=value,trial_num))
    interlayers_with_weights_islands <- rbind(interlayers_with_weights_islands,inter_fid)
  }
}

interlayer_edges_shuf_both<-interlayers_with_weights_islands 
#write.csv(interlayer_edges_shuf_both, "./csvs/Islands/Jac/interlayer_edges_shuf_both_islands_as_layers.csv",row.names = FALSE)

##---- create intraedges inverted versions and put weight ---------------------------------
dryad_intralayer_shuf_both <- read.csv("csvs/Islands/Jac/shuf_null_edge_list_both_islands_as_layers.csv") 
dryad_intralayer_shuf_both <- dryad_intralayer_shuf_both[, c(6,1,2,3,4,5)]

#----inverted versions
#pols
intralayer_inverted_shuf_both <- tibble(values= dryad_intralayer_shuf_both$layer_to, dryad_intralayer_shuf_both$node_to, 
                                          dryad_intralayer_shuf_both$layer_from, dryad_intralayer_shuf_both$node_from, 
                                          dryad_intralayer_shuf_both$weight, dryad_intralayer_shuf_both$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_both) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")

#--- weighted
## ----weighted intralayer edges--------------------------------------------------------------------------
#plants in from
tot_plant_shuf_both <- intralayer_inverted_shuf_both %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_both <- intralayer_inverted_shuf_both %>% left_join(tot_plant_shuf_both) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_both <- tot_plant_shuf_both[, c(3,1,2,4)]

#pols in from
tot_pol_shuf_both<- dryad_intralayer_shuf_both %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_both <- dryad_intralayer_shuf_both %>% left_join(tot_pol_shuf_both) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_both <- tot_pol_shuf_both[, c(3,1,2,4)]

## ----multilayer_extended_final--------------------------------------------------------------------------------------
edgelist_intralayer_shuf_both <- bind_rows(intralayer_weighted_shuf_both, intralayer_weighted_inverted_shuf_both)#intraedges
edgelist_intralayer_shuf_both<-edgelist_intralayer_shuf_both[,c(1,2,3,4,6,5)]#change order columns
colnames(edgelist_intralayer_shuf_both)[6]<-"trial_num"

dryad_interlayer_shuf_both<-interlayer_edges_shuf_both#interedges
dryad_interlayer_shuf_both<-interlayer_edges_shuf_both[,c(2,1,3,4,5,6)]#change order columns
dryad_interlayer_shuf_both$node_from<-as.integer(dryad_interlayer_shuf_both$node_from)
dryad_interlayer_shuf_both$node_to<-as.integer(dryad_interlayer_shuf_both$node_to)

dryad_edgelist_complete_shuf_both <- bind_rows(edgelist_intralayer_shuf_both, dryad_interlayer_shuf_both) #combine inter and intra
#write.csv(dryad_edgelist_complete_shuf_both, "./csvs/Islands/Jac/dryad_edgelist_complete_shuf_both_islands_as_layers.csv", row.names = FALSE)

##CALCULAR MODULARIDAD Y SIMILARIDAD PARA CADA UNO Y HACER LAS REGRESIONES (A PARTIR LINEA 827)

#ARMAR RMARKDOWN PARA MOSTRAR SHAI

