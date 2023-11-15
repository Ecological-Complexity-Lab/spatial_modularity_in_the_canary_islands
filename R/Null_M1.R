# ---- NULL MODEL M1: SHUFFLING SPECIES BETWEEN ISLANDS ------------------------------------------------------------------------

# this portion of the code shuffles species between layers in one of 3 versions:
# 1. shuffling plants among themselves
# 2. shuffling pollinators among themselves
# 3. shuffling plants among themselves and then pollinators among themselves
# the shuffled networks are then compared the empirical network to determine whether 
# species turnover is influencing the pattern found

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
library(ecodist)


##----get_data--------------------------------------------------------------------------------------------------------
setwd("D:/Trabajo/Papers/Canary_Island/spatial_modularity_in_the_canary_islands")
source("D:/Trabajo/Papers/Canary_Island/spatial_modularity_in_the_canary_islands/R/functions.R")

dryad_intralayer <- read.csv("./csvs_nuevo/intralayer_file.csv")

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
physical_nodes <- read.csv("./csvs_nuevo/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs_nuevo/layer_metadata_islands.csv")

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

#write.csv(shuf_null_edge_list, "./csvs_nuevo/shuf_null_edge_list_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix, "./csvs_nuevo/shuf_trial_matrix_islands_as_layers.csv", row.names = TRUE)
#shuf_null_edge_list <- read.csv("./csvs_nuevo/shuf_null_edge_list_islands_as_layers.csv")
shuf_trial_matrix <- read.csv("./csvs_nuevo/shuf_trial_matrix_islands_as_layers.csv", row.names = 1)

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
#write.csv(interlayer_edges_shuf_pols, "./csvs_nuevo/interlayer_edges_shuf_pols_islands_as_layers.csv",row.names = FALSE)

interlayer_edges_shuf_pols <- read.csv("./csvs_nuevo/interlayer_edges_shuf_pols_islands_as_layers.csv")
interlayer_edges_shuf_pols<-interlayer_edges_shuf_pols[,c(2,1,3,4,5,6)]#change order columns

#inverted version
interlayer_inverted <- tibble(values= interlayer_edges_shuf_pols$layer_to, interlayer_edges_shuf_pols$node_to, interlayer_edges_shuf_pols$layer_from, 
                              interlayer_edges_shuf_pols$node_from, interlayer_edges_shuf_pols$weight, interlayer_edges_shuf_pols$trial_num) #create an inverted copy for directed intralayers
colnames(interlayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight","trial_num")

#Create interedgelist
edgelist_interlayers_pols <- bind_rows(interlayer_edges_shuf_pols, interlayer_inverted) #combine inverted and non inverted versions of intra


##---- create intraedges inverted versions and put weight ---------------------------------
dryad_intralayer_shuf_pols <- read.csv("./csvs_nuevo/shuf_null_edge_list_islands_as_layers.csv") 
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

dryad_edgelist_complete_shuf_pols <- bind_rows(edgelist_intralayer_shuf_pols, edgelist_interlayers_pols) #combine inter and intra
#write.csv(dryad_edgelist_complete_shuf_pols, "./csvs_nuevo/dryad_edgelist_complete_shuf_pols_islands_as_layers.csv", row.names = FALSE)



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

#write.csv(shuf_null_edge_list_plants, "./csvs_nuevo/shuf_null_edge_list_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_plants, "./csvs_nuevo/shuf_trial_matrix_plants_islands_as_layers.csv", row.names = TRUE)
shuf_null_edge_list_plants <- read.csv("./csvs_nuevo/shuf_null_edge_list_plants_islands_as_layers.csv")
shuf_trial_matrix_plants <- read.csv("./csvs_nuevo/shuf_trial_matrix_plants_islands_as_layers.csv")

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
#write.csv(interlayer_edges_shuf_plants, "./csvs_nuevo/interlayer_edges_shuf_plants_islands_as_layers.csv",row.names = FALSE)

#interlayer_edges_shuf_plants <- read.csv("./csvs_nuevo/interlayer_edges_shuf_plants_islands_as_layers.csv")
interlayer_edges_shuf_plants<-interlayer_edges_shuf_plants[,c(2,1,3,4,5,6)]#change order columns

#inverted version
interlayer_inverted <- tibble(values= interlayer_edges_shuf_plants$layer_to, interlayer_edges_shuf_plants$node_to, interlayer_edges_shuf_plants$layer_from, 
                              interlayer_edges_shuf_plants$node_from, interlayer_edges_shuf_plants$weight, interlayer_edges_shuf_plants$trial_num) #create an inverted copy for directed intralayers
colnames(interlayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight","trial_num")

#Create interedgelist
edgelist_interlayers_plants <- bind_rows(interlayer_edges_shuf_plants, interlayer_inverted) #combine inverted and non inverted versions of intra


##---- create intraedges inverted versions and put weight ---------------------------------
dryad_intralayer_shuf_plants <- read.csv("./csvs_nuevo/shuf_null_edge_list_plants_islands_as_layers.csv") 
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

dryad_edgelist_complete_shuf_plants <- bind_rows(edgelist_intralayer_shuf_plants, edgelist_interlayers_plants) #combine inter and intra
#write.csv(dryad_edgelist_complete_shuf_plants, "./csvs_nuevo/dryad_edgelist_complete_shuf_plants_islands_as_layers.csv", row.names = FALSE)


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
#write.csv(shuf_null_edge_list_both, "./csvs_nuevo/shuf_null_edge_list_both_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_both, "./csvs_nuevo/shuf_trial_matrix_both_islands_as_layers.csv", row.names = TRUE)


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
#write.csv(interlayer_edges_shuf_both, "./csvs_nuevo/interlayer_edges_shuf_both_islands_as_layers.csv",row.names = FALSE)

#interlayer_edges_shuf_both <- read.csv("./csvs_nuevo/interlayer_edges_shuf_both_islands_as_layers.csv")
interlayer_edges_shuf_both <-interlayer_edges_shuf_both[,c(2,1,3,4,5,6)]#change order columns

#inverted version
interlayer_inverted <- tibble(values= interlayer_edges_shuf_both$layer_to, interlayer_edges_shuf_both$node_to, interlayer_edges_shuf_both$layer_from, 
                              interlayer_edges_shuf_both$node_from, interlayer_edges_shuf_both$weight, interlayer_edges_shuf_both$trial_num) #create an inverted copy for directed intralayers
colnames(interlayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight","trial_num")

#Create interedgelist
edgelist_interlayers_both <- bind_rows(interlayer_edges_shuf_both, interlayer_inverted) #combine inverted and non inverted versions of intra


##---- create intraedges inverted versions and put weight ---------------------------------
dryad_intralayer_shuf_both <- read.csv("./csvs_nuevo/shuf_null_edge_list_both_islands_as_layers.csv") 
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

dryad_edgelist_complete_shuf_both <- bind_rows(edgelist_intralayer_shuf_both, edgelist_interlayers_both) #combine inter and intra
#write.csv(dryad_edgelist_complete_shuf_both, "./csvs_nuevo/dryad_edgelist_complete_shuf_both_islands_as_layers.csv", row.names = FALSE)




## ---- Calculate number of modules and test distance decay in shuf networks --------------------------------------------------------

## ---- Number of modules ------------------------------------------------------------

# Input: An extended edge list.
dryad_edgelist_complete_shuf_pols <- read.csv("./csvs_nuevo/dryad_edgelist_complete_shuf_pols_islands_as_layers.csv")
colnames(dryad_edgelist_complete_shuf_pols)[6]<-"trial_number"

dryad_edgelist_complete_shuf_plants <- read.csv("./csvs_nuevo/dryad_edgelist_complete_shuf_plants_islands_as_layers.csv")
colnames(dryad_edgelist_complete_shuf_plants)[6]<-"trial_number"

dryad_edgelist_complete_shuf_both<- read.csv("./csvs_nuevo/dryad_edgelist_complete_shuf_both_islands_as_layers.csv")
colnames(dryad_edgelist_complete_shuf_both)[6]<-"trial_number"

layer_metadata <- read.csv("./csvs_nuevo/layer_metadata_islands.csv")
physical_nodes <- read.csv("./csvs_nuevo/physical_nodes_islands.csv")

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

#write.csv(dryad_multilayer_shuf_1000_pols_output, "./csvs_nuevo/dryad_multilayer_shuf_1000_pols_output_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_multilayer_shuf_1000_plants_output, "./csvs_nuevo/dryad_multilayer_shuf_1000_plants_output_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_multilayer_shuf_1000_both_output, "./csvs_nuevo/dryad_multilayer_shuf_1000_both_output_islands_as_layers.csv", row.names = FALSE)
dryad_multilayer_shuf_1000_pols_output <- read.csv("./csvs_nuevo/dryad_multilayer_shuf_1000_pols_output_islands_as_layers.csv")
dryad_multilayer_shuf_1000_plants_output <- read.csv("./csvs_nuevo/dryad_multilayer_shuf_1000_plants_output_islands_as_layers.csv")
dryad_multilayer_shuf_1000_both_output <- read.csv("./csvs_nuevo/dryad_multilayer_shuf_1000_both_output_islands_as_layers.csv")

##---- Distance decay of modules  ---------------------------------------------------------------
distances_with_ids <- read.csv("./csvs_nuevo/distances_with_ids_islands_as_layers.csv")

#pivot modules function for islands as layers
pivot_by_module_islands <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(data , id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s2<-na.omit(s2)
  s4 = t(s2) 
  s4 <- s4[-1,]
  colnames(s4) <- c(1,2,3,4,5,6,7)
  return(s4)
}

# this function calculates the Jaccard Similarity in modules between islands
module_distance_decay_islands_func <- function(multilayer_1000, 
                                               layers_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial) #to keep tab on how far along we are
    
    modules_for_similarity_shuf <- multilayer_1000  %>% filter(trial_num ==trial) #take only 1 trial
    
    #pivot modules
    module_pivoted_shuf <- pivot_by_module_islands(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    #create edge list with distances
    modules_edge_list_shuf <- NULL
    
    for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module_islands(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[k]
      if (is.null(modules_edge_list_shuf)) next
      modules_edge_list_shuf <- modules_edge_list_shuf %>% mutate(module = replace_na(module, current_module)) #add module number
    }
    
    edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
    #edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #we remove this line (which it removes NA), because we have a lot locations don't have modules in common
    
    
    for (i in 1:6){
      for (j in (i+1):7){
        modules_in_layer_from_shuf <- filter(modules_for_similarity_shuf, layer_id == i) %>% select(module) %>% unique() %>% unlist()
        modules_in_layer_to_shuf <- filter(modules_for_similarity_shuf, layer_id == j) %>% select(module) %>% unique() %>% unlist()
        #take all nodes in layer_from and all nodes in layer_to to check turnover
        int_both <- intersect(modules_in_layer_from_shuf, modules_in_layer_to_shuf) #how many nodes are found in both layers
        uni_both <- union(modules_in_layer_from_shuf, modules_in_layer_to_shuf)
        turnover <- length(int_both)/length(uni_both)
        module_layer_turnover_shuf <- rbind(module_layer_turnover_shuf, 
                                            tibble(layer_from = i, layer_to = j, turnover = turnover, trial = trial))
                                           
      }
    }
    layers_turnover_with_distnace <- edge_list_with_distances_shuf %>%
      merge(module_layer_turnover_shuf, by= c("layer_from", "layer_to")) #merge both versions
  }
  return(layers_turnover_with_distnace)
}

#function to create edge list per module function for islands as layers

edge_list_per_module_islands <- function(data,edge_list){
  #gets one row from a data frame and creates an edge list from it
  for (i in (1:6)){
    if (data[i]==0) next #only take layers where the module is present
    else {
      for (j in ((i+1):7)){
        if (data[j]==0) next #only take layers where the module is present
        else {
          edge_list <-rbind(edge_list,tibble(layer_from=i, layer_to=j, module=as.character(NA))) #create edge list of all the layer found in a module
        }
      }
    }
  }
  return(edge_list)
}


turnover_with_distance_pols <- NULL
turnover_with_distance_plants <- NULL
turnover_with_distance_both <- NULL
module_layer_turnover_shuf <- NULL
layers_turnover_with_distnace<-NULL


#shuff pols
all_edge_list_layer_combine_no_module_shuf_pols_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_pols_output,
                                                                                             turnover_with_distance_pols)


#write.csv(all_edge_list_layer_combine_no_module_shuf_pols_output, "./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_pols_output_islands.csv", row.names = FALSE)

#all_edge_list_layer_combine_no_module_shuf_pols_output <- read.csv("./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_pols_output_islands.csv")

#shuff plants
all_edge_list_layer_combine_no_module_shuf_plants_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_plants_output,
                                                                                               turnover_with_distance_plants)

#write.csv(all_edge_list_layer_combine_no_module_shuf_plants_output,  "./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_plants_output_islands.csv", row.names = FALSE)

#all_edge_list_layer_combine_no_module_shuf_plants_output <- read.csv("./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_plants_output_islands.csv")

#shuff both
all_edge_list_layer_combine_no_module_shuf_both_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_both_output,
                                                                                             turnover_with_distance_both)

#write.csv(all_edge_list_layer_combine_no_module_shuf_both_output, "./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_both_output_islands.csv", row.names = FALSE)

#all_edge_list_layer_combine_no_module_shuf_both_output <- read.csv("./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_both_output_islands.csv")


#---- create ave for jaccard islands -------------------------------------------------------------------------
#pols
ave_module_layer_turnover_shuf_pols <- all_edge_list_layer_combine_no_module_shuf_pols_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean(mean_distance)) %>% mutate(type="null_pollinators") #create mean and sd for each point


#write.csv(ave_module_layer_turnover_shuf_pols, "./csvs_nuevo/ave_module_layer_turnover_shuf_pols_islands.csv", row.names = FALSE)

#plants
ave_module_layer_turnover_shuf_plants <- all_edge_list_layer_combine_no_module_shuf_plants_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean(mean_distance)) %>% mutate(type="null_plants") #create mean and sd for each point


#write.csv(ave_module_layer_turnover_shuf_plants, "./csvs_nuevo/ave_module_layer_turnover_shuf_plants_islands.csv", row.names = FALSE)

#both
ave_module_layer_turnover_shuf_both <- all_edge_list_layer_combine_no_module_shuf_both_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean(mean_distance)) %>% mutate(type="null_both") #create mean and sd for each point

#write.csv(ave_module_layer_turnover_shuf_both,  "./csvs_nuevo/ave_module_layer_turnover_shuf_both_islands.csv",  row.names = FALSE)

#add empirical
islands_turnover_with_distnace_empirical <- read.csv("./csvs_nuevo/islands_turnover_with_distnace_empirical.csv")

empirical_turnover_for_modules_layers_shuf <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean_distance) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#combine empirical and null
jaccard_similarity_islands_empirical_and_null <- rbind(empirical_turnover_for_modules_layers_shuf, ave_module_layer_turnover_shuf_pols,
                                                       ave_module_layer_turnover_shuf_plants, ave_module_layer_turnover_shuf_both)

jaccard_similarity_layer_empirical_and_null_km <- jaccard_similarity_islands_empirical_and_null %>% 
  mutate(mean_dist_in_km = mean_distance/1000)

#write.csv(jaccard_similarity_layer_empirical_and_null_km,  "./csvs_nuevo/jaccard_similarity_layer_empirical_and_null_km_islands_m1.csv", row.names = FALSE)

#---- graphs for distance decay in modules shuf vs empirical--------------------------------
jaccard_similarity_layer_empirical_and_null_km <- read.csv("./csvs_nuevo/jaccard_similarity_layer_empirical_and_null_km_islands_m1.csv")

pdf('./graphs/M1_Modules_DD_Islands.pdf', 10, 6)
jaccard_similarity_layer_empirical_and_null_km %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+
  labs(x="Distance (Km)", y="Jaccard Similarity")+  
  scale_color_manual(name = "Null Model",  labels = c("E",expression("M"[1]^AP), expression("M"[1]^P),
                     expression("M"[1]^A)),values = c("#FB3B1E", "#15B7BC", 
                                                                                        "#72A323", "#BE75FA" )) +

  theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=15, color='black'),
        axis.title = element_text(size=17, color='black'),
        axis.line = element_blank(),
        legend.text.align = 0,
        legend.title =  element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11))

dev.off()

#------check if its significant for islands----------------------------------------------------------------
jaccard_similarity_layer_empirical_and_null_km <- read.csv("./csvs_nuevo/jaccard_similarity_layer_empirical_and_null_km_islands_m1.csv")

## check normality
pol<-jaccard_similarity_layer_empirical_and_null_km %>% filter(type=="null_pollinators")
pla<-jaccard_similarity_layer_empirical_and_null_km %>% filter(type=="null_plants")
both<-jaccard_similarity_layer_empirical_and_null_km %>% filter(type=="null_both")

shapiro.test(pol$ave)
shapiro.test(pla$ave)
shapiro.test(both$ave)

##models
m_pol<-MRM(ave ~ mean_dist_in_km,data=pol,nperm=9999 )
m_pol
m_pla<-MRM(ave ~ mean_dist_in_km,data=pla,nperm=9999 )
m_pla
m_both<-MRM(ave ~ mean_dist_in_km,data=both,nperm=9999 )
m_both


##correlation and r squared between jaccard and distance for each run ----------------------------------------------------------------------------

##correlation and r squared between jaccard and distance for each run ----------------------------------------------------------------------------

## -- pols
all_edge_list_layer_combine_no_module_shuf_pols_output <- read.csv("./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_pols_output_islands.csv")

iteration_correlation_pols <- NULL
iteration_correlation_data_pols <- all_edge_list_layer_combine_no_module_shuf_pols_output %>% subset(layer_from != layer_to) 

iteration_correlation_data_pols_km <- iteration_correlation_data_pols %>% 
  mutate(mean_dist_in_km = mean_distance/1000)%>% group_by(trial) %>% 
  distinct(layer_from,layer_to, .keep_all = T)


for (i in 1:1000){
  print(i)
  trial_pols = iteration_correlation_data_pols_km %>% filter(trial == i)
  iteration_correlation_new_pols <- cor.test(trial_pols$turnover, trial_pols$mean_dist_in_km, method = "pearson")
  mrm_val_pols <- MRM(turnover ~ mean_dist_in_km, data = trial_pols)
  iteration_correlation_pols <- rbind(iteration_correlation_pols, tibble(estimate = iteration_correlation_new_pols$estimate, 
                                                                         p_val = iteration_correlation_new_pols$p.value, 
                                                                         statistic = iteration_correlation_new_pols$statistic, 
                                                                         confidence_int_low = iteration_correlation_new_pols$conf.int[1],
                                                                         confidence_int_high = iteration_correlation_new_pols$conf.int[2],
                                                                         slope = mrm_val_pols$coef[2],
                                                                         intercept = mrm_val_pols$coef[1],
                                                                         rsquared = mrm_val_pols$r.squared[1],
                                                                         trial_num = i))
}


#write.csv(iteration_correlation_pols, "./csvs_nuevo/iteration_correlation_pols.csv", row.names = FALSE)
#iteration_correlation_pols <- read.csv("./csvs_nuevo/iteration_correlation_pols.csv")

p_rsquared_pols <- sum(iteration_correlation_pols$rsquared > 0.653)/1000
p_rsquared_pols ## -- distribution of rsquared and  empirical


#----plants
all_edge_list_layer_combine_no_module_shuf_plants_output <- read.csv("./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_plants_output_islands.csv")

iteration_correlation_plants <- NULL
iteration_correlation_data_plants <- all_edge_list_layer_combine_no_module_shuf_plants_output %>% subset(layer_from != layer_to) 

iteration_correlation_data_plants_km <- iteration_correlation_data_plants %>% 
  mutate(mean_dist_in_km = mean_distance/1000)  %>% group_by(trial) %>% 
  distinct(layer_from,layer_to, .keep_all = T)

for (i in 1:1000){
  print(i)
  trial_plants = iteration_correlation_data_plants_km %>% filter(trial == i)
  iteration_correlation_new_plants <- cor.test(trial_plants$turnover, trial_plants$mean_dist_in_km, method = "pearson")
  mrm_val_plants <- MRM(turnover ~ mean_dist_in_km, data = trial_plants)
  iteration_correlation_plants <- rbind(iteration_correlation_plants, tibble(estimate = iteration_correlation_new_plants$estimate, 
                                                                             p_val = iteration_correlation_new_plants$p.value, 
                                                                             statistic = iteration_correlation_new_plants$statistic, 
                                                                             confidence_int_low = iteration_correlation_new_plants$conf.int[1],
                                                                             confidence_int_high = iteration_correlation_new_plants$conf.int[2],
                                                                             slope = mrm_val_plants$coef[2],
                                                                             intercept = mrm_val_plants$coef[1],
                                                                             rsquared = mrm_val_plants$r.squared[1],
                                                                             trial_num = i))
}


#write.csv(iteration_correlation_plants, "./csvs_nuevo/iteration_correlation_plants.csv", row.names = FALSE)
#iteration_correlation_plants <- read.csv("./csvs_nuevo/iteration_correlation_plants.csv")

p_rsquared_plants <- sum(iteration_correlation_plants$rsquared > 0.653)/1000
p_rsquared_plants ## -- distribution of rsquared and add empirica


#---- both
all_edge_list_layer_combine_no_module_shuf_both_output <- read.csv("./csvs_nuevo/all_edge_list_layer_combine_no_module_shuf_both_output_islands.csv")

iteration_correlation_both <- NULL
iteration_correlation_data_both <- all_edge_list_layer_combine_no_module_shuf_both_output %>% subset(layer_from != layer_to) 

iteration_correlation_data_both_km <- iteration_correlation_data_both %>% 
  mutate(mean_dist_in_km = mean_distance/1000)%>% group_by(trial) %>% 
  distinct(layer_from,layer_to, .keep_all = T)


for (i in 1:1000){
  print(i)
  trial_both = iteration_correlation_data_both_km %>% filter(trial == i)
  iteration_correlation_new_both <- cor.test(trial_both$turnover, trial_both$mean_dist_in_km, method = "pearson")
  mrm_val_both <- MRM(turnover ~ mean_dist_in_km, data = trial_both)
  iteration_correlation_both <- rbind(iteration_correlation_both, tibble(estimate = iteration_correlation_new_both$estimate, 
                                                                         p_val = iteration_correlation_new_both$p.value, 
                                                                         statistic = iteration_correlation_new_both$statistic, 
                                                                         confidence_int_low = iteration_correlation_new_both$conf.int[1],
                                                                         confidence_int_high = iteration_correlation_new_both$conf.int[2],
                                                                         slope = mrm_val_both$coef[2],
                                                                         intercept = mrm_val_both$coef[1],
                                                                         rsquared = mrm_val_both$r.squared[1],
                                                                         trial_num = i))
}


write.csv(iteration_correlation_both, "./csvs_nuevo/iteration_correlation_both.csv", row.names = FALSE)
#iteration_correlation_both<- read.csv("./csvs_nuevo/iteration_correlation_both.csv")


p_rsquared_both <- sum(iteration_correlation_both$rsquared > 0.653)/1000
p_rsquared_both ##distribution of rsquared and empirical



## -- all 3 R squared in the same square
iteration_correlation_pols_M1 <- iteration_correlation_pols %>% mutate(type = "shuf_pollinators")
iteration_correlation_plants_M1 <- iteration_correlation_plants %>% mutate(type = "shuf_plants")
iteration_correlation_both_M1 <- iteration_correlation_both %>% mutate(type = "shuf_both")

rqsuares_M1_all <- rbind(iteration_correlation_pols_M1, iteration_correlation_plants_M1, iteration_correlation_both_M1)

#write.csv(rqsuares_M1_all, "./csvs_nuevo/rqsuares_M1_all.csv", row.names = FALSE)
#rqsuares_M1_all <- read.csv("./csvs_nuevo/rqsuares_M1_all.csv")


rqsuares_M1_all$type <- factor(rqsuares_M1_all$type, levels = c("shuf_plants","shuf_pollinators","shuf_both"))


pdf('./graphs/M1_r_squares_module_DD.pdf', 10, 6)
rqsuares_M1_all %>% 
  ggplot(aes(x = rsquared, fill = type))+ 
  geom_density(alpha = 0.5)+ 
  geom_vline(xintercept = 0.653, linetype = "dashed", color = "#FB3B1E")+ #line R squared empirical
  labs(x= expression("R"^2), y="Density")+  
  scale_fill_manual(name = "Null Model",  labels = c(expression("M"[1]^P),expression("M"[1]^A),
                                                     expression("M"[1]^AP)), values = c("#72A323","#A44CD3", "#15B7BC"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=15, color='black'),
        axis.title = element_text(size=17, color='black'),
        axis.line = element_blank(),
        legend.text.align = 0,
        legend.title =  element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11))
dev.off()
