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

#this portion of the code comes to check whether the bias in the original research in the form of species 8 and 9
#changed the results in a significant way

#---- finding which species are in all 14 layers ----------------------------------------------------------------------------------------
dryad_edgelist_complete_ids_del_all_14_from <- NULL #need to build the modules without species found in all 14 layers
dryad_edgelist_complete_ids_del_all_14 <- NULL

for(i in 1:max(dryad_edgelist_complete_ids$node_from)){ #for species in from
  current_node <- dryad_edgelist_complete_ids %>% filter(node_from == i) %>% select(node_from, layer_from) %>% unique()
  node <- dryad_edgelist_complete_ids %>% filter(node_from == i)
  if(nrow(current_node) == 14) next #delete all species found in all 14 layers
  else{
    dryad_edgelist_complete_ids_del_all_14_from <- rbind(dryad_edgelist_complete_ids_del_all_14_from, node) #size is 2876, 2 species deleted
  }
}

for(i in 1:max(dryad_edgelist_complete_ids_del_all_14_from$node_to)){
  current_node <- dryad_edgelist_complete_ids %>% filter(node_to == i) %>% select(node_to, layer_to) %>% unique()
  node <- dryad_edgelist_complete_ids %>% filter(node_to == i)
  if(nrow(current_node) == 14) next #delete all species found in all 14 layers
  else{
    dryad_edgelist_complete_ids_del_all_14 <- rbind(dryad_edgelist_complete_ids_del_all_14, node) #size is 2876, no more species deleted
  }
}

#del species deleted from the network from the physical nodes
which_missing <- anti_join(dryad_edgelist_complete_ids, dryad_edgelist_complete_ids_del_all_14) #nodes 8 and 9 were found in all 14 layers

physical_nodes_del_14 <- physical_nodes %>% subset(node_id != 8 & node_id != 9) #delete the missing values

dryad_edgelist_complete_ids_del_all_14_interactions <-subset(dryad_edgelist_complete_ids_del_all_14, !(node_from %in% which_missing$node_to)) #del all interactions 8 and 9 are a part of

## delete species who were only found in interaction with 8 or 9 from physical nodes
who_del <- subset(physical_nodes_del_14, !(node_id %in% dryad_edgelist_complete_ids_del_all_14_interactions$node_to)) #9 pol species who only interacted with 9 or 8
physical_nodes_del <- anti_join(physical_nodes_del_14, who_del) #only save the species who weren't interacting only with them

#---- ##remove those species and all partners -------------------------------------------------------------------------------------

dryad_intralayer_del <- dryad_intralayer %>% subset(node_from != "Euphorbia.balsamifera..f." & node_from != "Euphorbia.balsamifera..m.")
intralayer_inverted_del <- intralayer_inverted %>% subset(node_to != "Euphorbia.balsamifera..f." & node_to != "Euphorbia.balsamifera..m.")

#plants in from
tot_plant_del <- dryad_intralayer_del %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_del <- dryad_intralayer_del %>% left_join(tot_plant_del) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)

#pols in from
tot_pol_del <- intralayer_inverted_del %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_del <- intralayer_inverted_del %>% left_join(tot_pol_del) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight) 


all_species_all_layers_del <- rbind(tot_plant_del, tot_pol_del) %>% inner_join(physical_nodes_del, by= c("node_from" = "species")) %>% 
  inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>% subset(select = -c(layer_from, node_from, tot, type)) 
#change node names to ids and layer names to ids and remove unwanted columns

classic_layers_turnover_del <- NULL

for (i in (1:13)){
  for (j in ((i+1):14)){
    physical_nodes_in_layer_from <- filter(all_species_all_layers_del, layer_id == i) %>% select(node_id) %>% unlist()
    physical_nodes_in_layer_to <- filter(all_species_all_layers_del, (layer_id == j)) %>% select(node_id) %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
    int_both <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in both layers
    uni_both <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to)
    turnover <- length(int_both)/length(uni_both)
    classic_layers_turnover_del <- rbind(classic_layers_turnover_del, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

classic_layers_turnover_del <- classic_layers_turnover_del %>% unique()

classic_layers_turnover_with_distances_del <- right_join(classic_layers_turnover_del, distances_with_ids, by= c("layer_from", "layer_to"))
classic_layers_turnover_with_distances_del <- na.omit(classic_layers_turnover_with_distances_del) #remove NA and delete layer name

classic_layers_turnover_with_distances_del <- classic_layers_turnover_with_distances_del %>% mutate(distance_in_km=distance_in_meters/1000)

classic_layers_turnover_with_distances_del %>%
  ggplot(aes(x=distance_in_km, y=turnover))+ geom_point()+ theme_classic()+ stat_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in km", y="Jaccard Similarity") 

#---- modularity----------------------------------------------------------------------------------------------------

dryad_multilayer_del_14 <- create_multilayer_object(extended = dryad_edgelist_complete_ids_del_all_14, #taking edge list and returning multilayer network
                                                    nodes = physical_nodes_del, #make sure missing species are not found
                                                    layers = layer_metadata,
                                                    intra_output_extended = T)


#create modules for empirical network
modules_dryad_multilayer_del_14 <- run_infomap_multilayer(dryad_multilayer_del_14, 
                                                          infomap_executable = "../Infomap",
                                                          flow_model = 'directed',
                                                          relax = F, 
                                                          silent = T, 
                                                          trials = 100,
                                                          seed = 497294, 
                                                          temporal_network = F)

#---- similarity analysis------------------------------------------------------------------------------------------------
modules_for_similarity_num_del_14 <- modules_dryad_multilayer_del_14$modules %>% select(module, layer_id) %>% 
  unique() %>% group_by(module) %>% select(module) %>% unique()
modules_for_similarity_del_14 <- modules_dryad_multilayer_del_14$modules %>%
  filter(module %in% modules_for_similarity_num_del_14$module)

#pivot modules
module_pivoted_del_14 <- pivot_by_module(modules_for_similarity_del_14)


## check if modules change when they're deleted
modules_edge_list_del_14 <- NULL

for (i in (1:nrow(module_pivoted_del_14))){ #run the function for each row in the data frame
  modules_edge_list_del_14 <- edge_list_per_module(module_pivoted_del_14[i,], modules_edge_list_del_14) 
  current_module <- rownames(module_pivoted_del_14)[i]
  modules_edge_list_del_14 <- modules_edge_list_del_14 %>% mutate(module = replace_na(module, current_module)) #add module number
}

edge_list_with_distances_del_14 <- right_join(modules_edge_list_del_14, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances_del_14 <- na.omit(edge_list_with_distances_del_14) #remove NA and delete layer name

#arrange data to include coordinates and modules sizes
size_del_14 <- count(modules_dryad_multilayer_del_14$modules, module)  #create a data frame of all modules and how many nodes are in each (size of module)
module_data_del_14 <- merge(modules_dryad_multilayer_del_14$modules , size_del_14, by=c("module","module")) #merge size of module with all the other info about the modules
colnames(module_data_del_14)[7] <- "size_of_module" #rename column

module_data_with_loc_del_14 <- merge(module_data_del_14, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon_del_14 <- module_data_with_loc_del_14 %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon_del_14$count <- c(1)
##--------------------------------------------------------------------------------------------------------------------------------
#how many islands are within a module
modules_with_lat_lon_islands_del_14 <- module_data_with_loc_del_14 %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
modules_with_lat_lon_islands_del_14$layer_id[modules_with_lat_lon_islands_del_14$layer_id %in% old] <- new[match(modules_with_lat_lon_islands_del_14$layer_id, old)]
modules_with_lat_lon_islands_del_14$count <- c(1)

modules_with_lat_lon_islands_del_14 %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_id)))+ geom_bar(stat= "identity")+ theme_classic()+
  scale_x_continuous(breaks=seq(1,66,2))+ labs(y="number of physical nodes", x="module number")+ guides(fill=guide_legend(title="island\nnumber"))

#modules similarity pairwise distance between islands
edge_list_by_islands_del_14 <- edge_list_with_distances_del_14
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
edge_list_by_islands_del_14$layer_from[edge_list_by_islands_del_14$layer_from %in% old] <- new[match(edge_list_by_islands_del_14$layer_from, old)]
edge_list_by_islands_del_14$layer_to[edge_list_by_islands_del_14$layer_to %in% old] <- new[match(edge_list_by_islands_del_14$layer_to, old)]

#version with # of modules in layers
edge_list_by_islands_modules_del_14 <- edge_list_by_islands_del_14 %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_islands_modules_del_14$count <- c(1)
edge_list_by_islands_modules_del_14 <- edge_list_by_islands_modules_del_14 %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules) 

#version with correct average between layers
edge_list_by_islands_ave_del_14 <- edge_list_by_islands_del_14 %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_island_combine_del_14 <- edge_list_by_islands_ave_del_14 %>%
  merge(edge_list_by_islands_modules_del_14, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_island_combine_no_module_del_14 <- edge_list_island_combine_del_14 %>% select(-module) %>% unique() #have version where modules aren't present


#---- similarity for layers-----------------------------------------------------------
#version with # of modules in layers
edge_list_by_layers_modules_del_14 <- edge_list_with_distances_del_14 %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_layers_modules_del_14$count <- c(1)
edge_list_by_layers_modules_del_14 <- edge_list_by_layers_modules_del_14 %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules) 

#version with correct average between layers
edge_list_by_layers_ave_del_14 <- edge_list_with_distances_del_14 %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_layer_combine_del_14 <- edge_list_by_layers_ave_del_14 %>%
  merge(edge_list_by_layers_modules_del_14, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_layer_combine_no_module_del_14 <- edge_list_layer_combine_del_14 %>% select(-module) %>% unique() #have version where modules aren't present


module_island_turnover_del_14 <- NULL

island_list <- c("1","2","3","4east","4west","5","6")

for (i in island_list){
  for (j in island_list){
    print(i)
    modules_in_island_from <- filter(edge_list_by_islands_modules_del_14, layer_from == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_island_to <- filter(edge_list_by_islands_modules_del_14, layer_from == j) %>% select(module) %>% unique() %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
    int_both <- intersect(modules_in_island_from, modules_in_island_to) #how many nodes are found in both layers
    uni_both <- union(modules_in_island_from, modules_in_island_to)
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover_del_14 <- rbind(module_island_turnover_del_14, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

edge_list_by_islands_ave_del_14 <- edge_list_by_islands_del_14 %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

islands_turnover_with_distnace_empirical_del_14 <- edge_list_by_islands_ave_del_14 %>%
  merge(module_island_turnover_del_14, by= c("layer_from", "layer_to")) #merge both versions 


##same but with layers
module_layer_turnover_del_14 <- NULL

for (i in 1:14){
  for (j in 1:14){
    print(i)
    modules_in_layer_from <- filter(edge_list_by_layers_modules_del_14, layer_from == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_layer_to <- filter(edge_list_by_layers_modules_del_14, layer_from == j) %>% select(module) %>% unique() %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
    int_both <- intersect(modules_in_layer_from, modules_in_layer_to) #how many nodes are found in both layers
    uni_both <- union(modules_in_layer_from, modules_in_layer_to)
    turnover <- length(int_both)/length(uni_both)
    module_layer_turnover_del_14 <- rbind(module_layer_turnover_del_14, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

edge_list_by_layers_ave_del_14 <- edge_list_with_distances_del_14 %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

layers_turnover_with_distnace_empirical_del_14 <- edge_list_by_layers_ave_del_14 %>%
  merge(module_layer_turnover_del_14, by= c("layer_from", "layer_to")) #merge both versions 

#---- empirical-------------------------------------------------------------------------------------------------------------------
#add the empirical
empirical_turnover_for_module_layer_del_14 <- layers_turnover_with_distnace_empirical_del_14 %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null


empirical_turnover_for_module_layer_shuf_no_self_loop_del_14 <- empirical_turnover_for_module_layer_del_14 %>% subset(layer_from != layer_to) #for empirical only distance decay graph

empirical_turnover_for_module_layer_shuf_no_self_loop_km_del_14 <- empirical_turnover_for_module_layer_shuf_no_self_loop_del_14 %>%
  mutate(ave_dist_in_km = ave_dist/1000)

empirical_turnover_for_module_layer_shuf_no_self_loop_km_del_14 %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="distance in km", y="Jaccard Similarity")

#---- null models where species are shuffled between laters--------------------------------------------------------------------------------------------------------------------
#---- shuffling for pols
interlayer_edges_pols_shuf_del_14 <- NULL
interlayer_edges_plants_shuf_del_14 <- NULL
interlayer_edges_both_shuf_del_14 <- NULL

dryad_intralayer_del_14 <- dryad_intralayer %>% subset(node_from != "Euphorbia.balsamifera..f.", 
                                                       node_from != "Euphorbia.balsamifera..m.") #delete the two plant species from the data

intralayers_with_ids <- 
  dryad_intralayer_del_14 %>% 
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


intralayer_matrix <- intralayers_with_ids %>%  #turn edge list to matrix
  select(from, to, weight) %>%
  dcast(from ~ to, value.var = "weight", fill = 0) %>%
  column_to_rownames(var="from")

intralayer_matrix <- t(intralayer_matrix) #change plants and pol- shuffling pols 

shuf_trial_matrix_del <- NULL
shuf_null_edge_list_del <- NULL


for (i in 1:1000){ 
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
    if (wanted_value_layer_1 %in% rownames(trial)) next #make sure we dont already have the wanted combination in the layer
    if (wanted_value_layer_2 %in% rownames(trial)) next
    
    
    else{
      #if we dont have the combination:
      rownames(trial)[random_num_1] <- wanted_value_layer_1 #change the species to its new layer 
      rownames(trial)[random_num_2] <- wanted_value_layer_2 
    }
  }
  trial_with_shuf <- cbind(trial, i) #add trial number 
  shuf_trial_matrix_del <- rbind(shuf_trial_matrix_del, trial_with_shuf) #create big matrix of all matrices
  edge_list_version <- melt(as.matrix(trial)) %>% filter(value > 0) %>%
    select(from=Var1, to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
  edge_list_version$trial_number <- i #add trial number to the edge list
  shuf_null_edge_list_del <- rbind(shuf_null_edge_list_del, edge_list_version) #create mega edge list with all repetitions
}

shuf_null_edge_list_del$layer_from <- substr(shuf_null_edge_list_del$from, 1,2) #choose the layer portion of the mixed layer_node name
shuf_null_edge_list_del$node_from <- substr(shuf_null_edge_list_del$from, 5,8) #choose the node portion of the mixed layer_node name
shuf_null_edge_list_del$layer_to <- substr(shuf_null_edge_list_del$to, 1,2) #choose the layer portion of the mixed layer_node name
shuf_null_edge_list_del$node_to <- substr(shuf_null_edge_list_del$to, 5,8) #choose the node portion of the mixed layer_node name

shuf_null_edge_list_del <- shuf_null_edge_list_del %>% 
  select(layer_from, node_from, layer_to, node_to, weight, trial_number)

shuf_null_edge_list_del$layer_from <- as.numeric(shuf_null_edge_list_del$layer_from)
shuf_null_edge_list_del$layer_to <- as.numeric(shuf_null_edge_list_del$layer_to)
shuf_null_edge_list_del$node_from <- as.numeric(shuf_null_edge_list_del$node_from)
shuf_null_edge_list_del$node_to <- as.numeric(shuf_null_edge_list_del$node_to)

write.csv(shuf_null_edge_list_del, "./csvs/shuf_null_edge_list_del.csv", row.names = FALSE)
write.csv(shuf_trial_matrix_del, "./csvs/shuf_trial_matrix_del.csv", row.names = TRUE)
#shuf_null_edge_list <- read.csv("./csvs/shuf_null_edge_list.csv")
#shuf_trial_matrix <- read.csv("./csvs/shuf_trial_matrix.csv")

# 1-39 to
# 40-288 from

#interlayer edges
interlayer_edges_from_shuf_del <- shuf_null_edge_list_del %>% group_by(trial_number, node_from) %>%
  select(trial_number, layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3], 
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in


interlayer_edges_to_shuf_del <- shuf_null_edge_list_del %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_pols_shuf_del <- rbind(interlayer_edges_from_shuf_del, interlayer_edges_to_shuf_del) 


write.csv(interlayer_edges_pols_shuf_del, "./csvs/interlayer_edges_pols_shuf_del.csv", row.names = FALSE)

#---- shuffling for plants
intralayer_matrix_plants <- intralayer_matrix #change plants and pol- shuffling pols

shuf_trial_matrix_plants_del <- NULL
shuf_null_edge_list_plants_del <- NULL


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
  shuf_trial_matrix_plants_del <- rbind(shuf_trial_matrix_plants_del, trial_with_shuf_plants) #create big matrix of all matrices
  edge_list_version <- melt(as.matrix(trial)) %>% filter(value > 0) %>%
    select(from=Var1, to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
  edge_list_version$trial_number <- i #add trial number to the edge list
  shuf_null_edge_list_plants_del <- rbind(shuf_null_edge_list_plants_del, edge_list_version) #create mega edge list with all repetitions
}

shuf_null_edge_list_plants_del$layer_from <- substr(shuf_null_edge_list_plants_del$from, 1,2)
shuf_null_edge_list_plants_del$node_from <- substr(shuf_null_edge_list_plants_del$from, 5,8)
shuf_null_edge_list_plants_del$layer_to <- substr(shuf_null_edge_list_plants_del$to, 1,2)
shuf_null_edge_list_plants_del$node_to <- substr(shuf_null_edge_list_plants_del$to, 5,8)

shuf_null_edge_list_plants_del <- shuf_null_edge_list_plants_del %>% 
  select(layer_from, node_from, layer_to, node_to, weight, trial_number)

shuf_null_edge_list_plants_del$layer_from <- as.numeric(shuf_null_edge_list_plants_del$layer_from)
shuf_null_edge_list_plants_del$layer_to <- as.numeric(shuf_null_edge_list_plants_del$layer_to)
shuf_null_edge_list_plants_del$node_from <- as.numeric(shuf_null_edge_list_plants_del$node_from)
shuf_null_edge_list_plants_del$node_to <- as.numeric(shuf_null_edge_list_plants_del$node_to)

write.csv(shuf_null_edge_list_plants_del, "./csvs/shuf_null_edge_list_plants_del.csv", row.names = FALSE)
write.csv(shuf_trial_matrix_plants_del, "./csvs/shuf_trial_matrix_plants_del.csv", row.names = TRUE)
#shuf_null_edge_list_plants <- read.csv("./csvs/shuf_null_edge_list_plants.csv")
#shuf_trial_matrix_plants <- read.csv("./csvs/shuf_trial_matrix_plants.csv")

#interlayer edges
interlayer_edges_from_shuf_plants_del <- shuf_null_edge_list_plants_del %>% group_by(trial_number, node_from) %>% 
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3],
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in

interlayer_edges_to_shuf_plants_del <- shuf_null_edge_list_plants_del %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_plants_shuf_del <- rbind(interlayer_edges_from_shuf_plants_del, interlayer_edges_to_shuf_plants_del) 


write.csv(interlayer_edges_plants_shuf_del, "./csvs/interlayer_edges_plants_shuf_del.csv", row.names = FALSE)

#---- shuffling for both

intralayer_matrix_both <- shuf_trial_matrix_del #copying and now we'll change columns and not rows



shuf_trial_matrix_both_del <- NULL
shuf_null_edge_list_both_del <- NULL

for (i in 1:1000){ 
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
  shuf_trial_matrix_both_del <- rbind(shuf_trial_matrix_both_del, trial_with_shuf_both) #create big matrix of all matrices
  edge_list_version <- melt(as.matrix(trial)) %>% filter(value > 0) %>%
    select(from=Var1, to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
  edge_list_version$trial_number <- i #add trial number to the edge list
  shuf_null_edge_list_both_del <- rbind(shuf_null_edge_list_both_del, edge_list_version) #create mega edge list with all repetitions
}


shuf_null_edge_list_both_del$layer_from <- substr(shuf_null_edge_list_both_del$from, 1,2)
shuf_null_edge_list_both_del$node_from <- substr(shuf_null_edge_list_both_del$from, 5,8)
shuf_null_edge_list_both_del$layer_to <- substr(shuf_null_edge_list_both_del$to, 1,2)
shuf_null_edge_list_both_del$node_to <- substr(shuf_null_edge_list_both_del$to, 5,8)

shuf_null_edge_list_both_del <- shuf_null_edge_list_both_del %>% 
  select(layer_from, node_from, layer_to, node_to, weight, trial_number)

shuf_null_edge_list_both_del$layer_from <- as.numeric(shuf_null_edge_list_both_del$layer_from)
shuf_null_edge_list_both_del$layer_to <- as.numeric(shuf_null_edge_list_both_del$layer_to)
shuf_null_edge_list_both_del$node_from <- as.numeric(shuf_null_edge_list_both_del$node_from)
shuf_null_edge_list_both_del$node_to <- as.numeric(shuf_null_edge_list_both_del$node_to)

shuf_null_edge_list_both_del <- subset(shuf_null_edge_list_both_del, layer_to != "i") #delete every row where "i" is found in layer_to
shuf_null_edge_list_both_del <- subset(shuf_null_edge_list_both_del, layer_from != "i") #delete every row where "i" is found in layer_from 

write.csv(shuf_null_edge_list_both_del, "./csvs/shuf_null_edge_list_both_del.csv", row.names = FALSE)
write.csv(shuf_trial_matrix_both_del, "./csvs/shuf_trial_matrix_both_del.csv", row.names = TRUE)


#interlayer edges
interlayer_edges_from_shuf_both_del <- shuf_null_edge_list_both_del %>% group_by(trial_number, node_from) %>% 
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3],
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in


interlayer_edges_to_shuf_both_del <- shuf_null_edge_list_both_del %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_both_shuf_del <- rbind(interlayer_edges_from_shuf_both_del, interlayer_edges_to_shuf_both_del) 


write.csv(interlayer_edges_both_shuf_del, "./csvs/interlayer_edges_both_shuf_del.csv", row.names = FALSE)

#---- write results for all versions and HPC organization-------------------------------------------
write.csv(interlayer_edges_pols_shuf_del, "./HPC/interlayer_edges_pols_shuf_del.csv", row.names = FALSE) #create to run on HPC
write.csv(interlayer_edges_plants_shuf_del, "./HPC/interlayer_edges_plants_shuf_del.csv", row.names = FALSE) #create to run on HPC
write.csv(interlayer_edges_both_shuf_del, "./HPC/interlayer_edges_both_shuf_del.csv", row.names = FALSE) #create to run on HPC


## run of HPC and come back 

#both
files_both_del<- list.files("./HPC/csvs_both_del/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_both_del <- read_csv(files_both_del) %>% bind_rows() #create a long edge list with all the csvs

#pollinators
files_pollinators_del <- list.files("./HPC/csvs_pollinators_del/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_pol_del <- read_csv(files_pollinators_del) %>% bind_rows() #create a long edge list with all the csvs

#plants
files_plants_del <- list.files("./HPC/csvs_plants_del/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_plants_del <- read_csv(files_plants_del) %>% bind_rows() #create a long edge list with all the csvs

write.csv(my_merged_interlayer_shuf_both_del, "./csvs/my_merged_interlayer_shuf_both_del.csv", row.names = FALSE)
write.csv(my_merged_interlayer_shuf_pol_del, "./csvs/my_merged_interlayer_shuf_pol_del.csv", row.names = FALSE)
write.csv(my_merged_interlayer_shuf_plants_del, "./csvs/my_merged_interlayer_shuf_plants_del.csv", row.names = FALSE)

distances_with_weights_ids_del <- distances_with_weights %>%
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, layer_to=layer_id.y, weight)


#pols
interlayers_with_weights_shuf_pols_del <- my_merged_interlayer_shuf_pol_del %>% inner_join(distances_with_weights_ids_del, 
                                                                                           by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights_shuf_pols_del <- interlayers_with_weights_shuf_pols_del[!duplicated(interlayers_with_weights_shuf_pols_del[c(1,3,5,6)]),]

#plants
interlayers_with_weights_shuf_plants_del <- my_merged_interlayer_shuf_plants_del %>% inner_join(distances_with_weights_ids_del, 
                                                                                                by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights_shuf_plants_del <- interlayers_with_weights_shuf_plants_del[!duplicated(interlayers_with_weights_shuf_plants_del[c(1,3,5,6)]),]

#both
interlayers_with_weights_shuf_both_del <- my_merged_interlayer_shuf_both_del %>% inner_join(distances_with_weights_ids_del, 
                                                                                            by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights_shuf_both_del <- interlayers_with_weights_shuf_both_del[!duplicated(interlayers_with_weights_shuf_both_del[c(1,3,5,6)]),]

#write.csv(interlayers_with_weights_shuf_pols_del, "./csvs/interlayers_with_weights_shuf_pols_del.csv", row.names = FALSE)
#write.csv(interlayers_with_weights_shuf_plants_del, "./csvs/interlayers_with_weights_shuf_plants_del.csv", row.names = FALSE)
#write.csv(interlayers_with_weights_shuf_both_del, "./csvs/interlayers_with_weights_shuf_both_del.csv", row.names = FALSE)

#create inter and intra for the 1000 shuf trials
dryad_interlayer_shuf_pols_del <- read.csv("./csvs/interlayers_with_weights_shuf_pols_del.csv") #already has inverted
dryad_interlayer_shuf_plants_del <- read.csv("./csvs/interlayers_with_weights_shuf_plants_del.csv") #already has inverted
dryad_interlayer_shuf_both_del <- read.csv("./csvs/interlayers_with_weights_shuf_both_del.csv") #already has inverted


#pols
dryad_intralayer_shuf_pols_del <- read.csv("./csvs/shuf_null_edge_list_del.csv") 
dryad_intralayer_shuf_pols_del <- dryad_intralayer_shuf_pols_del[, c(6,1,2,3,4,5)]

#plants
dryad_intralayer_shuf_plants_del <- read.csv("./csvs/shuf_null_edge_list_plants_del.csv") 
dryad_intralayer_shuf_plants_del <- dryad_intralayer_shuf_plants_del[, c(6,1,2,3,4,5)]

#both
dryad_intralayer_shuf_both_del <- read.csv("./csvs/shuf_null_edge_list_both_del.csv") 
dryad_intralayer_shuf_both_del <- dryad_intralayer_shuf_both_del[, c(6,1,2,3,4,5)]

#create inverted versions
#pols
intralayer_inverted_shuf_pols_del <- tibble(values= dryad_intralayer_shuf_pols_del$layer_to, dryad_intralayer_shuf_pols_del$node_to, 
                                            dryad_intralayer_shuf_pols_del$layer_from, dryad_intralayer_shuf_pols_del$node_from, 
                                            dryad_intralayer_shuf_pols_del$weight, dryad_intralayer_shuf_pols_del$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_pols_del) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")


#plants
intralayer_inverted_shuf_plants_del <- tibble(values= dryad_intralayer_shuf_plants_del$layer_to, dryad_intralayer_shuf_plants_del$node_to, 
                                              dryad_intralayer_shuf_plants_del$layer_from, dryad_intralayer_shuf_plants_del$node_from, 
                                              dryad_intralayer_shuf_plants_del$weight, dryad_intralayer_shuf_plants_del$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_plants_del) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")

#both
intralayer_inverted_shuf_both_del <- tibble(values= dryad_intralayer_shuf_both_del$layer_to, dryad_intralayer_shuf_both_del$node_to, 
                                            dryad_intralayer_shuf_both_del$layer_from, dryad_intralayer_shuf_both_del$node_from, 
                                            dryad_intralayer_shuf_both_del$weight, dryad_intralayer_shuf_both_del$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_both_del) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")

# ----normalize intralayer edge weights--------------------------------------------------------------------------
##pols
#plants in from
tot_plant_shuf_pols_del <- intralayer_inverted_shuf_pols_del %>% #if it doesnt work need to change it to dryad_intralayer_shuf
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_pols_del <- intralayer_inverted_shuf_pols_del %>% left_join(tot_plant_shuf_pols_del) %>% 
  mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_pols_del <- tot_plant_shuf_pols_del[, c(3,1,2,4)]

#pols in from
tot_pol_shuf_pols_del <- dryad_intralayer_shuf_pols_del %>% #if it doesnt work need to change it to intralayer_inverted_shuf
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_pols_del <- dryad_intralayer_shuf_pols_del %>% left_join(tot_pol_shuf_pols_del) %>% 
  mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_pols_del <- tot_pol_shuf_pols_del[, c(3,1,2,4)]

##plants
#plants in from
tot_plant_shuf_plants_del <- intralayer_inverted_shuf_plants_del %>% #if it doesnt work need to change it to dryad_intralayer_shuf
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_plants_del <- intralayer_inverted_shuf_plants_del %>% left_join(tot_plant_shuf_plants_del) %>% 
  mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_plants_del <- tot_plant_shuf_plants_del[, c(3,1,2,4)]

#pols in from
tot_pol_shuf_plants_del <- dryad_intralayer_shuf_plants_del %>% #if it doesnt work need to change it to intralayer_inverted_shuf
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_plants_del <- dryad_intralayer_shuf_plants_del %>% left_join(tot_pol_shuf_plants_del) %>% 
  mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_plants_del <- tot_pol_shuf_plants_del[, c(3,1,2,4)]

##both
#plants in from
tot_plant_shuf_both_del <- intralayer_inverted_shuf_both_del %>% #if it doesnt work need to change it to dryad_intralayer_shuf
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_both_del <- intralayer_inverted_shuf_both_del %>% left_join(tot_plant_shuf_both_del) %>% 
  mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_both_del <- tot_plant_shuf_both_del[, c(3,1,2,4)]


#pols in from
tot_pol_shuf_both_del <- dryad_intralayer_shuf_both_del %>% #if it doesnt work need to change it to intralayer_inverted_shuf
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_both_del <- dryad_intralayer_shuf_both_del %>% left_join(tot_pol_shuf_both_del) %>% 
  mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_both_del <- tot_pol_shuf_both_del[, c(3,1,2,4)]

# ----multilayer_extended_final--------------------------------------------------------------------------------------
#edgelist_non_inverted_shuf <- bind_rows(intralayer_weighted_shuf, dryad_interlayer_shuf) #combine weighted version of intra with inter
#edgelist_inverted_shuf <- bind_rows(intralayer_weighted_inverted_shuf, interlayer_inverted_shuf) #combine inverted version of intra with inverted version of inter

edgelist_intralayer_shuf_pols_del <- bind_rows(intralayer_weighted_shuf_pols_del, intralayer_weighted_inverted_shuf_pols_del)
edgelist_intralayer_shuf_plants_del <- bind_rows(intralayer_weighted_shuf_plants_del, intralayer_weighted_inverted_shuf_plants_del)
edgelist_intralayer_shuf_both_del <- bind_rows(intralayer_weighted_shuf_both_del, intralayer_weighted_inverted_shuf_both_del)

#dryad_edgelist_complete_shuf <- bind_rows(edgelist_non_inverted_shuf, edgelist_inverted_shuf) #combine inverted and non inverted verions
dryad_edgelist_complete_shuf_pols_del <- bind_rows(edgelist_intralayer_shuf_pols_del, dryad_interlayer_shuf_pols_del) #combine inter and intra
dryad_edgelist_complete_shuf_plants_del <- bind_rows(edgelist_intralayer_shuf_plants_del, dryad_interlayer_shuf_plants_del) #combine inter and intra
dryad_edgelist_complete_shuf_both_del <- bind_rows(edgelist_intralayer_shuf_both_del, dryad_interlayer_shuf_both_del) #combine inter and intra


## ----classical distance decay for shuf networks --------------------------------------------------------------------
#pols
all_species_all_layers_all_trials_pols_del <- rbind(tot_plant_shuf_pols_del, tot_pol_shuf_pols_del) %>% subset(select = -c(tot)) %>% 
  rename(node_id = node_from, layer_id = layer_from)

#plants
all_species_all_layers_all_trials_plants_del <- rbind(tot_plant_shuf_plants_del, tot_pol_shuf_plants_del) %>% subset(select = -c(tot)) %>% 
  rename(node_id = node_from, layer_id = layer_from)

#both
all_species_all_layers_all_trials_both_del <- rbind(tot_plant_shuf_both_del, tot_pol_shuf_both_del) %>% subset(select = -c(tot)) %>% 
  rename(node_id = node_from, layer_id = layer_from)

#---- distance decay in species function---------------------------------------------------------

classic_layers_turnover_shuf_pols_del <- NULL
classic_layers_turnover_shuf_plants_del <- NULL
classic_layers_turnover_shuf_both_del <- NULL


#pols
classic_layers_turnover_shuf_output_pols_del <- distnace_decay_shuf(all_species_all_layers_all_trials_pols_del, 
                                                                    classic_layers_turnover_shuf_pols_del)

#plants
classic_layers_turnover_shuf_output_plants_del <- distnace_decay_shuf(all_species_all_layers_all_trials_plants_del, 
                                                                      classic_layers_turnover_shuf_plants_del)

#both
classic_layers_turnover_shuf_output_both_del <- distnace_decay_shuf(all_species_all_layers_all_trials_both_del, 
                                                                    classic_layers_turnover_shuf_both_del)

#pols
classic_layers_turnover_with_distances_shuf_pols_del <- right_join(classic_layers_turnover_shuf_output_pols_del, distances_with_ids, by= c("layer_from", "layer_to"))
classic_layers_turnover_with_distances_shuf_pols_del <- na.omit(classic_layers_turnover_with_distances_shuf_pols_del) #remove NA and delete layer name

#plants
classic_layers_turnover_with_distances_shuf_plants_del <- right_join(classic_layers_turnover_shuf_output_plants_del, distances_with_ids, by= c("layer_from", "layer_to"))
classic_layers_turnover_with_distances_shuf_plants_del <- na.omit(classic_layers_turnover_with_distances_shuf_plants_del) #remove NA and delete layer name

#both
classic_layers_turnover_with_distances_shuf_both_del <- right_join(classic_layers_turnover_shuf_output_both_del, distances_with_ids, by= c("layer_from", "layer_to"))
classic_layers_turnover_with_distances_shuf_both_del <- na.omit(classic_layers_turnover_with_distances_shuf_both_del) #remove NA and delete layer name


#write.csv(classic_layers_turnover_with_distances_shuf_pols_del, "./csvs/classic_layers_turnover_with_distances_shuf_pols_del.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_with_distances_shuf_plants_del, "./csvs/classic_layers_turnover_with_distances_shuf_plants_del.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_with_distances_shuf_both_del, "./csvs/classic_layers_turnover_with_distances_shuf_both_del.csv", row.names = FALSE)
classic_layers_turnover_with_distances_shuf_pols_del <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_pols_del.csv")
classic_layers_turnover_with_distances_shuf_plants_del <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_plants_del.csv")
classic_layers_turnover_with_distances_shuf_both_del <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_both_del.csv")

#----create an average for shuf with sd---------------------------------------------------------
#pols
ave_turnover_for_shuf_pols_del <- classic_layers_turnover_with_distances_shuf_pols_del %>% group_by(layer_from, layer_to, distance_in_meters) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_pollinators") #create mean and sd for each point

#plants
ave_turnover_for_shuf_plants_del <- classic_layers_turnover_with_distances_shuf_plants_del %>% group_by(layer_from, layer_to, distance_in_meters) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_plants") #create mean and sd for each point

#both
ave_turnover_for_shuf_both_del <- classic_layers_turnover_with_distances_shuf_both_del %>% group_by(layer_from, layer_to, distance_in_meters) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_both") #create mean and sd for each point

#add the emprical classical turnover
empirical_turnover_for_shuf_del <- classic_layers_turnover_with_distances_del %>% group_by(layer_from, layer_to, distance_in_meters) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#create graph
turnover_shuf_and_empirical_del <- rbind(empirical_turnover_for_shuf_del, ave_turnover_for_shuf_pols_del, ave_turnover_for_shuf_plants_del,
                                         ave_turnover_for_shuf_both_del)

turnover_shuf_and_empirical_del <- turnover_shuf_and_empirical_del %>% mutate(distance_in_km = distance_in_meters/1000)

turnover_shuf_and_empirical_del %>% ggplot(aes(x= distance_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13),legend.text = element_text(size = 13))+
  labs(x="distance in km", y="Jaccard Similarity")

#----making sure its significant-------------------------------------------------------
lm1_del = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical_del,turnover_shuf_and_empirical_del$type=="empirical")) #in empirical
lm2_del = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical_del,turnover_shuf_and_empirical_del$type=="null_pollinators")) #in null pols
lm3_del = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical_del,turnover_shuf_and_empirical_del$type=="null_plants")) #in null plants
lm4_del = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical_del,turnover_shuf_and_empirical_del$type=="null_both")) #in null plants

#get equations
lm1_equation_del <- paste("y=", coef(lm1_del)[[1]], "+", coef(lm1_del)[[2]], "*x")
lm2_equation_del <- paste("y=", coef(lm2_del)[[1]], "+", coef(lm2_del)[[2]], "*x")
lm3_equation_del <- paste("y=", coef(lm3_del)[[1]], "+", coef(lm3_del)[[2]], "*x")
lm4_equation_del <- paste("y=", coef(lm4_del)[[1]], "+", coef(lm4_del)[[2]], "*x")


b1_del <- summary(lm1_del)$coefficients[2,1]
se1_del <- summary(lm1_del)$coefficients[2,2]
b2_del <- summary(lm2_del)$coefficients[2,1]
se2_del <- summary(lm2_del)$coefficients[2,2]
b3_del <- summary(lm3_del)$coefficients[2,1]
se3_del <- summary(lm3_del)$coefficients[2,2]
b4_del <- summary(lm4_del)$coefficients[2,1]
se4_del <- summary(lm4_del)$coefficients[2,2]

p_value_pol_del = 2*pnorm(-abs(compare.coeff(b1_del,se1_del,b2_del,se2_del)))
p_value_pol_del

p_value_plants_del = 2*pnorm(-abs(compare.coeff(b1_del,se1_del,b3_del,se3_del)))
p_value_plants_del

p_value_both_del = 2*pnorm(-abs(compare.coeff(b1_del,se1_del,b4_del,se4_del)))
p_value_both_del

# ----multilayer_class-----------------------------------------------------------------------------------------------
# Input: An extended edge list.
dryad_edgelist_complete_shuf_pols_del <- dryad_edgelist_complete_shuf_pols_del[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_plants_del <- dryad_edgelist_complete_shuf_plants_del[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_both_del <- dryad_edgelist_complete_shuf_both_del[, c(1,2,3,4,6,5)] #make sure weight is #5

who_del_shuf <- subset(physical_nodes_del_14, !(node_id %in% dryad_edgelist_complete_shuf_pols_del$node_to)) #9 pol species who only interacted with 9 or 8
physical_nodes_del_shuf <- anti_join(physical_nodes_del_14, who_del_shuf) #only save the species who weren't interacting only with them


dryad_multilayer_shuf_1000_pols_del <- NULL
dryad_multilayer_shuf_1000_plants_del <- NULL
dryad_multilayer_shuf_1000_both_del <- NULL


#pols
dryad_multilayer_shuf_1000_pols_output_del <- modularity_for_shuf_del(dryad_edgelist_complete_shuf_pols_del, 
                                                                      dryad_multilayer_shuf_1000_pols_del)

#plants
dryad_multilayer_shuf_1000_plants_output_del <- modularity_for_shuf_del(dryad_edgelist_complete_shuf_plants_del, 
                                                                        dryad_multilayer_shuf_1000_plants_del)

#both
dryad_multilayer_shuf_1000_both_output_del <- modularity_for_shuf_del(dryad_edgelist_complete_shuf_both_del, 
                                                                      dryad_multilayer_shuf_1000_both_del)

write.csv(dryad_multilayer_shuf_1000_pols_output_del, "./csvs/dryad_multilayer_shuf_1000_pols_output_del.csv", row.names = FALSE)
write.csv(dryad_multilayer_shuf_1000_plants_output_del, "./csvs/dryad_multilayer_shuf_1000_plants_output_del.csv", row.names = FALSE)
write.csv(dryad_multilayer_shuf_1000_both_output_del, "./csvs/dryad_multilayer_shuf_1000_both_output_del.csv", row.names = FALSE)


#---- distance decay of modules shuf vs empirical---------------------------------------------------------------------------
all_edge_list_island_combine_no_module_shuf_pols_del <- NULL
all_edge_list_island_combine_no_module_shuf_plants_del <- NULL
all_edge_list_island_combine_no_module_shuf_both_del <- NULL
all_edge_list_layers_combine_no_module_shuf_del <- NULL


islands_turnover_with_distnace_pols_del <- NULL
islands_turnover_with_distnace_plants_del <- NULL
islands_turnover_with_distnace_both_del <- NULL

module_island_turnover_shuf_del <- NULL

dryad_multilayer_shuf_1000_pols_output_del <- as.data.frame(dryad_multilayer_shuf_1000_pols_output_del)
dryad_multilayer_shuf_1000_plants_output_del <- as.data.frame(dryad_multilayer_shuf_1000_plants_output_del) 
dryad_multilayer_shuf_1000_both_output_del <- as.data.frame(dryad_multilayer_shuf_1000_both_output_del) 

#---- distance decay in modules islands------------------------------------------------------------------------------------------------
## islands
#pols
all_edge_list_island_combine_no_module_shuf_pols_output_del <- module_distance_decay_func_del(dryad_multilayer_shuf_1000_pols_output_del,
                                                                                              islands_turnover_with_distnace_pols_del)


write.csv(all_edge_list_island_combine_no_module_shuf_pols_output_del, "./csvs/all_edge_list_island_combine_no_module_shuf_pols_output_del.csv", 
          row.names = FALSE)

all_edge_list_island_combine_no_module_shuf_pols_output_del <- read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_pols_output_del.csv")

#plants
all_edge_list_island_combine_no_module_shuf_plants_output_del <- module_distance_decay_func_del(dryad_multilayer_shuf_1000_plants_output_del,
                                                                                                islands_turnover_with_distnace_plants_del)

write.csv(all_edge_list_island_combine_no_module_shuf_plants_output_del, "./csvs/all_edge_list_island_combine_no_module_shuf_plants_output_del.csv", 
          row.names = FALSE)

all_edge_list_island_combine_no_module_shuf_plants_output_del <- read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_plants_output_del.csv")

#both
all_edge_list_island_combine_no_module_shuf_both_output_del <- module_distance_decay_func_del(dryad_multilayer_shuf_1000_both_output_del,
                                                                                              islands_turnover_with_distnace_both_del)

write.csv(all_edge_list_island_combine_no_module_shuf_both_output_del, "./csvs/all_edge_list_island_combine_no_module_shuf_both_output_del.csv", 
          row.names = FALSE)

all_edge_list_island_combine_no_module_shuf_both_output_del <- read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_both_output_del.csv")

#----distance decay in modules layers ---------------------------------------------------------------------------------------------
layers_turnover_with_distnace_pols_del <- NULL
layers_turnover_with_distnace_plants_del <- NULL
layers_turnover_with_distnace_both_del <- NULL

module_layer_turnover_shuf_del <- NULL


##layers
all_edge_list_layer_combine_no_module_shuf_pols_output_del <- module_distance_decay_layer_func(dryad_multilayer_shuf_1000_pols_output_del,
                                                                                               layers_turnover_with_distnace_pols_del)


write.csv(all_edge_list_layer_combine_no_module_shuf_pols_output_del, "./csvs/all_edge_list_layer_combine_no_module_shuf_pols_output_del.csv", 
          row.names = FALSE)

all_edge_list_layer_combine_no_module_shuf_pols_output_del <- read.csv("./csvs/all_edge_list_layer_combine_no_module_shuf_pols_output_del.csv")

#plants
all_edge_list_layer_combine_no_module_shuf_plants_output_del <- module_distance_decay_layer_func(dryad_multilayer_shuf_1000_plants_output_del,
                                                                                                 layers_turnover_with_distnace_plants_del)

write.csv(all_edge_list_layer_combine_no_module_shuf_plants_output_del, "./csvs/all_edge_list_layer_combine_no_module_shuf_plants_output_del.csv", 
          row.names = FALSE)

all_edge_list_layer_combine_no_module_shuf_plants_output_del <- read.csv("./csvs/all_edge_list_layer_combine_no_module_shuf_plants_output_del.csv")

#both
all_edge_list_layer_combine_no_module_shuf_both_output_del <- module_distance_decay_layer_func(dryad_multilayer_shuf_1000_both_output_del,
                                                                                               layers_turnover_with_distnace_both_del)

write.csv(all_edge_list_layer_combine_no_module_shuf_both_output_del, "./csvs/all_edge_list_layer_combine_no_module_shuf_both_output_del.csv", 
          row.names = FALSE)

all_edge_list_layer_combine_no_module_shuf_both_output_del <- read.csv("./csvs/all_edge_list_layer_combine_no_module_shuf_both_output_del.csv")

#---- create ave for jaccard islands ----------------------------------------------------------------------------------------
#pols
ave_module_island_turnover_shuf_pols_del <- all_edge_list_island_combine_no_module_shuf_pols_output_del %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_pollinators") #create mean and sd for each point

#plants
ave_module_island_turnover_shuf_plants_del <- all_edge_list_island_combine_no_module_shuf_plants_output_del %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_plants") #create mean and sd for each point

#both
ave_module_island_turnover_shuf_both_del <- all_edge_list_island_combine_no_module_shuf_both_output_del %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_both") #create mean and sd for each point

#add the empirical empirical
empirical_turnover_for_module_island_shuf_del <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

empirical_turnover_for_module_island_shuf_no_self_loop_del <- empirical_turnover_for_module_island_shuf_del %>% subset(layer_from != layer_to) #for empirical only distance decay graph

empirical_turnover_for_module_island_shuf_no_self_loop_km_del <- empirical_turnover_for_module_island_shuf_no_self_loop_del %>%
  mutate(ave_dist_in_km = ave_dist/1000)

#---- create ave for jaccard layers ----------------------------------------------------------------------------------------------
#pols
ave_module_layer_turnover_shuf_pols_del <- all_edge_list_layer_combine_no_module_shuf_pols_output_del %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_pollinators") #create mean and sd for each point

#plants
ave_module_layer_turnover_shuf_plants_del <- all_edge_list_layer_combine_no_module_shuf_plants_output_del %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_plants") #create mean and sd for each point

#both
ave_module_layer_turnover_shuf_both_del <- all_edge_list_layer_combine_no_module_shuf_both_output_del %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_both") #create mean and sd for each point

#add the empirical empirical
empirical_turnover_for_module_layer_shuf_del <- layers_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null


empirical_turnover_for_module_layer_shuf_no_self_loop_del <- empirical_turnover_for_module_layer_shuf_del %>% subset(layer_from != layer_to) #for empirical only distance decay graph

empirical_turnover_for_module_layer_shuf_no_self_loop_km_del <- empirical_turnover_for_module_layer_shuf_no_self_loop_del %>%
  mutate(ave_dist_in_km = ave_dist/1000)

#---- combine---------------------------------------------------------------------------------------------
#combine all islands
jaccard_similarity_empirical_and_null_del <- rbind(empirical_turnover_for_module_island_shuf_del, ave_module_island_turnover_shuf_pols_del,
                                                   ave_module_island_turnover_shuf_plants_del, ave_module_island_turnover_shuf_both_del)

jaccard_similarity_empirical_and_null_no_self_loop_del <- jaccard_similarity_empirical_and_null_del %>% subset(layer_from != layer_to)

jaccard_similarity_empirical_and_null_no_self_loop_km_del <- jaccard_similarity_empirical_and_null_no_self_loop_del %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#combine all layers
jaccard_similarity_layer_empirical_and_null_del <- rbind(empirical_turnover_for_module_layer_shuf_del, ave_module_layer_turnover_shuf_pols_del,
                                                         ave_module_layer_turnover_shuf_plants_del, ave_module_layer_turnover_shuf_both_del)

jaccard_similarity_layer_empirical_and_null_no_self_loop_del <- jaccard_similarity_layer_empirical_and_null_del %>% subset(layer_from != layer_to)

jaccard_similarity_layer_empirical_and_null_no_self_loop_km_del <- jaccard_similarity_layer_empirical_and_null_no_self_loop_del %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#---- graphs distance decay----------------------------------------------------
#no self loop island
#just emprical
empirical_turnover_for_module_island_shuf_no_self_loop_km_del %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")

#empirical and null
jaccard_similarity_empirical_and_null_no_self_loop_km_del %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")


#no self loop layer
#just emprical
empirical_turnover_for_module_layer_shuf_no_self_loop_km_del %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")

#empirical and null
jaccard_similarity_layer_empirical_and_null_no_self_loop_km_del %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")


#---- check if its significant-----------------------------------------
lm1_module_del = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop_del,
                                                jaccard_similarity_empirical_and_null_no_self_loop_del$type=="empirical")) #in empirical
lm2_module_del = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop_del,
                                                jaccard_similarity_empirical_and_null_no_self_loop_del$type=="null_pollinators")) #in null pols
lm3_module_del = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop_del,
                                                jaccard_similarity_empirical_and_null_no_self_loop_del$type=="null_plants")) #in null plants
lm4_module_del = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop_del,
                                                jaccard_similarity_empirical_and_null_no_self_loop_del$type=="null_both")) #in null both

#get equations
lm1_module_equation_del <- paste("y=", coef(lm1_module_del)[[1]], "+", coef(lm1_module_del)[[2]], "*x")
lm2_module_equation_del <- paste("y=", coef(lm2_module_del)[[1]], "+", coef(lm2_module_del)[[2]], "*x")
lm3_module_equation_del <- paste("y=", coef(lm3_module_del)[[1]], "+", coef(lm3_module_del)[[2]], "*x")
lm4_module_equation_del <- paste("y=", coef(lm4_module_del)[[1]], "+", coef(lm4_module_del)[[2]], "*x")

b1_module_del <- summary(lm1_module_del)$coefficients[2,1]
se1_module_del <- summary(lm1_module_del)$coefficients[2,2]
b2_module_del <- summary(lm2_module_del)$coefficients[2,1]
se2_module_del <- summary(lm2_module_del)$coefficients[2,2]
b3_module_del <- summary(lm3_module_del)$coefficients[2,1]
se3_module_del <- summary(lm3_module_del)$coefficients[2,2]
b4_module_del <- summary(lm4_module_del)$coefficients[2,1]
se4_module_del <- summary(lm4_module_del)$coefficients[2,2]

p_value_module_pols_del = 2*pnorm(-abs(compare.coeff(b1_module_del,se1_module_del,b2_module_del,se2_module_del)))
p_value_module_pols_del
p_value_module_plants_del = 2*pnorm(-abs(compare.coeff(b1_module_del,se1_module_del,b3_module_del,se3_module_del)))
p_value_module_plants_del
p_value_module_both_del = 2*pnorm(-abs(compare.coeff(b1_module_del,se1_module_del,b4_module_del,se4_module_del)))
p_value_module_both_del
