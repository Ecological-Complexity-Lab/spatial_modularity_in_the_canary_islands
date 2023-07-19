#This code turns the data into a multilayer network with 7 layers 
#(islands as layers). It then compares the distance decay in species and modules
#of the empirical network to null model 1 explained in the ms.

#LINES 12- 289: Create dataframe and exploratory analysis
#LINES 289-351: Distance decay in species - Empirical data
#LINES 351-482: Distance decay in species - Empirical data
#LINES 483- 1551: Null Model 1 (M1): Shuffling pollinator, plant and both species between islands and calculate distance decay.



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
library(reshape2)



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

#-- Total number of interactions per island
#Tot_int<- dryad_intralayer_islands_grouped %>% 
#  group_by(layer_from) %>% summarise(T_int = sum(sum_weight)) 

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

## ----interlayers as islands -----------------------------------------------------------------------
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

distances_normalized <- distances %>% filter(layer_from != layer_to) %>% #delete distances between sites in the same island
  group_by(layer_to, layer_from) %>% #group will contain 4 sites- site 1 and 1 of layer from and site 1 and 2 or layer to
  summarise(mean_distance = mean(distance_in_meters)) %>%unique() #use an average distance of the 4 sites in 2 different islands to determine the distance between the islands

distances_normalized <- distances_normalized[c("layer_from", "layer_to", "mean_distance")]

shortest_distance <- min(distances_normalized$mean_distance) #the shortest distance between islands

interlayer_weight <- function(d){
  #receives distance and normalizes it
  weight <- (1/log(d))/(1/log(shortest_distance))
  return(weight)
}

distances_with_weights <- distances_normalized %>% 
  mutate(weight = interlayer_weight(mean_distance)) %>% #add weight using the function
  subset(select = -mean_distance)

#write.csv(distances_with_weights, "./csvs/Islands/distances_with_weights_islands_as_layers.csv", row.names = FALSE)

interlayers_with_weights_islands <- interlayers_new_islands %>% inner_join(distances_with_weights, 
                                                           by = c("layer_from", "layer_to")) %>% unique()

#write.csv(interlayers_with_weights_islands, "./csvs/Islands/interlayers_with_weights_islands.csv", row.names = FALSE)


## ----multilayer_extended_final--------------------------------------------------------------------------------------
edgelist_intralayers_both <- bind_rows(intralayer_weighted, intralayer_weighted_inverted) #combine inverted and non inverted versions of intra

dryad_edgelist_complete <- bind_rows(edgelist_intralayers_both, interlayers_with_weights_islands) #combine weighted version of intra with inter

## ----node_metadata--------------------------------------------------------------------------------------------------                                            
pollinators <- sort(unique(intralayer_weighted$node_to)) #adding up only pol who haven't been added yet 
plants <- sort(unique(intralayer_weighted$node_from)) #adding up only plants who haven't been added yet
intersect(pollinators, plants) #making sure I don't have plants in pol or other way around
A <- length(pollinators) # Number of pollinators
P <- length(plants) # Number of plants
S <- A+P

island_names <- c("WesternSahara", #islands as layers
               "Fuerteventura",
               "GranCanaria",
               "TenerifeSouth",
               "TenerifeTeno",
               "Gomera",
               "Hierro")

# Create a table with node metadata
physical_nodes <- tibble(node_id=1:S, #1 till the last species
                         type=c(rep('plant',P),rep('pollinator',A)), #replicate the words P and A times
                         species=c(plants,pollinators)) #add species from plants and pollinators in accordance
layer_metadata <- tibble(layer_id=c(1:7), layer_name=island_names)  #give num to each layer

#write.csv(physical_nodes, "./csvs/Islands/physical_nodes_islands.csv", row.names = FALSE)
#write.csv(layer_metadata, "./csvs/Islands/layer_metadata_islands.csv", row.names = FALSE)


# Replace the node names with node_ids
dryad_edgelist_complete_ids <- 
  dryad_edgelist_complete %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)

#write.csv(dryad_edgelist_complete_ids, "./csvs/Islands/dryad_edgelist_complete_ids_islands.csv", row.names = FALSE)

#---- basic analysis for richness in island --------------------------------------------------------------------------------
tot_plant_ids <- tot_plant %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>%
  subset(select = -layer_from) %>% count(layer_id)
tot_plant_ids$type <- "plant"

tot_pol_ids <- tot_pol %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>%
  subset(select = -layer_from) %>% count(layer_id)
tot_pol_ids$type <- "pollinator"

richness_in_islands <- rbind(tot_plant_ids, tot_pol_ids)

richness_in_islands %>%
  ggplot(aes(x= layer_id, y=n, fill=type))+ geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ 
  scale_x_continuous(breaks=seq(1,7,1))+ labs(x= "Island ID", y= "Number of Species")+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())

## ----multilayer_class-----------------------------------------------------------------------------------------------
# Input: An extended edge list.
dryad_multilayer <- create_multilayer_object(extended = dryad_edgelist_complete_ids, #taking edge list and returning multilayer network
                                             nodes = physical_nodes,
                                             layers = layer_metadata,
                                             intra_output_extended = T)

# Input: intra non-extended edge lists and inter extended edge list
intra_nonextended <-
  dryad_edgelist_complete_ids %>% 
  filter(layer_from==layer_to) %>% #only intra
  dplyr::select(layer=layer_from, node_from, node_to, weight)
inter_extended <-
  dryad_edgelist_complete_ids %>% 
  filter(layer_from!=layer_to) #only inter

#write.csv(intra_nonextended, "./csvs/Islands/dryad_only_intralayer_edges_islands_as_layers.csv")
#write.csv(inter_extended, "./csvs/Islands/dryad_only_interlayer_edges_islands_as_layers.csv")

#create modules for empirical network
modules_dryad_multilayer <- modified_multi(dryad_multilayer, 
                                           infomap_executable = "Infomap",
                                           flow_model = 'directed',
                                           relax = F, 
                                           silent = T, 
                                           trials = 100,
                                           seed = 497294, 
                                           temporal_network = F)

#general info about modules
num_of_nodes_in_module <- modules_dryad_multilayer$modules %>% count(module) #num of nodes in module. biggest module has 44 nodes and smallest has 2
local_modules <- modules_dryad_multilayer$modules %>% select(module, layer_id)
num_of_layers_in_module <- local_modules %>% distinct() %>% count(module) #num of layers a module is found in
mean_num_of_layers_in_module <- mean(num_of_layers_in_module$n) #average number of layers (islands) a module is found in is 3.193

modules_in_network <- modules_dryad_multilayer$modules
#write.csv(modules_in_network, "./csvs/Islands/modules_in_network_islands_as_layers.csv",  row.names = FALSE)

#layer distribution in modulesn
num_of_layers_in_module <- num_of_layers_in_module %>% rename("number_of_layers" = "n")
num_of_modules_with_layers <- num_of_layers_in_module %>% group_by(number_of_layers) %>% count()

pdf('./graphs/modularity_analysis/number_of_layers_in_module.pdf', 10, 6)
num_of_modules_with_layers %>%
  ggplot(aes(x = number_of_layers, y = n))+ geom_bar(stat="identity")+ #stacked
  theme_classic()+ scale_x_continuous(breaks=seq(1,14,1))+ labs(x= "Number of Layers", y= "Number of Modules")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=12, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()






#---- interlayer and intralayer distribution-----------------------------------------------------------------------------------------------------
#interlayer distribution
inter_weight_distribution <- inter_extended %>% arrange(weight) %>% select(weight) 

#intralayer distribution
intra_weight_distribution <- intra_nonextended %>% arrange(weight) %>% select(weight) 

#inter and intra non directed together
inter_intra_non_directed <- data.frame(values= c(intra_nonextended$weight, inter_extended$weight), 
                                       group= c(rep("intra", nrow(intra_nonextended)), rep("inter", nrow(inter_extended))))


inter_intra_non_directed %>%
  ggplot(aes(x=values, fill=group))+ geom_histogram(position= "identity", alpha= 0.6, color= "black")+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())

#directed intralayer weights distribution
directed_weight_distribution_plants <- intralayer_weighted %>% arrange(weight) %>% select(weight)

mean(unlist(directed_weight_distribution_plants))
median(unlist(directed_weight_distribution_plants))

directed_weight_distribution_pols <- intralayer_weighted_inverted %>% arrange(weight) %>% select(weight)

mean(unlist(directed_weight_distribution_pols))
median(unlist(directed_weight_distribution_pols))

# new df with inter and intra in same df with grouping to create distribution
intra_inter_data_for_distibution <- data.frame(values= c(intralayer_weighted$weight, 
                                                         intralayer_weighted_inverted$weight, 
                                                         inter_extended$weight),
                                               group= c(rep("intra plants", nrow(intralayer_weighted)), 
                                                        rep("intra pollinators", nrow(intralayer_weighted_inverted)),
                                                        rep("inter", nrow(inter_extended))))

#write.csv(intra_inter_data_for_distibution, "./csvs/Islands/intra_inter_data_for_distibution_islands_as_layers.csv",  row.names = FALSE)

pdf('./graphs/Islands/intra and interlinks_islands.pdf', 10, 6)
intra_inter_data_for_distibution %>%
  ggplot(aes(x=values, fill=group))+ geom_histogram(position= "identity", alpha= 0.6, color= "black")+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()




# ---- DISTANCE DECAY IN SPECIES- EMPIRICAL DATA ------------------------------------------------------------------------
#distances data
distances_with_ids <- distances_normalized %>% left_join(layer_metadata, by= c("layer_from"="layer_name")) %>% 
  left_join(layer_metadata, by= c("layer_to"="layer_name")) %>% #add correct id to layer name
  select(mean_distance, layer_id.x, layer_id.y) #discard actual names of layers
names(distances_with_ids)[3] <- "layer_from" 
names(distances_with_ids)[4] <- "layer_to"
names(distances_with_ids)[1] <- "layer_name_to"

distances_with_ids <- distances_with_ids[c("layer_from", "layer_to", "mean_distance")]

#write.csv(distances_with_ids, "./csvs/Islands/distances_with_ids_islands_as_layers.csv", row.names = FALSE)
distances_with_ids <- read.csv("./csvs/Islands/distances_with_ids_islands_as_layers.csv")

# classic distance decay
physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")

all_species_all_layers <- rbind(tot_plant, tot_pol) %>% inner_join(physical_nodes, by= c("node_from" = "species")) %>% 
  inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>% subset(select = -c(layer_from, node_from, tot, type)) 
#change node names to ids and layer names to ids and remove unwanted columns

classic_layers_turnover <- NULL

for (i in (1:6)){
  for (j in ((i+1):7)){
    physical_nodes_in_layer_from <- filter(all_species_all_layers, layer_id == i) %>% select(node_id) %>% unlist()
    physical_nodes_in_layer_to <- filter(all_species_all_layers, (layer_id == j)) %>% select(node_id) %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
    int_both <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in both layers
    uni_both <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to)
    turnover <- length(int_both)/length(uni_both)
    classic_layers_turnover <- rbind(classic_layers_turnover, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

classic_layers_turnover <- classic_layers_turnover %>% unique()

classic_layers_turnover_with_distances <- right_join(classic_layers_turnover, distances_with_ids, by= c("layer_from", "layer_to"))
classic_layers_turnover_with_distances <- na.omit(classic_layers_turnover_with_distances) #remove NA and delete layer name

#write.csv(classic_layers_turnover_with_distances, "./csvs/Islands/classic_layers_turnover_with_distances_islands_as_layers.csv")
#classic_layers_turnover_with_distances <- read.csv("./csvs/Islands/classic_layers_turnover_with_distances_islands_as_layers.csv")

classic_layers_turnover_with_distances <- classic_layers_turnover_with_distances %>% 
  mutate(distance_in_km=mean_distance/1000) #turn to km

#distance decay in species empirical alone
pdf('./graphs/Islands/species distance decay alone_islands.pdf', 10, 6)
classic_layers_turnover_with_distances %>%
  ggplot(aes(x=distance_in_km, y=turnover))+ geom_point(color = "indianred2")+ theme_classic()+ 
  stat_smooth(method= "lm", se=F, color = "indianred2")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Km", y="Jaccard Similarity")+
  stat_cor(aes(label = after_stat(rr.label)), label.x = 400, label.y = c(0.36, 0.34, 0.32, 0.30))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()


# ---- DISTANCE DECAY IN MODULES - EMPIRICAL DATA ------------------------------------------------------------------------

#---- jaccard on islands---------------------------------------------------------------------------------------------------
modules<- read.csv("./csvs/Islands/modules_in_network_islands_as_layers.csv")

#write.csv(modules_in_network, "./csvs/Islands/modules_in_network_islands_as_layers.csv")
#similarity check 2 furthest apart
modules_for_similarity_num <- modules %>% select(module, layer_id) %>% 
  unique() %>% group_by(module) %>% select(module) %>% unique()
modules_for_similarity <- modules %>%
  filter(module %in% modules_for_similarity_num$module) #only save the modules that are found in 2 or more layers

#pivot modules function for islands as layers
pivot_by_module_islands <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(data, id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s3 = t(s2) 
  s3 <- s3[-1,]
  colnames(s3) <- c(1,2,3,4,5,6,7)
  return(s3)
}

module_pivoted <- pivot_by_module_islands(modules_for_similarity)

#write.csv(module_pivoted, "./csvs/Islands/module_pivoted_for_state_node_similarity_islands_as_layers.csv")
#module_pivoted <- read.csv("./csvs/Islands/module_pivoted_for_state_node_similarity_islands_as_layers.csv")

#edge list per module function for islands as layers
edge_list_per_module_islands <- function(data,edge_list){
  #gets one row from a data frame and creates an edge list from it
  for (i in (1:6)){
    if (data[i]==0) next #only take layers where the module is present
    else {
      for (j in ((i+1):7)){
        if (data[j]==0) next #only take layers where the module is present
        else {
          edge_list <- rbind(edge_list, tibble(layer_from=i, layer_to=j, module=as.character(NA))) #create edge list of all the layer found in a module
        }
      }
    }
  }
  return(edge_list)
}


modules_edge_list <- NULL

for (i in (1:nrow(module_pivoted))){ #run the function for each row in the data frame
  modules_edge_list <- edge_list_per_module_islands(module_pivoted[i,], modules_edge_list)
  current_module <- rownames(module_pivoted)[i]
  modules_edge_list <- modules_edge_list %>% mutate(module = replace_na(module, current_module)) #add module number
}

edge_list_with_distances <- right_join(modules_edge_list, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances <- na.omit(edge_list_with_distances) #remove NA and delete layer name

#arrange data to include coordinates and modules sizes
size <- count(modules_dryad_multilayer$modules, module)  #create a data frame of all modules and how many nodes are in each (size of module)
module_data <- merge(modules_dryad_multilayer$modules , size, by=c("module","module")) #merge size of module with all the other info about the modules
colnames(module_data)[7] <- "size_of_module" #rename column

#write.csv(module_data, "csvs/Islands/module_data_islands_as_layers.csv", row.names = FALSE)
#write.csv(size, "csvs/Islands/size_islands_as_layers.csv", row.names = FALSE)

lon_lat_data <- read_csv('./csvs/layers.csv') #create new data frame with just the layer data
lon_lat_data <- lon_lat_data %>% select(c("layer_id","lat","Lon")) %>% na.omit()  #only select layer id and coordinates

#layers as islands
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new_values <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)
lon_lat_data$layer_id[lon_lat_data$layer_id %in% old] <- new_values[match(lon_lat_data$layer_id, old)]

lon_lat_data <- lon_lat_data %>% unique() #delete duplicates caused by islands having 2 sites

#write.csv(lon_lat_data, "csvs/Islands/lon_lat_data_islands_as_layers.csv", row.names = FALSE)
#lon_lat_data <- read.csv("./csvs/Islands/lon_lat_data_islands_as_layers.csv")

module_data_with_loc <- merge(module_data, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon <- module_data_with_loc %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon$count <- c(1)

#write.csv(modules_with_lat_lon, "csvs/Islands/modules_with_lat_lon_islands_as_layers.csv", row.names = FALSE)
modules_with_lat_lon <- read.csv("csvs/Islands/modules_with_lat_lon_islands_as_layers.csv")

#version with # of modules in layers
edge_list_by_islands_modules <- edge_list_with_distances
edge_list_by_islands_modules$count <- c(1)
edge_list_by_islands_modules <- edge_list_by_islands_modules %>% group_by(layer_from, layer_to) %>% 
  mutate(number_of_modules= sum(count)) %>% #count how many modules are shared by 2 islands
  select(layer_from, layer_to, module, number_of_modules, mean_distance) 

#distance decay in modules
module_island_turnover <- NULL

for (i in (1:6)){
  for (j in ((1+i):7)){
    print(i)
    modules_in_island_from <- filter(modules_with_lat_lon, layer_id == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_island_to <- filter(modules_with_lat_lon, layer_id == j) %>% select(module) %>% unique() %>% unlist()
    #take all modules in one layer and all modules in another layer to check turnover
    int_both <- intersect(modules_in_island_from, modules_in_island_to) #how many modules are found in both layers
    uni_both <- union(modules_in_island_from, modules_in_island_to)
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover <- rbind(module_island_turnover, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

islands_turnover_with_distnace_empirical <- edge_list_by_islands_modules %>%
  merge(module_island_turnover, by= c("layer_from", "layer_to")) %>% 
  select(layer_from, layer_to, number_of_modules, mean_distance, turnover) %>% unique() #merge both versions 

islands_turnover_with_distnace_empirical <- islands_turnover_with_distnace_empirical %>% 
  mutate(distance_in_km=mean_distance/1000) #turn to km

#write.csv(islands_turnover_with_distnace_empirical, "csvs/Islands/islands_turnover_with_distnace_empirical.csv", row.names = FALSE)
#islands_turnover_with_distnace_empirical <- read.csv("csvs/Islands/islands_turnover_with_distnace_empirical.csv")

pdf('./graphs/Islands/Modules distance decay alone_islands.pdf', 10, 6)
islands_turnover_with_distnace_empirical %>%
  ggplot(aes(x=distance_in_km, y=turnover))+
  geom_point(color = "indianred2")+ 
  scale_x_continuous()+ stat_smooth(method= "lm", se=F, color = "indianred2")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ theme_bw()+
  stat_cor(aes(label = after_stat(rr.label)), label.x = 400, label.y = c(0.36, 0.34, 0.32, 0.30))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()



# ---- NULL MODEL M1: SHUFFLING SPECIES BETWEEN ISLANDS ------------------------------------------------------------------------

# this portion of the code shuffles the networks between layers in one of 3 versions:
# 1. shuffling plants among themselves
# 2. shuffling pollinators among themselves
# 3. shuffling plants among themselves and then pollinators among themselves
# the shuffled networks are then compared the empirical network to determine whether 
# certain network properties are random or not

#---- shuffle pollinators ---------------------------------------------------------------------------------------------
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

#write.csv(shuf_null_edge_list, "./csvs/shuf_null_edge_list_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix, "./csvs/shuf_trial_matrix_islands_as_layers.csv", row.names = TRUE)
#shuf_null_edge_list <- read.csv("./csvs/shuf_null_edge_list_islands_as_layers.csv")
#shuf_trial_matrix <- read.csv("./csvs/shuf_trial_matrix_islands_as_layers.csv")

# 1-39 to
# 40-288 from

#interlayer edges
interlayer_edges_from_shuf <- shuf_null_edge_list %>% group_by(trial_number, node_from) %>%
  select(trial_number, layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3], 
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7]) #all layers the species is found in


interlayer_edges_to_shuf <- shuf_null_edge_list %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_shuf <- rbind(interlayer_edges_from_shuf, interlayer_edges_to_shuf) 

#write.csv(interlayer_edges_shuf, "./csvs/Islands/interlayer_edges_shuf_islands_as_layers.csv",row.names = FALSE)

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

#write.csv(shuf_null_edge_list_plants, "./csvs/shuf_null_edge_list_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_plants, "./csvs/shuf_trial_matrix_plants_islands_as_layers.csv", row.names = TRUE)
#shuf_null_edge_list_plants <- read.csv("./csvs/shuf_null_edge_list_plants_islands_as_layers.csv")
#shuf_trial_matrix_plants <- read.csv("./csvs/shuf_trial_matrix_plants_islands_as_layers.csv")

#interlayer edges
interlayer_edges_from_shuf_plants <- shuf_null_edge_list_plants %>% group_by(trial_number, node_from) %>% 
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3],
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7]) #all layers the species is found in

interlayer_edges_to_shuf_plants <- shuf_null_edge_list_plants %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_shuf_plants <- rbind(interlayer_edges_from_shuf_plants, interlayer_edges_to_shuf_plants) 


#write.csv(interlayer_edges_shuf_plants, "./csvs/Islands/interlayer_edges_shuf_plants_islands_as_layers.csv",row.names = FALSE)

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
#write.csv(shuf_null_edge_list_both, "./csvs/Islands/shuf_null_edge_list_both_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_both, "./csvs/Islands/shuf_trial_matrix_both_islands_as_layers.csv", row.names = TRUE)


#interlayer edges
interlayer_edges_from_shuf_both <- shuf_null_edge_list_both %>% group_by(trial_number, node_from) %>% 
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3],
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7]) #all layers the species is found in


interlayer_edges_to_shuf_both <- shuf_null_edge_list_both %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_shuf_both <- rbind(interlayer_edges_from_shuf_both, interlayer_edges_to_shuf_both) 

#write.csv(interlayer_edges_shuf_both, "./csvs/Islands/interlayer_edges_shuf_both_islands_as_layers.csv",row.names = FALSE)

##---- run oh HPC and get results for analysis -----------------------------------------------------------------
#write.csv(interlayer_edges_shuf, "./HPC_Islands/shuf_between_layers_islands_as_layers/interlayer_edges_shuf_islands_as_layers.csv", row.names = FALSE) #create to run on HPC
#write.csv(interlayer_edges_shuf_plants, "./HPC_Islands/shuf_between_layers_islands_as_layers/interlayer_edges_shuf_plants_islands_as_layers.csv", row.names = FALSE) #create to run on HPC
#write.csv(interlayer_edges_shuf_both, "./HPC_Islands/shuf_between_layers_islands_as_layers/interlayer_edges_shuf_both_islands_as_layers.csv", row.names = FALSE) #create to run on HPC

#run on HPC and then come back with results
#each version has 3 code portions for the HPC:
# 1. HPC_network_x_shuffle.R (x being pollinators, plants or both)
# 2. 1_1000_x.sh (x being pollinators, plants or both)
# 3. i_x.sh (x being pollinators, plants or both)
# running 1_1000_x.sh manually in the cmd will make the other two run.
# all 1000 result csvs can be found in HPC/csvs_x (x being pollinators, plants or both) files

#both
files_both <- list.files("./HPC_Islands/shuf_between_layers_islands_as_layers/csvs_both/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_both <- read_csv(files_both) %>% bind_rows() #create a long edge list with all the csvs

#pollinators
files_pollinators <- list.files("./HPC_Islands/shuf_between_layers_islands_as_layers/csvs_pollinators/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_pol <- read_csv(files_pollinators) %>% bind_rows() #create a long edge list with all the csvs

#plants
files_plants <- list.files("./HPC_Islands/shuf_between_layers_islands_as_layers/csvs_plants/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_plants <- read_csv(files_plants) %>% bind_rows() #create a long edge list with all the csvs

#write.csv(my_merged_interlayer_shuf_both, "./csvs/my_merged_interlayer_shuf_both_islands_as_layers.csv", row.names = FALSE)
#write.csv(my_merged_interlayer_shuf_pol, "./csvs/my_merged_interlayer_shuf_pol_islands_as_layers.csv", row.names = FALSE)
#write.csv(my_merged_interlayer_shuf_plants, "./csvs/my_merged_interlayer_shuf_plants_islands_as_layers.csv", row.names = FALSE)

#---- interlayers with weights shuf version ------------------------------------------
distances_with_weights <- read.csv("./csvs/Islands/distances_with_weights_islands_as_layers.csv")
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")

distances_with_weights_ids <- distances_with_weights %>%
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, layer_to=layer_id.y, weight)

#write.csv(distances_with_weights_ids, "./csvs/Islands/distances_with_weights_ids_islands_as_layers.csv", row.names = FALSE)
#distances_with_weights_ids <- read.csv("./csvs/Islands/distances_with_weights_ids_islands_as_layers.csv")

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

#write.csv(interlayers_with_weights_shuf_pols, "./csvs/Islands/interlayer_shuf_file_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(interlayers_with_weights_shuf_plants, "./csvs/Islands/interlayer_shuf_file_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(interlayers_with_weights_shuf_both, "./csvs/Islands/interlayer_shuf_file_both_islands_as_layers.csv", row.names = FALSE)

#create inter and intra for the 1000 shuf trials
dryad_interlayer_shuf_pols <- read.csv("./csvs/Islands/interlayer_shuf_file_pols_islands_as_layers.csv") #already has inverted
dryad_interlayer_shuf_plants <- read.csv("./csvs/Islands/interlayer_shuf_file_plants_islands_as_layers.csv") #already has inverted
dryad_interlayer_shuf_both <- read.csv("./csvs/Islands/interlayer_shuf_file_both_islands_as_layers.csv") #already has inverted

#pols
#write.csv(shuf_null_edge_list_both, "./csvs/Islands/shuf_null_edge_list_both_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_both, "./csvs/Islands/shuf_trial_matrix_both_islands_as_layers.csv", row.names = TRUE)

dryad_intralayer_shuf_pols <- read.csv("csvs/Islands/shuf_null_edge_list_poll_islands_as_layers.csv") 
dryad_intralayer_shuf_pols <- dryad_intralayer_shuf_pols[, c(6,1,2,3,4,5)]

#plants
dryad_intralayer_shuf_plants <- read.csv("csvs/Islands/shuf_null_edge_list_plants_islands_as_layers.csv") 
dryad_intralayer_shuf_plants <- dryad_intralayer_shuf_plants[, c(6,1,2,3,4,5)]

#both
dryad_intralayer_shuf_both <- read.csv("csvs/Islands/shuf_null_edge_list_both_islands_as_layers.csv") 
dryad_intralayer_shuf_both <- dryad_intralayer_shuf_both[, c(6,1,2,3,4,5)]

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

#write.csv(dryad_edgelist_complete_shuf_pols, "./csvs/islands/dryad_edgelist_complete_shuf_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_edgelist_complete_shuf_plants, "./csvs/Islands/dryad_edgelist_complete_shuf_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_edgelist_complete_shuf_both, "./csvs/Islands/dryad_edgelist_complete_shuf_both_islands_as_layers.csv", row.names = FALSE)




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

#write.csv(all_species_all_layers_all_trials_pols, "./csvs/Islands/all_species_all_layers_all_trials_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(all_species_all_layers_all_trials_plants, "./csvs/Islands/all_species_all_layers_all_trials_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(all_species_all_layers_all_trials_both, "./csvs/Islands/all_species_all_layers_all_trials_both_islands_as_layers.csv", row.names = FALSE)
all_species_all_layers_all_trials_pols <- read.csv("./csvs/Islands/all_species_all_layers_all_trials_pols_islands_as_layers.csv") #already has inverted
all_species_all_layers_all_trials_plants <- read.csv("./csvs/Islands/all_species_all_layers_all_trials_plants_islands_as_layers.csv") #already has inverted
all_species_all_layers_all_trials_both<- read.csv("./csvs/Islands/all_species_all_layers_all_trials_both_islands_as_layers.csv") #already has inverted


#function for islands
islands_distnace_decay_shuf <- function(species_all_layers, output){
  for (k in 1:1000){ #1000 iterations
    current_trial <- species_all_layers %>% filter(trial_number == k)
    print(k) #to keep tab on how far along we are
    for (i in (1:6)){
      for (j in ((i+1):7)){
        physical_nodes_in_layer_from <- filter(current_trial, layer_id == i) %>% select(node_id) %>% unlist() #all species in layer from
        physical_nodes_in_layer_to <- filter(current_trial, (layer_id == j)) %>% select(node_id) %>% unlist() #all species in layer to
        int_both <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are common to both layers
        uni_both <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in both layers total
        turnover <- length(int_both)/length(uni_both)
        output <- rbind(output, tibble(trial_number = k, layer_from = i, layer_to = j, turnover = turnover))
      }
    }
  }
  
  output <- output %>% unique()
  
  return(output)
}

classic_layers_turnover_shuf_pols <- NULL
classic_layers_turnover_shuf_plants <- NULL
classic_layers_turnover_shuf_both <- NULL

#pols
classic_layers_turnover_shuf_output_pols <- islands_distnace_decay_shuf(all_species_all_layers_all_trials_pols, 
                                                                classic_layers_turnover_shuf_pols) 

#plants
classic_layers_turnover_shuf_output_plants <- islands_distnace_decay_shuf(all_species_all_layers_all_trials_plants, 
                                                                  classic_layers_turnover_shuf_plants)

#both
classic_layers_turnover_shuf_output_both <- islands_distnace_decay_shuf(all_species_all_layers_all_trials_both, 
                                                                classic_layers_turnover_shuf_both)


#write.csv(classic_layers_turnover_shuf_output_pols, "./csvs/Islands/classic_layers_turnover_shuf_output_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_shuf_output_plants, "./csvs/Islands/classic_layers_turnover_shuf_output_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_shuf_output_both, "./csvs/Islands/classic_layers_turnover_shuf_output_both_islands_as_layers.csv", row.names = FALSE)
classic_layers_turnover_shuf_output_pols <- read.csv("./csvs/Islands/classic_layers_turnover_shuf_output_pols_islands_as_layers.csv") #already has inverted
classic_layers_turnover_shuf_output_plants <- read.csv("./csvs/Islands/classic_layers_turnover_shuf_output_plants_islands_as_layers.csv") #already has inverted
classic_layers_turnover_shuf_output_both<- read.csv("./csvs/Islands/classic_layers_turnover_shuf_output_both_islands_as_layers.csv") #already has inverted

#pols
distances_with_ids <- read.csv("./csvs/Islands/distances_with_ids_islands_as_layers.csv")
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

#write.csv(classic_layers_turnover_with_distances_shuf_pols, "./csvs/Islands/classic_layers_turnover_with_distances_shuf_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_with_distances_shuf_plants, "./csvs/Islands/classic_layers_turnover_with_distances_shuf_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_with_distances_shuf_both, "./csvs/Islands/classic_layers_turnover_with_distances_shuf_both_islands_as_layers.csv", row.names = FALSE)
classic_layers_turnover_with_distances_shuf_pols <- read.csv("./csvs/Islands/classic_layers_turnover_with_distances_shuf_pols_islands_as_layers.csv")
classic_layers_turnover_with_distances_shuf_plants <- read.csv("./csvs/Islands/classic_layers_turnover_with_distances_shuf_plants_islands_as_layers.csv")
classic_layers_turnover_with_distances_shuf_both <- read.csv("./csvs/Islands/classic_layers_turnover_with_distances_shuf_both_islands_as_layers.csv")

#create an average for shuf with sd
#pols
ave_turnover_for_shuf_pols <- classic_layers_turnover_with_distances_shuf_pols %>% 
  group_by(layer_from, layer_to, mean_distance) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_pollinators") #create mean and sd for each point

#plants
ave_turnover_for_shuf_plants <- classic_layers_turnover_with_distances_shuf_plants %>% 
  group_by(layer_from, layer_to, mean_distance) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_plants") #create mean and sd for each point

#both
ave_turnover_for_shuf_both <- classic_layers_turnover_with_distances_shuf_both %>% 
  group_by(layer_from, layer_to, mean_distance) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null_both") #create mean and sd for each point

#add the emprical classical turnover
classic_layers_turnover_with_distances <- read.csv("./csvs/Islands/classic_layers_turnover_with_distances_islands_as_layers.csv")

empirical_turnover_for_shuf <- classic_layers_turnover_with_distances %>% 
  group_by(layer_from, layer_to, mean_distance) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#---- graphs for distance decay in species-----------------------------------------------------
turnover_shuf_and_empirical <- rbind(empirical_turnover_for_shuf, ave_turnover_for_shuf_pols, ave_turnover_for_shuf_plants,
                                     ave_turnover_for_shuf_both)

turnover_shuf_and_empirical <- turnover_shuf_and_empirical %>% mutate(distance_in_km = mean_distance/1000)

#write.csv(turnover_shuf_and_empirical, "./csvs/Islands/turnover_shuf_and_empirical_islands_as_layers.csv", row.names = FALSE)
turnover_shuf_and_empirical <- read.csv("./csvs/Islands/turnover_shuf_and_empirical_islands_as_layers.csv")

pdf('./graphs/Islands/M1_Species_DD_Islands.pdf', 10, 6)
turnover_shuf_and_empirical %>% ggplot(aes(x= distance_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  labs(x="Distance in Km", y="Jaccard Similarity")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank()) + stat_cor(aes(label = ..p.label..), label.x = 400)+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.36, 0.34, 0.32, 0.30))
dev.off()


#-------making sure its significant---------------------------------------------

lm1 = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical,
                                           turnover_shuf_and_empirical$type=="empirical")) #in empirical
lm2 = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical,
                                           turnover_shuf_and_empirical$type=="null_pollinators")) #in null pols
lm3 = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical,
                                           turnover_shuf_and_empirical$type=="null_plants")) #in null plants
lm4 = lm(ave ~ distance_in_km ,data=subset(turnover_shuf_and_empirical,
                                         turnover_shuf_and_empirical$type=="null_both")) #in null pols                                                   turnover_shuf_and_empirical$type=="empirical"))

b1 <- summary(lm1)$coefficients[2,1] #coef for empirical (slope)
se1 <- summary(lm1)$coefficients[2,2] #sd for empirical
se1
b2 <- summary(lm2)$coefficients[2,1] #coef for pols
se2 <- summary(lm2)$coefficients[2,2] #sd for pols

b3 <- summary(lm3)$coefficients[2,1] #coef for plants
se3 <- summary(lm3)$coefficients[2,2] #sd for plants
b3
se3
b4 <- summary(lm4)$coefficients[2,1] #coef for both
se4 <- summary(lm4)$coefficients[2,2] #sd for both


p_value_pol = 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
p_value_pol #result is 1.752491e-07

p_value_plants = 2*pnorm(-abs(compare.coeff(b1,se1,b3,se3)))
p_value_plants #result is 0.3198874

p_value_both = 2*pnorm(-abs(compare.coeff(b1,se1,b4,se4)))
p_value_both #result is 7.374046e-12





## ----distance decay in modules shuf networks --------------------------------------------------------------------

## ----multilayer_class-----------------------------------------------------------------------------------------------
# Input: An extended edge list.
dryad_edgelist_complete_shuf_pols <- dryad_edgelist_complete_shuf_pols[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_plants <- dryad_edgelist_complete_shuf_plants[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_both <- dryad_edgelist_complete_shuf_both[, c(1,2,3,4,6,5)] #make sure weight is #5

#write.csv(dryad_edgelist_complete_shuf_pols, "./csvs/Islands/dryad_edgelist_complete_shuf_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_edgelist_complete_shuf_plants, "./csvs/Islands/dryad_edgelist_complete_shuf_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_edgelist_complete_shuf_both, "./csvs/Islands/dryad_edgelist_complete_shuf_both_islands_as_layers.csv", row.names = FALSE)
dryad_edgelist_complete_shuf_pols <- read.csv("./csvs/Islands/dryad_edgelist_complete_shuf_pols_islands_as_layers.csv")
dryad_edgelist_complete_shuf_plants <- read.csv("./csvs/Islands/dryad_edgelist_complete_shuf_plants_islands_as_layers.csv")
dryad_edgelist_complete_shuf_both<- read.csv("./csvs/Islands/dryad_edgelist_complete_shuf_both_islands_as_layers.csv")

layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")
physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")

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

#write.csv(dryad_multilayer_shuf_1000_pols_output, "./csvs/Islands/dryad_multilayer_shuf_1000_pols_output_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_multilayer_shuf_1000_plants_output, "./csvs/Islands/dryad_multilayer_shuf_1000_plants_output_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_multilayer_shuf_1000_both_output, "./csvs/Islands/dryad_multilayer_shuf_1000_both_output_islands_as_layers.csv", row.names = FALSE)
dryad_multilayer_shuf_1000_pols_output <- read.csv("./csvs/Islands/dryad_multilayer_shuf_1000_pols_output_islands_as_layers.csv")
dryad_multilayer_shuf_1000_plants_output <- read.csv("./csvs/Islands/dryad_multilayer_shuf_1000_plants_output_islands_as_layers.csv")
dryad_multilayer_shuf_1000_both_output <- read.csv("./csvs/Islands/dryad_multilayer_shuf_1000_both_output_islands_as_layers.csv")


##---- distance decay of modules  ---------------------------------------------------------------
distances_with_ids <- read.csv("./csvs/Islands/distances_with_ids_islands_as_layers.csv")
physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")


#pivot modules function for islands as layers
pivot_by_module_islands <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(modules_for_similarity_shuf, id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s2<-na.omit(s2)
  s3 <-s2%>%
    select(where(~ any(. != 0)))
  s4 = t(s3) 
  s4 <- s4[-1,]
  colnames(s4) <- c(1,2,3,4,5,6,7)
  return(s4)
}



# this function calculates the Jaccard Similarity in modules between islands
module_distance_decay_islands_func <- function(multilayer_1000, 
                                             layers_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial) #to keep tab on how far along we are
    modules_for_similarity_shuf <- multilayer_1000 %>% filter(trial_num == trial) #take only 1 trial
    
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
    edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA
    

    
    
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
          edge_list <- rbind(edge_list, tibble(layer_from=i, layer_to=j, module=as.character(NA))) #create edge list of all the layer found in a module
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

#pols
all_edge_list_layer_combine_no_module_shuf_pols_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_pols_output,
                                                                                             turnover_with_distance_pols)


#write.csv(all_edge_list_layer_combine_no_module_shuf_pols_output, 
     #     "./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_pols_output_islands.csv", 
      #    row.names = FALSE)

#all_edge_list_layer_combine_no_module_shuf_pols_output <- read.csv("./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_pols_output_islands.csv")

#plants
all_edge_list_layer_combine_no_module_shuf_plants_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_plants_output,
                                                                                               turnover_with_distance_plants)


#write.csv(all_edge_list_layer_combine_no_module_shuf_plants_output, 
 #         "./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_plants_output_islands.csv", 
  #        row.names = FALSE)

#all_edge_list_layer_combine_no_module_shuf_plants_output <- read.csv("./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_plants_output_islands.csv")

#both
all_edge_list_layer_combine_no_module_shuf_both_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_both_output,
                                                                                             turnover_with_distance_both)

#write.csv(all_edge_list_layer_combine_no_module_shuf_both_output, 
#          "./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_both_output_islands.csv", 
 #         row.names = FALSE)

#all_edge_list_layer_combine_no_module_shuf_both_output <- read.csv("./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_both_output_islands.csv")


#---- create ave for jaccard islands -------------------------------------------------------------------------
#pols
ave_module_layer_turnover_shuf_pols <- all_edge_list_layer_combine_no_module_shuf_pols_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean(mean_distance)) %>% mutate(type="null_pollinators") #create mean and sd for each point


#write.csv(ave_module_layer_turnover_shuf_pols, 
       #   "./csvs/Islands/ave_module_layer_turnover_shuf_pols_islands.csv", 
      #   row.names = FALSE)

#plants
ave_module_layer_turnover_shuf_plants <- all_edge_list_layer_combine_no_module_shuf_plants_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean(mean_distance)) %>% mutate(type="null_plants") #create mean and sd for each point


#write.csv(ave_module_layer_turnover_shuf_plants, 
 #         "./csvs/Islands/ave_module_layer_turnover_shuf_plants_islands.csv", 
  #        row.names = FALSE)

#both
ave_module_layer_turnover_shuf_both <- all_edge_list_layer_combine_no_module_shuf_both_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean(mean_distance)) %>% mutate(type="null_both") #create mean and sd for each point

#write.csv(ave_module_layer_turnover_shuf_both, 
 #         "./csvs/Islands/ave_module_layer_turnover_shuf_both_islands.csv", 
  #        row.names = FALSE)

#add empirical
#islands_turnover_with_distnace_empirical <- read.csv("csvs/Islands/islands_turnover_with_distnace_empirical.csv")

empirical_turnover_for_modules_layers_shuf <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean_distance) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#combine empirical and null
jaccard_similarity_islands_empirical_and_null <- rbind(empirical_turnover_for_modules_layers_shuf, ave_module_layer_turnover_shuf_pols,
                                                     ave_module_layer_turnover_shuf_plants, ave_module_layer_turnover_shuf_both)

jaccard_similarity_layer_empirical_and_null_km <- jaccard_similarity_islands_empirical_and_null %>% 
  mutate(mean_dist_in_km = mean_distance/1000)

#write.csv(jaccard_similarity_layer_empirical_and_null_km, 
 #         "./csvs/Islands/jaccard_similarity_layer_empirical_and_null_km_islands_m1.csv", 
  #        row.names = FALSE)


#---- graphs for distance decay in modules shuf vs empirical--------------------------------
jaccard_similarity_layer_empirical_and_null_km <- read.csv("csvs/Islands/jaccard_similarity_layer_empirical_and_null_km_islands_m1.csv")

pdf('./graphs/Islands/M1_Modules_DD_Islands.pdf', 10, 6)
jaccard_similarity_layer_empirical_and_null_km %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ 
  stat_cor(aes(label = after_stat(rr.label)), label.x = 400, label.y = c(0.59, 0.56, 0.53, 0.50))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank()) + stat_cor(aes(label = ..p.label..), label.x = 400, label.y = c(0.76, 0.72, 0.68, 0.64))
  
dev.off()

#------check if its significant for islands----------------------------------------------------------------
lm1_module = lm(ave ~ mean_dist_in_km ,data=subset(jaccard_similarity_layer_empirical_and_null_km,
                                            jaccard_similarity_layer_empirical_and_null_km$type=="empirical")) #in empirical
lm2_module = lm(ave ~ mean_dist_in_km ,data=subset(jaccard_similarity_layer_empirical_and_null_km,
                                            jaccard_similarity_layer_empirical_and_null_km$type=="null_pollinators")) #in null pols
lm3_module = lm(ave ~ mean_dist_in_km ,data=subset(jaccard_similarity_layer_empirical_and_null_km,
                                            jaccard_similarity_layer_empirical_and_null_km$type=="null_plants")) #in null plants
lm4_module = lm(ave ~ mean_dist_in_km ,data=subset(jaccard_similarity_layer_empirical_and_null_km,
                                            jaccard_similarity_layer_empirical_and_null_km$type=="null_both")) #in null both
summary(lm1_module)
summary(lm2_module)
summary(lm3_module)
summary(lm4_module)

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
p_value_module_pols #value is 1.287614e-06
p_value_module_plants = 2*pnorm(-abs(compare.coeff(b1_module,se1_module,b3_module,se3_module)))
p_value_module_plants #value is 0.03441519
p_value_module_both = 2*pnorm(-abs(compare.coeff(b1_module,se1_module,b4_module,se4_module)))
p_value_module_both #value is 1.823412e-07

##correlation and r sqaured between jaccard and distance for each run ----------------------------------------------------------------------------
#pols
#all_edge_list_layer_combine_no_module_shuf_pols_output <- read.csv("./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_pols_output_islands.csv")
iteration_correlation_pols <- NULL
iteration_correlation_data_pols <- all_edge_list_layer_combine_no_module_shuf_pols_output %>% subset(layer_from != layer_to) 

iteration_correlation_data_pols_km <- iteration_correlation_data_pols %>% 
  mutate(mean_dist_in_km = mean_distance/1000)

for (i in 1:1000){
  print(i)
  trial_pols = iteration_correlation_data_pols_km %>% filter(trial == i)
  iteration_correlation_new_pols <- cor.test(trial_pols$turnover, trial_pols$mean_dist_in_km, method = "pearson")
  lm_val_pols <- lm(turnover ~ mean_dist_in_km, data = trial_pols)
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


#write.csv(iteration_correlation_pols, "./csvs/Islands/iteration_correlation_pols.csv", row.names = FALSE)
#iteration_correlation_pols <- read.csv("./csvs/Islands/iteration_correlation_pols.csv")

#correlation empirical
#classic_layers_turnover_with_distances <- read.csv("./csvs/Islands/classic_layers_turnover_with_distances_islands_as_layers.csv")
classic_layers_turnover_with_distances <- classic_layers_turnover_with_distances [,-1]
classic_layers_turnover_with_distances_km <- classic_layers_turnover_with_distances %>% 
  mutate(mean_dist_in_km = mean_distance/1000)

layer_turnover_with_distnace_empirical_no_loop <- classic_layers_turnover_with_distances_km %>% subset(layer_from != layer_to) 

correlation_empirical_data_pols <- cor.test(layer_turnover_with_distnace_empirical_no_loop$turnover, #ave is just the value of the turnover
                                            layer_turnover_with_distnace_empirical_no_loop$mean_dist_in_km, method = "pearson")

lm_val_empirical_pols <- lm(turnover ~ mean_dist_in_km, data = layer_turnover_with_distnace_empirical_no_loop)

correlation_empirical_pols <- tibble(estimate = correlation_empirical_data_pols$estimate, 
                                     p_val = correlation_empirical_data_pols$p.value, 
                                     statistic = correlation_empirical_data_pols$statistic, 
                                     confidence_int_low = correlation_empirical_data_pols$conf.int[1],
                                     confidence_int_high = correlation_empirical_data_pols$conf.int[2],
                                     slope = lm_val_empirical_pols$coefficients[2],
                                     intercept = lm_val_empirical_pols$coefficients[1], 
                                     rsquared = summary(lm_val_empirical_pols)$adj.r.squared)

#write.csv(correlation_empirical_pols, "./csvs/Islands/correlation_empirical_pols.csv", row.names = FALSE) #so it can be used for classical shuffling
#correlation_empirical_pols <- read.csv("./csvs/Islands/correlation_empirical_pols.csv")

##distribution of rsquared and add empirical
#pdf('./graphs/shuffle_between_layers/pols_r_squares_module_DD.pdf', 10, 6) SACAR SI NO SE USA
rsquared_pols <- iteration_correlation_pols %>% 
  ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#BE75FA", color = "#BE75FA", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_pols$rsquared, linetype = "dashed", color = "#F47069")+
  theme(axis.title=element_text(size=22))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
#dev.off()

rsquared_pols 
p_rsquared_pols <- sum(iteration_correlation_pols$rsquared > correlation_empirical_pols$rsquared)/1000
p_rsquared_pols

#----plants
#all_edge_list_layer_combine_no_module_shuf_plants_output <- read.csv("./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_plants_output_islands.csv")
iteration_correlation_plants <- NULL
iteration_correlation_data_plants <- all_edge_list_layer_combine_no_module_shuf_plants_output %>% subset(layer_from != layer_to) 

iteration_correlation_data_plants_km <- iteration_correlation_data_plants %>% 
  mutate(mean_dist_in_km = mean_distance/1000)

for (i in 1:1000){
  print(i)
  trial_plants = iteration_correlation_data_plants_km %>% filter(trial == i)
  iteration_correlation_new_plants <- cor.test(trial_plants$turnover, trial_plants$mean_dist_in_km, method = "pearson")
  lm_val_plants <- lm(turnover ~ mean_dist_in_km, data = trial_plants)
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

#write.csv(iteration_correlation_plants, "./csvs/Islands/iteration_correlation_plants.csv", row.names = FALSE)
#iteration_correlation_plants <- read.csv("./csvs/Islands/iteration_correlation_plants.csv")

##distribution of rsquared and add empirical
#pdf('./graphs/shuffle_between_layers/plants_r_squares_module_DD.pdf', 10, 6) SACAR SI NO SE USA
rsquared_plants <- iteration_correlation_plants %>% 
  ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#15B7BC", color = "#15B7BC", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_pols$rsquared, linetype = "dashed", color = "#F47069")+
  theme(axis.title=element_text(size=22))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
#dev.off()
rsquared_plants
p_rsquared_plants <- sum(iteration_correlation_plants$rsquared < correlation_empirical_pols$rsquared)/1000
p_rsquared_plants #0.997


#---- both
#all_edge_list_layer_combine_no_module_shuf_both_output <- read.csv("./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_both_output_islands.csv")
iteration_correlation_both <- NULL
iteration_correlation_data_both <- all_edge_list_layer_combine_no_module_shuf_both_output %>% subset(layer_from != layer_to) 

iteration_correlation_data_both_km <- iteration_correlation_data_both %>% 
  mutate(mean_dist_in_km = mean_distance/1000)


for (i in 1:1000){
  print(i)
  trial_both = iteration_correlation_data_both_km %>% filter(trial == i)
  iteration_correlation_new_both <- cor.test(trial_both$turnover, trial_both$mean_dist_in_km, method = "pearson")
  lm_val_both <- lm(turnover ~ mean_dist_in_km, data = trial_both)
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

#write.csv(iteration_correlation_both, "./csvs/Islands/iteration_correlation_both.csv", row.names = FALSE)
#iteration_correlation_both<- read.csv("./csvs/Islands/iteration_correlation_both.csv")

##distribution of rsquared and add empirical

#pdf('./graphs/shuffle_between_layers/both_r_squares_module_DD.pdf', 10, 6) SACAR SI NO SE USA
rsquared_both <- iteration_correlation_both %>% 
  ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#72A323", color = "#72A323", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_pols$rsquared, linetype = "dashed", color = "#F47069")+
  theme(axis.title=element_text(size=22))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
#dev.off()
rsquared_both
p_rsquared_both <- sum(iteration_correlation_both$rsquared > correlation_empirical_pols$rsquared)/1000
p_rsquared_both #0


#all 3 R squared in the same square
iteration_correlation_pols_M1 <- iteration_correlation_pols %>% mutate(type = "shuf_pollinators")
iteration_correlation_plants_M1 <- iteration_correlation_plants %>% mutate(type = "shuf_plants")
iteration_correlation_both_M1 <- iteration_correlation_both %>% mutate(type = "shuf_both")

rqsuares_M1_all <- rbind(iteration_correlation_pols_M1, iteration_correlation_plants_M1, iteration_correlation_both_M1)

#write.csv(rqsuares_M1_all, "./csvs/Islands/rqsuares_M1_all.csv", row.names = FALSE)
#rqsuares_M1_all <- read.csv("./csvs/Islands/rqsuares_M1_all.csv")


group_color <- c(shuf_pollinators = "#BE75FA", 
                 shuf_plants = "#15B7BC",
                 shuf_both = "#72A323")


pdf('./graphs/Islands/M1_r_squares_module_DD.pdf', 10, 6)
rqsuares_M1_all %>% 
  ggplot(aes(x = rsquared, fill = type))+ 
  geom_density(alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_pols$rsquared, linetype = "dashed", color = "#F47069")+
  theme(axis.title=element_text(size=22))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())+
  scale_fill_manual(values = group_color)
dev.off()

