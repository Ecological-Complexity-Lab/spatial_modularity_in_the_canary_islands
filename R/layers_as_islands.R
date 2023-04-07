
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
library(reshape2)

#this portion of code turns the data into a multilayer network with 7 layers 
#(islands as layers). It then compares the distance decay in species and modules
#of the empirical network to null models where plants, pollinators and both are
#shuffled.



##----get_data--------------------------------------------------------------------------------------------------------
#setwd("/Users/maya/Desktop/plant_pollinator_data/dryad_network")
#setwd("/Users/mayagoldstein/Desktop/project")
#setwd("/Users/golds/Desktop/spatial_modularity_in_the_canary_islands")
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
  summarise(sum_weight = sum(weight)) #turn sums of sites to sum of island

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

distances_normalized <- distances %>% filter(layer_from != layer_to) %>% #delete distances between sites in the same island
  group_by(layer_to, layer_from) %>% #group will contain 4 sites- site 1 and 1 of layer from and site 1 and 2 or layer to
  summarise(mean_distance = mean(distance_in_meters)) %>%unique() #use an average distance of the 4 sites in 2 different islands to determine the distance between the islands

distances_normalized <- distances_normalized[c("layer_from", "layer_to", "mean_distance")]

shortest_distance <- min(distances_normalized$mean_distance) #the shortest distance between islands

interlayer_weight <- function(d){
  #recieves distance and normalizes it
  weight <- (1/log(d))/(1/log(shortest_distance))
  return(weight)
}

distances_with_weights <- distances_normalized %>% 
  mutate(weight = interlayer_weight(mean_distance)) %>% #add weight using the function
  subset(select = -mean_distance)

#write.csv(distances_with_weights, "./csvs/distances_with_weights_islands_as_layers.csv", row.names = FALSE)

interlayers_with_weights_islands <- interlayers_new_islands %>% inner_join(distances_with_weights, 
                                                           by = c("layer_from", "layer_to")) %>% unique()

#write.csv(interlayers_with_weights_islands, "./csvs/interlayers_with_weights_islands.csv", row.names = FALSE)

#interlayer distribution
weight_distribution_islands_as_layers <- interlayers_with_weights_islands %>% arrange(weight) %>% select(weight) %>% 
  mutate(type = "islands")
weight_distribution_islands_as_layers %>%
  ggplot(aes(x=weight))+geom_density(fill = "#F47069", color = "#F47069", alpha = 0.4)+theme_classic()+ 
  geom_vline(xintercept = median(weight_distribution_islands_as_layers$weight), linetype = "dashed", color = "black", size = 1)+
  theme(axis.title=element_text(size=15))

mean_islands <- mean(unlist(weight_distribution_islands_as_layers$weight))
median(unlist(weight_distribution_islands_as_layers$weight))

#compare previous interayer edges to new interlayer edges
interlayer_extended <- read.csv("./csvs/dryad_only_interlayer_edges.csv") #edges based on sites

weight_distribution <- interlayer_extended %>% arrange(weight) %>% select(weight) %>%
  mutate(type = "sites")

mean_sites <- mean(unlist(weight_distribution$weight))

mean_both_versions <- data.frame(mean_weight = c(mean_islands, mean_sites),
                                 type = c("islands", "sites")) #create data frame of both means to add as vlines to graph

interlayers_both_versions <- rbind(weight_distribution, weight_distribution_islands_as_layers) #combine the two to compare interlayer edges distribution 

#graph comparing interlayer edges in both versions
interlayers_both_versions %>%
  ggplot(aes(x=weight, group = type))+ geom_density(aes(fill = type, color = type), alpha = 0.4)+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())+
  geom_vline(data = mean_both_versions, aes(xintercept = mean_weight, color = type))


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

#write.csv(dryad_edgelist_complete_ids, "./csvs/dryad_edgelist_complete_ids_islands.csv", row.names = FALSE)

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

#write.csv(intra_nonextended, "./csvs/dryad_only_intralayer_edges_islands_as_layers.csv")
#write.csv(inter_extended, "./csvs/dryad_only_interlayer_edges_islands_as_layers.csv")

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
#write.csv(modules_in_network, "./csvs/modules_in_network_islands_as_layers.csv")

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

intra_inter_data_for_distibution %>%
  ggplot(aes(x=values, fill=group))+ geom_histogram(position= "identity", alpha= 0.6, color= "black")+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())

# ---- distance decay in species in the network ------------------------------------------------------------------------
#distances of data
distances_with_ids <- distances_normalized %>% left_join(layer_metadata, by= c("layer_from"="layer_name")) %>% 
  left_join(layer_metadata, by= c("layer_to"="layer_name")) %>% #add correct id to layer name
  select(mean_distance, layer_id.x, layer_id.y) #discard actual names of layers
names(distances_with_ids)[3] <- "layer_from" 
names(distances_with_ids)[4] <- "layer_to"
names(distances_with_ids)[1] <- "layer_name_to"

distances_with_ids <- distances_with_ids[c("layer_from", "layer_to", "mean_distance")]

#write.csv(distances_with_ids, "./csvs/distances_with_ids_islands_as_layers.csv", row.names = FALSE)
#distances_with_ids <- read.csv("./csvs/distances_with_ids_islands_as_layers.csv")

# classic distnace decay
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

#write.csv(classic_layers_turnover_with_distances, "./csvs/classic_layers_turnover_with_distances_islands_as_layers.csv")
#classic_layers_turnover_with_distances <- read.csv("./csvs/classic_layers_turnover_with_distances_islands_as_layers.csv")

classic_layers_turnover_with_distances <- classic_layers_turnover_with_distances %>% 
  mutate(distance_in_km=mean_distance/1000) #turn to km

classic_layers_turnover_with_distances %>%
  ggplot(aes(x=distance_in_km, y=turnover))+ geom_point(color = "indianred2")+ theme_classic()+ 
  stat_smooth(method= "lm", se=F, color = "indianred2")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Km", y="Jaccard Similarity")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())

#---- jaccard on islands---------------------------------------------------------------------------------------------------
#similarity check 2 furthest apart
modules_for_similarity_num <- modules_dryad_multilayer$modules %>% select(module, layer_id) %>% 
  unique() %>% group_by(module) %>% select(module) %>% unique()
modules_for_similarity <- modules_dryad_multilayer$modules %>%
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

#write.csv(module_pivoted, "./csvs/module_pivoted_for_state_node_similarity_islands_as_layers.csv")

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

#write.csv(module_data, "csvs/module_data_islands_as_layers.csv", row.names = FALSE)
#write.csv(size, "csvs/size_islands_as_layers.csv", row.names = FALSE)


lon_lat_data <- read_csv('./csvs/layers.csv') #create new data frame with just the layer data
lon_lat_data <- lon_lat_data %>% select(c("layer_id","lat","Lon")) %>% na.omit()  #only select layer id and coordinates

#layers as islands
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new_values <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)
lon_lat_data$layer_id[lon_lat_data$layer_id %in% old] <- new_values[match(lon_lat_data$layer_id, old)]

lon_lat_data <- lon_lat_data %>% unique() #delete duplicates caused by islands having 2 sites

#write.csv(lon_lat_data, "csvs/lon_lat_data_islands_as_layers.csv", row.names = FALSE)
#lon_lat_data <- read.csv("./csvs/lon_lat_data_islands_as_layers.csv")

module_data_with_loc <- merge(module_data, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon <- module_data_with_loc %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon$count <- c(1)

modules_with_lat_lon %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_id)))+ geom_bar(stat= "identity")+ 
  scale_x_continuous(breaks=seq(1,88,5))+ labs(y="Number of Physical Nodes", x="Module Number")+
  guides(fill=guide_legend(title="Layer\nNumber"))+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())

#write.csv(modules_with_lat_lon, "csvs/modules_with_lat_lon_islands_as_layers.csv", row.names = FALSE)
#modules_with_lat_lon <- read.csv("csvs/modules_with_lat_lon_islands_as_layers.csv")

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


pdf('./graphs/layers_as_islands_and_not_sites/distance_decay_in_modules_empirical_islands_as_layers.pdf', 10, 6)
islands_turnover_with_distnace_empirical %>%
  ggplot(aes(x=distance_in_km, y=turnover))+
  geom_point(color = "indianred2")+ 
  scale_x_continuous()+ stat_smooth(method= "lm", se=F, color = "indianred2")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()

#write.csv(islands_turnover_with_distnace_empirical, "csvs/islands_turnover_with_distnace_empirical_islands_as_layers.csv", row.names = FALSE)

##---- shuf pols plants both-------------------------------------------------------------------------------------------
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

#write.csv(interlayer_edges_shuf, "./csvs/interlayer_edges_shuf_islands_as_layers.csv",row.names = FALSE)

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


#write.csv(interlayer_edges_shuf_plants, "./csvs/interlayer_edges_shuf_plants_islands_as_layers.csv",row.names = FALSE)

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
#write.csv(shuf_null_edge_list_both, "./csvs/shuf_null_edge_list_both_islands_as_layers.csv", row.names = FALSE)
#write.csv(shuf_trial_matrix_both, "./csvs/shuf_trial_matrix_both_islands_as_layers.csv", row.names = TRUE)


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

#write.csv(interlayer_edges_shuf_both, "./csvs/interlayer_edges_shuf_both_islands_as_layers.csv",row.names = FALSE)

##---- run oh HPC and get results for analysis -----------------------------------------------------------------
write.csv(interlayer_edges_shuf, "./HPC/shuf_between_layers_islands_as_layers/interlayer_edges_shuf_islands_as_layers.csv", row.names = FALSE) #create to run on HPC
write.csv(interlayer_edges_shuf_plants, "./HPC/shuf_between_layers_islands_as_layers/interlayer_edges_shuf_plants_islands_as_layers.csv", row.names = FALSE) #create to run on HPC
write.csv(interlayer_edges_shuf_both, "./HPC/shuf_between_layers_islands_as_layers/interlayer_edges_shuf_both_islands_as_layers.csv", row.names = FALSE) #create to run on HPC

#run on HPC and then come back with results
#each version has 3 code portions for the HPC:
# 1. HPC_network_x_shuffle.R (x being pollinators, plants or both)
# 2. 1_1000_x.sh (x being pollinators, plants or both)
# 3. i_x.sh (x being pollinators, plants or both)
# running 1_1000_x.sh manually in the cmd will make the other two run.
# all 1000 result csvs can be found in HPC/csvs_x (x being pollinators, plants or both) files

#both
files_both <- list.files("./HPC/shuf_between_layers_islands_as_layers/csvs_both/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_both <- read_csv(files_both) %>% bind_rows() #create a long edge list with all the csvs

#pollinators
files_pollinators <- list.files("./HPC/shuf_between_layers_islands_as_layers/csvs_pollinators/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_pol <- read_csv(files_pollinators) %>% bind_rows() #create a long edge list with all the csvs

#plants
files_plants <- list.files("./HPC/shuf_between_layers_islands_as_layers/csvs_plants/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_plants <- read_csv(files_plants) %>% bind_rows() #create a long edge list with all the csvs

#write.csv(my_merged_interlayer_shuf_both, "./csvs/my_merged_interlayer_shuf_both_islands_as_layers.csv", row.names = FALSE)
#write.csv(my_merged_interlayer_shuf_pol, "./csvs/my_merged_interlayer_shuf_pol_islands_as_layers.csv", row.names = FALSE)
#write.csv(my_merged_interlayer_shuf_plants, "./csvs/my_merged_interlayer_shuf_plants_islands_as_layers.csv", row.names = FALSE)

#---- interlayers with weights shuf version ------------------------------------------
distances_with_weights_ids <- distances_with_weights %>%
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, layer_to=layer_id.y, weight)

#write.csv(distances_with_weights_ids, "./csvs/distances_with_weights_ids_islands_as_layers.csv", row.names = FALSE)

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

#write.csv(interlayers_with_weights_shuf_pols, "./csvs/interlayer_shuf_file_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(interlayers_with_weights_shuf_plants, "./csvs/interlayer_shuf_file_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(interlayers_with_weights_shuf_both, "./csvs/interlayer_shuf_file_both_islands_as_layers.csv", row.names = FALSE)

#create inter and intra for the 1000 shuf trials
dryad_interlayer_shuf_pols <- read.csv("./csvs/interlayer_shuf_file_pols_islands_as_layers.csv") #already has inverted
dryad_interlayer_shuf_plants <- read.csv("./csvs/interlayer_shuf_file_plants_islands_as_layers.csv") #already has inverted
dryad_interlayer_shuf_both <- read.csv("./csvs/interlayer_shuf_file_both_islands_as_layers.csv") #already has inverted

#pols
dryad_intralayer_shuf_pols <- read.csv("csvs/shuf_null_edge_list_islands_as_layers.csv") 
dryad_intralayer_shuf_pols <- dryad_intralayer_shuf_pols[, c(6,1,2,3,4,5)]

#plants
dryad_intralayer_shuf_plants <- read.csv("csvs/shuf_null_edge_list_plants_islands_as_layers.csv") 
dryad_intralayer_shuf_plants <- dryad_intralayer_shuf_plants[, c(6,1,2,3,4,5)]

#both
dryad_intralayer_shuf_both <- read.csv("csvs/shuf_null_edge_list_both_islands_as_layers.csv") 
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

#write.csv(all_species_all_layers_all_trials_pols, "./csvs/all_species_all_layers_all_trials_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(all_species_all_layers_all_trials_plants, "./csvs/all_species_all_layers_all_trials_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(all_species_all_layers_all_trials_both, "./csvs/all_species_all_layers_all_trials_both_islands_as_layers.csv", row.names = FALSE)


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


#write.csv(classic_layers_turnover_shuf_output_pols, "./csvs/classic_layers_turnover_shuf_output_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_shuf_output_plants, "./csvs/classic_layers_turnover_shuf_output_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_shuf_output_both, "./csvs/classic_layers_turnover_shuf_output_both_islands_as_layers.csv", row.names = FALSE)


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

#write.csv(classic_layers_turnover_with_distances_shuf_pols, "./csvs/classic_layers_turnover_with_distances_shuf_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_with_distances_shuf_plants, "./csvs/classic_layers_turnover_with_distances_shuf_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(classic_layers_turnover_with_distances_shuf_both, "./csvs/classic_layers_turnover_with_distances_shuf_both_islands_as_layers.csv", row.names = FALSE)
classic_layers_turnover_with_distances_shuf_pols <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_pols_islands_as_layers.csv")
classic_layers_turnover_with_distances_shuf_plants <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_plants_islands_as_layers.csv")
classic_layers_turnover_with_distances_shuf_both <- read.csv("./csvs/classic_layers_turnover_with_distances_shuf_both_islands_as_layers.csv")

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
empirical_turnover_for_shuf <- classic_layers_turnover_with_distances %>% 
  group_by(layer_from, layer_to, mean_distance) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#---- graphs for distance decay in species-----------------------------------------------------
turnover_shuf_and_empirical <- rbind(empirical_turnover_for_shuf, ave_turnover_for_shuf_pols, ave_turnover_for_shuf_plants,
                                     ave_turnover_for_shuf_both)

turnover_shuf_and_empirical <- turnover_shuf_and_empirical %>% mutate(distance_in_km = mean_distance/1000)

#write.csv(turnover_shuf_and_empirical, "./csvs/turnover_shuf_and_empirical_islands_as_layers.csv", row.names = FALSE)

turnover_shuf_and_empirical %>% ggplot(aes(x= distance_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  labs(x="Distance in Km", y="Jaccard Similarity")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
  #+ stat_cor(aes(label = ..p.label..), label.x = 400)+
  #stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.36, 0.34, 0.32, 0.30))

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
p_value_pol #result is 1.752491e-07

p_value_plants = 2*pnorm(-abs(compare.coeff(b1,se1,b3,se3)))
p_value_plants #result is 0.3198874

p_value_both = 2*pnorm(-abs(compare.coeff(b1,se1,b4,se4)))
p_value_both #result is 7.374046e-12

## ----multilayer_class-----------------------------------------------------------------------------------------------
# Input: An extended edge list.
dryad_edgelist_complete_shuf_pols <- dryad_edgelist_complete_shuf_pols[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_plants <- dryad_edgelist_complete_shuf_plants[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_both <- dryad_edgelist_complete_shuf_both[, c(1,2,3,4,6,5)] #make sure weight is #5

#write.csv(dryad_edgelist_complete_shuf_pols, "./csvs/dryad_edgelist_complete_shuf_pols_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_edgelist_complete_shuf_plants, "./csvs/dryad_edgelist_complete_shuf_plants_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_edgelist_complete_shuf_both, "./csvs/dryad_edgelist_complete_shuf_both_islands_as_layers.csv", row.names = FALSE)

dryad_edgelist_complete_shuf_pols <- read.csv("./csvs/dryad_edgelist_complete_shuf_pols_islands_as_layers.csv")
dryad_edgelist_complete_shuf_plants <- read.csv("./csvs/dryad_edgelist_complete_shuf_plants_islands_as_layers.csv")
dryad_edgelist_complete_shuf_both <- read.csv("./csvs/dryad_edgelist_complete_shuf_both_islands_as_layers.csv")


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

#write.csv(dryad_multilayer_shuf_1000_pols_output, "./csvs/dryad_multilayer_shuf_1000_pols_output_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_multilayer_shuf_1000_plants_output, "./csvs/dryad_multilayer_shuf_1000_plants_output_islands_as_layers.csv", row.names = FALSE)
#write.csv(dryad_multilayer_shuf_1000_both_output, "./csvs/dryad_multilayer_shuf_1000_both_output_islands_as_layers.csv", row.names = FALSE)
dryad_multilayer_shuf_1000_pols_output <- read.csv("./csvs/dryad_multilayer_shuf_1000_pols_output_islands_as_layers.csv")
dryad_multilayer_shuf_1000_plants_output <- read.csv("./csvs/dryad_multilayer_shuf_1000_plants_output_islands_as_layers.csv")
dryad_multilayer_shuf_1000_both_output <- read.csv("./csvs/dryad_multilayer_shuf_1000_both_output_islands_as_layers.csv")


##---- distance decay of modules  ---------------------------------------------------------------
# this function calculates the Jaccard Similarity in modules between layers
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
    
    #version with # of modules in layers
    #edge_list_by_layers_modules_shuf <- edge_list_with_distances_shuf %>% group_by(layer_from, layer_to, module) %>%
    #  summarise(ave_distance= mean(mean_distance)) 
    #edge_list_by_layers_modules_shuf$count <- c(1)
    #edge_list_by_layers_modules_shuf <- edge_list_by_layers_modules_shuf %>% mutate(number_of_modules= sum(count)) %>%
    #  select(layer_from, layer_to, module, number_of_modules) 
    
    module_data_shuf <- merge(modules_for_similarity_shuf , size, by=c("module","module")) #merge size of module with all the other info about the modules
    colnames(module_data_shuf)[8] <- "size_of_module" #rename column

    
    module_data_with_loc_shuf<- merge(module_data_shuf, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates
    
    modules_with_lat_lon_shuf <- module_data_with_loc_shuf %>% 
      select(layer_id, module, lat, Lon, trial_num, size_of_module) %>% unique() #take only certain columns
    modules_with_lat_lon_shuf$count <- c(1)
    
    
    
    for (i in 1:6){
      for (j in (i+1):7){
        modules_in_layer_from_shuf <- filter(modules_with_lat_lon_shuf, layer_id == i) %>% select(module) %>% unique() %>% unlist() #modules in layer from
        modules_in_layer_to_shuf <- filter(modules_with_lat_lon_shuf, layer_id == j) %>% select(module) %>% unique() %>% unlist() #modules in layer to
        int_both <- intersect(modules_in_layer_from_shuf, modules_in_layer_to_shuf) #how many modules are common in both layers
        uni_both <- union(modules_in_layer_from_shuf, modules_in_layer_to_shuf) #how many modules are found in both layers in total
        turnover <- length(int_both)/length(uni_both)
        module_layer_turnover_shuf <- rbind(module_layer_turnover_shuf, tibble(layer_from= i, layer_to= j, turnover= turnover, trial = trial))
      }
    }
    
    #edge_list_by_layers_ave_shuf <- edge_list_with_distances_shuf %>% group_by(layer_from, layer_to) %>%
    #  summarise(ave_distance= mean(mean_distance)) %>% unique()
    
    
    layers_turnover_with_distnace <- edge_list_with_distances_shuf %>%
      merge(module_layer_turnover_shuf, by= c("layer_from", "layer_to")) #merge both versions

  }
  return(layers_turnover_with_distnace)
}

turnover_with_distance_pols <- NULL
turnover_with_distance_plants <- NULL
turnover_with_distance_both <- NULL

module_layer_turnover_shuf <- NULL

#pols
all_edge_list_layer_combine_no_module_shuf_pols_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_pols_output,
                                                                                             turnover_with_distance_pols)


write.csv(all_edge_list_layer_combine_no_module_shuf_pols_output, 
          "./csvs/all_edge_list_layer_combine_no_module_shuf_pols_output_islands.csv", 
          row.names = FALSE)

#plants
all_edge_list_layer_combine_no_module_shuf_plants_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_plants_output,
                                                                                               turnover_with_distance_plants)


write.csv(all_edge_list_layer_combine_no_module_shuf_plants_output, 
          "./csvs/all_edge_list_layer_combine_no_module_shuf_plants_output_islands.csv", 
          row.names = FALSE)

#both
all_edge_list_layer_combine_no_module_shuf_both_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_both_output,
                                                                                             turnover_with_distance_both)

write.csv(all_edge_list_layer_combine_no_module_shuf_both_output, 
          "./csvs/all_edge_list_layer_combine_no_module_shuf_both_output_islands.csv", 
          row.names = FALSE)

#---- create ave for jaccard islands -------------------------------------------------------------------------
#pols
ave_module_layer_turnover_shuf_pols <- all_edge_list_layer_combine_no_module_shuf_pols_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean_distance) %>% mutate(type="null_pollinators") #create mean and sd for each point

write.csv(ave_module_layer_turnover_shuf_pols, 
          "./csvs/ave_module_layer_turnover_shuf_pols_islands.csv", 
          row.names = FALSE)

#plants
ave_module_layer_turnover_shuf_plants <- all_edge_list_layer_combine_no_module_shuf_plants_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean_distance) %>% mutate(type="null_plants") #create mean and sd for each point

write.csv(ave_module_layer_turnover_shuf_plants, 
          "./csvs/ave_module_layer_turnover_shuf_plants_islands.csv", 
          row.names = FALSE)

#both
ave_module_layer_turnover_shuf_both <- all_edge_list_layer_combine_no_module_shuf_both_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean_distance) %>% mutate(type="null_both") #create mean and sd for each point

write.csv(ave_module_layer_turnover_shuf_both, 
          "./csvs/ave_module_layer_turnover_shuf_both_islands.csv", 
          row.names = FALSE)

#add empirical
empirical_turnover_for_modules_layers_shuf <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), mean_distance = mean_distance) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null


#combine empirical and null
jaccard_similarity_islands_empirical_and_null <- rbind(empirical_turnover_for_modules_layers_shuf, ave_module_layer_turnover_shuf_pols,
                                                     ave_module_layer_turnover_shuf_plants, ave_module_layer_turnover_shuf_both)

jaccard_similarity_layer_empirical_and_null_km <- jaccard_similarity_islands_empirical_and_null %>% 
  mutate(mean_dist_in_km = mean_distance/1000)

write.csv(jaccard_similarity_layer_empirical_and_null_km, 
          "./csvs/jaccard_similarity_layer_empirical_and_null_km_islands.csv", 
          row.names = FALSE)

#---- graphs for distance decay in modules shuf vs empirical--------------------------------
pdf('./graphs/layers_as_islands_and_not_sites/jaccard_similarity_layer_empirical_and_null_km_islands_as_layers.pdf', 10, 6)
jaccard_similarity_layer_empirical_and_null_km %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ 
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())#stat_cor(aes(label = ..p.label..), label.x = 400)+
  #stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.63, 0.60, 0.57, 0.54))
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

##---- classic shuffle within layers ---------------------------------------------

shuf_trial_matrix_classic <- NULL
shuf_null_edge_list_classic <- NULL

intralayers_with_ids_for_shuf <- 
  dryad_intralayer_islands_grouped %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, sum_weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, sum_weight)


for (i in 1:7){
  current_matrix <- intralayers_with_ids_for_shuf %>%
    filter(layer_from == layer_to) %>% filter(layer_from == i) %>% #take 1 layer at a time
    select(node_from, node_to, sum_weight) %>% #select interactions
    dcast(node_from ~ node_to, value.var = "sum_weight", fill = 0) %>%
    column_to_rownames(var="node_from") 
  
  current_matrix <- t(current_matrix) #put pols in rows
  
  print(i) #to keep tab on which layer we're on
  
  for (j in 1:1000){ #1000 iterations
    null <- vegan::nullmodel(current_matrix, method = 'r00_samp') #shuffle within layer
    shuffled_matrices <- simulate(null, nsim = 1)
    
    trial_with_shuf_num <- cbind(shuffled_matrices, j) #add trial number 
    shuf_trial_matrix_classic <- rbind(shuf_trial_matrix_classic, trial_with_shuf_num) #create big matrix of all matrices
    edge_list_version_classic <- melt(as.matrix(shuffled_matrices)) %>% filter(value > 0) %>%
      select(node_from=Var1, node_to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
    edge_list_version_classic$trial_number <- j #add trial number to the edge list
    edge_list_version_classic$layer_from <- i #intra so layer from and to are identical
    edge_list_version_classic$layer_to <- i
    shuf_null_edge_list_classic <- rbind(shuf_null_edge_list_classic, edge_list_version_classic) #create mega edge list with all repetitions
  }
}

shuf_null_edge_list_classic <- shuf_null_edge_list_classic[, c(5,1,6,2,3,4)] #change to regular order

#write.csv(shuf_null_edge_list_classic, "./csvs/shuf_null_edge_list_classic_islands_as_layers.csv", row.names = FALSE)

#interlayer edges
interlayer_edges_from_shuf_classic <- shuf_null_edge_list_classic %>% group_by(trial_number, node_from) %>%
  select(trial_number, layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3], 
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in


interlayer_edges_to_shuf_classic <- shuf_null_edge_list_classic %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind

interlayer_edges_shuf_classic <- rbind(interlayer_edges_from_shuf_classic, interlayer_edges_to_shuf_classic) 

#write.csv(interlayer_edges_shuf_classic, "./csvs/interlayer_edges_shuf_classic_islands_as_layers.csv",row.names = FALSE)

#---- run oh HPC and get results for analysis -----------------------------------------------------------------
write.csv(interlayer_edges_shuf_classic, "./HPC/shuf_within_layers_islands/interlayer_edges_shuf_classic_islands_as_layers.csv", row.names = FALSE) #create to run on HPC

#run on HPC and then come back with results

files_classic <- list.files("./HPC/shuf_within_layers_islands/csvs_classic/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_classic <- read_csv(files_classic) %>% bind_rows() #create a long edge list with all the csvs

