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

#this portion of code turns the data into a multilayer network with 14 layers (sites) and does
#modularity analysis


##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
dryad_intralayer <- read.csv("./csvs/intralayer_file.csv")
dryad_interlayer <- read.csv("./csvs/interlayer_file.csv") #already has inverted within

## ----multilayer_intra-----------------------------------------------------------------------------------------------
dryad_matrices <- NULL
for (layer in 1:14){
  d <- suppressMessages(read_excel('all_sites.xlsx', sheet = layer+2))
  web <- data.matrix(d[,2:ncol(d)])#convert variables to numeric 
  rownames(web) <- as.data.frame(d)[,1]
  web[is.na(web)] <- 0
  dryad_matrices[[layer]] <- web
}

names_dryad_matrices <- c("WesternSahara1", "WesternSahara2", "Fuerteventura1", "Fuerteventura2",
                          "GranCanaria1", "GranCanaria2", "TenerifeSouth1", "TenerifeSouth2",
                          "TenerifeTeno1", "TenerifeTeno2",  "Gomera1", "Gomera2", "Hierro1", "Hierro2")

# Layer dimensions
sapply(dryad_matrices, dim)

## ----dryad intralayer interlayer both ways-------------------------------------------------------------------------------------
intralayer_inverted <- tibble(values= dryad_intralayer$layer_to, dryad_intralayer$node_to, dryad_intralayer$layer_from, 
                                  dryad_intralayer$node_from, dryad_intralayer$weight) #create an inverted copy for directed intralayers
colnames(intralayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight")

## ---- normalize intralayer weights (for plants and pols)--------------------------------------------------------------------------
#plants in from
tot_plant <- dryad_intralayer %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted <- dryad_intralayer %>% left_join(tot_plant) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)


#pols in from
tot_pol <- intralayer_inverted %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted <- intralayer_inverted %>% left_join(tot_pol) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight) 

## ----multilayer_extended_final--------------------------------------------------------------------------------------
edgelist_intralayers_both <- bind_rows(intralayer_weighted, intralayer_weighted_inverted) #combine weighted version of intra with inter

dryad_edgelist_complete <- bind_rows(edgelist_intralayers_both, dryad_interlayer) #combine inverted and non inverted verions

## ----node_metadata--------------------------------------------------------------------------------------------------                                            
pollinators <- sort(unique(intralayer_weighted$node_to)) #adding up only pol who haven't been added yet 
plants <- sort(unique(intralayer_weighted$node_from)) #adding up only plants who haven't been added yet
intersect(pollinators, plants) #making sure I don't have plants in pol or other way around
A <- length(pollinators) # Number of pollinators
P <- length(plants) # Number of plants
S <- A+P

# Create a table with node metadata
physical_nodes <- tibble(node_id=1:S, #1 till the last species
                         type=c(rep('plant',P),rep('pollinator',A)), #replicate the words P and A times
                         species=c(plants,pollinators)) #add species from plants and pollinators in accordance
layer_metadata <- tibble(layer_id=c(1:14), layer_name=names_dryad_matrices)  #give num to each layer

#write.csv(layer_metadata, "./HPC/modularity/layer_metadata.csv", row.names = FALSE)
#write.csv(physical_nodes, "./HPC/modularity/physical_nodes.csv", row.names = FALSE)
#write.csv(layer_metadata, "./csvs/layer_metadata.csv", row.names = FALSE)
#write.csv(physical_nodes, "./csvs/physical_nodes.csv", row.names = FALSE)
#physical_nodes <- read.csv("./csvs/physical_nodes.csv")
#layer_metadata <-read.csv("./csvs/layer_metadata.csv")

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

#write.csv(dryad_edgelist_complete_ids, "./csvs/dryad_edgelist_complete_ids.csv", row.names = FALSE)

# Sanity checks:
dryad_edgelist_complete_ids %>% 
  filter(layer_from==layer_to) %>% #show only intra
  summarise_at(vars(node_from, node_to), .funs = c(mn=min,mx=max)) #apply funs of all nodes

# Number of intralayer edges
dryad_edgelist_complete_ids %>% 
  filter(layer_from==layer_to) %>% count()
nrow(intralayer_weighted)

# Number of interlayer edges
dryad_edgelist_complete_ids %>% 
  filter(layer_from!=layer_to) %>% count()
nrow(intralayer_weighted)

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

#write.csv(intra_nonextended, "./csvs/dryad_only_intralayer_edges.csv", row.names = FALSE)
#write.csv(inter_extended, "./csvs/dryad_only_interlayer_edges.csv", row.names = FALSE)
#inter_extended <- read.csv("./csvs/dryad_only_interlayer_edges.csv")
#intra_nonextended <- read.csv("./csvs/dryad_only_intralayer_edges.csv")

#create modules for empirical network
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")

modules_dryad_multilayer <- modified_multi(dryad_multilayer, 
                                             infomap_executable = "Infomap",
                                             flow_model = 'directed',
                                             relax = F, 
                                             silent = T, 
                                             trials = 100,
                                             seed = 497294, 
                                             temporal_network = F)

#write_csv(dryad_edgelist_complete_ids, './csvs/dryad_multilayer_edgelist.csv') 
#write_csv(dryad_multilayer$nodes, './csvs/dryad_multilayer_nodes.csv')
#write_csv(dryad_multilayer$layers, './csvs/dryad_multilayer_layers.csv')
#write_csv(modules_dryad_multilayer$modules, './csvs/modules_dryad_multilayer.csv')

modules <- read.csv('./csvs/modules_dryad_multilayer.csv')

#---- interlayer and intralayer distribution-----------------------------------------------------------------------------------------------------
intra_inter_data_for_distibution <- data.frame(values= c(intralayer_weighted$weight, 
                                                         intralayer_weighted_inverted$weight, 
                                                   inter_extended$weight),
                                         group= c(rep("intra plants", nrow(intralayer_weighted)), 
                                                  rep("intra pollinators", nrow(intralayer_weighted_inverted)),
                                                  rep("inter", nrow(inter_extended))))

pdf('./graphs/basic_analysis/intralayer_interlayer_directed_distribution.pdf', 10, 6)
intra_inter_data_for_distibution %>%
  ggplot(aes(x=values, fill=group))+ 
  geom_histogram(position= "identity", alpha= 0.6, color= "black")+ theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()



## ----distance decay in species--------------------------------------------------------------------------------------
#similarity check 2 furthest apart
modules_for_similarity_num <- modules_dryad_multilayer$modules %>% select(module, layer_id) %>% 
  unique() %>% group_by(module) %>% select(module) %>% unique()
modules_for_similarity <- modules_dryad_multilayer$modules %>%
  filter(module %in% modules_for_similarity_num$module) #only save the modules that are found in 2 or more layers

#pivot modules
module_pivoted <- pivot_by_module(modules_for_similarity)

#write.csv(module_pivoted, "./csvs/module_pivoted_for_state_node_similarity.csv")

#distances of data
distances <- read.csv("./csvs/distances_file.csv") #read the file that contains all geographical places with distances
distances_with_ids <- distances %>% left_join(layer_metadata, by= c("layer_from"="layer_name")) %>% 
  left_join(layer_metadata, by= c("layer_to"="layer_name")) %>% #add correct id to layer name
  select(distance_in_meters, layer_id.x, layer_id.y) #discard actual names of layers
names(distances_with_ids)[2] <- "layer_from" 
names(distances_with_ids)[3] <- "layer_to"

mean(distances_with_ids$distance_in_meters)
median(distances_with_ids$distance_in_meters)

#write.csv(distances_with_ids, "./csvs/distances_with_ids.csv", row.names = FALSE)
#distances_with_ids <- read.csv("./csvs/distances_with_ids.csv")

# classic distnace decay
all_species_all_layers <- rbind(tot_plant, tot_pol) %>% inner_join(physical_nodes, by= c("node_from" = "species")) %>% 
  inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>% subset(select = -c(layer_from, node_from, tot, type)) 
#change node names to ids and layer names to ids and remove unwanted columns
  
classic_layers_turnover <- NULL

for (i in (1:13)){
  for (j in ((i+1):14)){
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

#write.csv(classic_layers_turnover_with_distances, "./csvs/classic_layers_turnover_with_distances.csv")
#classic_layers_turnover_with_distances <- read.csv("./csvs/classic_layers_turnover_with_distances.csv")

classic_layers_turnover_with_distances <- classic_layers_turnover_with_distances%>% mutate(distance_in_km=distance_in_meters/1000)

pdf('./graphs/basic_analysis/distance_decay_species.pdf', 10, 6)
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



## ----distance decay in modules--------------------------------------------------------------------------------------
modules_edge_list <- NULL 

for (i in (1:nrow(module_pivoted))){ #run the function for each row (module) in the data frame
  modules_edge_list <- edge_list_per_module(module_pivoted[i,], modules_edge_list) 
  current_module <- rownames(module_pivoted)[i]
  modules_edge_list <- modules_edge_list %>% mutate(module = replace_na(module, current_module)) #add module number
}

#write.csv(modules_edge_list, "csvs/modules_edge_list.csv", row.names = FALSE)

edge_list_with_distances <- right_join(modules_edge_list, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances <- na.omit(edge_list_with_distances) #remove NA and delete layer name

#arrange data to include coordinates and modules sizes
size <- count(modules_dryad_multilayer$modules, module)  #create a data frame of all modules and how many nodes are in each (size of module)
module_data <- merge(modules_dryad_multilayer$modules , size, by=c("module","module")) #merge size of module with all the other info about the modules
colnames(module_data)[7] <- "size_of_module" #rename column

lon_lat_data <- read_csv('./csvs/layers.csv') #create new data frame with just the layer data
lon_lat_data <- lon_lat_data %>% select(c("layer_id","lat","Lon")) %>% na.omit()  #only select layer id and coordinates

#write.csv(lon_lat_data, "csvs/lon_lat_data.csv", row.names = FALSE)

module_data_with_loc <- merge(module_data, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon <- module_data_with_loc %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon$count <- c(1)

#write.csv(modules_with_lat_lon, "csvs/modules_with_lat_lon.csv", row.names = FALSE)

#---- jaccard on layers---------------------------------------------------------------------------------------------------
module_layer_turnover <- NULL

for (i in 1:14){
  for (j in (1+i):13){
    print(i)
    modules_in_layer_from <- filter(modules_with_lat_lon, layer_id == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_layer_to <- filter(modules_with_lat_lon, layer_id == j) %>% select(module) %>% unique() %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
    int_both <- intersect(modules_in_layer_from, modules_in_layer_to) #how many nodes are found in both layers
    uni_both <- union(modules_in_layer_from, modules_in_layer_to)
    turnover <- length(int_both)/length(uni_both)
    module_layer_turnover <- rbind(module_layer_turnover, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

#write.csv(module_layer_turnover, "csvs/module_layer_turnover.csv", row.names = FALSE)

layers_turnover_with_distnace_empirical <- edge_list_with_distances %>%
  merge(module_layer_turnover, by= c("layer_from", "layer_to")) #merge both versions 

layers_turnover_with_distnace_empirical <- layers_turnover_with_distnace_empirical %>% 
  mutate(distance_in_km=distance_in_meters/1000) %>% select(layer_from, layer_to, turnover, distance_in_km) %>%
  unique() #turn to km

#write.csv(layers_turnover_with_distnace_empirical, "csvs/layers_turnover_with_distnace_empirical.csv", row.names = FALSE)

pdf('./graphs/modularity_analysis/distance_decay_in_modules_sites.pdf', 10, 6)
layers_turnover_with_distnace_empirical %>%
  ggplot(aes(x=distance_in_km, y=turnover))+
  geom_point(color = "indianred2")+ scale_x_continuous()+theme_classic()+ 
  stat_smooth(method= "lm", se = F, color = "indianred2")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Km", y="Jaccard Similarity")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()


#-------making sure its significant---------------------------------------------
layers_turnover_with_distnace_empirical <- read.csv("./csvs/layers_turnover_with_distnace_empirical.csv")
lm1 = lm(turnover ~ distance_in_km ,data=layers_turnover_with_distnace_empirical) #in empirical
summary(lm1) 
#p < 0.001, R squared = 0.549


#KEEP IT FOR NOW (SEE WHAT SHAI SAYS)
#jaccard similarity on map 

#combine turnover data frame with coordinates
#empirical_turnover_for_module_island_shuf_no_self_loop_km <- read.csv("./csvs/empirical_turnover_for_module_island_shuf_no_self_loop_km.csv")
#islands_with_lon_lat_dif <- read.csv("./csvs/islands_with_lon_lat_dif.csv")

#join both data frames
#jaccard_similarity_on_map <- merge(empirical_turnover_for_module_island_shuf_no_self_loop_km, islands_with_lon_lat_dif, 
#                                   by = c("layer_from", "layer_to"))

#worldmap <- map_data("world")
#dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990))  
#dryad_map_jaccard <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+ #requires an API key now
#  geom_segment(aes(x= x, xend= xend, y= y, yend= yend, color = ave),
#               data= jaccard_similarity_on_map)+ scale_color_gradient(high="red",low="lightskyblue1")+
#  geom_point(aes(x=x, y=y),shape= 21 , fill= 'white', color= "black", data= lat_lon_nodes)+ 
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
#  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+ 
#  labs(x = "Longitude", y = "Latitude")+ labs(color = "Jaccard Similarity")

#dryad_map_jaccard

##---- edge list for layers and not islands-------------------------------------------------------------------------
#version with # of modules in layers
edge_list_by_layer_modules <- edge_list_with_distances %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #doesn't really change the distance as its layers
edge_list_by_layer_modules$count <- c(1)
edge_list_by_layer_modules <- edge_list_by_layer_modules %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules) 

#version with correct average between layers
edge_list_by_layers_ave <- edge_list_with_distances %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_layers_combine <- edge_list_by_layers_ave %>%
  merge(edge_list_by_layer_modules, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_layers_combine_no_module <- edge_list_layers_combine %>% select(-module) %>% unique() #have version where modules aren't present

#write.csv(edge_list_layers_combine_no_module, "./csvs/edge_list_layers_combine_no_module.csv")
