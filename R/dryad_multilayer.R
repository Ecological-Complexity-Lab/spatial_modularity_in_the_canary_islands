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
#getwd()
dryad_intralayer <- read.csv("./csvs/intralayer_file.csv")
#print(dryad_intralayer)
dryad_interlayer <- read.csv("./csvs/interlayer_file.csv") #already has inverted within
#print(dryad_interlayer)

## ----multilayer_intra-----------------------------------------------------------------------------------------------
dryad_matrices <- NULL
for (layer in 1:14){
  d <- suppressMessages(read_excel('all_sites.xlsx', sheet = layer+2))
  web <- data.matrix(d[,2:ncol(d)])
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

## ---- normalize intralayer weights--------------------------------------------------------------------------
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

# Distribution of edge weights
dryad_edgelist_complete %>% 
  mutate(edge_type=ifelse(layer_from==layer_to,'Intra','Inter')) %>% #if same species in to and from edge_type is intra else it's inter
  ggplot()+
  geom_histogram(aes(x=weight, fill=edge_type), alpha=0.6)+ #create histogram of edge_type and weight
  facet_wrap(~layer_from, scales = 'free')+ #create multi-panel plot aligning by layer_from
  scale_fill_manual(values = c('navy','plum'))+theme_classic() 

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


##---- basic analysis for layer and island --------------------------------------------------------------------------------
#richness in each layer
tot_plant_layer_ids <- tot_plant %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>%
  subset(select = -layer_from) %>% count(layer_id)
tot_plant_layer_ids$type <- "plant"

tot_pol_layer_ids <- tot_pol %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>%
  subset(select = -layer_from) %>% count(layer_id)
tot_pol_layer_ids$type <- "pollinator"

richness_in_layer <- rbind(tot_plant_layer_ids, tot_pol_layer_ids)

richness_in_layer %>%
  ggplot(aes(x= layer_id, y=n, fill=type))+ geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ 
  theme_classic()+ scale_x_continuous(breaks=seq(1,14,1))+ labs(x= "layer id", y= "number of species")

#richness in each island
richness_in_island <- richness_in_layer
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
richness_in_island$layer_id[richness_in_island$layer_id %in% old] <- new[match(richness_in_island$layer_id, old)]

richness_in_island <- richness_in_island %>% group_by(layer_id,type) %>% 
  summarise(num_of_species=sum(n)) #sum every two layers in an island by type

richness_in_island %>%
  ggplot(aes(x= layer_id, y=num_of_species, fill=type))+ geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ 
  theme_classic()+ labs(x= "island id", y= "number of species")


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

#write.csv(intra_nonextended, "./csvs/dryad_only_intralayer_edges.csv")
#write.csv(inter_extended, "./csvs/dryad_only_interlayer_edges.csv")

#create modules for empirical network
modules_dryad_multilayer <- run_infomap_multilayer(dryad_multilayer, 
                                             infomap_executable = "../Infomap",
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

#general info about modules
num_of_nodes_in_module <- modules_dryad_multilayer$modules %>% count(module) #num of nodes in module
local_modules <- modules_dryad_multilayer$modules %>% select(module, layer_id)
num_of_layers_in_module <- local_modules %>% distinct() %>% count(module) #num of layers a module is found in

num_of_layers_in_module %>%
  ggplot(aes(x= module, y=n))+ geom_bar(stat="identity")+ #stacked
  theme_classic()+ scale_x_continuous(breaks=seq(1,43,2))+ labs(x= "Module", y= "Number of Layers")

#----for each layer how many nodes are found in each module-----------------------------------------------------------

## modules
modules_dryad_multilayer$modules %>%
  group_by(layer_id, module) %>%
  summarise(n=n_distinct(node_id)) %>%
  filter(n>2) %>%
  ggplot(aes(layer_id, module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text()+ xlab("Layer Id") + ylab("Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ scale_x_continuous(breaks=seq(1,14,1))+
  scale_y_continuous(breaks=seq(1,43,2))+
  scale_fill_gradient(low= "lightskyblue1", high= "red")

#---- which physical node is found in which module---------------------------------------------------------------------

modules_dryad_multilayer_species_analysis <- modules_dryad_multilayer$modules

modules_dryad_multilayer_species_analysis <- modules_dryad_multilayer_species_analysis %>% select(node_id, layer_id, module) %>% unique() #take only certain columns
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
modules_dryad_multilayer_species_analysis$layer_id[modules_dryad_multilayer_species_analysis$layer_id %in% old] <- 
  new[match(modules_dryad_multilayer_species_analysis$layer_id, old)]
modules_dryad_multilayer_species_analysis <- modules_dryad_multilayer_species_analysis %>% unique() #delete repeats because of layer to island conversion

modules_dryad_multilayer_species_analysis_module <- modules_dryad_multilayer_species_analysis %>% group_by(node_id, module) %>%
  mutate(number_of_islands = n()) %>% select(node_id, module, number_of_islands) %>% unique() #count number of islands each physical node is found in

modules_dryad_multilayer_species_analysis_module %>% ggplot(aes(x = module, y = node_id, fill = number_of_islands))+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ scale_x_continuous(breaks=seq(1,43,2))+
  scale_y_continuous(breaks = seq(1,288,10))+
  geom_tile()+ theme_classic()+ scale_fill_distiller(palette = "RdPu")

#cutoff at minimum 2 islands
modules_dryad_multilayer_species_analysis_module_min_2 <- modules_dryad_multilayer_species_analysis_module %>% filter(number_of_islands > 1)

modules_dryad_multilayer_species_analysis_module_min_2 %>% ggplot(aes(x = module, y = node_id, fill = number_of_islands))+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ scale_x_continuous(breaks=seq(1,43,2))+
  scale_y_continuous(breaks = seq(1,288,10))+
  geom_tile()+ theme_classic()+ scale_fill_distiller(palette = "RdPu")

modules_dryad_multilayer_species_analysis_module_max_1 <- modules_dryad_multilayer_species_analysis_module %>% 
  filter(number_of_islands == 1) %>% ungroup() %>% select(node_id) %>% unique() #find out how many physical nodes are local to 1 island

# reorder y axis
modules_dryad_multilayer_species_analysis_module %>% ggplot(aes(x = module, y = reorder(node_id, number_of_islands), 
                                                                fill = number_of_islands))+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ scale_x_continuous(breaks=seq(1,43,2))+
  geom_tile()+ theme_classic()+ scale_fill_distiller(palette = "GnBu")

##as matrix
modules_dryad_multilayer_species_analysis_module_matrix <- dcast(modules_dryad_multilayer_species_analysis_module, 
                                                                 node_id ~ module, value.var = "number_of_islands")

modules_dryad_multilayer_species_analysis_module_matrix <- as.matrix(modules_dryad_multilayer_species_analysis_module_matrix[,-1])

heatmap(modules_dryad_multilayer_species_analysis_module_matrix, Colv = NA, Rowv = NA)

#---- interlayer and intralayer distribution-----------------------------------------------------------------------------------------------------
#interlayer distribution
weight_distribution <- inter_extended %>% arrange(weight) %>% select(weight) 
weight_distribution %>%
  ggplot(aes(x=weight))+geom_density(fill = "#F47069", color = "#F47069", alpha = 0.4)+theme_classic()+ 
  geom_vline(xintercept = median(weight_distribution$weight), linetype = "dashed", color = "black", size = 1)+
  theme(axis.title=element_text(size=15))

mean(unlist(weight_distribution))
median(unlist(weight_distribution))

#intralayer distribution
intra_weight_distribution <- intra_nonextended %>% arrange(weight) %>% select(weight) 
intra_weight_distribution %>%
  ggplot()+geom_histogram(aes(x=weight))+theme_classic()

mean(unlist(intra_weight_distribution))
median(unlist(intra_weight_distribution))

#inter and intra non directed together
inter_intra_non_directed <- data.frame(values= c(intra_nonextended$weight, inter_extended$weight), 
                                       group= c(rep("intra", nrow(intra_nonextended)), rep("inter", nrow(inter_extended))))

inter_intra_non_directed %>%
  ggplot(aes(x=values, fill=group))+ geom_histogram(position= "identity", alpha= 0.6, color= "black")+ theme_classic()

#directed intralayer weights distribution
directed_weight_distribution_plants <- intralayer_weighted %>% arrange(weight) %>% select(weight)
directed_weight_distribution_plants %>%
  ggplot()+ geom_histogram(aes(x=weight))+theme_classic()

mean(unlist(directed_weight_distribution_plants))
median(unlist(directed_weight_distribution_plants))

directed_weight_distribution_pols <- intralayer_weighted_inverted %>% arrange(weight) %>% select(weight)
directed_weight_distribution_pols %>%
  ggplot()+ geom_histogram(aes(x=weight))+theme_classic()

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
  ggplot(aes(x=values, fill=group))+ geom_histogram(position= "identity", alpha= 0.6, color= "black")+ theme_classic()

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

write.csv(classic_layers_turnover_with_distances, "./csvs/classic_layers_turnover_with_distances.csv")

classic_layers_turnover_with_distances <- classic_layers_turnover_with_distances%>% mutate(distance_in_km=distance_in_meters/1000)

classic_layers_turnover_with_distances %>%
  ggplot(aes(x=distance_in_km, y=turnover))+ geom_point(color = "indianred2")+ theme_classic()+ 
  stat_smooth(method= "lm", se=F, color = "indianred2")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400, color = "indianred2")+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.4), color = "indianred2") 



modules_edge_list <- NULL

for (i in (1:nrow(module_pivoted))){ #run the function for each row in the data frame
  modules_edge_list <- edge_list_per_module(module_pivoted[i,], modules_edge_list) 
  current_module <- rownames(module_pivoted)[i]
  modules_edge_list <- modules_edge_list %>% mutate(module = replace_na(module, current_module)) #add module number
}

#view(modules_edge_list)


edge_list_with_distances <- right_join(modules_edge_list, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances <- na.omit(edge_list_with_distances) #remove NA and delete layer name

#arrange data to include coordinates and modules sizes
size <- count(modules_dryad_multilayer$modules, module)  #create a data frame of all modules and how many nodes are in each (size of module)
module_data <- merge(modules_dryad_multilayer$modules , size, by=c("module","module")) #merge size of module with all the other info about the modules
colnames(module_data)[7] <- "size_of_module" #rename column

lon_lat_data <- read_csv('./csvs/layers.csv') #create new data frame with just the layer data
lon_lat_data <- lon_lat_data %>% select(c("layer_id","lat","Lon")) %>% na.omit()  #only select layer id and coordinates

module_data_with_loc <- merge(module_data, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon <- module_data_with_loc %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon$count <- c(1)

modules_with_lat_lon %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_id)))+ geom_bar(stat= "identity")+ theme_classic()+
  scale_x_continuous(breaks=seq(1,43,2))+ labs(y="number of physical nodes", x="module number")+
  guides(fill=guide_legend(title="layer\nnumber"))

#how many islands are within a module
modules_with_lat_lon_islands <- module_data_with_loc %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
modules_with_lat_lon_islands$layer_id[modules_with_lat_lon_islands$layer_id %in% old] <- new[match(modules_with_lat_lon_islands$layer_id, old)]
modules_with_lat_lon_islands$count <- c(1)

modules_with_lat_lon_islands %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_id)))+ geom_bar(stat= "identity")+ theme_classic()+
  scale_x_continuous(breaks=seq(1,43,2))+ labs(y="number of physical nodes", x="module number")+ 
  guides(fill=guide_legend(title="island\nnumber"))

#modules similarity pairwise distance between islands
edge_list_by_islands <- edge_list_with_distances
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
edge_list_by_islands$layer_from[edge_list_by_islands$layer_from %in% old] <- new[match(edge_list_by_islands$layer_from, old)]
edge_list_by_islands$layer_to[edge_list_by_islands$layer_to %in% old] <- new[match(edge_list_by_islands$layer_to, old)]

#version with # of modules in layers
edge_list_by_islands_modules <- edge_list_by_islands %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_islands_modules$count <- c(1)
edge_list_by_islands_modules <- edge_list_by_islands_modules %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules) 

#version with correct average between layers
edge_list_by_islands_ave <- edge_list_by_islands %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_island_combine <- edge_list_by_islands_ave %>%
  merge(edge_list_by_islands_modules, by= c("layer_from", "layer_to")) #merge both versions 
  
edge_list_island_combine_no_module <- edge_list_island_combine %>% select(-module) %>% unique() #have version where modules aren't present


#same but for layers
#version with # of modules in layers
edge_list_by_layers_modules <- edge_list_with_distances %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_layers_modules$count <- c(1)
edge_list_by_layers_modules <- edge_list_by_layers_modules %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules) 

#version with correct average between layers
edge_list_by_layers_ave <- edge_list_with_distances %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_layer_combine <- edge_list_by_layers_ave %>%
  merge(edge_list_by_layers_modules, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_layer_combine_no_module <- edge_list_layer_combine %>% select(-module) %>% unique() #have version where modules aren't present
#---- graphs modules in common -----------------------------------------------------------------------

# number of modules in common as func of distance between islands
edge_list_island_combine_no_module %>%
  ggplot(aes(x=ave_distance, y=number_of_modules))+
  geom_point()+ scale_x_continuous(breaks=seq(0,455736.67290,100000))+theme_classic()+ stat_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="number of modules in common")

# heatmap
edge_list_by_islands_heatmap <- edge_list_island_combine_no_module %>% select(layer_from, layer_to, number_of_modules) 
edge_list_by_islands_heatmap %>% 
  ggplot()+ geom_tile(aes(x=layer_from, y=layer_to, fill= number_of_modules))+ theme_classic()+
  scale_fill_gradient(low= "lightskyblue1", high= "red")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="island from", y="island to", fill= "number of modules \n in common")

# network graph
network_graph <- graph_from_data_frame(edge_list_by_islands_heatmap, directed= FALSE)
network_graph %>%
  plot(edge.width=edge_list_by_islands_heatmap$number_of_modules)

#network graph on world map
lat_lon_for_graph <- modules_with_lat_lon 
lat_lon_for_graph <- lat_lon_for_graph %>% select(layer_id, lat, Lon)

old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
lat_lon_for_graph$layer_id[lat_lon_for_graph$layer_id %in% old] <- new[match(lat_lon_for_graph$layer_id, old)]

islands_with_lat_lon <- edge_list_island_combine_no_module %>%
  inner_join(lat_lon_for_graph, by= c("layer_from" = "layer_id")) %>%  #add edge list to coordinates
  rename(x= Lon, y=lat) %>%
  inner_join(lat_lon_for_graph, by= c("layer_to" = "layer_id")) %>%
  rename(xend= Lon, yend= lat) %>%
  select(layer_from, layer_to, number_of_modules, y, yend, x, xend) %>% unique() 

lat_lon_nodes <- islands_with_lat_lon %>% group_by(layer_from) %>% arrange(desc(layer_to)) %>% slice(1) #save coordinates of islands

islands_with_lon_lat_dif <- islands_with_lat_lon %>% filter(layer_from != layer_to)

#write.csv(islands_with_lon_lat_dif, "./csvs/islands_with_lon_lat_dif.csv", row.names = FALSE)

worldmap <- map_data("world")
dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990))  
dryad_map <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+
  geom_segment(aes(x= x, xend= xend, y= y, yend= yend, color = number_of_modules),
             data= islands_with_lon_lat_dif)+ scale_color_gradient(high="red",low="lightskyblue1")+
  geom_point(aes(x=x, y=y),shape= 21 , fill= 'white', color= "black", data= lat_lon_nodes)+ 
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+ 
  labs(x = "Longitude", y = "Latitude")+ labs(color = "number of modules \n in common")

print(dryad_map)
print(top_10_scatterpies)
top_10_scatterpies_no_legend <- top_10_scatterpies + theme(legend.position = "none")

combine_map <- ggdraw()+ draw_plot(dryad_map)+
  draw_plot(top_10_scatterpies_no_legend, x = -0.07, y = 0.1, scale = 0.70)
print(combine_map) 


#---- jaccard on islands---------------------------------------------------------------------------------------------------
module_island_turnover <- NULL

island_list <- c("1","2","3","4east","4west","5","6")

for (i in island_list){
  for (j in island_list){
    print(i)
    modules_in_island_from <- filter(edge_list_by_islands_modules, layer_from == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_island_to <- filter(edge_list_by_islands_modules, layer_from == j) %>% select(module) %>% unique() %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
    int_both <- intersect(modules_in_island_from, modules_in_island_to) #how many nodes are found in both layers
    uni_both <- union(modules_in_island_from, modules_in_island_to)
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover <- rbind(module_island_turnover, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

edge_list_by_islands_ave <- edge_list_by_islands %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

islands_turnover_with_distnace_empirical <- edge_list_by_islands_ave %>%
  merge(module_island_turnover, by= c("layer_from", "layer_to")) #merge both versions 

##same but with layers
module_layer_turnover <- NULL

for (i in 1:14){
  for (j in 1:14){
    print(i)
    modules_in_layer_from <- filter(edge_list_by_layers_modules, layer_from == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_layer_to <- filter(edge_list_by_layers_modules, layer_from == j) %>% select(module) %>% unique() %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
    int_both <- intersect(modules_in_layer_from, modules_in_layer_to) #how many nodes are found in both layers
    uni_both <- union(modules_in_layer_from, modules_in_layer_to)
    turnover <- length(int_both)/length(uni_both)
    module_layer_turnover <- rbind(module_layer_turnover, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

edge_list_by_layers_ave <- edge_list_with_distances %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

layers_turnover_with_distnace_empirical <- edge_list_by_layers_ave %>%
  merge(module_layer_turnover, by= c("layer_from", "layer_to")) #merge both versions 


islands_turnover_with_distnace_empirical %>%
  ggplot(aes(x=ave_distance, y=turnover))+
  geom_point()+ scale_x_continuous(breaks=seq(0,455736.67290,100000))+theme_classic()+ stat_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Meters", y="Jaccard Similarity")

#jaccard similarity on map

#combine turnover data frame with coordinates
empirical_turnover_for_module_island_shuf_no_self_loop_km <- read.csv("./csvs/empirical_turnover_for_module_island_shuf_no_self_loop_km.csv")
islands_with_lon_lat_dif <- read.csv("./csvs/islands_with_lon_lat_dif.csv")

#join both data frames
jaccard_similarity_on_map <- merge(empirical_turnover_for_module_island_shuf_no_self_loop_km, islands_with_lon_lat_dif, 
                                   by = c("layer_from", "layer_to"))

worldmap <- map_data("world")
dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990))  
dryad_map_jaccard <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+
  geom_segment(aes(x= x, xend= xend, y= y, yend= yend, color = ave),
               data= jaccard_similarity_on_map)+ scale_color_gradient(high="red",low="lightskyblue1")+
  geom_point(aes(x=x, y=y),shape= 21 , fill= 'white', color= "black", data= lat_lon_nodes)+ 
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+ 
  labs(x = "Longitude", y = "Latitude")+ labs(color = "Jaccard Similarity")

dryad_map_jaccard
##---- edge list for layers and not islands-------------------------------------------------------------------------

#version with # of modules in layers
edge_list_by_layer_modules <- edge_list_with_distances %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
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


#---- layers furthest apart---------------------------------------------------------------
layers_furthest_apart <- edge_list_with_distances %>% group_by(module) %>% slice(which.max(distance_in_meters))  #get two furthest layers in each modul
layers_furthest_apart <- data.frame(lapply(layers_furthest_apart, as.numeric))
layers_furthest_apart_km <- layers_furthest_apart %>% mutate(distance_in_meters = replace(distance_in_meters, distance_in_meters>0, 
                                                                                          distance_in_meters/1000))
names(layers_furthest_apart_km)[4] <- "distance_in_km"

layers_furthest_apart %>%
  ggplot()+geom_histogram(aes(x=distance_in_meters))+theme_classic()+ scale_x_continuous(breaks=seq(0,455736.67290,50000))+
  labs(x= "distnace in meters", y= "number of modules")

layers_furthest_apart_km %>%
  ggplot()+geom_histogram(aes(x=distance_in_km))+theme_classic()+ scale_x_continuous(breaks=seq(0,457,50))

#distance decay species turnover with regards to modules
layers_turnover <- NULL

for (i in (1:nrow(module_pivoted))){
  focal_module <- filter(edge_list_with_distances, module == i) #look at one module at a time
  for (j in (1:nrow(focal_module))){
    module <- focal_module[j,]
    current_module <- module$module
    current_layer <- module$layer_from
    current_layer_to <- module$layer_to
    current_distance <- module$distance_in_meters
    physical_nodes_in_layer_from <- filter(modules_for_similarity, layer_id == current_layer) %>% select(node_id) %>%
      unlist() #take the nodes that are found in the module in layer_from
    physical_nodes_in_layer_to <- filter(modules_for_similarity, (layer_id == current_layer_to)) %>% select(node_id) %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
    int_both <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in both layers
    uni_both <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to)
    turnover <- length(int_both)/length(uni_both)
    layers_turnover <- rbind(layers_turnover, tibble(layer_from= current_layer, layer_to= current_layer_to, turnover= turnover, 
                                                     distance= current_distance))
  }
}

layers_turnover <- layers_turnover %>% unique()

ggplot(layers_turnover, aes(x=distance, y=turnover))+
  geom_point()+ scale_x_continuous(breaks=seq(0,455736.67290,100000))+theme_classic()+ stat_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="similarity")


#module size distribution
plants_and_pols <- modules_dryad_multilayer$modules %>% count(type) 
plants <- plants_and_pols[1,2]
pols <- plants_and_pols[2,2]
modules_count <- modules_dryad_multilayer$modules %>% count(module, type) #count how many plants and pollinators are found in each module
modules_count_not_proportion <- modules_count
for(i in (1:nrow(modules_count))){
  if (modules_count$type[i] == "plant"){
    modules_count$n[i] <- modules_count$n[i]/plants
  }
  else{
    modules_count$n[i] <- modules_count$n[i]/pols
  }
}

#distribution by species
modules_count_not_proportion %>% ggplot(aes(x=module, y=n, fill= type))+
  geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ theme_classic()+
  scale_x_continuous(breaks=seq(1,43,2))+ labs(x= "module number", y= "number of species")

#distribution by species and proportion
modules_count %>% ggplot(aes(x=as.numeric(module), y=as.numeric(n), fill= type))+
  geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ theme_classic()+
  scale_x_continuous(breaks=seq(1,43,2))+ labs(x= "module number", y= "proportion")

#distribution by module size
module_sizes <- modules_dryad_multilayer$modules %>% count(module)
mean(module_sizes$n)
median(module_sizes$n)

#number of modules who span over x layers
modules_over_i <- NULL
for (i in 1:max(modules_dryad_multilayer$modules$module)){
  modules_for_similarity_over_i_num <- modules_dryad_multilayer$modules %>% select(module, layer_id) %>%
    unique() %>% group_by(module) %>% filter(n()>=i) %>% select(module) %>% unique() #check which modules occur in i or more layers
  modules_for_similarity_over_i <- modules_dryad_multilayer$modules %>%
    filter(module %in% modules_for_similarity_over_i_num$module) %>% count(module)#only save the modules that are found in i or more layers
  modules_over_i <- rbind(modules_over_i, tibble(module= modules_for_similarity_over_i$module, i=i))
}

num_of_modules_over_i <- modules_over_i %>% count(i)
num_of_modules_over_i %>% ggplot(aes(x=i, y=n))+geom_bar(stat="identity")+theme_classic()+
  scale_x_continuous(breaks=seq(1,14,1))+ 
  labs(x="number of layers", y="number of modules found in x or more layers")


#number of module which are found in x layers
modules_in_i <- NULL
for (i in 1:max(modules_dryad_multilayer$modules$module)){
  modules_for_similarity_in_i_num <- modules_dryad_multilayer$modules %>% select(module, layer_id) %>%
    unique() %>% group_by(module) %>% filter(n()==i) %>% select(module) %>% unique() #check which modules occur in i or more layers
  modules_for_similarity_in_i <- modules_dryad_multilayer$modules %>%
    filter(module %in% modules_for_similarity_in_i_num$module) %>% count(module)#only save the modules that are found in i or more layers
  modules_in_i <- rbind(modules_in_i, tibble(module= modules_for_similarity_in_i$module, i=i))
}

num_of_modules_in_i <- modules_in_i %>% count(i)
num_of_modules_in_i %>% ggplot(aes(x=i, y=n))+geom_bar(stat="identity")+theme_classic()+
  scale_x_continuous(breaks=seq(1,14,1))+ 
  labs(x="number of layers", y="number of modules found in x layers")

#modules in x layers and size
layers_and_sizes <- right_join(modules_in_i, module_sizes, by= c("module"="module")) #combine number of modules with size of module

main_bar_plot <- ggplot(layers_and_sizes, aes(x=i, y=n, fill=factor(module)))+ 
  geom_bar(stat="identity", color="black", show.legend= FALSE)+ theme_classic()+
  scale_x_continuous(breaks=seq(1,14,1))+ 
  labs(x="number of layers", y="number of species")

print(main_bar_plot)

#num of layers a module is found in as a function of the distance
modules_layers_vs_distance <- right_join(layers_furthest_apart, modules_in_i, by= c("module"="module"))
modules_layers_vs_distance[is.na(modules_layers_vs_distance)] <- 0
modules_layers_vs_distance %>% 
  ggplot(aes(x=distance_in_meters, y=i))+ geom_point()+ scale_x_continuous(breaks=seq(0,455736.67290,100000))+ theme_classic()+
  scale_y_continuous(breaks=seq(1,14,1))+ labs(x="distance in meters", y="number of layers")+ geom_smooth(method= "lm")

#---- module partners-----------------------------------------------------------------------------------------------------------

#jaccard similarity
jaccard <- function(a,b){ #recieves two lists a and b
  intersection <- length(intersect(a,b)) #intersection in how many species are found both in a and b
  union <- length(union(a,b)) #union in how many are found in both a and b together
  return (intersection/union)
}

#finding all node's partners in modules across all layers (empirical)
#partners_empirical <- NULL
#partners_empirical_list <- NULL


#for (j in 1:length(modules_dryad_multilayer$modules$node_id)){
#  current_node <- modules_dryad_multilayer$modules$node_id[j] #current node i'm looking at as a view point
#  filter_by_node <- filter(modules_dryad_multilayer$modules, modules_dryad_multilayer$modules$node_id==current_node) #only rows where node is
#  modules_with_species <- unique(filter_by_node$module) #only take unique module numbers
#  module_partners <- filter(modules_dryad_multilayer$modules, modules_dryad_multilayer$modules$module==modules_with_species)
#  just_partners <- subset(module_partners, node_id != current_node) #remove the current id from the module
#  unique_partners <- unique(just_partners$node_id) #unique partners in a specific module 
#  partners_empirical <- rbind(partners_empirical, tibble(node_id=current_node, partners=unique_partners)) #all nodes and all their partners
#  partners_empirical <- distinct(partners_empirical)
#}

#partners_empirical <- partners_empirical %>%
#  merge(physical_nodes, by="node_id") %>%
#  subset(select= -species) 

#num_of_partners <- partners_empirical %>% count(node_id)
#num_of_partners <- num_of_partners %>% count(n)
#num_of_partners %>% ggplot(aes(x=n, y=nn))+geom_bar(stat="identity")+theme_classic()+
#  labs(x= "number of partners per node", y= "count")+ scale_x_continuous(breaks=seq(1,74,2))

#mean(num_of_partners$n)
#median(num_of_partners$n)

#write.csv(partners_empirical, "./csvs/module_partners_in_empirical.csv")

## ---- jaccard between every two furthest layer per module --------------------------------------------------------------------------
modules_for_furthest_jaccard <- modules_dryad_multilayer$modules %>% subset(select= -c(flow,type,species))
furthest_from <- modules_for_furthest_jaccard %>% inner_join(layers_furthest_apart, by = c("module"= "module", "layer_id"="layer_from")) %>%
  rename(layer_from=layer_id) %>% na.omit() #all species in the module in layer from
furthest_to <- modules_for_furthest_jaccard %>% right_join(layers_furthest_apart, by = c("module"= "module", "layer_id"="layer_to")) %>%
  rename(layer_to=layer_id) %>% na.omit() #all species in the module in layer to

furthest_jaccard <- NULL

for (i in (1:43)){
  layer_from <- filter(furthest_from, module == i) %>% select(layer_from) %>% unlist()
  layer_to <- filter(furthest_from, module == i) %>% select(layer_to) %>% unlist()
  distance <- filter(furthest_from, module == i) %>% select(distance_in_meters) %>% unlist()
   physical_nodes_in_layer_from <- filter(furthest_from, module == i) %>% select(node_id) %>% unlist()
   physical_nodes_in_layer_to <- filter(furthest_to, (module == i)) %>% select(node_id) %>% unlist()
    #take all nodes in layer_from and all nodes in layer_to to check turnover
   int_both <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in both layers
   uni_both <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to)
  turnover <- length(int_both)/length(uni_both)
  furthest_jaccard <- rbind(furthest_jaccard, tibble(layer_from= layer_from, layer_to= layer_to, module= i, 
                                                     distance_in_meters = distance, turnover= turnover))
}

furthest_jaccard <- furthest_jaccard %>% unique() 

furthest_jaccard %>%
  ggplot(aes(x= distance_in_meters, y= turnover))+ geom_point()+ theme_classic()+
  labs(x="distnace in meters between furthest layers a module is found in", y="similarity in module composition")

##------------------------------------------------------------------------------------------------------------------------------
#10% biggest modules
top_10_percent <- select(module_data_with_loc, c("module", "size_of_module")) %>% distinct %>% arrange(desc(size_of_module)) %>%
filter(quantile(size_of_module, probs=0.9) < size_of_module) #arrange all modules to see 90 percentile and filter by them
top_10_percent$size_of_module <- top_10_percent$module #don't need the sizes anymore
names(top_10_percent)[names(top_10_percent)=="size_of_module"] <- "group" #change name to group- only the top percentile
pie_chart_data <- top_10_percent %>% right_join(module_data_with_loc, by=c("module"="module")) #join to have coordinates as well
pie_chart_data$group <- as.character(pie_chart_data$group) 
pie_chart_data <- pie_chart_data %>% drop_na()

pie_chart_data <- lon_lat_data %>% right_join(pivot_by_country(pie_chart_data), by = c("layer_id"="layer_id")) 
pie_chart_data <- pie_chart_data %>% group_by(lat, Lon) %>% summarise(across(everything(), sum))

#view(pie_chart_data)

worldmap <- map_data("world")
dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990))  
 inset_map <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+ #create a terrain map based on dryad
   geom_scatterpie(aes(x=Lon, y=lat, group=layer_id, r=0.028), alpha=0.6, data= pie_chart_data, cols=colnames(pie_chart_data[,c(4:8)]))+ #create a pie chart for every location with the 90 precentile
   coord_fixed(xlim=c(-18.542076, -12.58351), ylim=c(26, 30.323990)) #set the coordinates to the location of dryad
 inset_map <- inset_map + guides(fill=guide_legend(title="module\nnumber"))
 
print(inset_map)

top_10_scatterpies <- ggplot()+ geom_scatterpie(aes(x=Lon, y=lat, group=layer_id, r=0.2), alpha=0.6, data= pie_chart_data, 
                                                cols=colnames(pie_chart_data[,c(4:8)]))+
  theme_void()+ coord_fixed()+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))

#10% smallest modules
bottom_15_percent <- select(module_data_with_loc, c("module", "size_of_module")) %>% distinct %>% arrange(desc(size_of_module)) %>%
  filter(quantile(size_of_module, probs=0.15) > size_of_module) #arrange all modules to see <25 percentile and filter by them
bottom_15_percent$size_of_module <- bottom_15_percent$module #don't need the sizes anymore
names(bottom_15_percent)[names(bottom_15_percent)=="size_of_module"] <- "group" #change name to group- only the top percentile
pie_chart_data_bottom <- bottom_15_percent %>% right_join(module_data_with_loc, by=c("module"="module")) #join to have coordinates as well
pie_chart_data_bottom$group <- as.character(pie_chart_data_bottom$group) 
pie_chart_data_bottom <- pie_chart_data_bottom %>% drop_na()

pie_chart_data_bottom <- lon_lat_data %>% right_join(pivot_by_country(pie_chart_data_bottom), by = c("layer_id"="layer_id"))
pie_chart_data_bottom["layer_id"][pie_chart_data_bottom["layer_id"] == 7] <- 666
pie_chart_data_bottom <- pie_chart_data_bottom %>% group_by(lat, Lon) %>% summarise(across(everything(), sum))

worldmap <- map_data("world")
dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990))  
inset_map_bottom <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+ #create a terrain map based on dryad
  geom_scatterpie(aes(x=Lon, y=lat, group=layer_id, r=0.028), alpha=0.6, data= pie_chart_data_bottom, cols=colnames(pie_chart_data_bottom[,c(4:10)]))+ #create a pie chart for every location with the 90 precentile
  coord_fixed(xlim=c(-18.542076, -12.58351), ylim=c(26, 30.323990)) #set the coordinates to the location of dryad
inset_map_bottom <- inset_map_bottom + guides(fill=guide_legend(title="module\nnumber"))

print(inset_map_bottom)

ggplot()+ geom_scatterpie(aes(x=Lon, y=lat, group=layer_id, r=0.2), alpha=0.6, data= pie_chart_data_bottom, 
                          cols=colnames(pie_chart_data_bottom[,c(4:10)]))+ theme_classic()+ coord_fixed() #create a pie chart for every location with the 90 precentile


#---- species ids and number of islands they're found in----------------------------------
plants_in_network <- tot_plant %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>% #output is layer, node id and type of species
  inner_join(physical_nodes, by = c("node_from" = "species")) %>% subset(select = -c(layer_from, tot, node_from))

old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
plants_in_network$layer_id[plants_in_network$layer_id %in% old] <- new[match(plants_in_network$layer_id, old)] #change form layer to island

pols_in_network <- tot_pol %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>% #output is layer, node id and type of species
  inner_join(physical_nodes, by = c("node_from" = "species")) %>% subset(select = -c(layer_from, tot, node_from))

pols_in_network$layer_id[pols_in_network$layer_id %in% old] <- new[match(pols_in_network$layer_id, old)] #change form layer to island

species_in_network <- rbind(plants_in_network, pols_in_network) %>% unique()
num_of_islands_per_species <- species_in_network %>% group_by(node_id) %>% mutate(number_of_islands = n())

num_of_islands_per_species %>% ggplot(aes(x=reorder(node_id, -number_of_islands), y=number_of_islands, fill= type))+
  geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ theme_classic()+
  labs(x= "Node Id", y= "Number of Islands")+ scale_x_discrete(labels = NULL, breaks = NULL)+ 
  scale_y_continuous(breaks=seq(1,7,1))+ theme(axis.title=element_text(size=15))+
  theme(axis.text.y=element_text(size=12))+ theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12))
  
num_of_islands_per_species %>% ggplot(aes(x = number_of_islands, fill = type, color = type))+ 
  geom_density(alpha = 0.4)+ theme_classic()+ scale_x_continuous(breaks=seq(1,7,1))+ 
  theme(axis.title=element_text(size=15))+ theme(axis.text.y=element_text(size=12))+ 
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12))+
  labs(x= "Number of Islands")

#average number of islands a plant or a pollinator is found in
ave_number_of_islands_pols <- num_of_islands_per_species %>% filter(type == "pollinator") 
mean(ave_number_of_islands_pols$number_of_islands)

ave_number_of_islands_plants <- num_of_islands_per_species %>% filter(type == "plant") 
mean(ave_number_of_islands_plants$number_of_islands)

