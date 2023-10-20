#This code turns the data into a multilayer network with 7 layers 
#(islands as layers) and test distance decay in modules composition
#of the empirical network.


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
library(extRC)

##----get_data--------------------------------------------------------------------------------------------------------
setwd("D:/Trabajo/Papers/Canary_Island/spatial_modularity_in_the_canary_islands")
source("D:/Trabajo/Papers/Canary_Island/spatial_modularity_in_the_canary_islands/R/functions.R")

## ---- Create edgelist of the empirical network using islands (locations) as layers

#dataframe of intralayer links for each site 
interactions_csv <- read.csv("./csvs_nuevo/interactions.csv")

dryad_intralayer <- NULL

dryad_intralayer <- interactions_csv %>% group_by(interaction_id) %>% #group 1 interaction at a time
  mutate(layer_from = value[2], node_from = value[1], layer_to = value[4], 
         node_to = value[3], weight = value[5]) %>% #create new columns based on data
  subset(select = -c(interaction_id, attribute, value)) %>% unique() #delete previous columns

#write.csv(dryad_intralayer, "./csvs_nuevo/intralayer_file.csv", row.names = FALSE)


##---- Aggregate from site to island level as layers (locations) -------------------------------------------------------------------------------
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

#Create intraedgelist
edgelist_intralayers_both <- bind_rows(intralayer_weighted, intralayer_weighted_inverted) #combine inverted and non inverted versions of intra


##---- create interedges according to jaccard similarity between partners --------------------------------------------------------
# keep only species that occur in 2 or more layers
co_occurrence<- edgelist_intralayers_both%>% 
  group_by(node_from) %>%
  mutate(num_layers_from=n_distinct(layer_from)) %>% 
  filter(num_layers_from>="2")

# a for loop that calculates all the interlayer edges based on jaccard index
interlayers_with_weights_islands <- NULL

for (i in unique(co_occurrence$node_from)) {
  print(i)
  partners_sp <- 
    co_occurrence %>%
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
      mutate(node_from=i, node_to =i) %>%
      select(c(node_from,layer_from=Var1, layer_to=Var2,node_to, weight=value))
    interlayers_with_weights_islands <- rbind(interlayers_with_weights_islands,inter_fid)
}

interlayers_with_weights_islands<-interlayers_with_weights_islands[,c(2,1,3,4,5)]#change order columns

#inverted version
interlayer_inverted <- tibble(values= interlayers_with_weights_islands$layer_to, interlayers_with_weights_islands$node_to, interlayers_with_weights_islands$layer_from, 
                              interlayers_with_weights_islands$node_from, interlayers_with_weights_islands$weight) #create an inverted copy for directed intralayers
colnames(interlayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight")

#Create interedgelist
edgelist_interlayers_both <- bind_rows(interlayers_with_weights_islands, interlayer_inverted) #combine inverted and non inverted versions of intra

## ----multilayer_extended_final--------------------------------------------------------------------------------------
dryad_edgelist_complete <- bind_rows(edgelist_intralayers_both, edgelist_interlayers_both) #combine weighted version of intra and interlayer links

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

#write.csv(physical_nodes, "./csvs_nuevo/physical_nodes_islands.csv", row.names = FALSE)
#write.csv(layer_metadata, "./csvs_nuevo/layer_metadata_islands.csv", row.names = FALSE)


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

#write.csv(dryad_edgelist_complete_ids, "./csvs_nuevo/dryad_edgelist_complete_ids_islands.csv", row.names = FALSE)


## ---- Multilayer_class to calculate modularity-----------------------------------------------------------------------------------------------
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


#write.csv(intra_nonextended, "./csvs_nuevo/dryad_only_intralayer_edges_islands_as_layers.csv")
#write.csv(inter_extended, "./csvs_nuevo/dryad_only_interlayer_edges_islands_as_layers.csv")

#calculate modules for empirical network
modules_dryad_multilayer <- modified_multi(dryad_multilayer, 
                                           infomap_executable = "Infomap",
                                           flow_model = 'directed',
                                           relax = F, 
                                           silent = T, 
                                           trials = 100,
                                           seed = 497294, 
                                           temporal_network = F)

modules<-modules_dryad_multilayer$modules
#write.csv(modules, "./csvs_nuevo/modules_in_network_islands_as_layers.csv")


# ---- DISTANCE DECAY IN MODULES - EMPIRICAL DATA ------------------------------------------------------------------------

## Calculate distances between islands (locations)-----
distances <- read.csv("./csvs/distances_file.csv")

distances$layer_from[distances$layer_from %in% old_names] <- 
  new_names[match(distances$layer_from, old_names)] #change to reflect layer = island

distances$layer_to[distances$layer_to %in% old_names] <- 
  new_names[match(distances$layer_to, old_names)] #change to reflect layer = island

distances_normalized <- distances %>% filter(layer_from != layer_to) %>% #delete distances between sites in the same island
  group_by(layer_to, layer_from) %>% #group will contain 4 sites- site 1 and 1 of layer from and site 1 and 2 or layer to
  summarise(mean_distance = mean(distance_in_meters)) %>%unique() #use an average distance of the 4 sites in 2 different islands to determine the distance between the islands

distances_normalized <- distances_normalized[c("layer_from", "layer_to", "mean_distance")]

#add id of locations
distances_with_ids <- distances_normalized %>% left_join(layer_metadata, by= c("layer_from"="layer_name")) %>% 
  left_join(layer_metadata, by= c("layer_to"="layer_name")) %>% #add correct id to layer name
  select(mean_distance, layer_id.x, layer_id.y) #discard actual names of layers
names(distances_with_ids)[3] <- "layer_from" 
names(distances_with_ids)[4] <- "layer_to"
names(distances_with_ids)[1] <- "layer_name_to"

distances_with_ids <- distances_with_ids[c("layer_from", "layer_to", "mean_distance")]

#write.csv(distances_with_ids, "./csvs_nuevo/distances_with_ids_islands_as_layers.csv", row.names = FALSE)

## Jaccard on islands (locations) -----
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

#write.csv(module_pivoted, "./csvs_nuevo/module_pivoted_for_state_node_similarity_islands_as_layers.csv")

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

#write.csv(module_data, "./csvs_nuevo/module_data_islands_as_layers.csv", row.names = FALSE)
#write.csv(size, "./csvs_nuevo/size_islands_as_layers.csv", row.names = FALSE)

lon_lat_data <- read_csv('./csvs_nuevo/layers.csv') #create new data frame with just the layer data
lon_lat_data <- lon_lat_data %>% select(c("layer_id","lat","Lon")) %>% na.omit()  #only select layer id and coordinates

#layers as islands
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new_values <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)
lon_lat_data$layer_id[lon_lat_data$layer_id %in% old] <- new_values[match(lon_lat_data$layer_id, old)]

lon_lat_data <- lon_lat_data %>% unique() #delete duplicates caused by islands having 2 sites

#write.csv(lon_lat_data, "./csvs_nuevo/lon_lat_data_islands_as_layers.csv", row.names = FALSE)
lon_lat_data <- read.csv("./csvs_nuevo/lon_lat_data_islands_as_layers.csv")

module_data_with_loc <- merge(module_data, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon <- module_data_with_loc %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon$count <- c(1)

#write.csv(modules_with_lat_lon, "./csvs_nuevo/modules_with_lat_lon_islands_as_layers.csv", row.names = FALSE)
#modules_with_lat_lon <- read.csv("./csvs_nuevo/modules_with_lat_lon_islands_as_layers.csv")

#version with # of modules in layers
edge_list_by_islands_modules <- edge_list_with_distances
edge_list_by_islands_modules$count <- c(1)
edge_list_by_islands_modules <- edge_list_by_islands_modules %>% group_by(layer_from, layer_to) %>% 
  mutate(number_of_modules= sum(count)) %>% #count how many modules are shared by 2 islands
  select(layer_from, layer_to, module, number_of_modules, mean_distance) 

## Distance decay in modules -----
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

#write.csv(islands_turnover_with_distnace_empirical,  "./csvs_nuevo/islands_turnover_with_distnace_empirical.csv",  row.names = FALSE)

## Statistical analysis -----
shapiro.test(islands_turnover_with_distnace_empirical$turnover)#normal
lm1_module = lm(turnover ~ distance_in_km ,data=islands_turnover_with_distnace_empirical) #in empirical
summary(lm1_module)



## Graph distance decay in modules of the empirical network--
#pdf('./graphs/Islands/Jac/Modules distance decay alone_islands.pdf', 10, 6)
#islands_turnover_with_distnace_empirical %>%
  #ggplot(aes(x=distance_in_km, y=turnover))+
  #geom_point(color = "#FB3B1E")+ 
  #scale_x_continuous()+ stat_smooth(method= "lm", se=F, color = "#FB3B1E")+
  #theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  #theme(axis.text.y=element_text(size=15))+
  #labs(x="Distance in Km", y="Jaccard Similarity")+ theme_bw()+
  #theme(panel.grid = element_blank(),
  #     panel.border = element_rect(color = "black",fill = NA,size = 1),
  #      panel.spacing = unit(0.5, "cm", data = NULL),
  #      axis.text = element_text(size=14, color='black'),
  #      axis.title = element_text(size=14, color='black'),
  #      axis.line = element_blank())
#dev.off()


