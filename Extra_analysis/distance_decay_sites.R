#This code tests for distance decay at site level

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
library(ecodist)

##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")

## ---- Create edgelist of the empirical network using islands (locations) as layers

#dataframe of intralayer links for each site 
interactions_csv <- read.csv("./csvs_nuevo/interactions.csv")

dryad_intralayer <- NULL

dryad_intralayer <- interactions_csv %>% group_by(interaction_id) %>% #group 1 interaction at a time
  mutate(layer_from = value[2], node_from = value[1], layer_to = value[4], 
         node_to = value[3], weight = value[5]) %>% #create new columns based on data
  subset(select = -c(interaction_id, attribute, value)) %>% unique() #delete previous columns

#write.csv(dryad_intralayer, "./csvs_nuevo/intralayer_file.csv", row.names = FALSE)


##---- Remove dataframe of mainland -------------------------------------------------------------------------------

# Remove dataframe of mainland
dryad_intralayer_sites <- dryad_intralayer %>% filter (!(layer_from == "WesternSahara1"|
                                                          layer_from == "WesternSahara2"))

dryad_intralayer_sites$weight<-as.numeric(dryad_intralayer_sites$weight)
## ----dryad intralayer interlayer both ways-------------------------------------------------------------------------------------
intralayer_inverted <- tibble(values= dryad_intralayer_sites$layer_to, dryad_intralayer_sites$node_to, dryad_intralayer_sites$layer_from, 
                              dryad_intralayer_sites$node_from, dryad_intralayer_sites$weight) #create an inverted copy for directed intralayers
colnames(intralayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight")

## ---- normalize intralayer weights--------------------------------------------------------------------------
#plants in from
tot_plant <- dryad_intralayer_sites %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted <- dryad_intralayer_sites %>% left_join(tot_plant) %>% mutate(rel_weight=weight/tot) %>% 
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
interlayers_with_weights_sites <- NULL

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
    interlayers_with_weights_sites <- rbind(interlayers_with_weights_sites,inter_fid)
}

interlayers_with_weights_sites<-interlayers_with_weights_sites[,c(2,1,3,4,5)]#change order columns


#inverted version
interlayer_inverted <- tibble(values= interlayers_with_weights_sites$layer_to, interlayers_with_weights_sites$node_to, interlayers_with_weights_sites$layer_from, 
                              interlayers_with_weights_sites$node_from, interlayers_with_weights_sites$weight) #create an inverted copy for directed intralayers
colnames(interlayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight")

#Create interedgelist
edgelist_interlayers_both <- bind_rows(interlayers_with_weights_sites, interlayer_inverted) #combine inverted and non inverted versions of intra

## ----multilayer_extended_final--------------------------------------------------------------------------------------
dryad_edgelist_complete <- bind_rows(edgelist_intralayers_both, edgelist_interlayers_both) #combine weighted version of intra and interlayer links

## ----node_metadata--------------------------------------------------------------------------------------------------                                            
pollinators <- sort(unique(intralayer_weighted$node_to)) #adding up only pol who haven't been added yet 
plants <- sort(unique(intralayer_weighted$node_from)) #adding up only plants who haven't been added yet
intersect(pollinators, plants) #making sure I don't have plants in pol or other way around
A <- length(pollinators) # Number of pollinators
P <- length(plants) # Number of plants
S <- A+P

sites_names <- unique(dryad_edgelist_complete$layer_from)

# Create a table with node metadata
physical_nodes <- tibble(node_id=1:S, #1 till the last species
                         type=c(rep('plant',P),rep('pollinator',A)), #replicate the words P and A times
                         species=c(plants,pollinators)) #add species from plants and pollinators in accordance
layer_metadata <- tibble(layer_id=c(1:12), layer_name= sites_names)  #give num to each layer

#write.csv(physical_nodes, "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/physical_nodes_sites.csv", row.names = FALSE)
#write.csv(layer_metadata, "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/layer_metadata_sites.csv", row.names = FALSE)


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

#write.csv(dryad_edgelist_complete_ids, "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/dryad_edgelist_complete_ids_sites.csv", 
                      #row.names = FALSE)


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


#write.csv(intra_nonextended, "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/dryad_only_intralayer_edges_sites_as_layers.csv")
#write.csv(inter_extended, "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/dryad_only_interlayer_edges_sites_as_layers.csv")

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
#write.csv(modules, "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/modules_in_network_sites_as_layers.csv")


# ---- DISTANCE DECAY IN MODULES - EMPIRICAL DATA AND SITE LEVEL ------------------------------------------------------------------------

## Calculate distances between islands (locations)-----
distances <- read.csv("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/distances_sites_file.csv", sep = ",") %>% 
  filter(!(layer_from == "WesternSahara1" |layer_to == "WesternSahara1"|
             layer_from == "WesternSahara2" |layer_to == "WesternSahara2")) %>%  # remove mainland
    rename("mean_distance" = "distance_in_meters")


#add id of locations
distances_with_ids <- distances%>% left_join(layer_metadata, by= c("layer_from"="layer_name")) %>% 
  left_join(layer_metadata, by= c("layer_to"="layer_name")) %>% #add correct id to layer name
  select(mean_distance, layer_id.x, layer_id.y) #discard actual names of layers
names(distances_with_ids)[2] <- "layer_from" 
names(distances_with_ids)[3] <- "layer_to"
names(distances_with_ids)[1] <- "mean_distance"

distances_with_ids <- distances_with_ids[c("layer_from", "layer_to", "mean_distance")]

distances_with_ids_inverted<- tibble(layer_from=distances_with_ids$layer_to, 
                                     layer_to = distances_with_ids$layer_from, 
                                     mean_distance = distances_with_ids$mean_distance)# to add all combination of layers from and to

distances_with_ids_final<-rbind(distances_with_ids,distances_with_ids_inverted) %>% unique()

#write.csv(distances_with_ids_final, "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/distances_with_ids_sites_as_layers.csv", row.names = FALSE)


## Jaccard on sites (locations) -----
#similarity check 2 furthest apart
modules_for_similarity_num <- modules %>% select(module, layer_id) %>% 
  unique() %>% group_by(module) %>% select(module) %>% unique()
modules_for_similarity <- modules %>%
  filter(module %in% modules_for_similarity_num$module) #only save the modules that are found in 2 or more layers


#pivot modules function for sites as layers
pivot_by_module_sites <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(modules_for_similarity, id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s3 = t(s2) 
  s3 <- s3[-1,]
  colnames(s3) <- c(1:12)
  return(s3)
}

module_pivoted <- pivot_by_module_sites(modules_for_similarity)

#write.csv(module_pivoted, "./csvs_nuevo/module_pivoted_for_state_node_similarity_justislands_as_layers.csv")

#edge list per module function for islands as layers
edge_list_per_module_sites <- function(data,edge_list){
  #gets one row from a data frame and creates an edge list from it
  for (i in (1:11)){
    if (data[i]==0) next #only take layers where the module is present
    else {
      for (j in ((i+1):12)){
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
  modules_edge_list <- edge_list_per_module_sites(module_pivoted[i,], modules_edge_list)
  current_module <- rownames(module_pivoted)[i]
  modules_edge_list <- modules_edge_list %>% mutate(module = replace_na(module, current_module)) #add module number
}
 

 
edge_list_with_distances <- left_join(modules_edge_list, distances_with_ids_final, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances <- na.omit(edge_list_with_distances) #remove NA and delete layer name


#arrange data to include coordinates and modules sizes
size <- count(modules_dryad_multilayer$modules, module)  #create a data frame of all modules and how many nodes are in each (size of module)
module_data <- merge(modules_dryad_multilayer$modules , size, by=c("module","module")) #merge size of module with all the other info about the modules
colnames(module_data)[7] <- "size_of_module" #rename column


lon_lat_data <- read_csv('./csvs_nuevo/layers.csv') %>% 
  filter(!(layer_name =="WesternSahara1"| layer_name =="WesternSahara2" )) #create new data frame with just the layer data

#cambiar ID segun el nuevo dataframe
lon_lat_data <- lon_lat_data %>% select(c("layer_id","lat","Lon")) %>% na.omit() %>%  #only select layer id and coordinates
                select(-layer_id) 
lon_lat_data<-tibble(lon_lat_data, layer_id = c(1:12))
module_data_with_loc <- merge(module_data, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon <- module_data_with_loc %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon$count <- c(1)

#version with # of modules in layers
edge_list_by_sites_modules <- edge_list_with_distances
edge_list_by_sites_modules$count <- c(1)
edge_list_by_sites_modules <- edge_list_by_sites_modules %>% group_by(layer_from, layer_to) %>% 
  mutate(number_of_modules= sum(count)) %>% #count how many modules are shared by 2 islands
  select(layer_from, layer_to, module, number_of_modules, mean_distance) 

## Distance decay in modules -----
module_sites_turnover <- NULL

for (i in (1:11)){
  for (j in ((1+i):12)){
    print(i)
    modules_in_site_from <- filter(modules_with_lat_lon, layer_id == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_site_to <- filter(modules_with_lat_lon, layer_id == j) %>% select(module) %>% unique() %>% unlist()
    #take all modules in one layer and all modules in another layer to check turnover
    int_both <- intersect(modules_in_site_from, modules_in_site_to) #how many modules are found in both layers
    uni_both <- union(modules_in_site_from, modules_in_site_to)
    
    # Calculate turnover
    turnover <- ifelse(length(uni_both) == 0, 0, length(int_both) / length(uni_both))
    
    module_sites_turnover <- rbind(module_sites_turnover, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

sites_turnover_with_distnace_empirical <- edge_list_by_sites_modules %>%
  merge(module_sites_turnover, by= c("layer_from", "layer_to")) %>% 
  select(layer_from, layer_to, number_of_modules, mean_distance, turnover) %>% unique() #merge both versions 

sites_turnover_with_distnace_empirical <- sites_turnover_with_distnace_empirical %>% 
  mutate(distance_in_km=mean_distance/1000) #turn to km

#write.csv(sites_turnover_with_distnace_empirical,  "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/sites_turnover_with_distnace_empirical.csv",  row.names = FALSE)

## Statistical analysis -----
distances_with_ids_final <- read.csv("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/distances_with_ids_sites_as_layers.csv", sep = ",") 
  
emp <- read.csv( "/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/Extra_analysis/sites_turnover_with_distnace_empirical.csv", sep =",") %>% 
        select(-number_of_modules,-mean_distance)

#add the layers that don't share any module
to_add<- distances_with_ids_final %>% filter((layer_from ==1 &layer_to ==11) |
                                               (layer_from ==1 &layer_to ==12)|(layer_from ==2 &layer_to ==11) |
                                               (layer_from ==2 &layer_to ==12)) 
to_add_final<- tibble(to_add,turnover = 0) %>% mutate(distance_in_km = mean_distance /1000) %>% select(-mean_distance)

#prepare final dataframe
emp_final<-rbind(emp,to_add_final)

shapiro.test(emp_final$turnover)#normal

set.seed(122)

m_emp_sites<-MRM(turnover ~ distance_in_km, data=emp_final, nperm=999 )
m_emp_sites

#pdf('./graphs/Modules_DD_sites.pdf', 10, 6)
emp_final %>% 
  ggplot(aes(x= distance_in_km, y= turnover, color = "fixed_color"))+
  geom_point(color = "#FB3B1E")+
  scale_color_manual(values = c("fixed_color" = "#FB3B1E"))+
  labs(x="Distance (Km)", y="Jaccard Similarity")+  
  theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=15, color='black'),
        axis.title = element_text(size=17, color='black'),
        axis.line = element_blank(),
        legend.text.align = 0,
        legend.title =  element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11))+
  guides(color = FALSE)

#dev.off()



