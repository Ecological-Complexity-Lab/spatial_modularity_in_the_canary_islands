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

##shows change in species composition in modules between layers while trying to remove the factor of turnover
##intersection of species in layer A module 1 and layer B module 1
##union of species in layer A module 1 and layer B
##removing factor of turnover is done by conditions requiring:
## 1. modules found in at least 3 layers
## 2. at least 4 nodes in common netween layer A module 1 and layer B

## ---- module sharing empirical-------------------------------------------------------------------------------------------------
#similarity check every two layers over 5
modules_for_similarity_over_5_num <- modules_dryad_multilayer$modules %>% select(module, layer_id) %>%
  unique() %>% group_by(module) %>% filter(n()>4) %>% select(module) %>% unique() #check which modules occur in 5 or more layers
modules_for_similarity_over_5 <- modules_dryad_multilayer$modules %>%
  filter(module %in% modules_for_similarity_over_5_num$module) #only save the modules that are found in 5 or more layers
view(modules_for_similarity_over_5)
#pivot module
module_pivoted_over_5 <- pivot_by_module(modules_for_similarity_over_5)
print(module_pivoted_over_5)

#write.csv(module_pivoted_over_5, "./csvs/module_pivoted_over_3_for_state_node_similarity.csv")

#module turnover across layers for modules found in >=5 layers
modules_edge_list_over_5 <- NULL

for (i in (1:nrow(module_pivoted_over_5))){ #run the function for each row in the data frame
  modules_edge_list_over_5 <- edge_list_per_module(module_pivoted_over_5[i,], modules_edge_list_over_5)
  current_module_over_5 <- rownames(module_pivoted_over_5)[i]
  modules_edge_list_over_5 <- modules_edge_list_over_5 %>% mutate(module = replace_na(module, current_module_over_5)) #add module number
}

edge_list_with_distances_over_5 <- right_join(modules_edge_list_over_5, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances_over_5 <- na.omit(edge_list_with_distances_over_5) #remove NA


#----intersection over intersection-----------------------------------------
module_sharing_data <- NULL

for (i in 1:30){
  focal_module <- filter(edge_list_with_distances_over_5, module == i) #look at one module at a time
  print(i)
  for (j in (1:nrow(focal_module))){
    module <- focal_module[j,]
    current_module <- module$module
    current_layer <- module$layer_from
    current_layer_to <- module$layer_to
    current_distance <- module$distance_in_meters
    physical_nodes_in_layer_from <- filter(modules_for_similarity_over_5, (module %in% current_module & layer_id %in% current_layer)) %>% 
      select(node_id)  #take the nodes that are found in the module in layer_from
    if (nrow(physical_nodes_in_layer_from)<4) next
    else {
      physical_nodes_in_layer_from <- physical_nodes_in_layer_from %>% unlist()
      physical_nodes_in_layer_to <- filter(modules_for_similarity_over_5, (layer_id %in% current_layer_to)) %>% select(node_id) %>% unlist()
      if (length(intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to))<4) next #only continue if more than 4 nodes ar in common
      else{ 
        S_OL <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in layer A in the module 
        #and also in layer B (but not necessarily in the module)
        layer_to_data <- filter(modules_for_similarity_over_5, (layer_id %in% current_layer_to)) 
        layer_to_and_also_module <- filter(layer_to_data, module %in% current_module) %>% select(node_id) %>% unlist()
        S_all <- (intersect(layer_to_and_also_module, S_OL)) #how many nodes are found in layer A in the module 
        #and also in layer B in the same module
        module_sharing <- length(S_all)/length(S_OL)
        module_sharing_data <- rbind(module_sharing_data, tibble(module= current_module, layer_from= current_layer,
                                                                 layer_to= current_layer_to, module_sharing= module_sharing,
                                                                 distance= current_distance))
      }
    }
  }
}

#view(module_sharing_data)
#write.csv(module_sharing_data, "./csvs/module_sharing_between_layers_over_3.csv")

#----intersection over union---------------------------------------------------------------
module_sharing_data_union <- NULL
union_two_layers <- NULL

for (i in 1:30){
  focal_module <- filter(edge_list_with_distances_over_5, module == i) #look at one module at a time
  for (j in (1:nrow(focal_module))){
    module <- focal_module[j,]
    current_module <- module$module
    current_layer <- module$layer_from
    current_layer_to <- module$layer_to
    current_distance <- module$distance_in_meters
    physical_nodes_in_layer_from <- filter(modules_for_similarity_over_5, (module %in% current_module & layer_id %in% current_layer)) %>% 
      select(node_id)  #take the nodes that are found in the module in layer_from
    if (nrow(physical_nodes_in_layer_from)<4) next #too strict
    else {
      physical_nodes_in_layer_from <- physical_nodes_in_layer_from %>% unlist()
      physical_nodes_in_layer_to <- filter(modules_for_similarity_over_5, (layer_id %in% current_layer_to)) %>% select(node_id) %>% unlist()
      physical_nodes_in_layer_to_and_module <- filter(modules_for_similarity_over_5, (module %in% current_module & layer_id %in% current_layer_to)) %>% 
        select(node_id) %>% unlist()
      if (length(intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to))<4) next #only continue if more than 5 nodes ar in common
      else{
        S_OL <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in layer A in the module 
        #and also in layer B (but not necessarily in the module)
        layer_to_data <- filter(modules_for_similarity_over_5, (layer_id %in% current_layer_to)) 
        layer_to_and_also_module <- filter(layer_to_data, module %in% current_module) %>% select(node_id) %>% unlist()
        S_union <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to_and_module)
        S_all <- (intersect(layer_to_and_also_module, S_OL)) #how many nodes are found in layer A in the module 
        #and also in layer B in the same module
        module_sharing <- length(S_all)/length(S_union)
        module_sharing_data_union <- rbind(module_sharing_data_union, tibble(module= current_module, layer_from= current_layer,
                                                                             layer_to= current_layer_to, module_sharing= module_sharing,
                                                                             distance= current_distance))
        union_two_layers <- rbind(union_two_layers, tibble(module= current_module,layer_from= current_layer, layer_to= current_layer_to, 
                                                           distance= current_distance, num_of_species= length(S_union)))
      }
    }
  }
}

#view(module_sharing_data_union)
#write.csv(module_sharing_data_union, "./csvs/module_sharing_between_layers_int_over_union_over_3.csv")
#write.csv(union_two_layers, "./csvs/speies_found_in_layers_with_module_sharing_over_3.csv")

#---- graphs for module sharing---------------------------------------------------------
#graph for module sharing between every two layers
ggplot(module_sharing_data, aes(x=distance, y=module_sharing, color=module))+
  geom_point()+ scale_x_continuous(breaks=seq(0,455736.67290,100000))+theme_classic()+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="module sharing")

#graph for module sharing between every two layers no distinction
ggplot(module_sharing_data, aes(x=distance, y=module_sharing))+
  geom_point()+ scale_x_continuous(breaks=seq(0,455736.67290,100000))+ theme_classic()+ geom_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="module sharing")

#number of species union two layers
ggplot(union_two_layers, aes(x=distance, y=num_of_species, color=module))+
  geom_point()+ scale_x_continuous()+ theme_classic()+ 
  theme(axis.title=element_text(size=15))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=12))+
  labs(x="distance in meters", y="number of species in the module in both layers")

ggplot(union_two_layers, aes(x=distance, y=num_of_species, color=module))+
  geom_point()+ scale_x_continuous()+ theme_classic()+ geom_smooth(se = FALSE, method = "lm")+
  theme(axis.title=element_text(size=15))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=12))+
  labs(x="distance in meters", y="number of species in the module in both layers")

#intersection over union
ggplot(module_sharing_data_union, aes(x=distance, y=module_sharing, color=module))+
  geom_point()+ scale_x_continuous()+theme_classic()+ geom_smooth()+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="module sharing")

ggplot(module_sharing_data_union, aes(x=distance, y=module_sharing, color=module))+
  geom_point()+ scale_x_continuous()+theme_classic()+ geom_smooth(se = FALSE, method = "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="module sharing")

#version with no distinction between modules
ggplot(module_sharing_data_union, aes(x=distance, y=module_sharing))+
  geom_point()+ scale_x_continuous()+ theme_classic()+ geom_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="module sharing")


#---- module sharing shuffled null models-----------------------------------------------------

module_sharing_all_shuf <- NULL
union_two_layers_all_shuf <- NULL

for (shuf in 1:max(shuffled_500$i)){ 
  current_shuf <- shuffled_500 %>% filter(i == shuf)
  print(shuf)
  #similarity check every two layers over 5
  modules_for_similarity_over_3_num_shuf <- current_shuf %>% select(module, layer_id) %>%
    unique() %>% group_by(module) %>% filter(n()>2) %>% select(module) %>% unique() #check which modules occur in 3 or more layers
  modules_for_similarity_over_3_shuf <- current_shuf %>%
    filter(module %in% modules_for_similarity_over_3_num_shuf$module) #only save the modules that are found in 3 or more layers
  
  #pivot module
  module_pivoted_over_3_shuf <- pivot_by_module(modules_for_similarity_over_3_shuf)
  module_pivoted_over_3_shuf <- module_pivoted_over_3_shuf[, -1]
  
  #view(module_pivoted_over_3_shuf)
  
  #module turnover across layers for modules found in >=5 layers
  modules_edge_list_over_3_shuf <- NULL
  
  for (i in (1:nrow(module_pivoted_over_3_shuf))){ #run the function for each row in the data frame
    modules_edge_list_over_3_shuf <- edge_list_per_module(module_pivoted_over_3_shuf[i,], modules_edge_list_over_3_shuf) 
  }
  
  edge_list_with_distances_over_3_shuf <- right_join(modules_edge_list_over_3_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
  edge_list_with_distances_over_3_shuf <- na.omit(edge_list_with_distances_over_3_shuf) #remove NA
  
  #module sharing null vs empirical
  module_sharing_data_union_shuf <- NULL
  union_two_layer_shuf <- NULL
  
  for (i in (1:nrow(module_pivoted_over_3_shuf))){
    focal_module <- filter(edge_list_with_distances_over_3_shuf, module == i) #look at one module at a time
    for (j in (1:nrow(focal_module))){
      module <- focal_module[j,]
      current_module <- module$module
      current_layer <- module$layer_from
      current_layer_to <- module$layer_to
      current_distance <- module$distance_in_meters
      physical_nodes_in_layer_from <- filter(modules_for_similarity_over_3_shuf, (module %in% current_module & layer_id %in% current_layer)) %>% 
        select(node_id)  #take the nodes that are found in the module in layer_from
      if (nrow(physical_nodes_in_layer_from)<4) next
      else {
        physical_nodes_in_layer_from <- physical_nodes_in_layer_from %>% unlist()
        physical_nodes_in_layer_to <- filter(modules_for_similarity_over_3_shuf, (layer_id %in% current_layer_to)) %>% select(node_id) %>% unlist()
        physical_nodes_in_layer_to_and_module <- filter(modules_for_similarity_over_3_shuf, (module %in% current_module & layer_id %in% current_layer_to)) %>% 
          select(node_id) %>% unlist()
        if (length(intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to))<4) next #only continue if more than 4 nodes ar in common
        else{
          S_OL <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in layer A in the module 
          #and also in layer B (but not necessarily in the module)
          layer_to_data <- filter(modules_for_similarity_over_3_shuf, (layer_id %in% current_layer_to)) 
          layer_to_and_also_module <- filter(layer_to_data, module %in% current_module) %>% select(node_id) %>% unlist()
          S_union <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to_and_module)
          S_all <- (intersect(layer_to_and_also_module, S_OL)) #how many nodes are found in layer A in the module 
          #and also in layer B in the same module
          module_sharing <- length(S_all)/length(S_union)
          module_sharing_data_union_shuf <- rbind(module_sharing_data_union_shuf, tibble(module= current_module, layer_from= current_layer,
                                                                                         layer_to= current_layer_to, module_sharing= module_sharing,
                                                                                         distance= current_distance))
          union_two_layer_shuf <- rbind(union_two_layer_shuf, tibble(module= current_module,layer_from= current_layer, layer_to= current_layer_to, 
                                                                     distance= current_distance, num_of_species= length(S_union)))
        }
      }
    }
  }
  module_sharing_all_shuf <- rbind(module_sharing_all_shuf, tibble(module_sharing_data_union_shuf, shuf))
  union_two_layers_all_shuf <- rbind(union_two_layers_all_shuf, tibble(union_two_layer_shuf, shuf))
}

ave_module_sharing_shuf <- module_sharing_all_shuf %>% group_by(module, layer_from, layer_to, distance) %>%
  summarise(ave=mean(module_sharing), sd=sd(module_sharing)) %>% mutate(type="null") #create mean and sd for each point

ave_module_sharing_empirical <- module_sharing_data_union %>% group_by(module, layer_from, layer_to, distance) %>%
  summarise(ave=mean(module_sharing), sd=sd(module_sharing)) %>% mutate(type="empirical") #create mean and sd for each point

module_sharing_null_and_empirical <- rbind(ave_module_sharing_empirical, ave_module_sharing_shuf)

#---- graph for module sharing between every two layers------------------------------------------------------
ggplot(module_sharing_null_and_empirical, aes(x=distance, y=ave, color=type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ 
  scale_x_continuous(breaks=seq(0,455736.67290,100000))+ theme_classic()+ geom_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="module sharing")



##---- hierarchical modularity module sharing---------------------------------------------------------------------------


#similarity check every two layers over 5
sub_modules_for_similarity_over_5_num <- modules_dryad_multilayer_multi_level_sub_modules %>% select(module_sub_module, layer_id) %>%
  unique() %>% group_by(module_sub_module) %>% filter(n()>4) %>% select(module_sub_module) %>% unique() #check which sub modules occur in 5 or more layers
sub_modules_for_similarity_over_5 <- modules_dryad_multilayer_multi_level_sub_modules %>%
  filter(module_sub_module %in% sub_modules_for_similarity_over_5_num$module_sub_module) #only save the sub modules that are found in 5 or more layers

#pivot module
sub_module_pivoted_over_5 <- pivot_by_sub_module(sub_modules_for_similarity_over_5) 

#write.csv(sub_module_pivoted_over_5, "./csvs/sub_module_pivoted_over_5.csv")

#module turnover across layers for modules found in >=5 layers
sub_modules_edge_list_over_5 <- NULL

for (i in (1:nrow(sub_module_pivoted_over_5))){ #run the function for each row in the data frame
  sub_modules_edge_list_over_5 <- edge_list_per_sub_module(sub_module_pivoted_over_5[i,], sub_modules_edge_list_over_5)
  current_sub_module_over_5 <- rownames(sub_module_pivoted_over_5)[i]
  sub_modules_edge_list_over_5 <- sub_modules_edge_list_over_5 %>% 
    mutate(module_sub_module = replace_na(module_sub_module, current_sub_module_over_5)) #add module number
}

edge_list_with_distances_over_5_sub_modules <- right_join(sub_modules_edge_list_over_5, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances_over_5_sub_modules <- na.omit(edge_list_with_distances_over_5_sub_modules) #remove NA

view(edge_list_with_distances_over_5_sub_modules)

#----intersection over union---------------------------------------------------------
sub_module_sharing_data_union <- NULL
union_two_layers_sub_modules <- NULL


for (i in 1:nrow(edge_list_with_distances_over_5_sub_modules)){
  focal_sub_module_number <- edge_list_with_distances_over_5_sub_modules[i,]$module_sub_module #get sub module number
  focal_sub_module <- filter(edge_list_with_distances_over_5_sub_modules, 
                             module_sub_module == focal_sub_module_number) #get only rows where the module_sub_module is the wanter value
  for (j in (1:nrow(focal_sub_module))){
    module_sub_module <- focal_sub_module[j,]
    current_sub_module <- module_sub_module$module_sub_module
    current_layer_sub_module <- module_sub_module$layer_from
    current_layer_to_sub_module <- module_sub_module$layer_to
    current_distance_sub_module <- module_sub_module$distance_in_meters
    physical_nodes_in_layer_from_sub_module <- filter(sub_modules_for_similarity_over_5, (module_sub_module %in% current_sub_module 
                                                                                          & layer_id %in% current_layer_sub_module)) %>% 
      select(node_id)  #take the nodes that are found in the module in layer_from
    if (nrow(physical_nodes_in_layer_from_sub_module)<4) next #too strict
    else {
      physical_nodes_in_layer_from_sub_module <- physical_nodes_in_layer_from_sub_module %>% unlist()
      physical_nodes_in_layer_to_sub_module <- filter(sub_modules_for_similarity_over_5, (layer_id %in% current_layer_to_sub_module)) %>% 
        select(node_id) %>% unlist()
      physical_nodes_in_layer_to_and_sub_module <- filter(sub_modules_for_similarity_over_5, (module_sub_module %in% current_sub_module & 
                                                                                                layer_id %in% current_layer_to_sub_module)) %>% 
        select(node_id) %>% unlist()
      if (length(intersect(physical_nodes_in_layer_from_sub_module, physical_nodes_in_layer_to_sub_module))<4) next #only continue if more than 5 nodes ar in common
      else{
        S_OL_sub_module <- intersect(physical_nodes_in_layer_from_sub_module, physical_nodes_in_layer_to_sub_module) #how many nodes are found in layer A in the module 
        #and also in layer B (but not necessarily in the module)
        layer_to_data_sub_module <- filter(sub_modules_for_similarity_over_5, (layer_id %in% current_layer_to_sub_module)) 
        layer_to_and_also_sub_module <- filter(layer_to_data_sub_module, module_sub_module %in% current_sub_module) %>% 
          select(node_id) %>% unlist()
        S_union_sub_module <- union(physical_nodes_in_layer_from_sub_module, physical_nodes_in_layer_to_and_sub_module)
        S_all_sub_module <- (intersect(layer_to_and_also_sub_module, S_OL_sub_module)) #how many nodes are found in layer A in the module 
        #and also in layer B in the same module
        sub_module_sharing <- length(S_all_sub_module)/length(S_union_sub_module)
        sub_module_sharing_data_union <- rbind(sub_module_sharing_data_union, tibble(module_sub_module= current_sub_module, 
                                                                                     layer_from= current_layer_sub_module,
                                                                                     layer_to= current_layer_to_sub_module, 
                                                                                     sub_module_sharing= sub_module_sharing,
                                                                                     distance= current_distance_sub_module))
        union_two_layers_sub_modules <- rbind(union_two_layers_sub_modules, tibble(module_sub_module= current_sub_module,
                                                                                   layer_from= current_layer_sub_module, 
                                                                                   layer_to= current_layer_to_sub_module, 
                                                                                   distance= current_distance_sub_module, 
                                                                                   num_of_species= length(S_union_sub_module)))
      }
    }
  }
}

#write.csv(sub_module_sharing_data_union, "./csvs/sub_module_sharing_data_union.csv")
#write.csv(union_two_layers_sub_modules, "./csvs/union_two_layers_sub_modules.csv")

sub_module_sharing_data_union <- sub_module_sharing_data_union %>% unique() #make sure there are no repeats

view(sub_module_sharing_data_union)

#----graphs module sharing hierarchial modularity-----------------------------------------------

ggplot(sub_module_sharing_data_union, aes(x=distance, y=sub_module_sharing, color=module_sub_module))+
  geom_point()+ scale_x_continuous()+theme_classic()+ geom_smooth()+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Meters", y="Sub Module Sharing")

ggplot(sub_module_sharing_data_union, aes(x=distance, y=sub_module_sharing, color=module_sub_module))+
  geom_point()+ scale_x_continuous()+theme_classic()+ geom_smooth(se = FALSE, method = "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Meters", y="Sub Module Sharing")

#version with no distinction between modules
ggplot(sub_module_sharing_data_union, aes(x=distance, y=sub_module_sharing))+
  geom_point()+ scale_x_continuous()+ theme_classic()+ geom_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="Distance in Meters", y="Sub Module Sharing")

