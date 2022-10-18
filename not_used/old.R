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


#---- turnover distance decay with shuffled networks for sanity check-------------------------------------

shuffled_500 <- read.csv("./csvs/modules_interlayer_suffled.csv")

all_layers_turnover <- NULL

for (shuf in 1:max(shuffled_500$i)){ 
  current_shuf <- shuffled_500 %>% filter(i == shuf)
  print(shuf)
  modules_for_similarity_num_shuf <- current_shuf %>% select(module, layer_id) %>% 
    unique() %>% group_by(module) %>% filter(n()>1) %>% select(module) %>% unique() #check which modules occur in 2 or more layers
  modules_for_similarity_shuffled <- current_shuf %>%
    filter(module %in% modules_for_similarity_num_shuf$module) #only save the modules that are found in 2 or more layers
  
  module_pivoted_shuf <- pivot_by_module(modules_for_similarity_shuffled)
  print(module_pivoted_shuf)
  module_pivoted_shuf <- module_pivoted_shuf[, -1]
  
  modules_edge_list_shuf <- NULL
  for (i in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the current_shuf frame
    modules_edge_list_shuf <- edge_list_per_module(module_pivoted_shuf[i,], modules_edge_list_shuf) 
  }
  
  
  edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
  edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA
  
  layers_turnover_shuf <- NULL
  
  for (i in (1:nrow(module_pivoted_shuf))){
    focal_module <- filter(edge_list_with_distances_shuf, module == i) #look at one module at a time
    for (j in (1:nrow(focal_module))){
      module <- focal_module[j,]
      current_module <- module$module
      current_layer <- module$layer_from
      current_layer_to <- module$layer_to
      current_distance <- module$distance_in_meters
      physical_nodes_in_layer_from <- filter(modules_for_similarity_shuffled, layer_id == current_layer) %>% select(node_id) %>%
        unlist() #take the nodes that are found in the module in layer_from
      physical_nodes_in_layer_to <- filter(modules_for_similarity_shuffled, (layer_id == current_layer_to)) %>% select(node_id) %>% unlist()
      #take all nodes in layer_from and all nodes in layer_to to check turnover
      int_both <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in both layers
      uni_both <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to)
      turnover <- length(int_both)/length(uni_both)
      layers_turnover_shuf <- rbind(layers_turnover_shuf, tibble(layer_from= current_layer, layer_to= current_layer_to, turnover= turnover, 
                                                                 distance= current_distance))
    }
  }
  all_layers_turnover <- rbind(all_layers_turnover, tibble(layers_turnover_shuf, shuf))
  all_layers_turnover <- all_layers_turnover %>% unique()
}

ave_turnover <- all_layers_turnover %>% group_by(layer_from, layer_to, distance) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="null") #create mean and sd for each point

#view(ave_turnover)

empirical_turnover <- layers_turnover %>% group_by(layer_from, layer_to, distance) %>%
  summarise(ave=mean(turnover), sd=sd(turnover)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#view(empirical_turnover)

turnover_null_and_empirical <- rbind(empirical_turnover, ave_turnover)

turnover_null_and_empirical %>% ggplot(aes(x= distance, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+
  labs(x="distance in meters", y="similarity")



## number of modules in two islands shuf vs empirical
edge_list_by_islands_shuf_500 <- NULL
edge_list_island_combine_no_module_shuf_500 <- NULL

for (shuf in 1:max(shuffled_500$i)){ 
  current_shuf <- shuffled_500 %>% filter(i == shuf)
  print(shuf)
  modules_for_similarity_num_shuf <- current_shuf %>% select(module, layer_id) %>% 
    unique() %>% group_by(module) %>% filter(n()>1) %>% select(module) %>% unique() #check which modules occur in 2 or more layers
  modules_for_similarity_shuffled <- current_shuf %>%
    filter(module %in% modules_for_similarity_num_shuf$module) #only save the modules that are found in 2 or more layers
  
  module_pivoted_shuf <- pivot_by_module(modules_for_similarity_shuffled)
  module_pivoted_shuf <- module_pivoted_shuf[, -1]
  
  modules_edge_list_shuf <- NULL
  for (i in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the current_shuf frame
    modules_edge_list_shuf <- edge_list_per_module(module_pivoted_shuf[i,], modules_edge_list_shuf) 
  }
  
  
  edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
  edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA
  
  
  edge_list_by_islands_shuf <- edge_list_with_distances_shuf
  old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
  #new <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)
  new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
  edge_list_by_islands_shuf$layer_from[edge_list_by_islands_shuf$layer_from %in% old] <- new[match(edge_list_by_islands_shuf$layer_from, old)]
  edge_list_by_islands_shuf$layer_to[edge_list_by_islands_shuf$layer_to %in% old] <- new[match(edge_list_by_islands_shuf$layer_to, old)]
  
  edge_list_by_islands_shuf_500 <- rbind(edge_list_by_islands_shuf_500, tibble(edge_list_by_islands_shuf, shuf))
  
  #version with # of modules in layers
  edge_list_by_islands_modules_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to, module) %>%
    summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
  edge_list_by_islands_modules_shuf$count <- c(1)
  edge_list_by_islands_modules_shuf <- edge_list_by_islands_modules_shuf %>% mutate(number_of_modules= sum(count)) %>%
    select(layer_from, layer_to, module, number_of_modules) 
  
  #version with correct average between layers
  edge_list_by_islands_ave_shuf <- edge_list_by_islands %>% group_by(layer_from, layer_to) %>%
    summarise(ave_distance= mean(distance_in_meters)) %>% unique()
  
  #combine
  edge_list_island_combine_shuf <- edge_list_by_islands_ave_shuf %>%
    merge(edge_list_by_islands_modules_shuf, by= c("layer_from", "layer_to")) #merge both versions 
  
  edge_list_island_combine_no_module_shuf <- edge_list_island_combine_shuf %>% select(-module) %>% unique() #have version where modules aren't present
  
  edge_list_island_combine_no_module_shuf_500 <- rbind(edge_list_island_combine_no_module_shuf_500, 
                                                       tibble(edge_list_island_combine_no_module_shuf, shuf))
}



edge_list_island_combine_no_module_empirical <- edge_list_island_combine_no_module
edge_list_island_combine_no_module_empirical$shuf <- NA

edge_list_island_combine_no_module_empirical <- edge_list_island_combine_no_module_empirical %>%
  slice(rep(1:n(), each = 100)) #duplicate 100 times to i can actually see it in the graph


edge_list_island_combine_no_module_shuf_empirical <- rbind(edge_list_island_combine_no_module_shuf_500, 
                                                           edge_list_island_combine_no_module_empirical)

edge_list_island_combine_no_module_shuf_empirical %>% 
  mutate(type=ifelse(is.na(shuf) ,'empirical','null')) %>% 
  group_by(layer_from, layer_to) %>%
  ggplot()+
  geom_histogram(aes(x=number_of_modules, fill= type), alpha=0.6)+
  facet_wrap(layer_from ~ layer_to , scales = 'free')+ #create multi-panel plot aligning by layer_from
  scale_fill_manual(values = c('navy','plum'))+ theme_classic()+ 
  labs(x= "number of modules", y= "count")+
  scale_x_continuous(breaks=seq(2,13,2))

#view(edge_list_island_combine_no_module_shuf_empirical)