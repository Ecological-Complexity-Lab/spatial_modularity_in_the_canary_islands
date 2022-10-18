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

#this portion of the code refers to the uniform analysis where the interlayer edges are all uniform
#and fixed to the median value. distance decay in modules of the empirical network is compared to the
#distance decay in modules shown in the uniform network

##---- fixed for empirical data non hierarchical---------------------------------------------------------------------------

interlayer_edges_change <- select(inter_extended, -weight) #create data frame where weight doesn't exist
interlayer_edges_change$weight <- 0.357602 #set fixed interlayer value to median of interlayer distribution

dryad_multilayer_fixed <- create_multilayer_object(intra = intra_nonextended, 
                                                   inter = interlayer_edges_change, #create multilayer with new fixed inter value
                                                   nodes = physical_nodes,
                                                   layers = layer_metadata,
                                                   intra_output_extended = T) 

modules_edge_change_fixed <- run_infomap_multilayer(dryad_multilayer_fixed,
                                                    infomap_executable = "../Infomap",
                                                    flow_model = 'directed',
                                                    relax = F, 
                                                    silent = T, 
                                                    trials = 100,
                                                    seed = 497294, 
                                                    temporal_network = F)

# fixed
modules_fixed <- modules_edge_change_fixed$modules #39 modules

## empirical
modules_dryad_multilayer_fixed_analysis <- modules_dryad_multilayer$modules #create version just for the analysis

## create edge list with distances ------------------------------------------------------------------------------------

#get module_sub_module for both fixed and empirical multi lvl
module_pivoted_fixed <- pivot_by_module(modules_fixed) #pivot for fixed
module_pivoted_empirical_fixed_analysis <- pivot_by_module(modules_dryad_multilayer_fixed_analysis) #pivot for empirical

#----create edge list with distances-----------------------------
#for fixed

modules_edge_list_fixed_analysis <- NULL


for (k in (1:nrow(module_pivoted_fixed))){ #run the function for each row in the data frame
  modules_edge_list_fixed_analysis <- edge_list_per_module(module_pivoted_fixed[k,], modules_edge_list_fixed_analysis) 
  current_module_fixed <- rownames(module_pivoted_fixed)[k]
  if (is.null(modules_edge_list_fixed_analysis)) next
  modules_edge_list_fixed_analysis <- modules_edge_list_fixed_analysis %>% mutate(module = replace_na(module, current_module_fixed)) #add module number
}


# for empirical
modules_edge_list_empirical_fixed_analysis <- NULL


for (k in (1:nrow(module_pivoted_empirical_fixed_analysis))){ #run the function for each row in the data frame
  modules_edge_list_empirical_fixed_analysis <- edge_list_per_module(module_pivoted_empirical_fixed_analysis[k,], modules_edge_list_empirical_fixed_analysis) 
  current_module_empirical <- rownames(module_pivoted_empirical_fixed_analysis)[k]
  if (is.null(modules_edge_list_empirical_fixed_analysis)) next
  modules_edge_list_empirical_fixed_analysis <- modules_edge_list_empirical_fixed_analysis %>% mutate(module = replace_na(module, current_module_empirical)) #add module number
}

#edge list by island for fixed------------------------------------------------------------------------------------------------------------

edge_list_with_distances_fixed <- right_join(modules_edge_list_fixed_analysis, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances_fixed <- na.omit(edge_list_with_distances_fixed) #remove NA and delete layer name

#modules similarity pairwise distance between islands
edge_list_by_islands_fixed <- edge_list_with_distances_fixed
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
edge_list_by_islands_fixed$layer_from[edge_list_by_islands_fixed$layer_from %in% old] <- new[match(edge_list_by_islands_fixed$layer_from, old)]
edge_list_by_islands_fixed$layer_to[edge_list_by_islands_fixed$layer_to %in% old] <- new[match(edge_list_by_islands_fixed$layer_to, old)]

#version with # of modules in layers
edge_list_by_islands_modules_fixed <- edge_list_by_islands_fixed %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_islands_modules_fixed$count <- c(1)
edge_list_by_islands_modules_fixed <- edge_list_by_islands_modules_fixed %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules)

#version with correct average between layers
edge_list_by_islands_ave_fixed <- edge_list_by_islands_fixed %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_island_combine_fixed <- edge_list_by_islands_ave_fixed %>%
  merge(edge_list_by_islands_modules_fixed, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_island_combine_no_module_fixed_analysis <- edge_list_island_combine_fixed %>% select(-module) %>% unique() #have version where modules aren't present


# total number of modules in each layer

module_island_turnover_fixed_analysis <- NULL

island_list <- c("1","2","3","4east","4west","5","6")

for (i in island_list){
  for (j in island_list){
    modules_in_island_from_fixed_analysis <- filter(edge_list_by_islands_modules_fixed, layer_from == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_island_to_fixed_analysis <- filter(edge_list_by_islands_modules_fixed, layer_from == j) %>% select(module) %>% unique() %>% unlist()
    #take all sub modules in layer_from and all sub modules in layer_to to check turnover
    int_both <- intersect(modules_in_island_from_fixed_analysis, modules_in_island_to_fixed_analysis) #how many sub modules are common in both layers
    uni_both <- union(modules_in_island_from_fixed_analysis, modules_in_island_to_fixed_analysis) #how many sub modules are found in both layers in total
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover_fixed_analysis <- rbind(module_island_turnover_fixed_analysis, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

module_island_turnover_fixed_analysis <- drop_na(module_island_turnover_fixed_analysis)

edge_list_by_islands_ave_fixed <- edge_list_by_islands_fixed %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

islands_turnover_with_distnace_fixed <- edge_list_by_islands_ave_fixed %>%
  merge(module_island_turnover_fixed_analysis, by= c("layer_from", "layer_to")) #merge both versions

#---- edge list by island empirical-------------------------------------------------------------------------------------------------------------

edge_list_with_distances_empirical <- right_join(modules_edge_list_empirical_fixed_analysis, 
                                                 distances_with_ids, by= c("layer_from", "layer_to"))
edge_list_with_distances_empirical <- na.omit(edge_list_with_distances_empirical) #remove NA 


#modules similarity pairwise distance between islands
edge_list_by_islands_empirical <- edge_list_with_distances_empirical
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
edge_list_by_islands_empirical$layer_from[edge_list_by_islands_empirical$layer_from %in% old] <- 
  new[match(edge_list_by_islands_empirical$layer_from, old)]
edge_list_by_islands_empirical$layer_to[edge_list_by_islands_empirical$layer_to %in% old] <- 
  new[match(edge_list_by_islands_empirical$layer_to, old)]

#version with # of modules in layers
edge_list_by_islands_modules_empirical <- edge_list_by_islands_empirical %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) 
edge_list_by_islands_modules_empirical$count <- c(1)
edge_list_by_islands_modules_empirical <- edge_list_by_islands_modules_empirical %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules)


#total number of modules in each layer

module_island_turnover_empirical <- NULL

island_list <- c("1","2","3","4east","4west","5","6")

for (i in island_list){
  for (j in island_list){
    modules_in_island_from_empirical_analysis <- filter(edge_list_by_islands_modules_empirical, layer_from == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_island_to_empirical_analysis <- filter(edge_list_by_islands_modules_empirical, layer_from == j) %>% select(module) %>% unique() %>% unlist()
    #take all sub modules in layer_from and all sub modules in layer_to to check turnover
    int_both <- intersect(modules_in_island_from_empirical_analysis, modules_in_island_to_empirical_analysis) #how many sub modules are found in both layers
    uni_both <- union(modules_in_island_from_empirical_analysis, modules_in_island_to_empirical_analysis)
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover_empirical <- rbind(module_island_turnover_empirical, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

module_island_turnover_empirical <- drop_na(module_island_turnover_empirical)

edge_list_by_islands_ave_empirical <- edge_list_by_islands_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

islands_turnover_with_distnace_empirical <- edge_list_by_islands_ave_empirical %>%
  merge(module_island_turnover_empirical, by= c("layer_from", "layer_to")) #merge both versions

#----prepare data for graph--------------------------------------------------------------------------------------------------------

ave_module_island_turnover_shuf_fixed <- islands_turnover_with_distnace_fixed %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_fixed") #create mean and sd for each point

## empirical
empirical_turnover_for_module_island_analysis <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#combine all islands
jaccard_similarity_empirical_and_fixed <- rbind(empirical_turnover_for_module_island_analysis, 
                                                ave_module_island_turnover_shuf_fixed)

jaccard_similarity_empirical_and_fixed_no_self_loop <- jaccard_similarity_empirical_and_fixed %>% subset(layer_from != layer_to)

jaccard_similarity_empirical_and_fixed_no_self_loop_km <- jaccard_similarity_empirical_and_fixed_no_self_loop %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#just empirical
jaccard_similarity_empirica_no_self_loop_km <- jaccard_similarity_empirical_and_fixed_no_self_loop_km %>%
  filter(type == "empirical")

##---- graphs-----------------------------------------------------------------------------------------------------------

jaccard_similarity_empirica_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")

#empirical and null
jaccard_similarity_empirical_and_fixed_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400)+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.65, 0.62))+ scale_color_manual(values = c("#F47069", "#c4ad06"))

#version with no trendline
jaccard_similarity_empirical_and_fixed_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ 
  geom_smooth(data = jaccard_similarity_empirica_no_self_loop_km, method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ scale_color_manual(values = c("#F47069", "#c4ad06"))
