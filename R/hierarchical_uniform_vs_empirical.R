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
library(ggtree)

#this portion of the code comes to compare the empirical network to the null model version where 
#interlayer edges are uniform.

#---- multilayer network with uniform interlayer edges-------------------------------------------------------------------------

interlayer_edges_change <- select(inter_extended, -weight) #create data frame where weight doesn't exist
interlayer_edges_change$weight <- 0.357602 #set fixed interlayer value to median of interlayer distribution

dryad_multilayer_multi_lvl_fixed <- create_multilayer_object(intra = intra_nonextended, 
                                                             inter = interlayer_edges_change, #create multilayer with new fixed inter value
                                                             nodes = physical_nodes,
                                                             layers = layer_metadata,
                                                             intra_output_extended = T) 

modules_edge_change <- multi_lvl_infomap(dryad_multilayer_multi_lvl_fixed, #run as multi lvl
                                         infomap_executable = "../Infomap",
                                         flow_model = 'directed',
                                         relax = F, 
                                         silent = T, 
                                         trials = 100,
                                         seed = 497294, 
                                         temporal_network = F)

# analysis will be done on sub modules so leave only rows where they are present
# fixed
modules_fixed_multi_lvl <- modules_edge_change$modules #39 modules and 4 sub modules. no sub sun modules.

modules_fixed_multi_lvl$module_sub_module <- paste(modules_fixed_multi_lvl$module, ".", modules_fixed_multi_lvl$sub_module) #create new column combining module and sub module
modules_fixed_multi_lvl_sub_module <- modules_fixed_multi_lvl %>% drop_na(sub_module) #remove every row where there's no sub module

#----empirical counterpart for the analysis------------------------------------------------------
#modules_dryad_multilayer_multi_level <- read.csv("csvs/modules_dryad_multilayer_multi_level.csv")


modules_dryad_multilayer_multi_level_fixed_analysis <- modules_dryad_multilayer_multi_level #create version just for the analysis

modules_dryad_multilayer_multi_level_fixed_analysis$module_sub_module <- 
  paste(modules_dryad_multilayer_multi_level_fixed_analysis$module, ".", modules_dryad_multilayer_multi_level_fixed_analysis$sub_module)
modules_dryad_multilayer_multi_level_fixed_analysis_sub_module <- modules_dryad_multilayer_multi_level_fixed_analysis %>% drop_na(sub_module)

#write.csv(modules_dryad_multilayer_multi_level_fixed_analysis_sub_module, 
#          "./csvs/modules_dryad_multilayer_multi_level_fixed_analysis_sub_module.csv", row.names = FALSE)

#---- create edge list with distances ------------------------------------------------------------------------------------

#get module_sub_module for both fixed and empirical multi lvl
module_pivoted_fixed_multi_lvl <- pivot_by_sub_module(modules_fixed_multi_lvl_sub_module) #pivot for fixed
module_pivoted_empirical_multi_lvl <- pivot_by_sub_module(modules_dryad_multilayer_multi_level_fixed_analysis_sub_module) #pivot for empirical


#create edge list with distances
#similarity between state nodes for modules found in 2 or more layers

#---- edge list for uniform--------------------------------------------------------------------------------------

modules_edge_list_fixed <- NULL

edge_list_per_sub_module_fixed <- function(data,edge_list){
  #gets one row from a data frame and creates an edge list from it
  for (i in (1:10)){
    if (data[i] == 0) next #only take layers where the module is present
    else {
      for (j in (i+1):11){
        if (data[j] == 0) next #only take layers where the module is present
        else {
          edge_list <- rbind(edge_list, tibble(layer_from=i, layer_to=j, module_sub_module=NA)) #create edge list of all the layer found in a module
        }
      }
    }
  }
  return(edge_list)
}

for (k in (1:nrow(module_pivoted_fixed_multi_lvl))){ #run the function for each row in the data frame
  modules_edge_list_fixed <- edge_list_per_sub_module_fixed(module_pivoted_fixed_multi_lvl[k,], modules_edge_list_fixed) 
  current_sub_module <- rownames(module_pivoted_fixed_multi_lvl)[k]
  if (is.null(modules_edge_list_fixed)) next
  modules_edge_list_fixed <- modules_edge_list_fixed %>% mutate(module_sub_module = replace_na(module_sub_module, current_sub_module)) #add module number
}

#---- edge list for empirical--------------------------------------------------------------------------------------------------

modules_edge_list_multi_lvl_empirical <- NULL


for (k in (1:nrow(module_pivoted_empirical_multi_lvl))){ #run the function for each row in the data frame
  modules_edge_list_multi_lvl_empirical <- edge_list_per_sub_module(module_pivoted_empirical_multi_lvl[k,], modules_edge_list_multi_lvl_empirical) 
  current_sub_module <- rownames(module_pivoted_empirical_multi_lvl)[k]
  if (is.null(modules_edge_list_multi_lvl_empirical)) next
  modules_edge_list_multi_lvl_empirical <- modules_edge_list_multi_lvl_empirical %>% mutate(module_sub_module = replace_na(module_sub_module, current_sub_module)) #add module number
}

#---- distance decay for islands-----------------------------------------------------------------------

edge_list_with_distances_multi_lvl_fixed <- right_join(modules_edge_list_fixed, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances_multi_lvl_fixed <- na.omit(edge_list_with_distances_multi_lvl_fixed) #remove NA and delete layer name

#modules similarity pairwise distance between islands
edge_list_by_islands_multi_lvl_fixed <- edge_list_with_distances_multi_lvl_fixed
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
edge_list_by_islands_multi_lvl_fixed$layer_from[edge_list_by_islands_multi_lvl_fixed$layer_from %in% old] <- new[match(edge_list_by_islands_multi_lvl_fixed$layer_from, old)]
edge_list_by_islands_multi_lvl_fixed$layer_to[edge_list_by_islands_multi_lvl_fixed$layer_to %in% old] <- new[match(edge_list_by_islands_multi_lvl_fixed$layer_to, old)]

#version with # of modules in layers
edge_list_by_islands_modules_multi_lvl_fixed <- edge_list_by_islands_multi_lvl_fixed %>% group_by(layer_from, layer_to, module_sub_module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_islands_modules_multi_lvl_fixed$count <- c(1)
edge_list_by_islands_modules_multi_lvl_fixed <- edge_list_by_islands_modules_multi_lvl_fixed %>% mutate(number_of_sub_modules= sum(count)) %>%
  select(layer_from, layer_to, module_sub_module, number_of_sub_modules)

#version with correct average between layers
edge_list_by_islands_ave_multi_lvl_fixed <- edge_list_by_islands_multi_lvl_fixed %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_island_combine_multi_lvl_fixed <- edge_list_by_islands_ave_multi_lvl_fixed %>%
  merge(edge_list_by_islands_modules_multi_lvl_fixed, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_island_combine_no_module_fixed <- edge_list_island_combine_multi_lvl_fixed %>% select(-module_sub_module) %>% unique() #have version where modules aren't present


#total number of modules in each layer
module_island_turnover_fixed <- NULL

island_list <- c("1","2","3","4east","4west","5","6")

for (i in island_list){
  for (j in island_list){
    modules_in_island_from_fixed <- filter(edge_list_by_islands_modules_multi_lvl_fixed, layer_from == i) %>% select(module_sub_module) %>% unique() %>% unlist()
    modules_in_island_to_fixed <- filter(edge_list_by_islands_modules_multi_lvl_fixed, layer_from == j) %>% select(module_sub_module) %>% unique() %>% unlist()
    #take all sub modules in layer_from and all sub modules in layer_to to check turnover
    int_both <- intersect(modules_in_island_from_fixed, modules_in_island_to_fixed) #how many sub modules are found in both layers
    uni_both <- union(modules_in_island_from_fixed, modules_in_island_to_fixed)
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover_fixed <- rbind(module_island_turnover_fixed, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

module_island_turnover_fixed <- drop_na(module_island_turnover_fixed)

edge_list_by_islands_ave_multi_lvl_fixed <- edge_list_by_islands_multi_lvl_fixed %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

islands_turnover_with_distnace_multi_lvl_fixed <- edge_list_by_islands_ave_multi_lvl_fixed %>%
  merge(module_island_turnover_fixed, by= c("layer_from", "layer_to")) #merge both versions

#---- distance decay for the empirical-----------------------------------------------------------------------------

edge_list_with_distances_multi_lvl_empirical <- right_join(modules_edge_list_multi_lvl_empirical, distances_with_ids, by= c("layer_from", "layer_to"))
edge_list_with_distances_multi_lvl_empirical <- na.omit(edge_list_with_distances_multi_lvl_empirical) #remove NA and delete layer name

#modules similarity pairwise distance between islands
edge_list_by_islands_multi_lvl_empirical <- edge_list_with_distances_multi_lvl_empirical
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
edge_list_by_islands_multi_lvl_empirical$layer_from[edge_list_by_islands_multi_lvl_empirical$layer_from %in% old] <- new[match(edge_list_by_islands_multi_lvl_empirical$layer_from, old)]
edge_list_by_islands_multi_lvl_empirical$layer_to[edge_list_by_islands_multi_lvl_empirical$layer_to %in% old] <- new[match(edge_list_by_islands_multi_lvl_empirical$layer_to, old)]

#version with # of modules in layers
edge_list_by_islands_modules_multi_lvl_empirical <- edge_list_by_islands_multi_lvl_empirical %>% group_by(layer_from, layer_to, module_sub_module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_islands_modules_multi_lvl_empirical$count <- c(1)
edge_list_by_islands_modules_multi_lvl_empirical <- edge_list_by_islands_modules_multi_lvl_empirical %>% mutate(number_of_sub_modules= sum(count)) %>%
  select(layer_from, layer_to, module_sub_module, number_of_sub_modules)

#version with correct average between layers
edge_list_by_islands_ave_multi_lvl_empirical <- edge_list_by_islands_multi_lvl_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_island_combine_multi_lvl_empirical <- edge_list_by_islands_ave_multi_lvl_empirical %>%
  merge(edge_list_by_islands_modules_multi_lvl_empirical, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_island_combine_no_module_empirical<- edge_list_island_combine_multi_lvl_empirical %>% select(-module_sub_module) %>% unique() #have version where modules aren't present



#total number of modules in each layer
module_island_turnover_multi_lvl_empirical <- NULL

island_list <- c("1","2","3","4east","4west","5","6")

for (i in island_list){
  for (j in island_list){
    modules_in_island_from_empirical <- filter(edge_list_by_islands_modules_multi_lvl_empirical, layer_from == i) %>% select(module_sub_module) %>% unique() %>% unlist()
    modules_in_island_to_empirical <- filter(edge_list_by_islands_modules_multi_lvl_empirical, layer_from == j) %>% select(module_sub_module) %>% unique() %>% unlist()
    #take all sub modules in layer_from and all sub modules in layer_to to check turnover
    int_both <- intersect(modules_in_island_from_empirical, modules_in_island_to_empirical) #how many sub modules are found in both layers
    uni_both <- union(modules_in_island_from_empirical, modules_in_island_to_empirical)
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover_multi_lvl_empirical <- rbind(module_island_turnover_multi_lvl_empirical, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

module_island_turnover_multi_lvl_empirical <- drop_na(module_island_turnover_multi_lvl_empirical)

edge_list_by_islands_ave_multi_lvl_empirical <- edge_list_by_islands_multi_lvl_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

islands_turnover_with_distnace_multi_lvl_empirical <- edge_list_by_islands_ave_multi_lvl_empirical %>%
  merge(module_island_turnover_multi_lvl_empirical, by= c("layer_from", "layer_to")) #merge both versions

#---- combining the data for graphs---------------------------------------------------------------------------------------------------
## fixed
ave_module_island_turnover_shuf_multi_lvl_fixed <- islands_turnover_with_distnace_multi_lvl_fixed %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_model") #create mean and sd for each point

## empirical
empirical_turnover_for_module_island_multi_lvl <- islands_turnover_with_distnace_multi_lvl_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#combine all islands
jaccard_similarity_empirical_and_fixed_multi_lvl <- rbind(empirical_turnover_for_module_island_multi_lvl, 
                                                          ave_module_island_turnover_shuf_multi_lvl_fixed)

jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop <- jaccard_similarity_empirical_and_fixed_multi_lvl %>% subset(layer_from != layer_to)

jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_km <- jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#just empirical
jaccard_similarity_empirica_multi_lvl_no_self_loop_km <- jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_km %>%
  filter(type == "empirical")

#---- graphs----------------------------------------------------------------------------------------------------

#no self loop island
#just emprical
jaccard_similarity_empirica_multi_lvl_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="distance in km", y="Jaccard Similarity")

#empirical and null
jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="distance in km", y="Jaccard Similarity")

#----statistical analysis-----------------------------------------------------------

lm1_sub_module_multi_lvl = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop,
                                                          jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop$type == "empirical")) #in empirical
lm2_sub_module_multi_lvl = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop,
                                                          jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop$type == "null_model")) #in null pols

#get equations
lm1_sub_module_multi_lvl_equation <- paste("y=", coef(lm1_sub_module_multi_lvl)[[1]], "+", coef(lm1_sub_module_multi_lvl)[[2]], "*x")
lm2_sub_module_multi_lvl_equation <- paste("y=", coef(lm2_sub_module_multi_lvl)[[1]], "+", coef(lm2_sub_module_multi_lvl)[[2]], "*x")
lm1_sub_module_multi_lvl_equation
lm2_sub_module_multi_lvl_equation

b1_sub_module_multi_lvl <- summary(lm1_sub_module_multi_lvl)$coefficients[2,1]
se1_sub_module_multi_lvl <- summary(lm1_sub_module_multi_lvl)$coefficients[2,2]
b2_sub_module_multi_lvl <- summary(lm2_sub_module_multi_lvl)$coefficients[2,1]
se2_sub_module_multi_lvl <- summary(lm2_sub_module_multi_lvl)$coefficients[2,2]

p_value_sub_module_multi_lvl = 2*pnorm(-abs(compare.coeff(b1_sub_module_multi_lvl,se1_sub_module_multi_lvl,
                                                          b2_sub_module_multi_lvl,se2_sub_module_multi_lvl)))
p_value_sub_module_multi_lvl

#----fixed vs empirical only sub modules where both plant and pols are found --------------------------------------------------------
#version where sub-modules made up of only plants or only pollinators are disregarded
modules_dryad_multilayer_multi_level_fixed_analysis_sub_module_hetero <- modules_dryad_multilayer_multi_level_fixed_analysis_sub_module %>%
  filter(!sub_module %in% c(41,40,39,38,37,36,35,34,32,30,27,22,21,20,19,15,10)) #delete all sub modules where there are only plants or only pols

#fixed version only has 4 sub modules and all of them have both plants and pols

#pivot
module_pivoted_empirical_multi_lvl_hetero <- pivot_by_sub_module(modules_dryad_multilayer_multi_level_fixed_analysis_sub_module_hetero) #pivot for empirical


# add distances
modules_edge_list_multi_lvl_empirical_hetero <- NULL

for (k in (1:nrow(module_pivoted_empirical_multi_lvl_hetero))){ #run the function for each row in the data frame
  modules_edge_list_multi_lvl_empirical_hetero <- edge_list_per_sub_module(module_pivoted_empirical_multi_lvl_hetero[k,], modules_edge_list_multi_lvl_empirical_hetero) 
  current_sub_module <- rownames(module_pivoted_empirical_multi_lvl_hetero)[k]
  if (is.null(modules_edge_list_multi_lvl_empirical_hetero)) next
  modules_edge_list_multi_lvl_empirical_hetero <- modules_edge_list_multi_lvl_empirical_hetero %>% mutate(module_sub_module = replace_na(module_sub_module, current_sub_module)) #add module number
}


edge_list_with_distances_multi_lvl_empirical_hetero <- right_join(modules_edge_list_multi_lvl_empirical_hetero, distances_with_ids, by= c("layer_from", "layer_to"))
edge_list_with_distances_multi_lvl_empirical_hetero <- na.omit(edge_list_with_distances_multi_lvl_empirical_hetero) #remove NA and delete layer name

#modules similarity pairwise distance between islands
edge_list_by_islands_multi_lvl_empirical_hetero <- edge_list_with_distances_multi_lvl_empirical_hetero
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
edge_list_by_islands_multi_lvl_empirical_hetero$layer_from[edge_list_by_islands_multi_lvl_empirical_hetero$layer_from %in% old] <- 
  new[match(edge_list_by_islands_multi_lvl_empirical_hetero$layer_from, old)]
edge_list_by_islands_multi_lvl_empirical_hetero$layer_to[edge_list_by_islands_multi_lvl_empirical_hetero$layer_to %in% old] <-
  new[match(edge_list_by_islands_multi_lvl_empirical_hetero$layer_to, old)]

#version with # of modules in layers
edge_list_by_islands_modules_multi_lvl_empirical_hetero <- edge_list_by_islands_multi_lvl_empirical_hetero %>% group_by(layer_from, layer_to, module_sub_module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_islands_modules_multi_lvl_empirical_hetero$count <- c(1)
edge_list_by_islands_modules_multi_lvl_empirical_hetero <- edge_list_by_islands_modules_multi_lvl_empirical_hetero %>% mutate(number_of_sub_modules= sum(count)) %>%
  select(layer_from, layer_to, module_sub_module, number_of_sub_modules)

#version with correct average between layers
edge_list_by_islands_ave_multi_lvl_empirical_hetero <- edge_list_by_islands_multi_lvl_empirical_hetero %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_island_combine_multi_lvl_empirical_hetero <- edge_list_by_islands_ave_multi_lvl_empirical_hetero %>%
  merge(edge_list_by_islands_modules_multi_lvl_empirical_hetero, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_island_combine_no_module_empirical_hetero <- edge_list_island_combine_multi_lvl_empirical_hetero %>% select(-module_sub_module) %>% unique() #have version where modules aren't present



#total number of modules in each layer
module_island_turnover_multi_lvl_empirical_hetero <- NULL

island_list <- c("1","2","3","4east","4west","5","6")

for (i in island_list){
  for (j in island_list){
    modules_in_island_from_empirical_hetero <- filter(edge_list_by_islands_modules_multi_lvl_empirical_hetero, layer_from == i) %>% 
      select(module_sub_module) %>% unique() %>% unlist()
    modules_in_island_to_empirical_hetero <- filter(edge_list_by_islands_modules_multi_lvl_empirical_hetero, layer_from == j) %>% 
      select(module_sub_module) %>% unique() %>% unlist()
    #take all sub modules in layer_from and all sub modules in layer_to to check turnover
    int_both <- intersect(modules_in_island_from_empirical_hetero, modules_in_island_to_empirical_hetero) #how many sub modules are found in both layers
    uni_both <- union(modules_in_island_from_empirical_hetero, modules_in_island_to_empirical_hetero)
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover_multi_lvl_empirical_hetero <- rbind(module_island_turnover_multi_lvl_empirical_hetero, 
                                                               tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

module_island_turnover_multi_lvl_empirical_hetero <- drop_na(module_island_turnover_multi_lvl_empirical_hetero)

edge_list_by_islands_ave_multi_lvl_empirical_hetero <- edge_list_by_islands_multi_lvl_empirical_hetero %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

islands_turnover_with_distnace_multi_lvl_empirical_hetero <- edge_list_by_islands_ave_multi_lvl_empirical_hetero %>%
  merge(module_island_turnover_multi_lvl_empirical_hetero, by= c("layer_from", "layer_to")) #merge both versions


## empirical
empirical_turnover_for_module_island_multi_lvl_hetero <- islands_turnover_with_distnace_multi_lvl_empirical_hetero %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% 
  mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#combine all islands
jaccard_similarity_empirical_and_fixed_multi_lvl_hetero <- rbind(empirical_turnover_for_module_island_multi_lvl_hetero, 
                                                                 ave_module_island_turnover_shuf_multi_lvl_fixed)

jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_hetero <- jaccard_similarity_empirical_and_fixed_multi_lvl_hetero %>% 
  subset(layer_from != layer_to)

jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_km_hetero <- jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_hetero %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#just empirical
jaccard_similarity_empirica_multi_lvl_no_self_loop_km_hetero <- jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_km_hetero %>%
  filter(type == "empirical")

#----graphs---------------------------------------------------------------------------------------------------
#can't really compare the two as the null model only has 4 sub-modules
#no self loop island
#just emprical
jaccard_similarity_empirica_multi_lvl_no_self_loop_km_hetero %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")

#empirical and null
jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_km_hetero %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")

#---- statistical analysis------------------------------------------------------------
#check if its significant
lm1_sub_module_multi_lvl_hetero = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_hetero,
                                                                 jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_hetero$type == "empirical")) #in empirical
lm2_sub_module_multi_lvl_hetero = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_hetero,
                                                                 jaccard_similarity_empirical_and_fixed_multi_lvl_no_self_loop_hetero$type == "null_model")) #in null pols

#get equations
lm1_sub_module_multi_lvl_equation_hetero <- paste("y=", coef(lm1_sub_module_multi_lvl_hetero)[[1]], "+", coef(lm1_sub_module_multi_lvl_hetero)[[2]], "*x")
lm2_sub_module_multi_lvl_equation_hetero <- paste("y=", coef(lm2_sub_module_multi_lvl_hetero)[[1]], "+", coef(lm2_sub_module_multi_lvl_hetero)[[2]], "*x")
lm1_sub_module_multi_lvl_equation_hetero
lm2_sub_module_multi_lvl_equation_hetero

b1_sub_module_multi_lvl_hetero <- summary(lm1_sub_module_multi_lvl_hetero)$coefficients[2,1]
se1_sub_module_multi_lvl_hetero <- summary(lm1_sub_module_multi_lvl_hetero)$coefficients[2,2]
b2_sub_module_multi_lvl_hetero <- summary(lm2_sub_module_multi_lvl_hetero)$coefficients[2,1]
se2_sub_module_multi_lvl_hetero <- summary(lm2_sub_module_multi_lvl_hetero)$coefficients[2,2]

p_value_sub_module_multi_lvl_hetero = 2*pnorm(-abs(compare.coeff(b1_sub_module_multi_lvl_hetero,se1_sub_module_multi_lvl_hetero,
                                                                 b2_sub_module_multi_lvl_hetero,se2_sub_module_multi_lvl_hetero)))
p_value_sub_module_multi_lvl_hetero




#---- dendrogram for fixed---------------------------------------------------------------------------------------------------

multi_lvl_to_dendogram_fixed <- modules_fixed_multi_lvl %>% select(-module_sub_module)
multi_lvl_to_dendogram_fixed$state_node <- paste(multi_lvl_to_dendogram_fixed$node_id, "_", multi_lvl_to_dendogram_fixed$layer_id) #create a column for state nodes

modules_for_dendrogram <- multi_lvl_to_dendogram_fixed$module %>% unlist() %>% unique()

network2module_fixed <- data.frame(network = "network",
                                   module = modules_for_dendrogram) #create 1st split in network to 2 modules


#network level 
edges_network2module_level_fixed <- network2module_fixed %>% select(network, module) %>% 
  unique %>% rename(from = network, to = module) #create edge list from module to sub module

edges_network2module_level_fixed$to <- paste(edges_network2module_level_fixed$to, "_", "m") #add m to describe module


#top level
edges_module2sub_level_fixed <- multi_lvl_to_dendogram_fixed %>% select(module, sub_module) %>% 
  unique %>% rename(from = module, to = sub_module) #create edge list from module to sub module

edges_module2sub_level_fixed <- drop_na(edges_module2sub_level_fixed) #delete all occasions where there's no sub sub module

edges_module2sub_level_fixed$from <- paste(edges_module2sub_level_fixed$from, "_", "m") #add m to describe module
edges_module2sub_level_fixed$to <- paste(edges_module2sub_level_fixed$to, "_", "sm") #add sm to describe sub module

edges_module2sub_level_fixed$to <- paste("3." , edges_module2sub_level_fixed$to) #add module id to sub module

#no sub sub modules

#bottom level option 1
edges_sub2state_level_fixed <- multi_lvl_to_dendogram_fixed %>% select(sub_module, state_node) %>% 
  unique %>% rename(from = sub_module, to = state_node) #create edge list from sub module to state node

edges_sub2state_level_fixed <- drop_na(edges_sub2state_level_fixed) #delete all occasions where there's no sub sub module

edges_sub2state_level_fixed$from <- paste(edges_sub2state_level_fixed$from, "_", "sm") #add sm to describe sub module

edges_sub2state_level_fixed$from <- paste("3." , edges_sub2state_level_fixed$from) #distinguish between sub modules with the same id under different modules


#bottom level option 2
edges_module2state_level_fixed <- multi_lvl_to_dendogram_fixed %>% select(module, state_node) %>% 
  unique %>% rename(from = module, to = state_node) #create edge list from sub module to state node

edges_module2state_level_fixed <- edges_module2state_level_fixed[!(edges_module2state_level_fixed$to %in% edges_sub2state_level_fixed$to),] ## if there is a sub module 
#need to delete the corresponding edge in the module to state node

edges_module2state_level_fixed$from <- paste(edges_module2state_level_fixed$from, "_", "m") #add sm to describe sub module



## create edge list combined for dendrogram
edge_list_for_dendrogram_fixed <- rbind(edges_network2module_level_fixed, edges_module2sub_level_fixed, edges_sub2state_level_fixed,
                                        edges_module2state_level_fixed)

no_state_nodes_fixed <- rbind(edges_network2module_level_fixed, edges_module2sub_level_fixed) #version with no state nodes

#---- dendrogram graphs-------------------------------------------------------------------------------

## graphs
pre_graph_fixed <- graph_from_data_frame(edge_list_for_dendrogram_fixed)

pre_graph_no_state_nodes_fixed <- graph_from_data_frame(no_state_nodes_fixed)


ggtree(pre_graph_no_state_nodes_fixed, layout = "circular") + geom_tiplab()

ggtree(pre_graph_fixed, layout = "circular") 
