########### --------- EXTRA ANALYSIS --------

## -- 1. SENSITIVTY ANALYSIS NULL MODEL 3 

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


setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")

## -- Example to show Shai

## -- The results are different when I change manually the interedges weight and use relax rates. For example, if I fix
#interedges weights and relax rate of 0.9, the results are completely different.

physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")

inter_extended <- read.csv("./csvs/Islands/dryad_only_interlayer_edges_islands_as_layers.csv")
intra_nonextended <- read.csv("./csvs/Islands/dryad_only_intralayer_edges_islands_as_layers.csv")

inter_extended <- inter_extended [,-1]
intra_nonextended <- intra_nonextended [,-1]

#Interedges weight
inter_extended$weight <- 0.9

man<-create_multilayer_object(intra = intra_nonextended,
                              inter = inter_extended,
                              nodes = physical_nodes, #nodes are always the same. we're not deleting nodes.
                              layers = layer_metadata, #layers are always the same. we're not deleting layers.
                              intra_output_extended = T)


man_mod <- run_infomap_multilayer(man, relax = F, silent = T, trials = 100, seed = 497294, multilayer_relax_rate = r, multilayer_relax_limit_up = 1, multilayer_relax_limit_down = 0, temporal_network = F)

#Relax rate
r=0.9
rel<-create_multilayer_object(intra = intra_nonextended,
                              nodes = physical_nodes, #nodes are always the same. we're not deleting nodes.
                              layers = layer_metadata, #layers are always the same. we're not deleting layers.
                              intra_output_extended = F)


rel_mod <- run_infomap_multilayer(rel, relax = T, silent = T, trials = 100, seed = 497294, multilayer_relax_rate = r, multilayer_relax_limit_up = 1, multilayer_relax_limit_down = 0, temporal_network = F)




##---- fixed for empirical data ---------------------------------------------------------------------------
physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")

inter_extended <- read.csv("./csvs/Islands/dryad_only_interlayer_edges_islands_as_layers.csv")
intra_nonextended <- read.csv("./csvs/Islands/dryad_only_intralayer_edges_islands_as_layers.csv")

inter_extended <- inter_extended [,-1]
intra_nonextended <- intra_nonextended [,-1]

## - Set fixed interlayer value to median of interlayer distribution
interlayer_edges_change <- select(inter_extended, -weight) #create data frame where weight doesn't exist

##Create scenarios where we change the interedges weight from 0 to 1 
k<-seq(0.1,1,by=0.1)

Sen_list_classic = NULL

for(t in k){
  output<- interlayer_edges_change %>% 
    mutate(weight = t)#set fixed interlayer value to t
  
  Sen_list_classic<- rbind(Sen_list_classic,output)
}

  
# ----multilayer_class-----------------------------------------------------------------------------------------------
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")
physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")

# Input: An extended edge list.
dryad_edgelist_complete_fixed <- Sen_list_classic
dryad_multilayer_complete_fixed <- NULL

#Function to calculate modularity according to each interedge weight
modularity_for_fixed <- function(edge_list, output){
  for(trial_num in k){
    current_trial_edgelist <- edge_list%>% filter(weight == trial_num) #take 1 trial at a time to create multilayer
    dryad_multilayer_shuf_trial <- create_multilayer_object(intra = intra_nonextended,
                                                            inter = current_trial_edgelist,
                                                            nodes = physical_nodes, #nodes are always the same. we're not deleting nodes.
                                                            layers = layer_metadata, #layers are always the same. we're not deleting layers.
                                                            intra_output_extended = T)
    
    
    #create modules for empirical network
    modules_dryad_multilayer_shuf_1000 <- modified_multi(dryad_multilayer_shuf_trial, 
                                                         infomap_executable = "Infomap",
                                                         flow_model = 'directed',
                                                         relax = F, 
                                                         silent = T, 
                                                         trials = 100,
                                                         seed = 497294, #always the same seed
                                                         temporal_network = F)
    
    
    output <- rbind(output, tibble(modules_dryad_multilayer_shuf_1000$modules, trial_num)) #
  }  
  return(output)
}

#calculate modularity according to each interedges weight
dryad_multilayer_complete_fixed_output <- modularity_for_fixed(dryad_edgelist_complete_fixed, 
                                                               dryad_multilayer_complete_fixed)

dryad_multilayer_complete_fixed_output <- dryad_multilayer_complete_fixed_output %>% drop_na() 


##---- distance decay of modules  ---------------------------------------------------------------
distances_with_ids <- read.csv("./csvs/Islands/distances_with_ids_islands_as_layers.csv")

# this function calculates the Jaccard Similarity in modules between islands (modified for fixed values)
module_distance_decay_islands_func <- function(multilayer_1000, 
                                               layers_turnover_with_distnace){
  
  for (trial in k){
    print(trial) #to keep tab on how far along we are
    modules_for_similarity_shuf <- multilayer_1000  %>% filter(trial_num == trial) #take only 1 trial

    #pivot modules
    module_pivoted_shuf <- pivot_by_module_islands(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    #create edge list with distances
    modules_edge_list_shuf <- NULL
    
    for (u in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module_islands(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[u]
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



#pivot modules function for islands as layers
pivot_by_module_islands <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(data, id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s2<-na.omit(s2)
  s3 <-s2%>%
    select(where(~ any(. != 0)))
  s4 = t(s3) 
  s4 <- s4[-1,]
  colnames(s4) <- c(1,2,3,4,5,6,7)
  return(s4)
}#delete all NAs that are in the physical nodes but not the network (and were not assigned to modules)

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

turnover_with_distance_classic <- NULL
layers_turnover_with_distnace<-NULL
module_layer_turnover_shuf <- NULL


all_edge_list_islands_combine_no_module_fixed_output <- module_distance_decay_islands_func(dryad_multilayer_complete_fixed_output,
                                                                                                  turnover_with_distance_classic)

#---- create ave for jaccard layers-------------------------------------------------------------------------
ave_module_islands_turnover_fixed <- all_edge_list_islands_combine_no_module_fixed_output  %>% 
  group_by(trial,layer_from, layer_to) %>% 
  summarise(ave=mean(turnover), ave_dist=mean(mean_distance)) %>% select("layer_from","layer_to","ave","ave_dist",
                                                                        "trial")
ave_module_islands_turnover_fixed$trial<-as.character(ave_module_islands_turnover_fixed$trial)

## empirical
islands_turnover_with_distnace_empirical <- read.csv("csvs/Islands/islands_turnover_with_distnace_empirical.csv")

empirical_turnover_for_modules_layers <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), ave_dist = mean_distance) %>% mutate(trial="empirical")

#combine all layers
jaccard_similarity_empirical_and_fixed <- rbind(empirical_turnover_for_modules_layers, 
                                                ave_module_islands_turnover_fixed)

jaccard_similarity_empirical_and_fixed_km <- jaccard_similarity_empirical_and_fixed  %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#write.csv(jaccard_similarity_empirical_and_fixed_km, "./csvs/Islands/sensitivity_nullmodel4.csv", row.names = FALSE) #so it can be used for classical shuffling

#Plot
jaccard_similarity_empirical_and_fixed_km <- read.csv("csvs/Islands/sensitivity_nullmodel4.csv")

pdf('./graphs/Islands/sensitivity_M4.pdf', 10, 6)
jaccard_similarity_empirical_and_fixed_km %>% 
  ggplot(aes(x= ave_dist_in_km , y= ave, group= trial, color= trial))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  scale_color_discrete(name = "Interedges weight")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank()) 
dev.off()

#Linear regressions
l_emp = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                       jaccard_similarity_empirical_and_fixed_km$trial=="empirical")) 
l_1 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                       jaccard_similarity_empirical_and_fixed_km$trial=="0.1")) 
l_2 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                       jaccard_similarity_empirical_and_fixed_km$trial=="0.2")) 
l_3 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                       jaccard_similarity_empirical_and_fixed_km$trial=="0.3")) 
l_4 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                       jaccard_similarity_empirical_and_fixed_km$trial=="0.4"))
l_5 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                     jaccard_similarity_empirical_and_fixed_km$trial=="0.5")) 
l_6 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                     jaccard_similarity_empirical_and_fixed_km$trial=="0.6")) 
l_7 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                     jaccard_similarity_empirical_and_fixed_km$trial=="0.7")) 
l_8 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                     jaccard_similarity_empirical_and_fixed_km$trial=="0.8"))
l_9 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                     jaccard_similarity_empirical_and_fixed_km$trial=="0.9")) 
l_10 = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_km,
                                     jaccard_similarity_empirical_and_fixed_km$trial=="1"))

summary(l_emp)
summary(l_1)
summary(l_2)
summary(l_3)
summary(l_4)
summary(l_5)
summary(l_6)
summary(l_7)
summary(l_8)
summary(l_9)
summary(l_10)


#------#WE CHECK THE SENSITIVITY ANALYSIS REMOVING THE WEIGHT OF INTERLAYEREDGES AND CHANGING THE RELAX RATES

  ##---- fixed for empirical data ---------------------------------------------------------------------------

physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")

inter_extended <- read.csv("./csvs/Islands/dryad_only_interlayer_edges_islands_as_layers.csv")
intra_nonextended <- read.csv("./csvs/Islands/dryad_only_intralayer_edges_islands_as_layers.csv")


inter_extended <- inter_extended [,-1]
intra_nonextended <- intra_nonextended [,-1]

#Create multilayer object changing the relax rates to simulate the sensitivity analysis


dryad_multilayer_relax <- create_multilayer_object(intra = intra_nonextended,
                                                        inter = inter_extended,
                                                        nodes = physical_nodes, #nodes are always the same. we're not deleting nodes.
                                                        layers = layer_metadata, #layers are always the same. we're not deleting layers.
                                                        intra_output_extended = F)

# Ignore interlayer edges
dryad_multilayer_relax$inter <- NULL

#Run Infomap and return the network's modular structure at increasing relax-rates.
relaxrate_modules <- NULL
for (r in seq(0.1,1, by = 0.1)){
  print(r)
  modules_relax_rate <- run_infomap_multilayer(dryad_multilayer_relax, relax = T, silent = T, trials = 100, seed = 497294, multilayer_relax_rate = r, multilayer_relax_limit_up = 1, multilayer_relax_limit_down = 0, temporal_network = F)
  modules_relax_rate$modules$relax_rate <- r
  relaxrate_modules <- rbind(relaxrate_modules, modules_relax_rate$modules)
}




##---- distance decay of modules  ---------------------------------------------------------------
distances_with_ids <- read.csv("./csvs/Islands/distances_with_ids_islands_as_layers.csv")

# this function calculates the Jaccard Similarity in modules between islands (modified for fixed values)
module_distance_decay_islands_func <- function(multilayer_1000, 
                                               layers_turnover_with_distnace){
  k<-seq(0.1,1,by=0.1)
  for (trial in k){
    print(trial) #to keep tab on how far along we are
    modules_for_similarity_shuf <-  multilayer_1000 %>% filter(relax_rate == trial) #take only 1 trial
    
    #pivot modules
    module_pivoted_shuf <- pivot_by_module_islands(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    #create edge list with distances
    modules_edge_list_shuf <- NULL
    
    for (t in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module_islands(module_pivoted_shuf[t,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[t]
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



#pivot modules function for islands as layers
pivot_by_module_islands <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(data, id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s2<-na.omit(s2)
  s3 <-s2%>%
    select(where(~ any(. != 0)))
  s4 = t(s3) 
  s4 <- s4[-1,]
  colnames(s4) <- c(1,2,3,4,5,6,7)
  return(s4)
}#delete all NAs that are in the physical nodes but not the network (and were not assigned to modules)

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

turnover_with_distance_classic <- NULL
layers_turnover_with_distnace<-NULL
module_layer_turnover_shuf <- NULL


all_edge_list_islands_combine_no_module_relax_output<- module_distance_decay_islands_func(relaxrate_modules,
                                                                                           turnover_with_distance_classic)
#---- create ave for jaccard layers-------------------------------------------------------------------------
ave_module_islands_turnover_relax<- all_edge_list_islands_combine_no_module_relax_output %>% 
  group_by(trial,layer_from, layer_to) %>% 
  summarise(ave=mean(turnover), ave_dist=mean(mean_distance)) %>% select("layer_from","layer_to","ave","ave_dist",
                                                                         "trial")
ave_module_islands_turnover_relax$trial<-as.character(ave_module_islands_turnover_relax$trial)

## empirical
islands_turnover_with_distnace_empirical <- read.csv("csvs/Islands/islands_turnover_with_distnace_empirical.csv")

empirical_turnover_for_modules_layers <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), ave_dist = mean_distance) %>% mutate(trial="empirical")


#combine all layers
jaccard_similarity_empirical_and_relax <- rbind(empirical_turnover_for_modules_layers, 
                                                ave_module_islands_turnover_relax)

jaccard_similarity_empirical_and_relax_km <- jaccard_similarity_empirical_and_relax  %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#write.csv(jaccard_similarity_empirical_and_relax_km, "./csvs/Islands/relax_nullmodel4.csv", row.names = FALSE) #so it can be used for classical shuffling

#Plot
jaccard_similarity_empirical_and_relax_km <- read.csv("csvs/Islands/relax_nullmodel4.csv")

pdf('./graphs/Islands/relaxrate_M4.pdf', 10, 6)
jaccard_similarity_empirical_and_relax_km %>% 
  ggplot(aes(x= ave_dist_in_km , y= ave, group= trial, color= trial))+
  geom_point()+ theme_classic()+ geom_smooth(method= "lm", se=F)+
 scale_color_discrete(name = "Relax rate")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank()) 
dev.off()


