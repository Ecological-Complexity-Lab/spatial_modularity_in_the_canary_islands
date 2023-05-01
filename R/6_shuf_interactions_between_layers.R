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

#this portions of code shuffles interactions between layers

#setwd("/Users/maya/Desktop/plant_pollinator_data/dryad_network")
#getwd()
#---- number of layers each species is found in------------------------------------------------------
#plants
dryad_edgelist_complete_ids <- read.csv("./csvs/dryad_edgelist_complete_ids.csv")

plant_half <- dryad_edgelist_complete_ids %>% filter(layer_from == layer_to) %>% #take only intralayers
  filter(node_from <= 39) %>% select(layer_from, node_from) #one directed for this purpose, only plants

#old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
#revised <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)
#plant_half$layer_from[plant_half$layer_from %in% old] <- revised[match(plant_half$layer_from, old)] #change layers to islands

plant_half <- plant_half %>% unique() #remove doubles created due to aggragation
plant_half$count <- 1 
plant_half <- plant_half %>% group_by(node_from) %>% mutate(number_of_layers = sum(count)) #count number of layers 

#pollinators
pollinator_half <- dryad_edgelist_complete_ids %>% filter(layer_from == layer_to) %>% #take only intralayers
  filter(node_to > 39) %>% select(layer_from, node_to) #one directed for this purpose, only plants

#pollinator_half$layer_from[pollinator_half$layer_from %in% old] <- revised[match(pollinator_half$layer_from, old)] #change layers to islands

pollinator_half <- pollinator_half %>% unique() #remove doubles created due to aggregation
pollinator_half$count <- 1 
pollinator_half <- pollinator_half %>% group_by(node_to) %>% mutate(number_of_layers = sum(count)) #count number of layers 

#---- pairs of plant-pollinator------------------------------------------------------------------------

#only plants found in >1 sites
plant_half_for_pairs <- plant_half %>% filter(number_of_layers > 1)

#only pollinators found in >1 sites
pollinator_half_for_pairs <- pollinator_half %>% filter(number_of_layers > 1)

#----number of times species interacted--------------------------------------------------------------
dryad_edgelist_one_sided <- dryad_edgelist_complete_ids %>% filter(node_from <= 39) #make sure we don't have double the results

interaction_count <- dryad_edgelist_one_sided %>% filter(layer_from == layer_to) #only intra as we need just pollination interactions to change

#old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
#revised <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)
#interaction_count$layer_from[interaction_count$layer_from %in% old] <- revised[match(interaction_count$layer_from, old)] #change layers to islands
#interaction_count$layer_to[interaction_count$layer_to %in% old] <- revised[match(interaction_count$layer_to, old)] #change layers to islands

interaction_count <- interaction_count %>% select(layer_from, node_from, layer_to, node_to) %>% unique() #remove doubles created due to aggregation
interaction_count$count <- 1 #add count we can sum later
interaction_count <- aggregate(count ~ node_from + node_to, interaction_count, FUN = sum) #count number of times nodes interacted
interaction_count <- rename(interaction_count, interactions = count) #rename to have distinct name

#----number of times species co-occurred--------------------------------------------------------------
co_occurrence_count <- NULL

co_occurrence <- dryad_edgelist_one_sided %>% filter(layer_from == layer_to) #only intra as we need just pollination interactions to change

#co_occurrence$layer_from[co_occurrence$layer_from %in% old] <- revised[match(co_occurrence$layer_from, old)] #change layers to islands
#co_occurrence$layer_to[co_occurrence$layer_to %in% old] <- revised[match(co_occurrence$layer_to, old)] #change layers to islands

co_occurrence <- co_occurrence %>% select(layer_from, node_from, layer_to, node_to) %>% unique() #remove doubles created due to aggregation

for (i in 1:39){
  current_plant_interactions <- co_occurrence %>% filter(node_from == i) #focal point is plants
  unique_partners <- current_plant_interactions %>% select(node_from, node_to) %>% unique() 
  for (j in 1:nrow(unique_partners)){ #find all unique partners to not have doubles
    current_interaction <- unique_partners[j,] 
    current_plant <- current_interaction$node_from #who's current plant
    current_pollinator <- current_interaction$node_to #who's current pollinator
    
    layers_for_plant <- co_occurrence %>% filter(node_from == current_plant) %>%
      select(layer_from) %>% unique() %>% unlist() #find all islands current plant is found in
    layers_for_pollinator <- co_occurrence %>% filter(node_to == current_pollinator) %>% 
      select(layer_from) %>% unique() %>% unlist() #find all islands current pollinator is found in
    
    intersects <- intersect(layers_for_plant, layers_for_pollinator) #find which islands have both species
    intersect_value <- length(intersects) #count number of islands that have both species
    co_occurrence_count <- rbind(co_occurrence_count, tibble(node_from = current_plant, 
                                                             node_to = current_pollinator, 
                                                             co_occurrences = intersect_value))
  }
}

write.csv(co_occurrence_count, "./csvs/co_occurrence_count.csv", row.names = FALSE)

#---- combine interactions and co-occurrences--------------------------------------------
interactions_co_occurences <- merge(co_occurrence_count, interaction_count, by = c("node_from", "node_to")) #combine two data frames
#co_occurrences_several_islands <- interactions_co_occurences %>% filter(co_occurrences > 1) #only save rows where co-occur

#maybe use sample to randomly select value 
#need to take every pair
#need to sample their interactions number of times equal to number of times their interact
#need to redistribute said interactions
#for (i in 1:1000){ #1000 iterations
#  for (j in 1:nrow(co_occurrences_several_islands)){ #every pair
#    current_pair <- co_occurrences_several_islands[j,] #one pair at a time
    
#    pollinator_layers_pair <- pollinator_half_for_pairs %>% 
#      filter(node_to == current_pair$node_to) %>% select(layer_from) #find layer ids for pols
#    plant_layers_pair <- plant_half_for_pairs %>%
#      filter(node_from == current_pair$node_from) %>% select(layer_from) #find layer ids for plants
    
#    common_layers_pair <- inner_join(pollinator_layers_pair, plant_layers_pair) #nodes and common layers for them
    
#    layers_list <- common_layers_pair$layer_from #list of layers common to pair
    
#    for(k in 1:current_pair$interactions){
#    random_layer <- sample(layers_list, 1, replace = FALSE) #pick random common layer id- pick matrix
    #take interaction from existing data
#    }
#  }
#}


#get common (co-occurrence) layer identities
layers_for_pairs <- NULL

  for (j in 1:nrow(interactions_co_occurences)){ #every pair
    current_pair <- interactions_co_occurences[j,] #one pair at a time
    
    pollinator_layers_pair <- pollinator_half %>% 
      filter(node_to == current_pair$node_to) %>% select(layer_from) #find layer ids for pols
    plant_layers_pair <- plant_half %>%
      filter(node_from == current_pair$node_from) %>% select(layer_from) #find layer ids for plants
    
    common_layers_pair <- inner_join(pollinator_layers_pair, plant_layers_pair) #nodes and common layers for them
    layers_for_pairs <- rbind(layers_for_pairs, common_layers_pair) #create data frame with layers identities
  }

#write.csv(layers_for_pairs, "./csvs/layers_for_pairs.csv", row.names = FALSE)

intralayer_edges_ids <- 
  dryad_intralayer %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)

#intralayer_edges_ids$layer_from[intralayer_edges_ids$layer_from %in% old] <- revised[match(intralayer_edges_ids$layer_from, old)] #change layers to islands
#intralayer_edges_ids$layer_to[intralayer_edges_ids$layer_to %in% old] <- revised[match(intralayer_edges_ids$layer_to, old)] #change layers to islands

#---- shuffle the interactions between layers---------------------------------------------------------

pair_layer <- left_join(layers_for_pairs, interactions_co_occurences, by = c("node_from", "node_to")) #join layer ids with the pair 

pairs_with_weights <- left_join(pair_layer, intralayer_edges_ids) %>%
  mutate(weight = replace(weight, is.na(weight), 0)) %>%
  select(node_from, node_to, layer_from, co_occurrences, interactions, weight) #join pairs and layers with weights


interactions_shuf <- NULL

for (i in 1:1000){ #1000 iterations
  print(i)
  for(j in 1:nrow(interactions_co_occurences)){
    current_node_from <- interactions_co_occurences[j,]$node_from
    current_node_to <- interactions_co_occurences[j,]$node_to
    current_pair_data <- pairs_with_weights %>% 
      filter(node_from == current_node_from, node_to == current_node_to) #take all layers a pair is found in
    if (nrow(current_pair_data) == 1){
      interactions_shuf <- rbind(interactions_shuf, 
                                 tibble(current_pair_data, trial_num = i)) 
    } #if the pair is only found in 1 layer just add it and continue
    else{ #if not in just 1 layer also shuffle before adding
      shuffling_interactions <- transform(current_pair_data, layer_from = sample(layer_from)) #shuffle the layer ids
      interactions_shuf <- rbind(interactions_shuf, 
                                 tibble(shuffling_interactions, trial_num = i)) #create new data frame where interactions are shuffled
    } 
    }
}


write.csv(interactions_shuf, "./csvs/interactions_shuf.csv", row.names = FALSE)
#interactions_shuf <- read.csv("./csvs/interactions_shuf.csv")
interactions_shuf$layer_to <- interactions_shuf$layer_from #add layer_to which is identical to layer_from as we're only dealing with intralayers
interactions_shuf_intralayers <- interactions_shuf %>% select(layer_from, node_from, layer_to, node_to, weight, trial_num) #organize
interactions_shuf_intralayers <- interactions_shuf_intralayers %>% filter(weight != 0) #delete all 0 weights

#---- create inverted intralayer version-------------------------------------------------------------------
interactions_shuf_intralayers_inverted <- tibble(values= interactions_shuf_intralayers$layer_to, interactions_shuf_intralayers$node_to, 
                                                 interactions_shuf_intralayers$layer_from, interactions_shuf_intralayers$node_from, 
                                                 interactions_shuf_intralayers$weight, interactions_shuf_intralayers$trial_num) #create an inverted copy for directed intralayers
colnames(interactions_shuf_intralayers_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_num")


#---- normalize intralayer weights--------------------------------------------------------------------------
#plants in from
tot_plant_shuf <- interactions_shuf_intralayers %>% 
  group_by(layer_from,node_from, trial_num) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf <- interactions_shuf_intralayers %>% left_join(tot_plant_shuf) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)

intralayer_weighted_shuf <- intralayer_weighted_shuf[, c(1,2,3,4,6,5)] #make sure weight is col number 5

#pols in from
tot_pol_shuf <- interactions_shuf_intralayers_inverted %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf <- interactions_shuf_intralayers_inverted %>% left_join(tot_pol_shuf) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight) 

intralayer_weighted_inverted_shuf <- intralayer_weighted_inverted_shuf[, c(1,2,3,4,6,5)] #make sure weight is col number 5

intralayer_edges_shuf <- rbind(intralayer_weighted_shuf, intralayer_weighted_inverted_shuf) #combine inverted and non-inverted versions

write.csv(intralayer_edges_shuf, "./csvs/intralayer_edges_shuf.csv") #create file to run interlayers on HPC
#intralayer_edges_shuf <- read.csv("./csvs/intralayer_edges_shuf.csv")
#---- interlayer edges--------------------------------------------------------------------------------------
interlayer_interactions_from_shuf <- intralayer_edges_shuf %>% group_by(trial_num, node_from) %>%
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3], 
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in


interlayer_interactions_to_shuf <- intralayer_edges_shuf %>% group_by(trial_num, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind

interlayer_edges_interactions_shuf <- rbind(interlayer_interactions_from_shuf, interlayer_interactions_to_shuf) 

write.csv(interlayer_edges_interactions_shuf, "./HPC/shuf_interactions/interlayer_edges_interactions_shuf.csv")
interlayer_edges_interactions_shuf <- read.csv("./HPC/shuf_interactions/interlayer_edges_interactions_shuf.csv")
#---- run on HPC and come back
#there are 3 code portions for the HPC:
# 1. HPC_network_interactions_shuffle.R 
# 2. 1_1000_interactions.sh 
# 3. i_interactions.sh 
# running 1_1000_interactions.sh manually in the cmd will make the other two run.
# all 1000 result csvs can be found in HPC/csvs_interactions files

files_interactions <- list.files("./HPC/shuf_interactions/csvs_interactions/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_interactions <- read_csv(files_interactions) %>% bind_rows() #create a long edge list with all the csvs

#write.csv(my_merged_interlayer_shuf_interactions, "./csvs/my_merged_interlayer_shuf_interactions.csv", row.names = FALSE)

#---- interlayers with weights shuf version ------------------------------------------
distances_with_weights_ids <- read.csv("./csvs/distances_with_weights_ids.csv") #need to read this to run next part


interlayers_with_weights_shuf_interactions <- my_merged_interlayer_shuf_interactions %>% inner_join(distances_with_weights_ids, 
                                                                                   by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights_shuf_interactions <- 
  interlayers_with_weights_shuf_interactions[!duplicated(interlayers_with_weights_shuf_interactions[c(1,3,5,6)]),]

#write.csv(interlayers_with_weights_shuf_interactions, "./csvs/interlayers_with_weights_shuf_interactions.csv", row.names = FALSE)

## ----multilayer_extended_final--------------------------------------------------------------------------------------

names(intralayer_edges_shuf)[6] <- "trial_number" #change name of column so inter and intra correspond

dryad_edgelist_complete_shuf_interactions <- bind_rows(intralayer_edges_shuf, 
                                                       interlayers_with_weights_shuf_interactions) #combine inter and intra

## ----multilayer_class-----------------------------------------------------------------------------------------------
physical_nodes <- read.csv("./csvs/physical_nodes.csv")
layer_metadata <- read.csv("./csvs/layer_metadata.csv")

dryad_multilayer_shuf_1000_interactions <- NULL


dryad_multilayer_shuf_1000_interactions_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_interactions, 
                                                                      dryad_multilayer_shuf_1000_interactions)

write.csv(dryad_multilayer_shuf_1000_interactions_output, "./csvs/dryad_multilayer_shuf_1000_interactions_output.csv", row.names = FALSE)

##---- distance decay of modules in islands shuf vs empirical----------------------------------------------------
all_edge_list_island_combine_no_module_shuf_interations <- NULL
islands_turnover_with_distnace_interactions <- NULL
module_island_turnover_interactions_shuf <- NULL
islands_turnover_with_distnaces_interactions <- NULL


#---- outputs islands----------------------------------------------------------------------------
#distances_with_ids <- read.csv("./csvs/distances_with_ids.csv")

all_edge_list_island_combine_no_module_shuf_interactions_output <- module_distance_decay_func(dryad_multilayer_shuf_1000_interactions_output,
                                                                                      islands_turnover_with_distnace_interactions,
                                                                                      islands_turnover_with_distnaces_interactions)


write.csv(all_edge_list_island_combine_no_module_shuf_interactions_output, 
          "./csvs/all_edge_list_island_combine_no_module_shuf_interactions_output.csv", 
          row.names = FALSE)

#all_edge_list_island_combine_no_module_shuf_interactions_output <- 
#read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_interactions_output.csv")

#---- create ave for jaccard islands -------------------------------------------------------------------------
ave_module_island_turnover_shuf_interactions <- all_edge_list_island_combine_no_module_shuf_interactions_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_model") #create mean and sd for each point

#add the empirical 
islands_turnover_with_distnace_empirical <- read.csv("./csvs/islands_turnover_with_distnace_empirical.csv")

empirical_turnover_for_module_island_shuf <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

empirical_turnover_for_module_island_interactions_shuf_no_self_loop <- empirical_turnover_for_module_island_shuf %>% 
  subset(layer_from != layer_to) #for empirical only distance decay graph

empirical_turnover_for_module_island_interactions_shuf_no_self_loop_km <- empirical_turnover_for_module_island_interactions_shuf_no_self_loop %>%
  mutate(ave_dist_in_km = ave_dist/1000)

write.csv(empirical_turnover_for_module_island_interactions_shuf_no_self_loop,
          "./csvs/empirical_turnover_for_module_island_interactions_shuf_no_self_loop.csv", row.names = FALSE)

#---- combine for jaccard analysis------------------------------------------------------------------------------------------------------
#combine all islands
jaccard_similarity_empirical_and_null_interactions <- rbind(empirical_turnover_for_module_island_shuf, 
                                                            ave_module_island_turnover_shuf_interactions)

jaccard_similarity_empirical_and_null_interactions_no_self_loop <- jaccard_similarity_empirical_and_null_interactions %>% 
  subset(layer_from != layer_to)

jaccard_similarity_empirical_and_null_interactions_no_self_loop_km <- jaccard_similarity_empirical_and_null_interactions_no_self_loop %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

#---- graphs for distance decay in modules shuf vs empirical--------------------------------
#island
#empirical and null
jaccard_similarity_empirical_and_null_interactions_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400)+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.63, 0.60, 0.57, 0.54))+
  scale_color_manual(values = c("#F47069", "#0033cc"))

#version without trendline
jaccard_similarity_empirical_and_null_interactions_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ 
  geom_smooth(data = empirical_turnover_for_module_island_shuf_no_self_loop_km, method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+
  scale_color_manual(values = c("#F47069", "#0033cc"))

#---- statistical analysis--------------------------------------------------
#---- linear regression
lm1_module_interactions = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_interactions_no_self_loop,
                                                         jaccard_similarity_empirical_and_null_interactions_no_self_loop$type=="empirical")) #in empirical
lm2_module_interactions = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_interactions_no_self_loop,
                                                         jaccard_similarity_empirical_and_null_interactions_no_self_loop$type=="null_model")) #in null model


#get equations
lm1_module_equation_interactions <- paste("y=", coef(lm1_module_interactions)[[1]], "+", coef(lm1_module_interactions)[[2]], "*x")
lm2_module_equation_interactions <- paste("y=", coef(lm2_module_interactions)[[1]], "+", coef(lm2_module_interactions)[[2]], "*x")
lm1_module_equation_interactions
lm2_module_equation_interactions

b1_module_interactions <- summary(lm1_module_interactions)$coefficients[2,1]
se1_module_interactions <- summary(lm1_module_interactions)$coefficients[2,2]
b2_module_interactions <- summary(lm2_module_interactions)$coefficients[2,1]
se2_module_interactions <- summary(lm2_module_interactions)$coefficients[2,2]


p_value_module_interactions = 2*pnorm(-abs(compare.coeff(b1_module_interactions,se1_module_interactions,
                                                         b2_module_interactions,se2_module_interactions)))
p_value_module_interactions

#---- linear regression and correlation

iteration_correlation_interactions <- NULL
iteration_correlation_data_interactions <- all_edge_list_island_combine_no_module_shuf_interactions_output %>% subset(layer_from != layer_to) 

for (i in 1:1000){
  trial_interactions = iteration_correlation_data_interactions %>% filter(trial_num == i)
  iteration_correlation_new_interactions <- cor.test(trial_interactions$turnover, trial_interactions$ave_distance, method = "pearson")
  lm_val_interactions <- lm(turnover ~ ave_distance, data = trial_interactions)
  iteration_correlation_interactions <- rbind(iteration_correlation_interactions, tibble(estimate = iteration_correlation_new_interactions$estimate, 
                                                                               p_val = iteration_correlation_new_interactions$p.value, 
                                                                               statistic = iteration_correlation_new_interactions$statistic, 
                                                                               confidence_int_low = iteration_correlation_new_interactions$conf.int[1],
                                                                               confidence_int_high = iteration_correlation_new_interactions$conf.int[2],
                                                                               slope = lm_val_interactions$coefficients[2],
                                                                               intercept = lm_val_interactions$coefficients[1],
                                                                               rsquared = summary(lm_val_interactions)$adj.r.squared,
                                                                               trial_num = i))
}

write.csv(iteration_correlation_interactions, "./csvs/iteration_correlation_interactions.csv", row.names = FALSE)
#correlation_empirical_interactions <- read.csv("./csvs/correlation_empirical_pols.csv")

##distribution of rsquared and add empirical
iteration_correlation_interactions %>% ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#0033cc", color = "#0033cc", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_interactions$rsquared, linetype = "dashed", color = "#F47069") 

p_rsquared_interactions <- sum(iteration_correlation_interactions$rsquared > correlation_empirical_interactions$rsquared)/1000
p_rsquared_interactions

##compare both parts of beta diversity
iteration_correlation_pols <- read.csv("./csvs/iteration_correlation_pols.csv")

interactions_for_beta_diversity <- iteration_correlation_interactions
interactions_for_beta_diversity$type <- "shuffling_interactions" 

species_for_beta_diversity <- iteration_correlation_pols
species_for_beta_diversity$type <- "shuffling_species"

correlation_beta_diversity <- rbind(interactions_for_beta_diversity, species_for_beta_diversity)

correlation_beta_diversity %>% ggplot(aes(x = rsquared, color = type, fill = type))+ #overlay shuffling of interactions and species
  geom_density(alpha = 0.4)+ scale_color_manual(values = c("#0033cc", "#BE75FA"))+ 
  scale_fill_manual(values = c("#0033cc", "#BE75FA"))+
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_interactions$rsquared, linetype = "dashed", color = "#F47069") 


