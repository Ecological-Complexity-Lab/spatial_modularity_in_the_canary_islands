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
library(taxize)
library(taxizedb)

# This portion of code comes to prove that the interlayer
#shuffling that breaks the distance decay in species is
#the reason for the break in distance decay in modules.

#---- shuffle within layers---------------------------------------------------------
# all the csvs needed to run first portion of code:
dryad_intralayer <- read.csv("./csvs/intralayer_file.csv")
physical_nodes <- read.csv("./csvs/physical_nodes.csv")
layer_metadata <- read.csv("./csvs/layer_metadata.csv")


shuf_trial_matrix_classic <- NULL
shuf_null_edge_list_classic <- NULL

intralayers_with_ids_for_shuf <- 
  dryad_intralayer %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)


for (i in 1:14){
  current_matrix <- intralayers_with_ids_for_shuf %>%
    filter(layer_from == layer_to) %>% filter(layer_from == i) %>% #take 1 layer at a time
    select(node_from, node_to, weight) %>% #select interactions
    dcast(node_from ~ node_to, value.var = "weight", fill = 0) %>%
    column_to_rownames(var="node_from") 
  
  current_matrix <- t(current_matrix) #put pols in rows
  
  print(i) #to keep tab on which layer we're on
  
  for (j in 1:1000){ #1000 iterations
    null <- vegan::nullmodel(current_matrix, method = 'r00_samp') #shuffle within layer
    shuffled_matrices <- simulate(null, nsim = 1)
    
    trial_with_shuf_num <- cbind(shuffled_matrices, j) #add trial number 
    shuf_trial_matrix_classic <- rbind(shuf_trial_matrix_classic, trial_with_shuf_num) #create big matrix of all matrices
    edge_list_version_classic <- melt(as.matrix(shuffled_matrices)) %>% filter(value > 0) %>%
      select(node_from=Var1, node_to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
    edge_list_version_classic$trial_number <- j #add trial number to the edge list
    edge_list_version_classic$layer_from <- i #intra so layer from and to are identical
    edge_list_version_classic$layer_to <- i
    shuf_null_edge_list_classic <- rbind(shuf_null_edge_list_classic, edge_list_version_classic) #create mega edge list with all repetitions
  }
}

shuf_null_edge_list_classic <- shuf_null_edge_list_classic[, c(5,1,6,2,3,4)] #change to regular order

#write.csv(shuf_null_edge_list_classic, "./csvs/shuf_null_edge_list_classic.csv", row.names = FALSE)

#interlayer edges
interlayer_edges_from_shuf_classic <- shuf_null_edge_list_classic %>% group_by(trial_number, node_from) %>%
  select(trial_number, layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3], 
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7], loc8 = layer_from[8], loc9 = layer_from[9], 
         loc10 = layer_from[10], loc11 = layer_from[11], loc12 = layer_from[12],
         loc13 = layer_from[13], loc14 = layer_from[14]) #all layers the species is found in


interlayer_edges_to_shuf_classic <- shuf_null_edge_list_classic %>% group_by(trial_number, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7], loc8 = layer_to[8], loc9 = layer_to[9], 
         loc10 = layer_to[10], loc11 = layer_to[11], loc12 = layer_to[12],
         loc13 = layer_to[13], loc14 = layer_to[14]) %>% #all layers the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind


interlayer_edges_shuf_classic <- rbind(interlayer_edges_from_shuf_classic, interlayer_edges_to_shuf_classic) 

#write.csv(interlayer_edges_shuf_classic, "./csvs/interlayer_edges_shuf_classic.csv",row.names = FALSE)

#---- run oh HPC and get results for analysis -----------------------------------------------------------------
write.csv(interlayer_edges_shuf_classic, "./HPC/shuf_within_layers/interlayer_edges_shuf_classic.csv", row.names = FALSE) #create to run on HPC

#run on HPC and then come back with results

files_classic <- list.files("./HPC/shuf_within_layers/csvs_classic/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_classic <- read_csv(files_classic) %>% bind_rows() #create a long edge list with all the csvs

#write.csv(my_merged_interlayer_shuf_classic, "./csvs/my_merged_interlayer_shuf_classic.csv", row.names = FALSE) #create to run on HPC

interlayers_with_weights_shuf_classic <- my_merged_interlayer_shuf_classic %>% inner_join(distances_with_weights_ids, 
                                                                                          by = c("layer_from", "layer_to")) %>% unique()

interlayers_with_weights_shuf_classic <- interlayers_with_weights_shuf_classic[!duplicated(interlayers_with_weights_shuf_classic[c(1,3,5,6)]),]

#write.csv(interlayers_with_weights_shuf_classic, "./csvs/interlayers_with_weights_shuf_classic.csv", row.names = FALSE)

#create inter and intra for the 1000 shuf trials
dryad_interlayer_shuf_classic <- read.csv("./csvs/interlayers_with_weights_shuf_classic.csv") #already has inverted 


dryad_intralayer_shuf_classic <- read.csv("csvs/shuf_null_edge_list_classic.csv") 
dryad_intralayer_shuf_classic <- dryad_intralayer_shuf_classic[, c(6,1,2,3,4,5)]

#create inverted versions
intralayer_inverted_shuf_classic <- tibble(values= dryad_intralayer_shuf_classic$layer_to, dryad_intralayer_shuf_classic$node_to, 
                                           dryad_intralayer_shuf_classic$layer_from, dryad_intralayer_shuf_classic$node_from, 
                                           dryad_intralayer_shuf_classic$weight, dryad_intralayer_shuf_classic$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_classic) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")

# ----normalize weight--------------------------------------------------------------------------
#plants in from
tot_plant_shuf_classic <- intralayer_inverted_shuf_classic %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_classic <- intralayer_inverted_shuf_classic %>% left_join(tot_plant_shuf_classic) %>% 
  mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_classic <- tot_plant_shuf_classic[, c(3,1,2,4)]

#pols in from
tot_pol_shuf_classic <- dryad_intralayer_shuf_classic %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_classic <- dryad_intralayer_shuf_classic %>% left_join(tot_pol_shuf_classic) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_classic <- tot_pol_shuf_classic[, c(3,1,2,4)]

# ----multilayer_extended_final--------------------------------------------------------------------------------------
edgelist_intralayer_shuf_classic <- bind_rows(intralayer_weighted_shuf_classic, intralayer_weighted_inverted_shuf_classic)

dryad_edgelist_complete_shuf_classic <- bind_rows(edgelist_intralayer_shuf_classic, dryad_interlayer_shuf_classic) #combine inter and intra

# ----multilayer_class-----------------------------------------------------------------------------------------------
# Input: An extended edge list.
dryad_edgelist_complete_shuf_classic <- dryad_edgelist_complete_shuf_classic[, c(1,2,3,4,6,5)] #make sure weight is #5
dryad_edgelist_complete_shuf_classic <- as.data.frame(dryad_edgelist_complete_shuf_classic) 

dryad_multilayer_shuf_1000_classic <- NULL

modularity_for_shuf_classic <- function(edge_list, output){
  for(trial_num in 1:1000){
    print(trial_num) #keep tab on which trial is running

    current_trial_edgelist <- edge_list %>% filter(trial_number == trial_num) #take 1 trial at a time to create multilayer
    
    #check which species were deleted
    which_missing_classic_shuf <- anti_join(dryad_edgelist_complete_ids, current_trial_edgelist) #interactions found in the empirical network but not in the trial because of the shuffle
    
    candidates_to <- which_missing_classic_shuf %>% select(node_to) %>% distinct() #take all node_to that might not be in the trial
    candidates_from <- which_missing_classic_shuf %>% select(node_from) %>% distinct() #take all node_from that might not be in the trial
    
    missing_to_in_shuf <- subset(candidates_to, !(node_to %in% current_trial_edgelist$node_to)) #check which ones are not found in the trial
    missing_from_in_shuf <- subset(candidates_from, !(node_from %in% current_trial_edgelist$node_from)) #check which ones are not found in the trial
    
    physical_nodes_del_classic <- subset(physical_nodes, node_id %in% c(missing_from_in_shuf$node_from, missing_to_in_shuf$node_to)) #the unwanted values
    
    
    candidate_partners <- subset(which_missing_classic_shuf, node_from %in% physical_nodes_del_classic$node_id) #all candidate partners
    candidate_partners <- subset(candidate_partners, !(node_from == node_to)) %>% select(node_to) %>% distinct() #only stay with candidate
    eliminate_partners_only_of_del <- test %>% filter(node_to %in% candidate_partners$node_to) %>% subset(!(node_from == node_to)) #take candidate partners and delete all occations of interlayers
    who_to_eliminate <- candidate_partners %>% subset(!(node_to %in% eliminate_partners_only_of_del$node_from)) %>% select(node_to)
    
    physical_nodes_del_classic_partners <- anti_join(physical_nodes_del_classic, who_to_eliminate) #only save the species who weren't interacting only with them
    #need to delete all interactions the deleted species were a part of
    
    dryad_multilayer_shuf_trial <- create_multilayer_object(extended = current_trial_edgelist, #taking edge list and returning multilayer network
                                                            nodes = physical_nodes_del_classic_partners,
                                                            layers = layer_metadata,
                                                            intra_output_extended = T)
    
    
    #create modules for empirical network
    modules_dryad_multilayer_shuf_1000 <- run_infomap_multilayer(dryad_multilayer_shuf_trial, 
                                                                infomap_executable = "../Infomap",
                                                                flow_model = 'directed',
                                                                relax = F, 
                                                                silent = T, 
                                                                trials = 100,
                                                                seed = 497294, 
                                                                temporal_network = F)
    
    
    
    output <- rbind(output, tibble(modules_dryad_multilayer_shuf_1000$modules, trial_num)) #save modules for 1000 null models
  }  
  return(output)
}

dryad_multilayer_shuf_1000_classic_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_classic, 
                                                                 dryad_multilayer_shuf_1000_classic)

dryad_multilayer_shuf_1000_classic_output <- dryad_multilayer_shuf_1000_classic_output %>% drop_na() 
#delete all NAs that are in the physical nodes but not the network (and were not assigned to modules)

#write.csv(dryad_multilayer_shuf_1000_classic_output, "./csvs/dryad_multilayer_shuf_1000_classic_output.csv", row.names = FALSE)

#---- distance decay of modules shuf vs empirical----------------------
islands_turnover_with_distnace_classic <- NULL

module_island_turnover_shuf <- NULL

all_edge_list_island_combine_no_module_shuf_classic <- module_distance_decay_func(dryad_multilayer_shuf_1000_classic_output,
                                                                                  islands_turnover_with_distnace_classic)

write.csv(all_edge_list_island_combine_no_module_shuf_classic, "./csvs/all_edge_list_island_combine_no_module_shuf_classic.csv", 
          row.names = FALSE)

#all_edge_list_island_combine_no_module_shuf_classic <- read.csv("./csvs/all_edge_list_island_combine_no_module_shuf_classic.csv")

ave_module_island_turnover_shuf_classic <- all_edge_list_island_combine_no_module_shuf_classic %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_model_within") #create mean and sd for each point

#add the empirical empirical
empirical_turnover_for_module_island_shuf <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

empirical_turnover_for_module_island_shuf_no_self_loop_classic <- empirical_turnover_for_module_island_shuf %>% 
  subset(layer_from != layer_to) #distance decay graph for empirical 

empirical_turnover_for_module_island_shuf_no_self_loop_km_classic <- empirical_turnover_for_module_island_shuf_no_self_loop_classic %>%
  mutate(ave_dist_in_km = ave_dist/1000)

##---- combine empirical and shuffle 

jaccard_similarity_empirical_and_null_model <- rbind(empirical_turnover_for_module_island_shuf, ave_module_island_turnover_shuf_classic)

jaccard_similarity_empirical_and_null_no_self_loop_classic <- jaccard_similarity_empirical_and_null_model %>% subset(layer_from != layer_to)

jaccard_similarity_empirical_and_null_no_self_loop_km_classic <- jaccard_similarity_empirical_and_null_no_self_loop_classic %>% 
  mutate(ave_dist_in_km = ave_dist/1000)

just_empirical_within <- jaccard_similarity_empirical_and_null_no_self_loop_km_classic %>% filter(type == "empirical")

#no self loop
#empirical and null
jaccard_similarity_empirical_and_null_no_self_loop_km_classic %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400)+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.57, 0.54))+ scale_color_manual(values = c("#F47069", "#c4067c"))

#no trendline
jaccard_similarity_empirical_and_null_no_self_loop_km_classic %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ 
  geom_smooth(data = just_empirical_within, method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor( data = just_empirical_within, aes(label = ..p.label..), label.x = 400)+
  stat_cor(data = just_empirical_within, aes(label = ..rr.label..), label.x = 400, label.y = 0.6)+ scale_color_manual(values = c("#F47069", "#c4067c"))

#---- statistical analysis--------------------------------------------------
#---- linear regression
lm1_module_classic = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop_classic,
                                                    jaccard_similarity_empirical_and_null_no_self_loop_classic$type=="Empirical")) #in empirical
lm2_module_classic = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_null_no_self_loop_classic,
                                                    jaccard_similarity_empirical_and_null_no_self_loop_classic$type=="Null_Model")) #in null model


#get equations
lm1_module_equation_classic <- paste("y=", coef(lm1_module_classic)[[1]], "+", coef(lm1_module_classic)[[2]], "*x")
lm2_module_equation_classic <- paste("y=", coef(lm2_module_classic)[[1]], "+", coef(lm2_module_classic)[[2]], "*x")
lm1_module_equation_classic
lm2_module_equation_classic

b1_module_classic <- summary(lm1_module_classic)$coefficients[2,1]
se1_module_classic <- summary(lm1_module_classic)$coefficients[2,2]
b2_module_classic <- summary(lm2_module_classic)$coefficients[2,1]
se2_module_classic <- summary(lm2_module_classic)$coefficients[2,2]


p_value_module_classic = 2*pnorm(-abs(compare.coeff(b1_module_classic,se1_module_classic,b2_module_classic,se2_module_classic)))
p_value_module_classic


#---- linear regression and correlation

iteration_correlation_classic <- NULL
iteration_correlation_data_classic <- all_edge_list_island_combine_no_module_shuf_classic %>% subset(layer_from != layer_to) 

for (i in 1:1000){
  trial_classic = iteration_correlation_data_classic %>% filter(trial == i)
  iteration_correlation_new_classic <- cor.test(trial_classic$turnover, trial_classic$ave_distance, method = "pearson")
  lm_val_classic <- lm(turnover ~ ave_distance, data = trial_classic)
  iteration_correlation_classic <- rbind(iteration_correlation_classic, tibble(estimate = iteration_correlation_new_classic$estimate, 
                                                                         p_val = iteration_correlation_new_classic$p.value, 
                                                                         statistic = iteration_correlation_new_classic$statistic, 
                                                                         confidence_int_low = iteration_correlation_new_classic$conf.int[1],
                                                                         confidence_int_high = iteration_correlation_new_classic$conf.int[2],
                                                                         slope = lm_val_classic$coefficients[2],
                                                                         intercept = lm_val_classic$coefficients[1],
                                                                         rsquared = summary(lm_val_classic)$adj.r.squared,
                                                                         trial_num = i))
}

correlation_empirical_classic <- read.csv("./csvs/correlation_empirical_pols.csv")

##distribution of rsquared and add empirical
iteration_correlation_classic %>% ggplot(aes(x = rsquared))+ 
  geom_density(fill = "#c4067c", color = "#c4067c", alpha = 0.4)+ 
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_classic$rsquared, linetype = "dashed", color = "#F47069") 

p_rsquared_classic <- sum(iteration_correlation_classic$rsquared > correlation_empirical_classic$rsquared)/1000
p_rsquared_classic

##distribution of slope and add empirical
iteration_correlation_classic %>% ggplot(aes(x = slope))+ 
  geom_density(fill = "#c4067c", color = "#c4067c", alpha = 0.4)+ 
  theme_classic()+ labs(x = "slope")+
  geom_vline(xintercept = correlation_empirical_classic$slope, linetype = "dashed", color = "#F47069") 

p_slope_classic <- sum(iteration_correlation_classic$slope < correlation_empirical_classic$slope)/1000
p_slope_classic
