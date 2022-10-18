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

#null model with fixed interlayer edges-----------------------------------------------------------------------------------------------
#compare species composition, number of modules and which species are found in which modules

#I_or <- NULL #I between empirical and fixed
#I_or_layer = NULL #I between empirical and fixed for each layer
#num_of_modules <- NULL
#sanity_check_N <-NULL
#M_A <- modules_dryad_multilayer$modules %>% slice(which.max(module)) %>% select(module) #how many modules are in the empirical network
#num_of_modules <- rbind(num_of_modules, tibble(fixed_interlayer_value=NA ,number_of_modules=M_A))
#partner_comparison <- NULL

#for (val in seq(0,1,0.05)){ #create loop changing the module_edge_change from 0 to 1 using 0.1 jumps 
#  print(val)
#  partners_null_model <- NULL
#  interlayer_edges_change <- select(inter_extended, -weight) #create data frame where weight doesn't exist
#  weight <- data.frame() #create data frame of just weight
#  for (i in 1:length(inter_extended$layer_from)){
#    weight <- rbind(weight, val) #add fixed weight 
#  }
#  interlayer_edges_change <- cbind(interlayer_edges_change, weight) #bind the nodes and layers with the new fixed weight
#  colnames(interlayer_edges_change)[5] <- "weight"

#  dryad_multilayer_null <- create_multilayer_object(intra = intra_nonextended, 
#                                                 inter = interlayer_edges_change, #create multilayer with new fixed inter value
#                                                 nodes = physical_nodes,
#                                                 layers = layer_metadata,
#                                                 intra_output_extended = T) 

#  modules_edge_change <- run_infomap_multilayer(dryad_multilayer_null, 
#                                             infomap_executable = "../Infomap",
#                                             flow_model = 'directed',
#                                             relax = F, 
#                                             silent = T, 
#                                             trials = 100,
#                                             seed = 497294, 
#                                             temporal_network = F,remove_auxilary_files = F)

#compare mutual information between empirical and fixed
#A <- modules_dryad_multilayer$modules %>% select(module, node_id, layer_id)  #A should always stay the same
#B <- modules_edge_change$modules %>% select(module, node_id, layer_id) #B should change each round

#N <- inner_join(A,B,by=c('node_id','layer_id')) %>% 
#group_by(module.y) %>% 
#select(module.x) %>% table()
#I_or <- rbind(I_or, tibble(fixed_interlayer_value=val, mutual_information=NMI(N), layer=NA))

#sanity_check_N <- rbind(sanity_check_N, tibble(N=N, fixed_interlayer_value=val))

#for (i in (1:max(modules_dryad_multilayer$modules$layer_id))){ #do action for each layer on its own
#  N_layer <- inner_join(A,B,by=c('node_id','layer_id')) %>% filter(layer_id==i) %>% #each round take only one layer
#  group_by(module.y) %>%
#  select(module.x) %>% table()
#  I_or_layer <- rbind(I_or_layer, tibble(fixed_interlayer_value=val, mutual_information=NMI(N_layer), layer=i))
#}

#I_or_layer <- rbind(I_or_layer, I_or) #combine MI of all layers and each layer on its own

#compare number of modules between empirical and fixed
#M_B <- modules_edge_change$modules %>% slice(which.max(module)) %>% select(module) #how many modules are in the current fixed edge module
#num_of_modules <- rbind(num_of_modules, tibble(fixed_interlayer_value=val, number_of_modules=M_B))

#compare partners in modules
#for (j in 1:length(modules_edge_change$modules$node_id)){
#  current_node_nm <- modules_edge_change$modules$node_id[j] #current node I'm looking at as a view point
#  filter_by_node_nm <- filter(modules_edge_change$modules, modules_edge_change$modules$node_id==current_node_nm) #only rows where node is
#  modules_with_species_nm <- unique(filter_by_node_nm$module) #only take unique module numbers
#  module_partners_nm <- filter(modules_edge_change$modules, modules_edge_change$modules$module==modules_with_species_nm)
#  just_partners_nm <- subset(module_partners_nm, node_id != current_node_nm) #remove the current id from the module
#  unique_partners_nm <- unique(just_partners_nm$node_id) #unique partners in a specific module 
#  partners_null_model <- rbind(partners_null_model, tibble(node_id=current_node_nm, partners=unique_partners_nm)) #all nodes and all their partners
#  partners_null_model <- distinct(partners_null_model)
#}

#partners_null_model <- partners_null_model %>%
#  merge(physical_nodes, by="node_id") %>%
#  subset(select= -species)

#write_csv(partners_null_model, "./csvs/module_partners_with_fixed_interlayers.csv")

#jaccard_val <- NULL #restart every run

#for (i in 1:max(partners_empirical$node_id)){ #for all the nodes
#  filter_empirical <- filter(partners_empirical, partners_empirical$node_id==i) #filter only by node number i in empirical
#  filter_null_model <- filter(partners_null_model, partners_null_model$node_id==i) #filter only by node number i in null model
#  jaccard_similarity <- jaccard(filter_empirical$partners, filter_null_model$partners) #jaccard on node i in current val
#  jaccard_val <- rbind(jaccard_val, tibble(node=i, similarity=jaccard_similarity, type=filter_empirical$type)) #contains node and similarity
#}

#partner_comparison <- rbind(partner_comparison, tibble(fixed_interlayer_value=val, similarity=jaccard_val))
#partner_comparison <- distinct(partner_comparison) #only 1 occurence of each row

#}

#create graph for species composition
#ggplot(I_or_layer, aes(x=fixed_interlayer_value, y=mutual_information, group= layer, colour= as.factor(layer)))+geom_point()+geom_line()+
#  scale_x_continuous(breaks=seq(0,1,0.1))+ylim(0,1)+theme_classic()+
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
#  theme(axis.text.y=element_text(size=15))+
#  labs(x="fixed interlayer value", y="mutual information")+
#  theme(plot.title=element_text(hjust=0.5)) + guides(fill=guide_legend(title="layer"))

#again but just for entire network
#ggplot(I_or, aes(x=fixed_interlayer_value, y=mutual_information))+geom_point()+geom_line()+
#  scale_x_continuous(breaks=seq(0,1,0.1))+ylim(0,1)+theme_classic()+
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
#  theme(axis.text.y=element_text(size=15))+
#  labs(x="fixed interlayer value", y="mutual information")+
#  theme(plot.title=element_text(hjust=0.5)) + guides(fill=guide_legend(title="layer"))

#create graph for number of modules
#num_of_modules$number_of_modules <- unlist(num_of_modules$number_of_modules)
#threshold <- num_of_modules$number_of_modules

#num_of_modules %>%
#  ggplot(aes(x=fixed_interlayer_value, y=number_of_modules))+
#  geom_histogram(stat="identity")+
#  geom_hline(aes(yintercept=threshold[1]), color="red")+
#  scale_x_continuous(breaks=seq(0,1,0.1))+theme_classic()+
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
#  theme(axis.text.y=element_text(size=15))+
#  geom_text(x=0.2, y=threshold[1]+2, label="empirical network", color= "red", size=5, family= "Tahoma")+
#  labs(x="fixed interlayer value (null model)", y="number of modules")


#create graph for partner similarity
#partner_comparison$fixed_interlayer_value <- as.factor(partner_comparison$fixed_interlayer_value)

#partner_comparison %>%
#  ggplot(aes(x=fixed_interlayer_value, y=similarity$similarity, fill=similarity$type))+geom_boxplot()+
#  scale_x_discrete(breaks=seq(0,1,0.1))+theme_classic()+
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
#  theme(axis.text.y=element_text(size=15))+
#  labs(x="fixed interlayer value", y="similarity", fill="species type")


#write.csv(I_or, "./csvs/mutual_information_fixed_vs_empirical.csv")
#write.csv(num_of_modules, "./csvs/umber_of_modules_fixed_vs_empirical.csv")
#write.csv(partner_comparison, "./csvs/partner_similarity_fixed_vs_empirical.csv")
#write.csv(modules_dryad_multilayer$modules, "./csvs/modules_empirical.csv")
#write.csv(modules_edge_change$modules, "./csvs/modules_interlayer_fixed.csv")

#shuffling interlayer edges--------------------------------------------------------------------------------------------------------
#Num_of_modules <- NULL
#M_C <- modules_dryad_multilayer$modules %>% slice(which.max(module)) %>% select(module) #how many modules are in the empirical network
#Num_of_modules <- rbind(Num_of_modules, tibble(shuffle_trial=NA ,number_of_modules=M_C))
#I_es <- NULL #I between empirical and shuffled
#partner_comparison_shuf <- NULL
#shuffled_trial <- NULL
#shuffled_500 <- NULL
#all_interlayers_500 <- NULL

#for(i in (1:500)){

#shuffle distances -------------------------------------------------------------------------------
#  just_distances <- distances %>% subset(select = distance_in_meters) #take only distances
#  distances_shuf <- distances %>% subset(select = -distance_in_meters) #remove distances and stay with locations

#  partners_null_model_shuf <- NULL

#  print(i)
#  interlayer_shuffle <- just_distances[sample(1:nrow(just_distances)),] #shuffle distances
#  distance_combine_shuf <- cbind(distances_shuf, interlayer_shuffle) %>% rename(distance_in_meters=interlayer_shuffle) #stick shuffled distances to locations


#  distances_with_weights_shuf <- distance_combine_shuf %>% 
#    mutate(weight = interlayer_weight(distance_in_meters)) %>% #change distance to weight based on log function
#    subset(select = -distance_in_meters) #delete the distance and stay only with weights


#  interlayers_with_weights_shuf <- interlayers_new %>% inner_join(distances_with_weights_shuf, #combine weights and layers with nodes
#                                                             by = c("layer_from", "layer_to"))



#  intralayer_inverted_shuf <- tibble(values= dryad_intralayer$layer_to, dryad_intralayer$node_to, dryad_intralayer$layer_from, 
#                                dryad_intralayer$node_from, dryad_intralayer$weight) #create an inverted copy for directed intralayers
#  colnames(intralayer_inverted_shuf) <- c("layer_from", "node_from", "layer_to", "node_to", "weight") 



#  interlayer_inverted_shuf <- tibble(values= interlayers_with_weights_shuf$layer_to, interlayers_with_weights_shuf$node_to,
#                                interlayers_with_weights_shuf$layer_from, 
#                                interlayers_with_weights_shuf$node_from) #create an inverted copy for directed interlayers
#  colnames(interlayer_inverted_shuf) <- c("layer_from", "node_from", "layer_to", "node_to") 
#  interlayer_inverted_shuf <- interlayer_inverted_shuf %>% inner_join(distances_with_weights_shuf, by = c("layer_from", "layer_to"))


#plants on from
#  tot_plant_shuf <- dryad_intralayer %>% 
#    group_by(layer_from,node_from) %>% 
#    summarise(tot=sum(weight))
#  intralayer_weighted_shuf <- dryad_intralayer %>% left_join(tot_plant) %>% mutate(rel_weight=weight/tot) %>% 
#    select(-weight,-tot) %>% rename(weight=rel_weight)


#pols in from
#  tot_pol_shuf <- intralayer_inverted %>% 
#    group_by(layer_from,node_from) %>% 
#    summarise(tot=sum(weight))
#  intralayer_weighted_inverted_shuf <- intralayer_inverted %>% left_join(tot_pol) %>% mutate(rel_weight=weight/tot) %>% 
#    select(-weight,-tot) %>% rename(weight=rel_weight)

#  edgelist_non_inverted_shuf <- bind_rows(intralayer_weighted_shuf, interlayers_with_weights_shuf) #combine weighted version of intra with inter
#  edgelist_inverted_shuf <- bind_rows(intralayer_weighted_inverted_shuf, interlayer_inverted_shuf) #combine inverted version of intra with inverted version of inter

#  all_interlayers <- rbind(interlayers_with_weights_shuf, interlayer_inverted_shuf)
#  all_interlayers_500 <- rbind(all_interlayers_500, all_interlayers)

#  dryad_edgelist_complete_shuf <- bind_rows(edgelist_non_inverted_shuf, edgelist_inverted_shuf) #combine inverted and non inverted verions

#view(dryad_edgelist_complete_shuf)
#view(distances_with_weights_shuf)

#  dryad_edgelist_complete_ids_shuf <- 
#    dryad_edgelist_complete_shuf %>% 
#    left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
#    left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
#    dplyr::select(-node_from, -node_to) %>% #choose said columns
#    dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
#    left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
#    left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
#    dplyr::select(-layer_from, -layer_to) %>% 
#    dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)


#  intra_nonextended_shuf <-
#    dryad_edgelist_complete_ids_shuf %>% 
#    filter(layer_from==layer_to) %>% #only intra
#    dplyr::select(layer=layer_from, node_from, node_to, weight)
#  inter_extended_shuf <-
#    dryad_edgelist_complete_ids_shuf %>% 
#    filter(layer_from!=layer_to) #only inter

#run shuffle --------------------------------------------------------------------------------------------

#  shuffled_multilayer <- create_multilayer_object(intra = intra_nonextended_shuf,     
#                                                  inter = inter_extended_shuf,  #inter is new every time
#                                                  nodes = physical_nodes,
#                                                  layers = layer_metadata,
#                                                  intra_output_extended = T)

#  shuffled_modules <- run_infomap_multilayer(shuffled_multilayer, 
#                                             infomap_executable = "../Infomap",
#                                             flow_model = 'directed',
#                                             relax = F, 
#                                             silent = T, 
#                                             trials = 100,
#                                             seed = 497294, 
#                                             temporal_network = F,remove_auxilary_files = F)
#  C <- modules_dryad_multilayer$modules %>% select(module, node_id, layer_id)  #C should always stay the same
#  D <- shuffled_modules$modules %>% select(module, node_id, layer_id) #D should change each round


#  J <- inner_join(C,D,by=c('node_id','layer_id')) %>% 
#    group_by(module.y) %>% 
#    select(module.x) %>% table()
#  I_es <- rbind(I_es, tibble(i=i, mutual_information=NMI(J)))

#  M_D <-shuffled_modules$modules %>% slice(which.max(module)) %>% select(module) #how many modules are in the current shuffle trial
#  Num_of_modules <- rbind(Num_of_modules, tibble(shuffle_trial=i, number_of_modules=M_D))

#  shuffled_trial <- rbind(shuffled_trial, tibble(shuffle_trial= i, module= shuffled_modules$modules))

#  shuffled_500 <- rbind(shuffled_500, tibble(shuffled_modules$modules, i))

#compare partners in modules
#  for (j in 1:length(shuffled_modules$modules$node_id)){
#    current_node_shuf <- shuffled_modules$modules$node_id[j] #current node i'm looking at as a view point
#    filter_by_node_shuf <- filter(shuffled_modules$modules, shuffled_modules$modules$node_id==current_node_shuf) #only rows where node is
#    modules_with_species_shuf <- unique(filter_by_node_shuf$module) #only take unique module numbers
#    module_partners_shuf <- filter(shuffled_modules$modules, shuffled_modules$modules$module==modules_with_species_shuf)
#    just_partners_shuf <- subset(module_partners_shuf, node_id != current_node_shuf) #remove the current id from the module
#    unique_partners_shuf <- unique(just_partners_shuf$node_id) #unique partners in a specific module 
#    partners_null_model_shuf <- rbind(partners_null_model_shuf, tibble(node_id=current_node_shuf, partners=unique_partners_shuf)) #all nodes and all their partners
#    partners_null_model_shuf <- distinct(partners_null_model_shuf)
#  }

#  partners_null_model_shuf <- partners_null_model_shuf %>%
#    merge(physical_nodes, by="node_id") %>%
#    subset(select= -species)

#write_csv(partners_null_model, "./csvs/module_partners_with_fixed_interlayers.csv")

#  jaccard_val_shuf <- NULL #restart every run

#  for (k in 1:max(partners_empirical$node_id)){ #for all the nodes
#    filter_empirical <- filter(partners_empirical, partners_empirical$node_id==k) #filter only by node number k in empirical
#    filter_null_model_shuf <- filter(partners_null_model_shuf, partners_null_model_shuf$node_id==k) #filter only by node number k in null model
#    jaccard_similarity <- jaccard(filter_empirical$partners, filter_null_model_shuf$partners) #jaccard on node k in current val
#    jaccard_val_shuf <- rbind(jaccard_val_shuf, tibble(node=k, similarity=jaccard_similarity, type=filter_empirical$type)) #contains node and similarity
#  }

#  partner_comparison_shuf <- rbind(partner_comparison_shuf, tibble(trial_number=i, similarity=jaccard_val_shuf))
#  partner_comparison_shuf <- distinct(partner_comparison_shuf)

#}

#create graph for species composition 
#ggplot(I_es, aes(x=mutual_information))+geom_histogram()+
#  xlim(0.97,1)+ theme_classic()+
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
#  theme(axis.text.y=element_text(size=15))+
#  labs(x="mutual information", y="count")

#create graph for number of modules
#Num_of_modules$number_of_modules <- unlist(Num_of_modules$number_of_modules)
#threshold <- unlist(num_of_modules$number_of_modules)

#Num_of_modules %>%
#  ggplot(aes(x=number_of_modules))+
#  geom_histogram()+
#  geom_vline(aes(xintercept=threshold[1]), color="red")+
#  theme_classic()+
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
#  theme(axis.text.y=element_text(size=15))+
#  geom_text(y=400, x=87, label="empirical network", color= "red", size=5, family= "Tahoma")+
#  labs(x="number of modules", y="count")


#create graph for partner similarity
#partner_comparison_shuf %>%
#  ggplot(aes(x=similarity$similarity))+geom_histogram()+
#  theme_classic()+xlim(0,1)+
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
#  theme(axis.text.y=element_text(size=15))+
#  labs(x="similarity", y="count") 

#density plot for interlayer edges
#interlayer_500 <- all_interlayers_500 %>% cbind(type = "null")
#dryad_interlayers <- dryad_interlayers %>% cbind(type = "empirical")
#interlayer_density <- rbind(dryad_interlayers, interlayer_500)

#interlayer_density %>%
#  ggplot(aes(weight, fill = type, colour = type))+ geom_density(alpha = 0.4)+ theme_classic()

#write.csv(I_es, "./csvs/mutual_information_shuffled_vs_empirical.csv")
#write.csv(Num_of_modules, "./csvs/number_of_modules_shuffled_vs_empirical.csv")
#write.csv(partner_comparison_shuf, "./csvs/partner_similarity_shuffled_vs_empirical.csv")
#write.csv(shuffled_500, "./csvs/modules_interlayer_suffled.csv")


##put empirical networks on world map--------------------------------------------------------------------------------

#pick 75 percentile and put on map
#top_25_percent <- select(module_data_with_loc, c("module", "size_of_module")) %>% distinct %>% arrange(desc(size_of_module)) #arrange all modules to see biggest 5
#top_25_percent <- top_25_percent[1:5,] #only take top 5 modules
#top_25_percent$size_of_module <- top_25_percent$module #don't need the sizes anymore
#names(top_25_percent)[names(top_25_percent)=="size_of_module"] <- "group" #change name to group- only the top 5 groups
#pie_chart_data <- top_25_percent %>% right_join(module_data_with_loc, by=c("module"="module")) #join to have coordinates as well
#pie_chart_data$group <- as.character(pie_chart_data$group) %>% replace_na('other') #remove NA and replace with other- anything that's not top 5 modules
