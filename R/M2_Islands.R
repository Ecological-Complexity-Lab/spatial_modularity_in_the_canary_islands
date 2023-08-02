
# ---- NULL MODEL M2: SHUFFLING INTERACTION BETWEEN ISLANDS ------------------------------------------------------------------------

# this portion of the code shuffles interactions between layers. We randomly shuffled interactions of each pair of species between all 
#the layers in which they co-occur. For example, if a plant and a pollinator co-occurred in layers 1 and 3 but interacted only in layer 1,
#they would still co-occur in the same layers but may interact in layer 3 after shuffling. Hence, null model M3 shuffles intralayer links 
#between layers but conserves the interlayer links. The shuffled networks are then compared the empirical network to determine whether 
# turnover of interactions is influencing the pattern found.

#Also, we asses the contribution of both beta diversity components (turnover in species and turnover in interactions)in the patter found.


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
library(ggpubr)

##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")
physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")
dryad_edgelist_complete_ids <- read.csv("./csvs/Islands/dryad_edgelist_complete_ids_islands.csv")

#---- number of layers each species is found in------------------------------------------------------
#plants
plant_half <- dryad_edgelist_complete_ids %>% filter(layer_from == layer_to) %>% #take only intralayers
  filter(node_from <= 39) %>% select(layer_from, node_from) #one directed for this purpose, only plants

plant_half <- plant_half %>% unique() #remove doubles created due to aggregation
plant_half$count <- 1 
plant_half <- plant_half %>% group_by(node_from) %>% mutate(number_of_layers = sum(count)) #count number of layers 

#pollinators
pollinator_half <- dryad_edgelist_complete_ids %>% filter(layer_from == layer_to) %>% #take only intralayers
  filter(node_to > 39) %>% select(layer_from, node_to) #one directed for this purpose, only plants

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
interaction_count <- interaction_count %>% select(layer_from, node_from, layer_to, node_to) %>% unique() #remove doubles created due to aggregation
interaction_count$count <- 1 #add count we can sum later
interaction_count <- aggregate(count ~ node_from + node_to, interaction_count, FUN = sum) #count number of times nodes interacted
interaction_count <- rename(interaction_count, interactions = count) #rename to have distinct name

#----number of times species co-occurred--------------------------------------------------------------
co_occurrence_count <- NULL
co_occurrence <- dryad_edgelist_one_sided %>% filter(layer_from == layer_to) #only intra as we need just pollination interactions to change
co_occurrence <- co_occurrence %>% select(layer_from, node_from, layer_to, node_to) %>% unique() #remove doubles created due to aggregation

for (i in 1:39){
  current_plant_interactions <- co_occurrence %>% filter(node_from == i) #focal point is plants
  unique_partners <- current_plant_interactions %>% select(node_from, node_to) %>% unique() 
  for (j in 1:nrow(unique_partners)){ #find all unique partners to not have doubles
    current_interaction <- unique_partners[j,] 
    current_plant <- current_interaction$node_from #who's current plant
    current_pollinator <- current_interaction$node_to #who's current pollinator
    
    layers_for_plant <- co_occurrence %>% filter(node_from == current_plant) %>%
      select(layer_from) %>% unique() %>% unlist() #find all sites current plant is found in
    layers_for_pollinator <- co_occurrence %>% filter(node_to == current_pollinator) %>% 
      select(layer_from) %>% unique() %>% unlist() #find all sites current pollinator is found in
    
    intersects <- intersect(layers_for_plant, layers_for_pollinator) #find which sites have both species
    intersect_value <- length(intersects) #count number of sites that have both species
    co_occurrence_count <- rbind(co_occurrence_count, tibble(node_from = current_plant, 
                                                             node_to = current_pollinator, 
                                                             co_occurrences = intersect_value))
  }
}

#write.csv(co_occurrence_count, "./csvs/Islands/co_occurrence_count.csv", row.names = FALSE)

#---- combine interactions and co-occurrences--------------------------------------------
interactions_co_occurences <- merge(co_occurrence_count, interaction_count, by = c("node_from", "node_to")) #combine two data frames

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

#write.csv(layers_for_pairs, "./csvs/Islands/layers_for_pairs.csv", row.names = FALSE)


##---- Calculate the raw weight of intraedges using layers as islands -------------------------------------------------------------------------------
#physical_nodes <- read.csv("./csvs/Islands/physical_nodes.csv")
#layer_metadata <-read.csv("./csvs/Islands/layer_metadata.csv")
dryad_intralayer <- read.csv("./csvs/intralayer_file.csv")

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

intralayer_edges_ids <- # Replace the node names with node_ids
  dryad_intralayer_islands_grouped %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, sum_weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight = sum_weight)

prueba1<- dryad_edgelist_complete_ids %>% filter(layer_from == layer_to) %>% filter (layer_from ==1)  #recheck porque tengo mas interacciones
#en dryad (por ejemolo sp1 en capa 1)
  
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

#write.csv(interactions_shuf, "./csvs/Islands/interactions_shuf.csv", row.names = FALSE)
interactions_shuf <- read.csv("./csvs/Islands/interactions_shuf.csv")
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

#write.csv(intralayer_edges_shuf, "./csvs/Islands/intralayer_edges_shuf.csv", row.names = FALSE) #create file to run interlayers on HPC
#intralayer_edges_shuf <- read.csv("./csvs/Islands/intralayer_edges_shuf.csv")

#---- interlayer edges-------------------------------------------------------------------------------------- # see if should copy the script of m2
interlayer_interactions_from_shuf <- intralayer_edges_shuf %>% group_by(trial_num, node_from) %>%
  select(layer_from, node_from) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_from[1], loc2 = layer_from[2], loc3 = layer_from[3], 
         loc4 = layer_from[4], loc5 = layer_from[5], loc6 = layer_from[6],
         loc7 = layer_from[7]) #all islands the species is found in

interlayer_interactions_to_shuf <- intralayer_edges_shuf %>% group_by(trial_num, node_to) %>% 
  select(layer_to, node_to) %>% unique() %>% #group by species and find only locations
  mutate(loc1 = layer_to[1], loc2 = layer_to[2], loc3 = layer_to[3],
         loc4 = layer_to[4], loc5 = layer_to[5], loc6 = layer_to[6],
         loc7 = layer_to[7])  %>% #all islands the species is found in
  dplyr::rename(layer_from = layer_to, node_from = node_to) #make sure they look the same for rbind

interlayer_edges_interactions_shuf <- rbind(interlayer_interactions_from_shuf, interlayer_interactions_to_shuf) 

#write.csv(interlayer_edges_interactions_shuf, "./HPC_Islands/shuf_interactions_layers_islands/interlayer_edges_interactions_shuf.csv",  row.names = FALSE)
#interlayer_edges_interactions_shuf <- read.csv("./HPC_Islands/shuf_interactions_layers_islands/interlayer_edges_interactions_shuf.csv")


#---- run on HPC and come back
#there are 3 code portions for the HPC:
# 1. HPC_network_interactions_shuffle.R 
# 2. 1_1000_interactions.sh 
# 3. i_interactions.sh 
# running 1_1000_interactions.sh manually in the cmd will make the other two run.
# all 1000 result csvs can be found in /HPC_Islands/shuf_interactions_layers_islands/csvs_interactions files

files_interactions <- list.files("./HPC_Islands/shuf_interactions_layers_islands/csvs_interactions/", pattern = ".csv$", recursive = TRUE, full.names = TRUE)
my_merged_interlayer_shuf_interactions <- read_csv(files_interactions) %>% bind_rows() #create a long edge list with all the csvs

#write.csv(my_merged_interlayer_shuf_interactions, "./csvs/Islands/my_merged_interlayer_shuf_interactions.csv", row.names = FALSE)

#---- interlayers with weights shuf version ------------------------------------------
distances_with_weights_ids <- read.csv("./csvs/Islands/distances_with_weights_ids_islands_as_layers.csv")

interlayers_with_weights_shuf_interactions <- my_merged_interlayer_shuf_interactions %>% inner_join(distances_with_weights_ids, 
                                                                                                    by = c("layer_from", "layer_to")) %>% unique()
interlayers_with_weights_shuf_interactions <- 
  interlayers_with_weights_shuf_interactions[!duplicated(interlayers_with_weights_shuf_interactions[c(1,3,5,6)]),]

#write.csv(interlayers_with_weights_shuf_interactions, "./csvs/Islands/interlayers_with_weights_shuf_interactions.csv", row.names = FALSE)


## ----multilayer_extended_final--------------------------------------------------------------------------------------
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")
physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")
#intralayer_edges_shuf <- read.csv("./csvs/Islands/intralayer_edges_shuf.csv")

# Input: An extended edge list.
names(intralayer_edges_shuf)[6] <- "trial_number" #change name of column so inter and intra correspond

dryad_edgelist_complete_shuf_interactions <- bind_rows(intralayer_edges_shuf, 
                                                       interlayers_with_weights_shuf_interactions) #combine inter and intra

#write.csv(dryad_edgelist_complete_shuf_interactions, "./csvs/Islands/dryad_edgelist_complete_shuf_interactions_islands_as_layers.csv", row.names = FALSE)
#dryad_edgelist_complete_shuf_interactions <- read.csv("./csvs/Islands/dryad_edgelist_complete_shuf_interactions_islands_as_layers.csv")

#calculate modularity for each simulation
dryad_multilayer_shuf_1000_interactions <- NULL
dryad_multilayer_shuf_1000_interactions_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_interactions, 
                                                                      dryad_multilayer_shuf_1000_interactions)

#write.csv(dryad_multilayer_shuf_1000_interactions_output, "./csvs/Islands/dryad_multilayer_shuf_1000_interactions_output.csv", row.names = FALSE)
dryad_multilayer_shuf_1000_interactions_output <- read.csv("./csvs/Islands/dryad_multilayer_shuf_1000_interactions_output.csv")

##---- distance decay of modules  ---------------------------------------------------------------
distances_with_ids <- read.csv("./csvs/Islands/distances_with_ids_islands_as_layers.csv")
physical_nodes <- read.csv("./csvs/Islands/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs/Islands/layer_metadata_islands.csv")

#pivot modules function for islands as layers
pivot_by_module_islands <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(modules_for_similarity_shuf, id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s2<-na.omit(s2)
  s3 <-s2%>%
    select(where(~ any(. != 0)))
  s4 = t(s3) 
  s4 <- s4[-1,]
  colnames(s4) <- c(1,2,3,4,5,6,7)
  return(s4)
}



# this function calculates the Jaccard Similarity in modules between islands
module_distance_decay_islands_func <- function(multilayer_1000, 
                                               layers_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial) #to keep tab on how far along we are
    
    modules_for_similarity_shuf <- multilayer_1000 %>% filter(trial_num == trial) #take only 1 trial
    
    #pivot modules
    module_pivoted_shuf <- pivot_by_module_islands(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    #create edge list with distances
    modules_edge_list_shuf <- NULL
    
    for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module_islands(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[k]
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

turnover_with_distnaces_interactions <- NULL
layers_turnover_with_distnace<-NULL
module_layer_turnover_shuf <- NULL

all_edge_list_layer_combine_no_module_shuf_interactions_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_interactions_output,
                                                                                                     turnover_with_distnaces_interactions)


#write.csv(all_edge_list_layer_combine_no_module_shuf_interactions_output, 
 #    "./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_interactions_output_islands_as_layers.csv", 
  #  row.names = FALSE)

#all_edge_list_layer_combine_no_module_shuf_interactions_output <- read.csv("./csvs/Islands/all_edge_list_layer_combine_no_module_shuf_interactions_output_islands_as_layers.csv")

#---- create ave for jaccard layers-------------------------------------------------------------------------
ave_module_layer_turnover_shuf_interactions <- all_edge_list_layer_combine_no_module_shuf_interactions_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(mean_distance)) %>% mutate(type="null_model") #create mean and sd for each point

#add empirical
islands_turnover_with_distnace_empirical <- read.csv("csvs/Islands/islands_turnover_with_distnace_empirical.csv")

empirical_turnover_for_modules_layers_shuf <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), ave_dist = mean_distance) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#---- combine for jaccard analysis------------------------------------------------------------------------------------------------------
#combine all layers
jaccard_similarity_empirical_and_null_interactions <- rbind(empirical_turnover_for_modules_layers_shuf, 
                                                            ave_module_layer_turnover_shuf_interactions)

jaccard_similarity_islands_empirical_and_null <- rbind(empirical_turnover_for_modules_layers_shuf,ave_module_layer_turnover_shuf_interactions)

jaccard_similarity_layer_empirical_and_null_km <- jaccard_similarity_islands_empirical_and_null %>% 
  mutate(mean_dist_in_km = ave_dist/1000)

#write.csv(jaccard_similarity_layer_empirical_and_null_km,
 #    "./csvs/Islands/jaccard_similarity_layer_empirical_and_null_km_islands.csv", row.names = FALSE)


#---- graphs for distance decay in modules shuf vs empirical--------------------------------
jaccard_similarity_layer_empirical_and_null_km <- read.csv("./csvs/Islands/jaccard_similarity_layer_empirical_and_null_km_islands.csv") #need to read this to run next part

pdf('./graphs/Islands/M2_Modules_DD_Islands.pdf', 10, 6)
jaccard_similarity_layer_empirical_and_null_km %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ 
  theme_classic()+ geom_smooth(method= "lm", se=F)+
   scale_color_manual (name = "Null Model", labels = c("E",expression("M"[2])),
                                                       values = c("#FB3B1E",  "#E6AB02" ))+
  labs(x="Distance (Km)", y="Jaccard Similarity")+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(color = "black",fill = NA,size = 1),
          panel.spacing = unit(0.5, "cm", data = NULL),
          axis.text = element_text(size=15, color='black'),
          axis.title = element_text(size=17, color='black'),
          axis.line = element_blank(),
          legend.text.align = 0,
          legend.title =  element_text(size = 13, color = "black"),
          legend.text = element_text(size = 11))

dev.off()

#---- statistical analysis--------------------------------------------------

emp<-jaccard_similarity_layer_empirical_and_null_km %>% filter(type=="empirical")
null<-jaccard_similarity_layer_empirical_and_null_km %>% filter(type=="null_model")


shapiro.test(emp$ave)
shapiro.test(null$ave)#not normal

#---- linear regression
lm1_module_interactions = lm(ave ~ ave_dist ,data=emp) #in empirical
lm2_module_interactions = lm(ave ~ ave_dist ,data=null) #in empirical

#glm2_module_interactions = glm(ave ~ ave_dist, family = Gamma,data=null) #in null model
#summary(glm2_module_interactions)
#with(summary(glm2_module_interactions), 1 - deviance/null.deviance) #Rsquared for glm

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

##correlation and r sqaured between jaccard and distance for each run ----------------------------------------------------------------------------

#---- linear regression and correlation
iteration_correlation_interactions <- NULL
iteration_correlation_data_interactions <- all_edge_list_layer_combine_no_module_shuf_interactions_output %>% subset(layer_from != layer_to) 

iteration_correlation_data_interactions_km <- iteration_correlation_data_interactions %>% 
  mutate(mean_dist_in_km = mean_distance/1000)


for (i in 1:1000){
  trial_interactions = iteration_correlation_data_interactions_km %>% filter(trial == i)
  iteration_correlation_new_interactions <- cor.test(trial_interactions$turnover, trial_interactions$mean_dist_in_km, method = "pearson")
  lm_val_interactions <- lm(turnover ~ mean_dist_in_km, data = trial_interactions)
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

#write.csv(iteration_correlation_interactions, "./csvs/Islands/iteration_correlation_interactions_islands.csv", row.names = FALSE)

#correlation empirical
#classic_layers_turnover_with_distances <- read.csv("./csvs/Islands/classic_layers_turnover_with_distances_islands_as_layers.csv")
classic_layers_turnover_with_distances <- classic_layers_turnover_with_distances [,-1]
classic_layers_turnover_with_distances_km <- classic_layers_turnover_with_distances %>% 
  mutate(mean_dist_in_km = mean_distance/1000)

layer_turnover_with_distnace_empirical_no_loop <- classic_layers_turnover_with_distances_km %>% subset(layer_from != layer_to) 

correlation_empirical_data_interactions <- cor.test(layer_turnover_with_distnace_empirical_no_loop$turnover, #ave is just the value of the turnover
                                            layer_turnover_with_distnace_empirical_no_loop$mean_dist_in_km, method = "pearson")

lm_val_empirical_interactions <- lm(turnover ~ mean_dist_in_km, data = layer_turnover_with_distnace_empirical_no_loop)

correlation_empirical_interactions <- tibble(estimate = correlation_empirical_data_interactions$estimate, 
                                     p_val = correlation_empirical_data_interactions$p.value, 
                                     statistic = correlation_empirical_data_interactions$statistic, 
                                     confidence_int_low = correlation_empirical_data_interactions$conf.int[1],
                                     confidence_int_high = correlation_empirical_data_interactions$conf.int[2],
                                     slope = lm_val_empirical_interactions$coefficients[2],
                                     intercept = lm_val_empirical_interactions$coefficients[1], 
                                     rsquared = summary(lm_val_empirical_interactions)$adj.r.squared)

#write.csv(correlation_empirical_interactions, "./csvs/Islands/correlation_empirical_interactions.csv", row.names = FALSE) #so it can be used for classical shuffling
#correlation_empirical_interactions <- read.csv("./csvs/Islands/correlation_empirical_interactions.csv")

##distribution of rsquared and add empirical
correlation_empirical_interactions <- read.csv("./csvs/Islands/correlation_empirical_interactions.csv")
iteration_correlation_interactions <- read.csv("./csvs/Islands/iteration_correlation_interactions_islands.csv")
iteration_correlation_interactions2<-iteration_correlation_interactions %>% mutate(Type = "null_int")

pdf('./graphs/Islands/M2_r_squares_module_DD.pdf', 10, 6)
iteration_correlation_interactions2 %>% ggplot(aes(x = rsquared, fill= Type))+ 
  geom_density(alpha = 0.6)+ 
  geom_vline(xintercept = correlation_empirical_interactions$rsquared, linetype = "dashed", color = "#FB3B1E") +
  labs(x= expression("R"^2), y="Density")+  
  scale_fill_manual(name = "Null Model",  label = expression("M"^2), values= "#E6AB02")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=15, color='black'),
        axis.title = element_text(size=17, color='black'),
        axis.line = element_blank(),
        legend.text.align = 0,
        legend.title =  element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11))

dev.off()

p_rsquared_interactions <- sum(iteration_correlation_interactions$rsquared > correlation_empirical_interactions$rsquared)/1000
p_rsquared_interactions

##compare both parts of beta diversity (shuffling pollinator species )
iteration_correlation_interactions <- read.csv("./csvs/Islands/iteration_correlation_interactions_islands.csv")
iteration_correlation_pols <- read.csv("./csvs/Islands/iteration_correlation_pols.csv")

interactions_for_beta_diversity <- iteration_correlation_interactions
interactions_for_beta_diversity$type <- "shuffling_interactions" 

species_for_beta_diversity <- iteration_correlation_pols
species_for_beta_diversity$type <- "shuffling_species"

correlation_beta_diversity <- rbind(interactions_for_beta_diversity, species_for_beta_diversity)


pdf('./graphs/Islands/Betadiversity_module_DD.pdf', 10, 6)
correlation_beta_diversity %>% ggplot(aes(x = rsquared, color = type, fill = type))+ #overlay shuffling of interactions and species
  geom_density(alpha = 0.6)+ scale_color_manual(values = c("#E6AB02", "#A44CD3"))+ 
  scale_fill_manual(values = c("#E6AB02", "#A44CD3"))+
  theme_classic()+ labs(x = "R squared")+
  geom_vline(xintercept = correlation_empirical_interactions$rsquared, linetype = "dashed", color = "#FB3B1E") 

dev.off()

#non-parametric test to compare two groups
wilcox.test(iteration_correlation_pols$rsquared,iteration_correlation_interactions$rsquared,) # where y and x are numeric
 




