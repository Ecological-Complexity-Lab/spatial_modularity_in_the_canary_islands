
# ---- NULL MODEL M3: SHUFFLING INTERACTION WITHIN ISLANDS ------------------------------------------------------------------------

# this portion of the code shuffles interactions within layers.
# the shuffled networks are then compared the empirical network to determine whether 
# local factors influence the pattern found.


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
setwd("D:/Trabajo/Papers/Canary_Island/spatial_modularity_in_the_canary_islands")
source("D:/Trabajo/Papers/Canary_Island/spatial_modularity_in_the_canary_islands/R/functions.R")

physical_nodes <- read.csv("./csvs_nuevo/physical_nodes_islands.csv")
layer_metadata <- read.csv("./csvs_nuevo/layer_metadata_islands.csv")
dryad_intralayer <- read.csv("./csvs_nuevo/intralayer_file.csv")


##---- Preparing data - layers as islands -------------------------------------------------------------------------------
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
dryad_intralayer_islands_grouped2 <- dryad_intralayer_islands %>% 
  group_by(layer_from, node_from, layer_to, node_to) %>% 
  summarise(sum_weight = sum(weight)) #turn sums of sites to sum of island


intralayers_with_ids_for_shuf <- 
  dryad_intralayer_islands_grouped2 %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, sum_weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, sum_weight)

#write.csv(intralayers_with_ids_for_shuf , "./csvs_nuevo/intralayers_with_ids_for_shuf.csv", row.names = FALSE)


##---- classic shuffle within layers ---------------------------------------------
shuf_trial_matrix_classic <- NULL
shuf_null_edge_list_classic <- NULL

for (i in 1:7){
  current_matrix <- intralayers_with_ids_for_shuf %>%
    filter(layer_from == layer_to) %>% filter(layer_from == i) %>% #take 1 layer at a time
    select(node_from, node_to, sum_weight) %>% #select interactions
    dcast(node_from ~ node_to, value.var = "sum_weight", fill = 0) %>%
    column_to_rownames(var="node_from") 
  
  current_matrix <- t(current_matrix) #put pols in rows
  
  print(i) #to keep tab on which layer we're on
  
  
  for (j in 1:1000){ #1000 iterations
    null <- vegan::nullmodel(current_matrix, method = 'r00_samp') #shuffle within island
    shuffled_matrices <- simulate(null, nsim = 1)
    
    trial_with_shuf_num <- cbind(shuffled_matrices, j) #add trial number 
    shuf_trial_matrix_classic <- rbind(shuf_trial_matrix_classic, trial_with_shuf_num) #create big matrix of all matrices
    edge_list_version_classic <- melt(shuffled_matrices) %>% filter(value > 0) %>%
      select(node_from=Var1, node_to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
    edge_list_version_classic$trial_number <- j #add trial number to the edge list
    edge_list_version_classic$layer_from <- i #intra so layer from and to are identical
    edge_list_version_classic$layer_to <- i
    shuf_null_edge_list_classic <- rbind(shuf_null_edge_list_classic, edge_list_version_classic) #create mega edge list with all repetitions
  }
}


#write.csv(shuf_null_edge_list_classic, "./csvs_nuevo/shuf_null_edge_list_classic_islands_as_layers.csv", row.names = FALSE)


##---- create interedges  --------------------------------------------------------
# keep only species that occur in 2 or more layers
co_occurrence_from_shuf<- shuf_null_edge_list_classic  %>% 
  group_by(trial_number,node_from) %>%
  mutate(num_layers_from=n_distinct(layer_from)) %>% 
  filter(num_layers_from>="2")

# a for loop that calculates all the interlayer edges based on jaccard index
interlayers_with_weights_islands <- NULL

for (trial in 1:1000){
  print(trial) #to keep tab on how far along we are
  
  co_occurrence_from_shuf2 <- co_occurrence_from_shuf %>% filter(trial_number == trial) #take only 1 trial
  
  for (i in unique(co_occurrence_from_shuf2$node_from)) {
    print(i)
    partners_sp <- 
      co_occurrence_from_shuf2 %>% 
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
      mutate(node_from=i, node_to =i, trial_num=trial) %>%
      select(c(node_from,layer_from=Var1, layer_to=Var2,node_to, weight=value,trial_num))
    interlayers_with_weights_islands <- rbind(interlayers_with_weights_islands,inter_fid)
  }
}
#write.csv(interlayers_with_weights_islands, "./csvs_nuevo/interlayer_edges_shuf_classic_islands_as_layers.csv", row.names = FALSE)

interlayer_edges_shuf_classic <- read.csv("./csvs_nuevo/interlayer_edges_shuf_classic_islands_as_layers.csv") 

interlayer_edges_shuf_classic<-interlayer_edges_shuf_classic[,c(2,1,3,4,5,6)]#change order columns

#inverted version
interlayer_inverted <- tibble(values= interlayer_edges_shuf_classic$layer_to, interlayer_edges_shuf_classic$node_to, interlayer_edges_shuf_classic$layer_from, 
                              interlayer_edges_shuf_classic$node_from, interlayer_edges_shuf_classic$weight,interlayer_edges_shuf_classic$trial_num) #create an inverted copy for directed intralayers
colnames(interlayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight","trial_num")

#Create interedgelist
edgelist_interlayers_both <- bind_rows(interlayer_edges_shuf_classic, interlayer_inverted) #combine inverted and non inverted versions of intra


##---- create intraedges inverted versions and put weight ---------------------------------
dryad_intralayer_shuf_classic <- read.csv("./csvs_nuevo/shuf_null_edge_list_classic_islands_as_layers.csv") 
dryad_intralayer_shuf_classic<- dryad_intralayer_shuf_classic[, c(5,1,6,2,3,4)]

#----inverted versions
#pols
intralayer_inverted_shuf_classic <- tibble(values= dryad_intralayer_shuf_classic$layer_to, dryad_intralayer_shuf_classic$node_to, 
                                           dryad_intralayer_shuf_classic$layer_from, dryad_intralayer_shuf_classic$node_from, 
                                           dryad_intralayer_shuf_classic$weight, dryad_intralayer_shuf_classic$trial_number) #create an inverted copy for directed intralayers
colnames(intralayer_inverted_shuf_classic) <- c("layer_from", "node_from", "layer_to", "node_to", "weight", "trial_number")

#--- weighted
## ----weighted intralayer edges--------------------------------------------------------------------------
#plants in from
tot_plant_shuf_classic <- intralayer_inverted_shuf_classic %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_shuf_classic <- intralayer_inverted_shuf_classic %>% left_join(tot_plant_shuf_classic) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_plant_shuf_classic <- tot_plant_shuf_classic[, c(3,1,2,4)]

#pols in from
tot_pol_shuf_classic <- dryad_intralayer_shuf_classic %>% 
  group_by(layer_from,node_from, trial_number) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted_shuf_classic <- dryad_intralayer_shuf_classic %>% left_join(tot_pol_shuf_classic) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)
tot_pol_shuf_classic <- tot_pol_shuf_classic[, c(3,1,2,4)]

## ----multilayer_extended_final--------------------------------------------------------------------------------------
edgelist_intralayer_shuf_classic <- bind_rows(intralayer_weighted_shuf_classic, intralayer_weighted_inverted_shuf_classic)#intraedges
edgelist_intralayer_shuf_classic<-edgelist_intralayer_shuf_classic[,c(1,2,3,4,6,5)]#change order columns
colnames(edgelist_intralayer_shuf_classic)[6]<-"trial_num"


dryad_edgelist_complete_shuf_classic <- bind_rows(edgelist_intralayer_shuf_classic, edgelist_interlayers_both) #combine inter and intra
#write.csv(dryad_edgelist_complete_shuf_classic, "./csvs_nuevo/dryad_edgelist_complete_shuf_classic_islands_as_layers.csv", row.names = FALSE)

# ----multilayer_class-----------------------------------------------------------------------------------------------
layer_metadata <- read.csv("./csvs_nuevo/layer_metadata_islands.csv")
physical_nodes <- read.csv("./csvs_nuevo/physical_nodes_islands.csv")

## ---- Calculate number of modules and test distance decay in shuf networks --------------------------------------------------------

## ---- Number of modules ------------------------------------------------------------
dryad_edgelist_complete_shuf_classic <- as.data.frame(dryad_edgelist_complete_shuf_classic) 
colnames(dryad_edgelist_complete_shuf_classic)[6]<-"trial_number"
dryad_multilayer_shuf_1000_classic <- NULL

#calculate modularity for each simulation
dryad_multilayer_shuf_1000_classic_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_classic, 
                                                                 dryad_multilayer_shuf_1000_classic)

dryad_multilayer_shuf_1000_classic_output <- dryad_multilayer_shuf_1000_classic_output %>% drop_na() 
#delete all NAs that are in the physical nodes but not the network (and were not assigned to modules)

#write.csv(dryad_multilayer_shuf_1000_classic_output, "./csvs_nuevo/dryad_multilayer_shuf_1000_classic_output_islands.csv", row.names = FALSE)

##---- Distance decay of modules  ---------------------------------------------------------------
distances_with_ids <- read.csv("./csvs_nuevo/distances_with_ids_islands_as_layers.csv")
dryad_multilayer_shuf_1000_classic_output<- read.csv("./csvs_nuevo/dryad_multilayer_shuf_1000_classic_output_islands.csv")

#pivot modules function for islands as layers
pivot_by_module_islands <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(data , id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s2<-na.omit(s2)
  s4 = t(s2) 
  s4 <- s4[-1,]
  colnames(s4) <- c(1,2,3,4,5,6,7)
  return(s4)
}

# this function calculates the Jaccard Similarity in modules between islands
module_distance_decay_islands_func <- function(multilayer_1000, 
                                               layers_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial) #to keep tab on how far along we are
    
    modules_for_similarity_shuf <- multilayer_1000  %>% filter(trial_num ==trial) #take only 1 trial
    
    #pivot modules
    module_pivoted_shuf <- pivot_by_module_islands(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    # sanity check
    module_pivoted_shuf[module_pivoted_shuf > 0] <-  1
    counted <- sum(rowSums(module_pivoted_shuf) > 1) # if this is larger then 0 then there is similarity
    counted
          
    #create edge list with distances
    modules_edge_list_shuf <- NULL
    
    if (counted == 0) {
      modules_edge_list_shuf <- tibble(layer_from = int(), layer_to = int(), module = character())
    }
    else{
      for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
        modules_edge_list_shuf <- edge_list_per_module_islands(module_pivoted_shuf[k,], modules_edge_list_shuf, k) 
        current_module <- rownames(module_pivoted_shuf)[k]
        if (is.null(modules_edge_list_shuf)) next
        modules_edge_list_shuf <- modules_edge_list_shuf %>% mutate(module = replace_na(module, current_module)) #add module number
      }
    }

    edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
    #edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #we remove this line (which it removes NA), because we have a lot locations don't have modules in common
    
    
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
edge_list_per_module_islands <- function(data,edge_list, k){
  #gets one row from a data frame and creates an edge list from it
  
  for (i in (1:6)){
    if (data[i]==0) next #only take layers where the module is present
    else {
      for (j in ((i+1):7)){
        if (data[j]==0) next #only take layers where the module is present
        else {
          edge_list <-rbind(edge_list,tibble(layer_from=i, layer_to=j, module=as.character(NA))) #create edge list of all the layer found in a module
        }
      }
    }
  }
  return(edge_list)
}

turnover_with_distance_classic <- NULL
layers_turnover_with_distnace<-NULL
module_layer_turnover_shuf <- NULL


all_edge_list_islands_combine_no_module_shuf_classic_output <- module_distance_decay_islands_func(dryad_multilayer_shuf_1000_classic_output,
                                                                                                  turnover_with_distance_classic)

#write.csv(all_edge_list_islands_combine_no_module_shuf_classic_output, "./csvs_nuevo/all_edge_list_islands_combine_no_module_shuf_classic.csv", row.names = FALSE)

#---- create ave for jaccard layers-------------------------------------------------------------------------
all_edge_list_islands_combine_no_module_shuf_classic_output <- read.csv("./csvs_nuevo/all_edge_list_islands_combine_no_module_shuf_classic.csv")

ave_module_islands_turnover_shuf_classic <- all_edge_list_islands_combine_no_module_shuf_classic_output %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(mean_distance)) %>% mutate(type="null_model_within") #create mean and sd for each point



#add empirical
islands_turnover_with_distnace_empirical <- read.csv("./csvs_nuevo/islands_turnover_with_distnace_empirical.csv")

empirical_turnover_for_modules_layers_shuf <- islands_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave = mean(turnover), sd = sd(turnover), ave_dist = mean_distance) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#---- combine for jaccard analysis------------------------------------------------------------------------------------------------------
#combine all layers
jaccard_similarity_empirical_and_null_classic <- rbind(empirical_turnover_for_modules_layers_shuf, 
                                                       ave_module_islands_turnover_shuf_classic)

jaccard_similarity_layer_empirical_and_null_km_classic <- jaccard_similarity_empirical_and_null_classic %>% 
  mutate(mean_dist_in_km = ave_dist/1000)

#write.csv(jaccard_similarity_layer_empirical_and_null_km_classic,"./csvs_nuevo/jaccard_similarity_layer_empirical_and_null_km_classic_islands.csv", row.names = FALSE)


#---- graphs for distance decay in modules shuf vs empirical--------------------------------
#jaccard_similarity_layer_empirical_and_null_km_classic <- read.csv("./csvs_nuevo/jaccard_similarity_layer_empirical_and_null_km_classic_islands.csv") #need to read this to run next part

pdf('./graphs/M3_Modules_DD_Islands.pdf', 10, 6)
jaccard_similarity_layer_empirical_and_null_km_classic %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  scale_color_manual (name = "Null Model", labels = c("E",expression("M"[3])),
                      values = c("#FB3B1E","#FA86F2"))+
  labs(x="Distance (Km)", y="Jaccard Similarity")+  #stat_cor(aes(label = ..rr.label..))+
  
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
emp<-jaccard_similarity_layer_empirical_and_null_km_classic %>% filter(type=="empirical")
null<-jaccard_similarity_layer_empirical_and_null_km_classic %>% filter(type=="null_model_within")

shapiro.test(emp$ave)
shapiro.test(null$ave)#normal

#---- linear regression
lm1_module_classic = lm(ave ~ mean_dist_in_km ,data=emp) #in empirical
lm2_module_classic = lm(ave ~ mean_dist_in_km, data=null) #in null model
summary(lm1_module_classic)
summary(lm2_module_classic)

#----correlation and r sqaured between jaccard and distance for each run
all_edge_list_islands_combine_no_module_shuf_classic <- read.csv("./csvs_nuevo/all_edge_list_islands_combine_no_module_shuf_classic.csv")

iteration_correlation_classic <- NULL
iteration_correlation_data_classic <- all_edge_list_islands_combine_no_module_shuf_classic %>% subset(layer_from != layer_to) 
iteration_correlation_data_classic_km <- iteration_correlation_data_classic %>% 
  mutate(mean_dist_in_km = mean_distance/1000)


for (i in 1:1000){
  trial_classic = iteration_correlation_data_classic_km %>% filter(trial == i)
  iteration_correlation_new_classic <- cor.test(trial_classic$turnover, trial_classic$mean_dist_in_km, method = "pearson")
  lm_val_classic <- lm(turnover ~ mean_dist_in_km, data = trial_classic)
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

iteration_correlation_classic
#write.csv(iteration_correlation_classic, "./csvs_nuevo/iteration_correlation_classic_islands.csv", row.names = FALSE)


##distribution of rsquared and empirical
correlation_empirical_classic <- read.csv("./csvs_nuevo/correlation_empirical.csv")
iteration_correlation_classic <- read.csv("./csvs_nuevo/iteration_correlation_classic_islands.csv")
iteration_correlation_classic2<-iteration_correlation_classic %>% mutate(Type = "null_class")

pdf('./graphs/M3_r_squares_module_DD.pdf', 10, 6)
iteration_correlation_classic2 %>% ggplot(aes(x = rsquared, fill= Type))+ 
  geom_density(alpha = 0.6)+ 
  geom_vline(xintercept = correlation_empirical_classic$rsquared, linetype = "dashed", color = "#FB3B1E") +
  labs(x= expression("R"^2), y="Density")+  
  scale_fill_manual(name = "Null Model",  label = expression("M"[3]), values= "#FA86F2")+
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

p_rsquared_classic <- sum(iteration_correlation_classic$rsquared > correlation_empirical_classic$rsquared)/1000
p_rsquared_classic

