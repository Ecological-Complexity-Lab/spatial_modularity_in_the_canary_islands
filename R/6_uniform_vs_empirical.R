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

#this portion of the code refers to the uniform analysis where the interlayer edges are all uniform
#and fixed to the median value. distance decay in modules of the empirical network is compared to the
#distance decay in modules shown in the uniform network

setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")

##---- fixed for empirical data ---------------------------------------------------------------------------
inter_extended <- read.csv("./csvs/dryad_only_interlayer_edges.csv")
physical_nodes <- read.csv("./csvs/physical_nodes.csv")
layer_metadata <-read.csv("./csvs/layer_metadata.csv")
intra_nonextended <- read.csv("./csvs/dryad_only_intralayer_edges.csv")

interlayer_edges_change <- select(inter_extended, -weight) #create data frame where weight doesn't exist
interlayer_edges_change$weight <- 0.357602 #set fixed interlayer value to median of interlayer distribution

dryad_multilayer_fixed <- create_multilayer_object(intra = intra_nonextended, 
                                                   inter = interlayer_edges_change, #create multilayer with new fixed inter value
                                                   nodes = physical_nodes,
                                                   layers = layer_metadata,
                                                   intra_output_extended = T) 

modules_edge_change_fixed <- modified_multi(dryad_multilayer_fixed,
                                                    infomap_executable = "Infomap",
                                                    flow_model = 'directed',
                                                    relax = F, 
                                                    silent = T, 
                                                    trials = 100,
                                                    seed = 497294, 
                                                    temporal_network = F)

# fixed
modules_fixed <- modules_edge_change_fixed$modules #35 modules

## empirical
modules <- read.csv('./csvs/modules_dryad_multilayer.csv')

modules_dryad_multilayer_fixed_analysis <- modules #create version just for the analysis

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

#edge list by layer for fixed------------------------------------------------------------------------------------------------------------
distances_with_ids <- read.csv("./csvs/distances_with_ids.csv")
edge_list_with_distances_fixed <- right_join(modules_edge_list_fixed_analysis, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
edge_list_with_distances_fixed <- na.omit(edge_list_with_distances_fixed) #remove NA and delete layer name

#version with # of modules in layers
edge_list_by_layer_modules_fixed <- edge_list_with_distances_fixed %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) 
edge_list_by_layer_modules_fixed$count <- c(1)
edge_list_by_layer_modules_fixed <- edge_list_by_layer_modules_fixed %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules)

#version with correct average between layers
edge_list_by_layer_ave_fixed <- edge_list_with_distances_fixed %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_layer_combine_fixed <- edge_list_by_layer_ave_fixed %>%
  merge(edge_list_by_layer_modules_fixed, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_layer_combine_no_module_fixed_analysis <- edge_list_layer_combine_fixed %>% select(-module) %>% unique() #have version where modules aren't present

#module data
size <- count(modules_edge_change_fixed$modules, module)  #create a data frame of all modules and how many nodes are in each (size of module)
module_data <- merge(modules_edge_change_fixed$modules , size, by=c("module","module")) #merge size of module with all the other info about the modules
colnames(module_data)[7] <- "size_of_module" #rename column

lon_lat_data <- read_csv('./csvs/layers.csv') #create new data frame with just the layer data
lon_lat_data <- lon_lat_data %>% select(c("layer_id","lat","Lon")) %>% na.omit()  #only select layer id and coordinates

module_data_with_loc <- merge(module_data, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon <- module_data_with_loc %>% select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon$count <- c(1)

#write.csv(modules_with_lat_lon, "csvs/modules_with_lat_lon_fixed.csv", row.names = FALSE)


# total number of modules in each layer
module_layer_turnover_fixed_analysis <- NULL


for (i in 1:14){
  for (j in 1:14){
    modules_in_layer_from_fixed_analysis <- filter(modules_with_lat_lon, layer_id == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_layer_to_fixed_analysis <- filter(modules_with_lat_lon, layer_id == j) %>% select(module) %>% unique() %>% unlist()
    #take all sub modules in layer_from and all sub modules in layer_to to check turnover
    int_both <- intersect(modules_in_layer_from_fixed_analysis, modules_in_layer_to_fixed_analysis) #how many sub modules are common in both layers
    uni_both <- union(modules_in_layer_from_fixed_analysis, modules_in_layer_to_fixed_analysis) #how many sub modules are found in both layers in total
    turnover <- length(int_both)/length(uni_both)
    module_layer_turnover_fixed_analysis <- rbind(module_layer_turnover_fixed_analysis, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

module_layer_turnover_fixed_analysis <- drop_na(module_layer_turnover_fixed_analysis)


layer_turnover_with_distnace_fixed <- edge_list_by_layer_ave_fixed %>%
  merge(module_layer_turnover_fixed_analysis, by= c("layer_from", "layer_to")) #merge both versions

#---- edge list by layer empirical-------------------------------------------------------------------------------------------------------------
edge_list_with_distances_empirical <- right_join(modules_edge_list_empirical_fixed_analysis, 
                                                 distances_with_ids, by= c("layer_from", "layer_to"))
edge_list_with_distances_empirical <- na.omit(edge_list_with_distances_empirical) #remove NA 

#version with # of modules in layers
edge_list_by_layer_modules_empirical <- edge_list_with_distances_empirical %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) 
edge_list_by_layer_modules_empirical$count <- c(1)
edge_list_by_layer_modules_empirical <- edge_list_by_layer_modules_empirical %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules)

#total number of modules in each layer
modules_with_lat_lon_empirical <- read.csv("csvs/modules_with_lat_lon.csv")

module_layer_turnover_empirical <- NULL


for (i in 1:14){
  for (j in 1:14){
    modules_in_layer_from_empirical_analysis <- filter(modules_with_lat_lon_empirical, layer_id == i) %>% select(module) %>% unique() %>% unlist()
    modules_in_layer_to_empirical_analysis <- filter(modules_with_lat_lon_empirical, layer_id == j) %>% select(module) %>% unique() %>% unlist()
    #take all sub modules in layer_from and all sub modules in layer_to to check turnover
    int_both <- intersect(modules_in_layer_from_empirical_analysis, modules_in_layer_to_empirical_analysis) #how many sub modules are found in both layers
    uni_both <- union(modules_in_layer_from_empirical_analysis, modules_in_layer_to_empirical_analysis)
    turnover <- length(int_both)/length(uni_both)
    module_layer_turnover_empirical <- rbind(module_layer_turnover_empirical, tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

module_layer_turnover_empirical <- drop_na(module_layer_turnover_empirical)

edge_list_by_layer_ave_empirical <- edge_list_with_distances_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

layer_turnover_with_distnace_empirical <- edge_list_by_layer_ave_empirical %>%
  merge(module_layer_turnover_empirical, by= c("layer_from", "layer_to")) #merge both versions

#----prepare data for graph--------------------------------------------------------------------------------------------------------
ave_module_layer_turnover_shuf_fixed <- layer_turnover_with_distnace_fixed %>% 
  group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="null_fixed") #create mean and sd for each point

## empirical
empirical_turnover_for_module_layer_analysis <- layer_turnover_with_distnace_empirical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#combine all layers
jaccard_similarity_empirical_and_fixed <- rbind(empirical_turnover_for_module_layer_analysis, 
                                                ave_module_layer_turnover_shuf_fixed)

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
pdf('./graphs/uniform_interlayer_values/M4_distance_decay_in_modules.pdf', 10, 6)
jaccard_similarity_empirical_and_fixed_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ scale_color_manual(values = c("#F47069", "#c4ad06"))+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
#+ stat_cor(aes(label = ..p.label..), label.x = 400)+
  #stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.65, 0.62))
dev.off()

#write.csv(jaccard_similarity_empirical_and_fixed_no_self_loop_km, "csvs/jaccard_similarity_empirical_and_fixed_no_self_loop_km.csv", row.names = FALSE)

#version with no trendline
#jaccard_similarity_empirical_and_fixed_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave, group= type, color= type))+
 # geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ 
#  geom_smooth(data = jaccard_similarity_empirica_no_self_loop_km, method= "lm", se=F)+
#  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
 # theme(axis.text.y=element_text(size=15))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  #labs(x="Distance in Km", y="Jaccard Similarity")+ scale_color_manual(values = c("#F47069", "#c4ad06"))

#------check if its significant for layers----------------------------------------------------------------
lm1_module = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_no_self_loop_km,
                                            jaccard_similarity_empirical_and_fixed_no_self_loop_km$type=="empirical")) #in empirical
lm2_module = lm(ave ~ ave_dist ,data=subset(jaccard_similarity_empirical_and_fixed_no_self_loop_km,
                                            jaccard_similarity_empirical_and_fixed_no_self_loop_km$type=="null_fixed"))

#get equations
lm1_module_equation <- paste("y=", coef(lm1_module)[[1]], "+", coef(lm1_module)[[2]], "*x")
lm2_module_equation <- paste("y=", coef(lm2_module)[[1]], "+", coef(lm2_module)[[2]], "*x")

b1_module <- summary(lm1_module)$coefficients[2,1]
se1_module <- summary(lm1_module)$coefficients[2,2]
b2_module <- summary(lm2_module)$coefficients[2,1]
se2_module <- summary(lm2_module)$coefficients[2,2]

p_value_module = 2*pnorm(-abs(compare.coeff(b1_module,se1_module,
                                                    b2_module,se2_module)))
p_value_module #same pattern

##correlation and r sqaured between jaccard and distance for each run ----------------------------------------------------------------------------
iteration_correlation_data_fixed <- layer_turnover_with_distnace_fixed %>% subset(layer_from != layer_to) 


  trial_fixed = iteration_correlation_data_fixed 
  iteration_correlation_new_fixed <- cor.test(trial_fixed$turnover, trial_fixed$ave_distance, method = "pearson")
  lm_val_fixed <- lm(turnover ~ ave_distance, data = trial_fixed)
  iteration_correlation_fixed <- tibble(estimate = iteration_correlation_new_fixed$estimate, 
                                                                         p_val = iteration_correlation_new_fixed$p.value, 
                                                                         statistic = iteration_correlation_new_fixed$statistic, 
                                                                         confidence_int_low = iteration_correlation_new_fixed$conf.int[1],
                                                                         confidence_int_high = iteration_correlation_new_fixed$conf.int[2],
                                                                         slope = lm_val_fixed$coefficients[2],
                                                                         intercept = lm_val_fixed$coefficients[1],
                                                                         rsquared = summary(lm_val_fixed)$adj.r.squared)



#write.csv(iteration_correlation_fixed, "./csvs/iteration_correlation_fixed.csv", row.names = FALSE)
#iteration_correlation_fixed <- read.csv("./csvs/iteration_correlation_fixed.csv")

#correlation empirical
correlation_empirical <- read.csv("./csvs/correlation_empirical_pols.csv")

# R2 null_fixed = 0.5078
# R2 empirical =0.549

#pdf('./graphs/uniform_interlayer_values/M4_iteration_correlation.pdf', 10, 6)
#rsquared_fixed <- iteration_correlation_fixed %>% 
#  ggplot(aes(x = rsquared))+ 
#  geom_bar(fill = "#c4ad06", color = "#c4ad06", alpha = 0.4)+ 
#  theme_classic()+ labs(x = "R squared")+
#  geom_vline(xintercept = correlation_empirical$rsquared, linetype = "dashed", color = "#F47069")+
#  theme(axis.title=element_text(size=22))+
 # theme(panel.grid = element_blank(),
  #    panel.border = element_rect(color = "black",fill = NA,size = 1),
#        panel.spacing = unit(0.5, "cm", data = NULL),
#        axis.text = element_text(size=14, color='black'),
#        axis.title = element_text(size=14, color='black'),
#        axis.line = element_blank())
#dev.off()
#rsquared_fixed


