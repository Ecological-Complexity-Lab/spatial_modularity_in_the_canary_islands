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
library(ggtree)
library(ggpubr)

setwd("/Users/maya/Desktop/spatial_modularity_in_the_canary_islands")
getwd()

## This portion deletes all species found in the mainland so we can better aseess only the species in the islands

#---- read all csvs needed for analysis -----------------------------------------------------------------------------
dryad_edgelist_complete_ids <- read.csv("csvs/dryad_edgelist_complete_ids.csv")
physical_nodes <- read.csv("csvs/physical_nodes.csv")
layer_metadata <- read.csv("csvs/layer_metadata.csv")
distances_with_ids <- read.csv("csvs/distances_with_ids.csv")
lon_lat_data <- read.csv("csvs/lon_lat_data.csv")


#---- delete mainland and create hierarchical multilayer network --------------------------------------------------------------
only_islands_edgelist <- dryad_edgelist_complete_ids %>% filter(!layer_from %in% c(1,2) & !layer_to %in% c(1,2)) #filter out when layers are mainland

dryad_multilayer_islands <- create_multilayer_object(extended = only_islands_edgelist, #taking edge list and returning multilayer network
                                             nodes = physical_nodes,
                                             layers = layer_metadata,
                                             intra_output_extended = T)

output_multilayer_islands <- multi_lvl_infomap(dryad_multilayer_islands, 
                                       infomap_executable = "Infomap",
                                       flow_model = 'directed',
                                       relax = F, 
                                       silent = T, 
                                       trials = 100,
                                       seed = 497294, 
                                       temporal_network = F)

modules_dryad_multilayer_multi_level_islands <- output_multilayer_islands$modules #has module, sub_module and sub_sub_module but only for sub_module 3 and 17
modules_dryad_multilayer_multi_level_islands <- modules_dryad_multilayer_multi_level_islands %>% 
  drop_na(flow) #remove all rows where species were only in layers 1,2

#this version has 41 modules and fewer sun-modules. only modules 3 and 16 have sub_modules

#write.csv(modules_dryad_multilayer_multi_level_islands, "csvs/modules_dryad_multilayer_multi_level_islands.csv", row.names = FALSE)

##---- dendrogram--------------------------------------------------------------
## hierarchical as dendogram
hierarchical_to_dendogram <- modules_dryad_multilayer_multi_level_islands 
hierarchical_to_dendogram$state_node <- paste(hierarchical_to_dendogram$node_id, "_", hierarchical_to_dendogram$layer_id) #create a column for state nodes

network2module_islands <- data.frame(network = "network",
                             module = seq(1,41,1)) #create 1st split in network to 2 modules


#---- network level 
edges_network2module_level_islands <- network2module_islands %>% select(network, module) %>% 
  unique %>% rename(from = network, to = module) #create edge list from module to sub module

edges_network2module_level_islands$to <- paste(edges_network2module_level_islands$to, "_", "m") #add m to describe module


#---- top level
edges_module2sub_module <- hierarchical_to_dendogram %>% select(module, sub_module) %>% 
  unique %>% rename(from = module, to = sub_module) #create edge list from module to sub module

edges_module2sub_module$from <- paste(edges_module2sub_module$from, "_", "m") #add m to describe module
edges_module2sub_module$to <- paste(edges_module2sub_module$to, "_", "sm") #add sm to describe sub module

module_3 <- edges_module2sub_module %>% filter(from == "3 _ m") 
module_3$to <- paste("3." , module_3$to) #distinguish between sub modules with the same id under different modules

module_16 <- edges_module2sub_module %>% filter(from == "16 _ m") 
module_16$to <- paste("16." , module_16$to) #distinguish between sub modules with the same id under different modules

edges_module2sub_module <- rbind(module_3, module_16) #create new version where sub modules have indication to which module they belong to

# no sub_sub_modules

#---- bottom level option 1 (sm to state)
edges_sub_module2state_level <- hierarchical_to_dendogram %>% select(module, sub_module, state_node) %>% 
  unique %>% rename(from = sub_module, to = state_node) #create edge list from sub module to state node

edges_sub_module2state_level$from <- paste(edges_sub_module2state_level$from, "_", "sm") #add sm to describe sub module

module_3_state <- edges_sub_module2state_level %>% filter(module == 3) 
module_3_state$from <- paste("3." , module_3_state$from) #distinguish between sub modules with the same id under different modules

module_16_state <- edges_sub_module2state_level %>% filter(module == 16) 
module_16_state$from <- paste("16." , module_16_state$from) #distinguish between sub modules with the same id under different modules

edges_sub_module2state_level <- rbind(module_3_state, module_16_state) %>% select(from, to) #create new version where sub modules have indication to which module they belong to


#---- bottom level option 2 (m to state)
edges_module2state_level_islands <- hierarchical_to_dendogram %>% select(module, state_node) %>% 
  unique %>% rename(from = module, to = state_node) #create edge list from sub module to state node

edges_module2state_level_islands <- edges_module2state_level_islands[!(edges_module2state_level_islands$to %in% 
                                                                         edges_sub_module2state_level$to),] ## if there is a sub module 
#need to delete the corresponding edge in the module to state node

edges_module2state_level_islands$from <- paste(edges_module2state_level_islands$from, "_", "m") #add m to describe module



## create edge list combined for dendrogram
edge_list_for_dendrogram_islands <- rbind(edges_network2module_level_islands, edges_module2sub_module, edges_sub_module2state_level,
                                          edges_module2state_level_islands)

no_state_nodes_islands <- rbind(edges_network2module_level_islands, edges_module2sub_module) #version with no state nodes

#---- dendrogram graphs-------------------------------------------------------------------------------

## graphs
pre_graph_islands <- graph_from_data_frame(edge_list_for_dendrogram_islands)

pre_graph_no_state_nodes_islands <- graph_from_data_frame(no_state_nodes_islands)


ggtree(pre_graph_no_state_nodes_islands, layout = "circular") + geom_tiplab()

ggtree(pre_graph_islands, layout = "circular") 



#---- functions suited for only islands -------------------------------------------------------------------
pivot_by_module_islands <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(data, id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s3 = t(s2) 
  s3 <- s3[-1,]
  colnames(s3) <- c(3,4,5,6,7,8,9,10,11,12,13,14) #version without mainland
  return(s3)
}

edge_list_per_module_islands <- function(data,edge_list){
  #gets one row from a data frame and creates an edge list from it
  for (i in (1:11)){
    if (data[i]==0) next #only take layers where the module is present
    else {
      for (j in ((i+1):12)){
        if (data[j]==0) next #only take layers where the module is present
        else {
          edge_list <- rbind(edge_list, tibble(layer_from=i+2, layer_to=j+2, module=NA)) #create edge list of all the layer found in a module
        } # +2 is to make sure we get the layer id (id in the column) and not the number of the column
      }
    }
  }
  return(edge_list)
}

#---- create edge list with distances for modules-----------------------------------------------------------

pivoted_by_module_empirical_hierarchical_islands <- pivot_by_module_islands(modules_dryad_multilayer_multi_level_islands) #pivot modules for map

modules_edge_list_hierarchical_islands <- NULL

for (k in (1:nrow(pivoted_by_module_empirical_hierarchical_islands))){ #run the function for each row in the data frame
  modules_edge_list_hierarchical_islands <- edge_list_per_module_islands(pivoted_by_module_empirical_hierarchical_islands[k,], 
                                                                 modules_edge_list_hierarchical_islands) 
  current_module <- rownames(pivoted_by_module_empirical_hierarchical_islands)[k]
  if (is.null(modules_edge_list_hierarchical_islands)) next
  modules_edge_list_hierarchical_islands <- modules_edge_list_hierarchical_islands %>% mutate(module = replace_na(module, current_module)) #add module number
}

hierarchical_edge_list_with_distances_islands <- right_join(modules_edge_list_hierarchical_islands, distances_with_ids, 
                                                            by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
hierarchical_edge_list_with_distances_islands <- na.omit(hierarchical_edge_list_with_distances_islands) #remove NA and delete layer name


##---- island map 5 biggest sub modules-----------------------------------------------------------------------------------

#modules_dryad_multilayer_multi_level_islands ## edge lit with layer from, layer to, module sub module and distance

#arrange data to include coordinates and modules sizes
module_size_islands <- count(modules_dryad_multilayer_multi_level_islands, module)  #create a data frame of all modules and how many nodes are in each (size of module)
module_data_islands <- merge(modules_dryad_multilayer_multi_level_islands , module_size_islands, 
                         by=c("module","module")) #merge size of module with all the other info about the modules
colnames(module_data_islands)[9] <- "size_of_module" #rename column

module_data_with_loc_islands <- merge(module_data_islands, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a module
modules_with_lat_lon_islands <- module_data_with_loc_islands %>% 
  select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
modules_with_lat_lon_islands$count <- c(1)

#number of physical nodes in module
modules_with_lat_lon_islands %>% 
  ggplot(aes(x= module, y= count ,fill= factor(layer_id)))+ geom_bar(stat= "identity")+ theme_classic()+
  labs(y="Number of Physical Nodes", x="Module Number")+ scale_x_continuous(breaks = seq(1,41,2))+
  guides(fill=guide_legend(title="Layer\nNumber"))

#number of state nodes in module
modules_with_lat_lon_state_nodes_islands <- modules_with_lat_lon_islands %>% select(module, size_of_module) %>% unique()

modules_with_lat_lon_state_nodes_islands %>%
  ggplot(aes(x= module, y= size_of_module))+ geom_bar(stat= "identity")+ theme_classic()+
  labs(y="Number of State Nodes", x="Module Number")+ scale_x_continuous(breaks = seq(1,41,2))

#how many islands are within a module
modules_with_lat_lon_only_islands <- module_data_with_loc_islands %>% 
  select(layer_id, module, lat, Lon, size_of_module) %>% unique() #take only certain columns
old <- c(3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("2","2","3","3","4east","4east","4west","4west","5","5","6","6")
modules_with_lat_lon_only_islands$layer_id[modules_with_lat_lon_only_islands$layer_id %in% old] <- 
  new[match(modules_with_lat_lon_only_islands$layer_id, old)] #change layer number to island number
modules_with_lat_lon_only_islands$count <- c(1)

modules_with_lat_lon_only_islands %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_id)))+ geom_bar(stat= "identity")+ theme_classic()+
  scale_x_continuous(breaks = seq(1,41,2))+ 
  labs(y="Number of Physical Nodes", x="Module Number")+ 
  guides(fill=guide_legend(title="Island\nNumber"))

#---- distance decay for the empirical-----------------------------------------------------------------------------
edge_list_with_distances_hierarchical_islands <- right_join(modules_edge_list_hierarchical_islands, 
                                                            distances_with_ids, by= c("layer_from", "layer_to"))
edge_list_with_distances_hierarchical_islands <- na.omit(edge_list_with_distances_hierarchical_islands) #remove NA and delete layer name

#modules similarity pairwise distance between islands
edge_list_by_islands_hierarchical_islands <- edge_list_with_distances_hierarchical_islands
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
edge_list_by_islands_hierarchical_islands$layer_from[edge_list_by_islands_hierarchical_islands$layer_from %in% old] <- 
  new[match(edge_list_by_islands_hierarchical_islands$layer_from, old)]
edge_list_by_islands_hierarchical_islands$layer_to[edge_list_by_islands_hierarchical_islands$layer_to %in% old] <- 
  new[match(edge_list_by_islands_hierarchical_islands$layer_to, old)]

#version with # of modules in layers
edge_list_by_islands_modules_hierarchical <- edge_list_by_islands_hierarchical_islands %>% group_by(layer_from, layer_to, module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_islands_modules_hierarchical$count <- c(1)
edge_list_by_islands_modules_hierarchical <- edge_list_by_islands_modules_hierarchical %>% mutate(number_of_modules= sum(count)) %>%
  select(layer_from, layer_to, module, number_of_modules)

#version with correct average between layers
edge_list_by_islands_ave_hierarchical <- edge_list_by_islands_hierarchical_islands %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

#combine
edge_list_island_combine_hierarchical <- edge_list_by_islands_ave_hierarchical %>%
  merge(edge_list_by_islands_modules_hierarchical, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_island_combine_no_module_hierarchical<- edge_list_island_combine_hierarchical %>% select(-module) %>% unique() #have version where modules aren't present

#total number of modules in each layer
module_island_turnover_hierarchical_islands <- NULL

only_island_list <- c("2","3","4east","4west","5","6")

for (i in only_island_list){
  for (j in only_island_list){
    modules_in_island_from_hierarchical <- filter(edge_list_by_islands_modules_hierarchical, layer_from == i) %>% 
      select(module) %>% unique() %>% unlist()
    modules_in_island_to_hierarchical <- filter(edge_list_by_islands_modules_hierarchical, layer_from == j) %>% 
      select(module) %>% unique() %>% unlist()
    #take all sub modules in layer_from and all sub modules in layer_to to check turnover
    int_both <- intersect(modules_in_island_from_hierarchical, modules_in_island_to_hierarchical) #how many sub modules are found in both layers
    uni_both <- union(modules_in_island_from_hierarchical, modules_in_island_to_hierarchical)
    turnover <- length(int_both)/length(uni_both)
    module_island_turnover_hierarchical_islands <- rbind(module_island_turnover_hierarchical_islands, 
                                                         tibble(layer_from= i, layer_to= j, turnover= turnover))
  }
}

module_island_turnover_hierarchical_islands <- drop_na(module_island_turnover_hierarchical_islands)

edge_list_by_islands_ave_hierarchical <- edge_list_by_islands_hierarchical_islands %>% group_by(layer_from, layer_to) %>%
  summarise(ave_distance= mean(distance_in_meters)) %>% unique()

islands_turnover_with_distnace_hierarchical <- edge_list_by_islands_ave_hierarchical %>%
  merge(module_island_turnover_hierarchical_islands, by= c("layer_from", "layer_to")) #merge both versions

#---- combining the data for graphs---------------------------------------------------------------------------------------------------
## empirical
empirical_turnover_for_module_hierarchical_islands <- islands_turnover_with_distnace_hierarchical %>% group_by(layer_from, layer_to) %>%
  summarise(ave=mean(turnover), sd=sd(turnover), ave_dist=mean(ave_distance)) %>% mutate(type="empirical") #make sure sd is 0 cause its the empirical and not null

#just empirical
jaccard_similarity_hierarchical_no_self_loop_km <- empirical_turnover_for_module_hierarchical_islands %>% 
  subset(layer_from != layer_to) %>% mutate(ave_dist_in_km = ave_dist/1000)

#---- graphs----------------------------------------------------------------------------------------------------
#just emprical
jaccard_similarity_hierarchical_no_self_loop_km %>% ggplot(aes(x= ave_dist_in_km, y= ave))+
  geom_point(color = "#1EA784")+ theme_classic()+ geom_smooth(method= "lm", se=F, color = "#1EA784")+
  theme(axis.title=element_text(size=22))+theme(axis.text.x=element_text(size=15))+
  theme(axis.text.y=element_text(size=15))+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))+
  labs(x="Distance in Km", y="Jaccard Similarity")+ stat_cor(aes(label = ..p.label..), label.x = 400, color = "#1EA784")+
  stat_cor(aes(label = ..rr.label..), label.x = 400, label.y = c(0.65), color = "#1EA784")


