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

#this portion of code creates a hierarchial modularity version of the network, alongside modularity analysis.

#---- hierarchical modularity is run-----------------------------------------------------------------------------------------------------


output_multilayer <- multi_lvl_infomap(dryad_multilayer, 
                                       infomap_executable = "../Infomap",
                                       flow_model = 'directed',
                                       relax = F, 
                                       silent = T, 
                                       trials = 100,
                                       seed = 497294, 
                                       temporal_network = F)

modules_dryad_multilayer_multi_level <- output_multilayer$modules #has module, sub_module and sub_sub_module but only for sub_module 3 and 17


#write.csv(modules_dryad_multilayer_multi_level, "csvs/modules_dryad_multilayer_multi_level.csv")
#modules_dryad_multilayer_multi_level <- read.csv("csvs/modules_dryad_multilayer_multi_level.csv")

##---- for each layer how many nodes are found in each module-------------------------------------------------------

modules_dryad_multilayer_multi_level_graph <- modules_dryad_multilayer_multi_level
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("Western Sahara","Western Sahara", #change layers to islands
         "Fuerteventura","Fuerteventura",
         "GranCanaria","GranCanaria",
         "Tenerife South","Tenerife South",
         "Tenerife Teno","Tenerife Teno",
         "Gomera","Gomera",
         "Hierro","Hierro")
modules_dryad_multilayer_multi_level_graph$layer_id[modules_dryad_multilayer_multi_level_graph$layer_id %in% old] <- 
  new[match(modules_dryad_multilayer_multi_level_graph$layer_id, old)]



##---- graphs hierarchical modularity levels number of species in different levels--------------------------------------------

## modules
modules_dryad_multilayer_multi_level_graph %>%
  group_by(layer_id, module) %>%
  summarise(n=n_distinct(node_id)) %>%
  ggplot(aes(layer_id, module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text(color = "white")+ xlab("Island Name") + ylab("Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text())+ 
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_y_continuous(breaks=seq(1,2,1))+
  scale_fill_viridis()

## sub_modules module 1
modules_dryad_multilayer_multi_level_graph %>%
  filter(module == 1) %>%
  group_by(layer_id, sub_module) %>%
  summarise(n=n_distinct(node_id)) %>%
  filter(n>2) %>%
  ggplot(aes(layer_id, sub_module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text(color = "white")+ xlab("Island Name") + ylab("Sub Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_y_continuous(breaks=seq(1,35,3))+
  scale_fill_viridis()

## sub_modules module 2
modules_dryad_multilayer_multi_level_graph %>%
  filter(module == 2) %>%
  group_by(layer_id, sub_module) %>%
  summarise(n=n_distinct(node_id)) %>%
  filter(n>2) %>%
  ggplot(aes(layer_id, sub_module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text(color = "white")+ xlab("Island Name") + ylab("Sub Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ 
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_y_continuous(breaks=seq(1,35,1))+
  scale_fill_viridis()

## sub_sub_modules
modules_dryad_multilayer_multi_level_graph %>%
  filter(module == 1) %>%
  filter(sub_module == 3 | sub_module == 17) %>%
  group_by(layer_id, sub_sub_module) %>%
  summarise(n=n_distinct(node_id)) %>%
  ggplot(aes(layer_id, sub_sub_module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text(color = "white")+ xlab("Island Name") + ylab("Sub Sub Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ 
  scale_x_discrete(guide = guide_axis(angle = 45))+
  scale_y_continuous(breaks=seq(1,35,1))+
  scale_fill_viridis()


####
##---- dendrogram--------------------------------------------------------------
####

## multi lvl as dendogram
multi_lvl_to_dendogram <- modules_dryad_multilayer_multi_level 
multi_lvl_to_dendogram$state_node <- paste(multi_lvl_to_dendogram$node_id, "_", multi_lvl_to_dendogram$layer_id) #create a column for state nodes

network2module <- data.frame(network = c("network", "network"),
                             module = c(1,2)) #create 1st split in network to 2 modules


#----network level 
edges_network2module_level <- network2module %>% select(network, module) %>% 
  unique %>% rename(from = network, to = module) #create edge list from module to sub module

edges_network2module_level$to <- paste(edges_network2module_level$to, "_", "m") #add m to describe module


#----top level
edges_module2sub_level <- multi_lvl_to_dendogram %>% select(module, sub_module) %>% 
  unique %>% rename(from = module, to = sub_module) #create edge list from module to sub module

edges_module2sub_level$from <- paste(edges_module2sub_level$from, "_", "m") #add m to describe module
edges_module2sub_level$to <- paste(edges_module2sub_level$to, "_", "sm") #add sm to describe sub module

module_1 <- edges_module2sub_level %>% filter(from == "1 _ m") 
module_1$to <- paste("1." , module_1$to) #distinguish between sub modules with the same id under different modules

module_2 <- edges_module2sub_level %>% filter(from == "2 _ m") 
module_2$to <- paste("2." , module_2$to) #distinguish between sub modules with the same id under different modules

edges_module2sub_level <- rbind(module_1, module_2) #create new version where sub modules have indication to which module they belong to


#----middle level
edges_sub2subsub_level <- multi_lvl_to_dendogram %>% select(sub_module, sub_sub_module) %>% 
  unique %>% rename(from = sub_module, to = sub_sub_module) #create edge list from sub module to sub sub module

edges_sub2subsub_level <- drop_na(edges_sub2subsub_level) #delete all occations where there's no sub sub module

edges_sub2subsub_level$from <- paste("1." , edges_sub2subsub_level$from, "_", "sm") #add sm to describe sub module
edges_sub2subsub_level$to <- paste(edges_sub2subsub_level$to, "_", "ssm") #add ssm to describe sub sub module

sub_3 <- edges_sub2subsub_level %>% filter(from == "1. 3 _ sm")
sub_3$to <- paste("1.3." , sub_3$to)

sub_17 <- edges_sub2subsub_level %>% filter(from == "1. 17 _ sm")
sub_17$to <- paste("1.17." , sub_17$to)

edges_sub2subsub_level <- rbind(sub_3, sub_17) #create new version where sub modules have indication to which module they belong to


#----bottom level option 1 (ssm)
edges_subsub2state_level <- multi_lvl_to_dendogram %>% select(sub_module, sub_sub_module, state_node) %>% 
  unique %>% rename(from = sub_sub_module, to = state_node) #create edge list from sub sub module to state node

edges_subsub2state_level <- drop_na(edges_subsub2state_level) ## clean up to make sure there are no NAs

edges_subsub2state_level$from <- paste(edges_subsub2state_level$from, "_", "ssm") #add ssm to describe sub sub module

sub_3_state <- edges_subsub2state_level %>% filter(sub_module == 3)
sub_3_state$from <- paste("1.3." , sub_3_state$from)

sub_17_state <- edges_subsub2state_level %>% filter(sub_module == 17)
sub_17_state$from <- paste("1.17." , sub_17_state$from)

edges_subsub2state_level <- rbind(sub_3_state, sub_17_state) %>% select(from, to)


#----bottom level option 2 (sm)
edges_sub2state_level <- multi_lvl_to_dendogram %>% select(module, sub_module, state_node) %>% 
  unique %>% rename(from = sub_module, to = state_node) #create edge list from sub module to state node

edges_sub2state_level <- edges_sub2state_level[!(edges_sub2state_level$to %in% edges_subsub2state_level$to),] ## if there is a sub sub module 
#need to delete the corresponding edge in the sub module to state node

edges_sub2state_level$from <- paste(edges_sub2state_level$from, "_", "sm") #add sm to describe sub module

module_1_state <- edges_sub2state_level %>% filter(module == 1) 
module_1_state$from <- paste("1." , module_1_state$from) #distinguish between sub modules with the same id under different modules

module_2_state <- edges_sub2state_level %>% filter(module == 2) 
module_2_state$from <- paste("2." , module_2_state$from) #distinguish between sub modules with the same id under different modules

edges_sub2state_level <- rbind(module_1_state, module_2_state) %>% select(from, to) #create new version where sub modules have indication to which module they belong to


#---- create edge list combined for dendrogram---------------------------------------------------------------
edge_list_for_dendrogram <- rbind(edges_network2module_level, edges_module2sub_level, edges_sub2subsub_level,
                                  edges_subsub2state_level, edges_sub2state_level)

no_state_nodes <- rbind(edges_network2module_level, edges_module2sub_level, edges_sub2subsub_level) #version with no state nodes


#---- graphs dendrogram--------------------------------------------------------------------------------------------------
pre_graph <- graph_from_data_frame(edge_list_for_dendrogram)

pre_graph_no_state_nodes <- graph_from_data_frame(no_state_nodes,)

ggtree(pre_graph_no_state_nodes, layout = "circular") +geom_tiplab() + 
  geom_hilight(node = c(1:45), fill = "aquamarine3") + geom_hilight(node = c(46:54), fill = "thistle3")

ggtree(pre_graph, layout = "circular") 



#---- create edge list with distances for modules-----------------------------------------------------------

#modules_dryad_multilayer_multi_level_fixed_analysis_sub_module <- 
#  read.csv("csvs/modules_dryad_multilayer_multi_level_fixed_analysis_sub_module.csv")

pivoted_by_module_empirical_multi_lvl <- pivot_by_module(modules_dryad_multilayer_multi_level_fixed_analysis_sub_module) #pivot modules for map

modules_edge_list_multi_lvl <- NULL

for (k in (1:nrow(pivoted_by_module_empirical_multi_lvl))){ #run the function for each row in the data frame
  modules_edge_list_multi_lvl <- edge_list_per_module(pivoted_by_module_empirical_multi_lvl[k,], modules_edge_list_multi_lvl) 
  current_module <- rownames(pivoted_by_module_empirical_multi_lvl)[k]
  if (is.null(modules_edge_list_multi_lvl)) next
  modules_edge_list_multi_lvl <- modules_edge_list_multi_lvl %>% mutate(module = replace_na(module, current_module)) #add module number
}

multi_lvl_edge_list_with_distances <- right_join(modules_edge_list_multi_lvl, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
multi_lvl_edge_list_with_distances <- na.omit(multi_lvl_edge_list_with_distances) #remove NA and delete layer name

#----arrange data to include coordinates and modules sizes--------------------------------------------------------------


module_from_multi_lvl <- multi_lvl_edge_list_with_distances %>% select(layer_from, module) %>% unique()
colnames(module_from_multi_lvl)[1] <- "layer_id" #rename column

module_to_multi_lvl <- multi_lvl_edge_list_with_distances %>% select(layer_to, module) %>% unique()
colnames(module_to_multi_lvl)[1] <- "layer_id" #rename column

module_multi_lvl <- rbind(module_from_multi_lvl, module_to_multi_lvl) %>% unique()
module_multi_lvl_with_loc <- merge(module_multi_lvl, lon_lat_data, by= c("layer_id","layer_id"))

module_1_polygon <- module_multi_lvl_with_loc %>% filter(module == 1)
module_2_polygon <- module_multi_lvl_with_loc %>% filter(module == 2)

#edge_list_with_distances_multi_lvl_empirical ## module with sub module

worldmap <- map_data("world")
dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990))  
dryad_map_multi_lvl <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+
  geom_polygon(data = module_1_polygon, aes(x = Lon, y = lat), alpha = 0.3, color = "darkseagreen3", fill = "darkseagreen2")+
  geom_point(data = module_1_polygon, aes(x = Lon, y = lat), color = "darkseagreen3")+
  geom_point(data = module_2_polygon, aes(x = Lon, y = lat), color = "lightpink3")+
  theme(axis.title=element_text(size=14))+theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  labs(x = "Longitude", y = "Latitude")

dryad_map_multi_lvl

##---- map but with sub modules-----------------------------------------------------------------------------------------


sub_module_from_multi_lvl <- edge_list_with_distances_multi_lvl_empirical %>% select(layer_from, module_sub_module) %>% unique()
colnames(sub_module_from_multi_lvl)[1] <- "layer_id" #rename column

sub_module_to_multi_lvl <- edge_list_with_distances_multi_lvl_empirical %>% select(layer_to, module_sub_module) %>% unique()
colnames(sub_module_to_multi_lvl)[1] <- "layer_id" #rename column

sub_module_multi_lvl <- rbind(sub_module_from_multi_lvl, sub_module_to_multi_lvl) %>% unique()
sub_module_multi_lvl_with_loc <- merge(sub_module_multi_lvl, lon_lat_data, by= c("layer_id","layer_id")) %>% 
  group_by(module_sub_module) %>% filter(n()>5)

dryad_map_multi_lvl_sub_modules <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+
  geom_polygon(data = sub_module_multi_lvl_with_loc, aes(x = Lon, y = lat, color = module_sub_module, fill = module_sub_module), 
               alpha = 0.1)+
  theme(axis.title=element_text(size=14))+theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+
  labs(x = "Longitude", y = "Latitude")



##---- island map 5 biggest sub modules-----------------------------------------------------------------------------------


#edge_list_with_distances_multi_lvl_empirical ## edge lit with layer from, layer to, module sub module and distance

modules_dryad_multilayer_multi_level_sub_modules <- modules_dryad_multilayer_multi_level
modules_dryad_multilayer_multi_level_sub_modules$module_sub_module <- paste(modules_dryad_multilayer_multi_level_sub_modules$module, ".", 
                                                                            modules_dryad_multilayer_multi_level_sub_modules$sub_module)

#arrange data to include coordinates and modules sizes
sub_module_size <- count(modules_dryad_multilayer_multi_level_sub_modules, module_sub_module)  #create a data frame of all modules and how many nodes are in each (size of module)
sub_module_data <- merge(modules_dryad_multilayer_multi_level_sub_modules , sub_module_size, 
                         by=c("module_sub_module","module_sub_module")) #merge size of module with all the other info about the modules
colnames(sub_module_data)[10] <- "size_of_sub_module" #rename column

sub_module_data_with_loc <- merge(sub_module_data, lon_lat_data, by= c("layer_id","layer_id")) #merge modules with module size with the coordinates

#how many layers are within a sub module
sub_modules_with_lat_lon <- sub_module_data_with_loc %>% 
  select(layer_id, module_sub_module, lat, Lon, size_of_sub_module) %>% unique() #take only certain columns
sub_modules_with_lat_lon$count <- c(1)

#number of physical nodes in sub module
sub_modules_with_lat_lon %>% 
  ggplot(aes(x= module_sub_module, y= count ,fill= factor(layer_id)))+ geom_bar(stat= "identity")+ theme_classic()+
  labs(y="Number of Physical Nodes", x="Sub Module Number")+ scale_x_discrete(guide = guide_axis(angle = 90))+
  guides(fill=guide_legend(title="Layer\nNumber"))

#number of state nodes in sub module
sub_modules_with_lat_lon_state_nodes <- sub_modules_with_lat_lon %>% select(module_sub_module, size_of_sub_module) %>% unique()

sub_modules_with_lat_lon_state_nodes %>%
  ggplot(aes(x= module_sub_module, y= size_of_sub_module))+ geom_bar(stat= "identity")+ theme_classic()+
  labs(y="Number of State Nodes", x="Sub Module Number")+ scale_x_discrete(guide = guide_axis(angle = 90))

#how many islands are within a module
sub_modules_with_lat_lon_islands <- sub_module_data_with_loc %>% 
  select(layer_id, module_sub_module, lat, Lon, size_of_sub_module) %>% unique() #take only certain columns
old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
#new <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
sub_modules_with_lat_lon_islands$layer_id[sub_modules_with_lat_lon_islands$layer_id %in% old] <- 
  new[match(sub_modules_with_lat_lon_islands$layer_id, old)] #change layer number to island number
sub_modules_with_lat_lon_islands$count <- c(1)

sub_modules_with_lat_lon_islands %>% 
  ggplot(aes(x=module_sub_module, y= count ,fill= factor(layer_id)))+ geom_bar(stat= "identity")+ theme_classic()+
  scale_x_discrete(guide = guide_axis(angle = 90))+ 
  labs(y="Number of Physical Nodes", x="Sub Module Number")+ 
  guides(fill=guide_legend(title="Island\nNumber"))

#modules similarity pairwise distance between islands
#edge_list_by_islands_multi_lvl_empirical #edge list with distances by island

#----version with number of modules in layers-----------------------------------------------------------
edge_list_by_islands_sub_modules <- edge_list_by_islands_multi_lvl_empirical %>% group_by(layer_from, layer_to, module_sub_module) %>%
  summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
edge_list_by_islands_sub_modules$count <- c(1)
edge_list_by_islands_sub_modules <- edge_list_by_islands_sub_modules %>% mutate(number_of_sub_modules= sum(count)) %>%
  select(layer_from, layer_to, module_sub_module, number_of_sub_modules) 

#version with correct average between layers
#edge_list_by_islands_ave_multi_lvl_empirical #edge list with correct distances between layers

#combine
edge_list_island_combine_multi_lvl_empirical_sub_modules <- edge_list_by_islands_ave_multi_lvl_empirical %>%
  merge(edge_list_by_islands_sub_modules, by= c("layer_from", "layer_to")) #merge both versions 

edge_list_island_combine_no_sub_module <- edge_list_island_combine_multi_lvl_empirical_sub_modules %>% 
  select(-module_sub_module) %>% unique() #have version where modules aren't present

##create on map
lat_lon_for_graph_sub_modules <- sub_modules_with_lat_lon 
lat_lon_for_graph_sub_modules <- lat_lon_for_graph_sub_modules %>% select(layer_id, lat, Lon)

old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
lat_lon_for_graph_sub_modules$layer_id[lat_lon_for_graph_sub_modules$layer_id %in% old] <- 
  new[match(lat_lon_for_graph_sub_modules$layer_id, old)]

islands_with_lat_lon_sub_module <- edge_list_island_combine_no_sub_module %>%
  inner_join(lat_lon_for_graph_sub_modules, by= c("layer_from" = "layer_id")) %>%  #add edge list to coordinates
  rename(x= Lon, y=lat) %>%
  inner_join(lat_lon_for_graph_sub_modules, by= c("layer_to" = "layer_id")) %>%
  rename(xend= Lon, yend= lat) %>%
  select(layer_from, layer_to, number_of_sub_modules, y, yend, x, xend) %>% unique() 

lat_lon_nodes_sub_module <- islands_with_lat_lon_sub_module %>% 
  group_by(layer_from) %>% arrange(desc(layer_to)) %>% slice(1) #save coordinates of islands

islands_with_lon_lat_dif_sub_modules <- islands_with_lat_lon_sub_module %>% filter(layer_from != layer_to)

worldmap <- map_data("world")
dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990))  
dryad_map_sub_modules <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+
  geom_segment(aes(x= x, xend= xend, y= y, yend= yend, color = number_of_sub_modules),
               data= islands_with_lon_lat_dif_sub_modules)+ scale_color_gradient(high="red",low="lightskyblue1")+
  geom_point(aes(x=x, y=y),shape= 21 , fill= 'white', color= "black", data= lat_lon_nodes_sub_module)+ 
  theme(axis.title=element_text(size=15))+theme(axis.text.x=element_text(size=13))+
  theme(axis.text.y=element_text(size=13))+ theme(legend.title = element_text(size = 10), legend.text = element_text(size = 10))+ 
  labs(x = "Longitude", y = "Latitude")+ labs(color = "Number of Sub Modules \n in Common")

#dryad_map_sub_modules

#---- biggest sub modules in the network-------------------------------------------------------------
top_10_percent_sub_modules <- select(sub_module_data_with_loc, c("module_sub_module", "size_of_sub_module")) %>% distinct %>% 
  arrange(desc(size_of_sub_module)) %>%
  filter(quantile(size_of_sub_module, probs=0.90) < size_of_sub_module) #arrange all modules to see 90 percentile and filter by them

names(top_10_percent_sub_modules)[names(top_10_percent_sub_modules)=="size_of_sub_module"] <- "group" #change name to group- only the top percentile
pie_chart_data_sub_module <- top_10_percent_sub_modules %>% left_join(sub_module_data_with_loc, 
                                                                      by=c("module_sub_module" = "module_sub_module")) #join to have coordinates as well


pie_chart_data_sub_module_output <- lon_lat_data %>% right_join(pivot_by_island_sub_module(pie_chart_data_sub_module), 
                                                                by = c("layer_id" = "layer_id")) 

pie_chart_data_sub_module_output$layer_id <- as.character(pie_chart_data_sub_module_output$layer_id) #change layer id to not summarise it

pie_chart_data_sub_module_output <- pie_chart_data_sub_module_output %>% group_by(lat, Lon) %>% 
  summarise(across(where(is.numeric), sum)) #group by island

pie_chart_data_sub_module_output$layer_id <- c("1", "2", "3", "4east" , "4west" , "5", "6") #add island name


worldmap <- map_data("world")
dryad_location <- make_bbox(lon= c(-18.542076, -12.58351), lat= c(26, 30.323990))  
inset_map_sub_modules <- get_map(location=dryad_location, zoom=10, maptype="terrain") %>% ggmap()+ #create a terrain map based on dryad
  geom_scatterpie(aes(x=Lon, y=lat, group=layer_id, r=0.028), alpha=0.6, data= pie_chart_data_sub_module_output, 
                  cols=colnames(pie_chart_data_sub_module_output[,c(3:7)]))+ #create a pie chart for every location with the 90 percentile
  coord_fixed(xlim=c(-18.542076, -12.58351), ylim=c(26, 30.323990)) #set the coordinates to the location of dryad
inset_map_sub_modules <- inset_map_sub_modules + guides(fill=guide_legend(title="Sub Module\nNumber"))

#just scatterpies
top_10_scatterpies_sub_modules <- ggplot()+ geom_scatterpie(aes(x=Lon, y=lat, group=layer_id, r=0.2), alpha=0.6, 
                                                            data= pie_chart_data_sub_module_output, 
                                                            cols=colnames(pie_chart_data_sub_module_output[,c(3:7)]))+
  theme_void()+ coord_fixed()+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))

top_10_scatterpies_no_legend_sub_module <- top_10_scatterpies_sub_modules + theme(legend.position = "none")

top_10_scatterpies_no_legend_sub_module


#----number of islands a sub-module is found in statistics-------------------------------------------------------------------

module_pivoted_empirical_multi_lvl_by_island <- module_pivoted_empirical_multi_lvl %>% as.data.frame()

sub_module_and_island <- NULL

for(i in 1:nrow(module_pivoted_empirical_multi_lvl_by_island)){ 
  current_row <- module_pivoted_empirical_multi_lvl_by_island[i,]
  for(j in 1:14){
    current_layer <- current_row[,j]
    if (current_layer > 0){
      sub_module_and_island <- rbind(sub_module_and_island, tibble(sub_module = row.names(current_row), 
                                                                   layer = colnames(current_row[j]))) 
    } #create data frame of sub-modules with layers they are found in
  }
}

old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("Western Sahara","Western Sahara", #change layers to islands
         "Fuerteventura","Fuerteventura",
         "GranCanaria","GranCanaria",
         "Tenerife South","Tenerife South",
         "Tenerife Teno","Tenerife Teno",
         "Gomera","Gomera",
         "Hierro","Hierro")
sub_module_and_island$layer[sub_module_and_island$layer %in% old] <- new[match(sub_module_and_island$layer, old)] #replace id with name

sub_module_and_island_unique <- sub_module_and_island %>% unique() #merge where two names are identical (two consecutive ids)
sub_module_and_island_unique_count <- sub_module_and_island_unique %>% group_by(sub_module) %>% count() #count # of islands sub-module was in

mean(sub_module_and_island_unique_count$n)
