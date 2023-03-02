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


dryad_edgelist_complete_ids <- read.csv("./csvs/dryad_edgelist_complete_ids.csv")
physical_nodes <- read.csv("./csvs/physical_nodes.csv")
layer_metadata <- read.csv("./csvs/layer_metadata.csv")

#---- analysis only on islands--------------------------------------
dryad_edgelist_complete_ids_islands <- dryad_edgelist_complete_ids %>%
  subset(layer_from > 2) #only islands with no mainland

# Input: An extended edge list.
islands_multilayer <- create_multilayer_object(extended = dryad_edgelist_complete_ids_islands, #taking edge list and returning multilayer network
                                               nodes = physical_nodes,
                                               layers = layer_metadata,
                                               intra_output_extended = T)

#create modules for empirical network
modules_islands_multilayer <- modified_multi(islands_multilayer, 
                                                     infomap_executable = "Infomap",
                                                     flow_model = 'directed',
                                                     relax = F, 
                                                     silent = T, 
                                                     trials = 100,
                                                     seed = 497294, 
                                                     temporal_network = F)

islands_network <- modules_islands_multilayer$modules %>% na.omit() #delete NAs as some nodes found in physical_nodes are local to the mainland

#write_csv(islands_network, './csvs/islands_network.csv')


modules_in_islands <- islands_network

old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
modules_in_islands$layer_id[modules_in_islands$layer_id %in% old] <- 
  new[match(modules_in_islands$layer_id, old)]

#number of module which are found in x islands
modules_in_only_islands <- NULL
for (i in 1:6){
  modules_for_similarity_in_only_islands_num <- modules_in_islands %>% select(module, layer_id) %>%
    unique() %>% group_by(module) %>% filter(n()==i) %>% select(module) %>% unique() #check which modules occur in i or more layers
  modules_for_similarity_in_only_islands <- modules_in_islands %>%
    filter(module %in% modules_for_similarity_in_only_islands_num$module) %>% count(module)#only save the modules that are found in i or more layers
  modules_in_only_islands <- rbind(modules_in_only_islands, tibble(module= modules_for_similarity_in_only_islands$module, i=i))
}

num_of_modules_in_only_islands <- modules_in_only_islands %>% count(i)
num_of_modules_in_only_islands %>% ggplot(aes(x=i, y=n))+geom_bar(stat="identity")+theme_classic()+
  scale_x_continuous(breaks=seq(1,6,1))+ scale_y_continuous(breaks=seq(1,10,1))+
  labs(x="number of islands", y="Number of Modules")

#number of modules found in each island
islands_names_network <- islands_network %>% select(module,layer_id) #create to replace with island names later

island_names <- c("Western Sahara","Western Sahara", #change layers to island names
                  "Fuerteventura","Fuerteventura",
                  "GranCanaria","GranCanaria",
                  "Tenerife South","Tenerife South",
                  "Tenerife Teno","Tenerife Teno",
                  "Gomera","Gomera",
                  "Hierro","Hierro")

islands_names_network$layer_id[islands_names_network$layer_id %in% old] <- 
  island_names[match(islands_names_network$layer_id, old)]

islands_names_network <- islands_names_network %>% unique() #delete doubles caused by layers turning to islands

modules_per_island_no_mainland <- islands_names_network %>% group_by(layer_id) %>% count() #number of modules found per island

modules_per_island_no_mainland %>% ggplot(aes(x=layer_id, y=n))+geom_bar(stat="identity")+theme_classic()+
  scale_x_discrete()+ 
  labs(x="Island Name", y="Number of Modules in Island")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.55, hjust = 0.4))





