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

#this portion of the code comes to check whether species found in all 12 layers affect the modularity
#when it comes to hierarchical modulairty

#----hierarchical modularity on species not found in all 14 layers--------------------------------------------------------


output_multilayer_del_14 <- multi_lvl_infomap(dryad_multilayer_del_14, 
                                              infomap_executable = "../Infomap",
                                              flow_model = 'directed',
                                              relax = F, 
                                              silent = T, 
                                              trials = 100,
                                              seed = 497294, 
                                              temporal_network = F)

modules_dryad_multilayer_multi_level_del_14 <- output_multilayer_del_14$modules #has module, sub_module and sub_sub_module but only for sub_module 20 and 47

#write.csv(modules_dryad_multilayer_multi_level_del_14, "csvs/modules_dryad_multilayer_multi_level_del_14.csv")
#modules_dryad_multilayer_multi_level_del_14 <- read.csv("csvs/modules_dryad_multilayer_multi_level_del_14.csv")

#---- for each layer how many nodes are found in each module-----------------------------------------------

## modules
modules_dryad_multilayer_multi_level_del_14 %>%
  group_by(layer_id, module) %>%
  summarise(n=n_distinct(node_id)) %>%
  ggplot(aes(layer_id, module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text()+ xlab("Layer Id") + ylab("Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ scale_x_continuous(breaks=seq(1,14,1))+
  scale_y_continuous(breaks=seq(1,2,1))+
  scale_fill_gradient(low= "lightskyblue1", high= "red")

## sub_modules module 1
modules_dryad_multilayer_multi_level_del_14 %>%
  filter(module == 1) %>%
  group_by(layer_id, sub_module) %>%
  summarise(n=n_distinct(node_id)) %>%
  filter(n>2) %>%
  ggplot(aes(layer_id, sub_module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text()+ xlab("Layer Id") + ylab("Sub Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ scale_x_continuous(breaks=seq(1,14,1))+
  scale_y_continuous(breaks=seq(1,64,3))+
  scale_fill_gradient(low= "lightskyblue1", high= "red")

## sub_modules module 2
modules_dryad_multilayer_multi_level_del_14 %>%
  filter(module == 2) %>%
  group_by(layer_id, sub_module) %>%
  summarise(n=n_distinct(node_id)) %>%
  filter(n>1) %>%
  ggplot(aes(layer_id, sub_module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text()+ xlab("Layer Id") + ylab("Sub Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ scale_x_continuous(breaks=seq(1,14,1))+
  scale_y_continuous(breaks=seq(1,35,1))+
  scale_fill_gradient(low= "lightskyblue1", high= "red")

## sub_sub_modules
modules_dryad_multilayer_multi_level_del_14 %>%
  filter(module == 1) %>%
  filter(sub_module == 20 | sub_module == 47) %>%
  group_by(layer_id, sub_sub_module) %>%
  summarise(n=n_distinct(node_id)) %>%
  ggplot(aes(layer_id, sub_sub_module, fill=n, label=n))+geom_tile()+ theme_classic()+
  geom_text()+ xlab("Layer Id") + ylab("Sub Sub Module Number")+
  theme(axis.text = element_text(size=13),
        axis.text.x = element_text(hjust=1))+ scale_x_continuous(breaks=seq(1,14,1))+
  scale_y_continuous(breaks=seq(1,35,1))+
  scale_fill_gradient(low= "lightskyblue1", high= "red")

#----create shuffled versions-------------------------------------------------------------------------------------------
dryad_multilayer_shuf_1000_pols_multi_lvl_del_14 <- NULL
dryad_multilayer_shuf_1000_plants_multi_lvl_del_14 <- NULL
dryad_multilayer_shuf_1000_both_multi_lvl_del_14 <- NULL

#pols
dryad_multilayer_shuf_1000_pols_output_multi_lvl_del_14 <- modularity_for_shuf_multi_lvl(dryad_edgelist_complete_shuf_pols_del, 
                                                                                         dryad_multilayer_shuf_1000_pols_multi_lvl_del_14)

#plants
dryad_multilayer_shuf_1000_plants_output_multi_lvl_del_14 <- modularity_for_shuf_multi_lvl(dryad_edgelist_complete_shuf_plants_del, 
                                                                                           dryad_multilayer_shuf_1000_plants_multi_lvl_del_14)


#both
dryad_multilayer_shuf_1000_both_output_multi_lvl_del_14 <- modularity_for_shuf_multi_lvl(dryad_edgelist_complete_shuf_both_del, 
                                                                                         dryad_multilayer_shuf_1000_both_multi_lvl_del_14)

write.csv(dryad_multilayer_shuf_1000_pols_output_multi_lvl_del_14, "./csvs/dryad_multilayer_shuf_1000_pols_output_multi_lvl_del_14.csv", row.names = FALSE)
write.csv(dryad_multilayer_shuf_1000_plants_output_multi_lvl_del_14, "./csvs/dryad_multilayer_shuf_1000_plants_output_multi_lvl_del_14.csv", row.names = FALSE)
write.csv(dryad_multilayer_shuf_1000_both_output_multi_lvl_del_14, "./csvs/dryad_multilayer_shuf_1000_both_output_multi_lvl_del_14.csv", row.names = FALSE)

