####-----Modularity plots

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
library(sf)
library(remotes)
library(tmap)

##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")

##----upload data frame islands -------
modules_dryad_multilayer <- read.csv("./csvs/Islands/modules_in_network_islands_as_layers.csv") 

#general info about modules
num_of_nodes_in_module <- modules_dryad_multilayer %>% count(module) #num of nodes in module. biggest module has 44 nodes and smallest has 2
local_modules <- modules_dryad_multilayer %>% select(module, layer_id)
num_of_layers_in_module <- local_modules %>% distinct() %>% count(module) #num of layers a module is found in
mean_num_of_layers_in_module <- mean(num_of_layers_in_module$n) #average number of layers (islands) a module is found in is 3.193

#number of module which are found in x islands
modules_with_lat_lon <- read.csv("csvs/Islands/modules_with_lat_lon_islands_as_layers.csv") 

modules_in_i_islands <- NULL
for (i in 1:7){
  modules_for_similarity_in_i_islands_num <- modules_with_lat_lon %>% select(module, layer_id) %>%
    unique() %>% group_by(module) %>% filter(n()==i) %>% select(module) %>% unique() #check which modules occur in i or more layers
  modules_for_similarity_in_i_islands <- modules_with_lat_lon %>%
    filter(module %in% modules_for_similarity_in_i_islands_num$module) %>% count(module)#only save the modules that are found in i or more layers
  modules_in_i_islands <- rbind(modules_in_i_islands, tibble(module= modules_for_similarity_in_i_islands$module, i=i))
}

num_of_modules_in_i_islands <- modules_in_i_islands %>% count(i)
num_of_modules_in_i_islands %>% ggplot(aes(x=i, y=n))+geom_bar(stat="identity")+theme_classic()+
  scale_x_continuous(breaks=seq(1,7,1))+ 
  labs(x="Number of Spatial Locations", y="Number of Modules")

#number of modules found in each island with names
modules_with_lat_lon$layer_id<- as.character(modules_with_lat_lon$layer_id)

modules_with_lat_lon_id <- modules_with_lat_lon %>% 
  mutate(layer_name= case_when(layer_id == '1' ~ 'WesternSahara',
                                layer_id == '2' ~ 'Fuerteventura',
                                layer_id == '3' ~ 'GranCanaria',
                                layer_id == '4' ~ 'TenerifeSouth',
                                layer_id == '5' ~ 'TenerifeTeno',
                                layer_id == '6' ~ 'Gomera',
                                layer_id == '7' ~ 'Hierro'))

modules_with_lat_lon_id <- modules_with_lat_lon_id %>% unique() #delete doubles caused by layers turning to islands

modules_per_island <- modules_with_lat_lon_id %>% group_by(layer_name) %>% count() #number of modules found per island

modules_per_island %>% ggplot(aes(x=layer_name, y=n))+geom_bar(stat="identity")+theme_classic()+
  scale_x_discrete()+ 
  labs(x="Island Name", y="Number of Modules in Island")+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.55, hjust = 0.4))


#how many islands are within a module
islands_in_modules<-modules_with_lat_lon_id %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_name)))+ geom_bar(stat= "identity")+ 
  scale_x_continuous(breaks=seq(1,88,5))+ labs(y="Number of Physical Nodes", x="Module Number")+
  guides(fill=guide_legend(title="Layer\nNumber"))+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())

modules_in_islands <- modules_with_lat_lon_id %>% 
  ggplot(aes(x=layer_name, y= module))+ geom_bar(stat= "identity")+ 
  scale_x_continuous(breaks=seq(1,7,5))+ labs(y="Number of Physical Nodes", x="Module Number")+
  guides(fill=guide_legend(title="Layer\nNumber"))+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())


#-----------Final Plots -------------

#Distribution of modules in islands
#modules_with_lat_lon <- read.csv("csvs/Islands/modules_with_lat_lon_islands_as_layers.csv") 

modules_with_lat_lon$layer_id<- as.character(modules_with_lat_lon$layer_id)

modules_with_lat_lon_id <- modules_with_lat_lon %>% 
  mutate(layer_name= case_when(layer_id == '1' ~ 'WesternSahara',
                               layer_id == '2' ~ 'Fuerteventura',
                               layer_id == '3' ~ 'GranCanaria',
                               layer_id == '4' ~ 'TenerifeSouth',
                               layer_id == '5' ~ 'TenerifeTeno',
                               layer_id == '6' ~ 'Gomera',
                               layer_id == '7' ~ 'Hierro'))

modules_with_lat_lon_id <- modules_with_lat_lon_id %>% unique() #delete doubles caused by layers turning to islands

pdf('./graphs/Islands/Distribution_modules_islands.pdf', 10, 6)
modules_with_lat_lon_id %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_name)))+ geom_bar(stat= "identity")+ theme_classic()+
  scale_x_continuous(breaks=seq(1,88,2))+ labs(y="Number of islands", x="Module number")+
  guides(fill=guide_legend(title="Island"))
dev.off()

### ---Module size distribution
plants_and_pols <- modules_dryad_multilayer%>% count(type) 
plants <- plants_and_pols[1,2]
pols <- plants_and_pols[2,2]
modules_count <- modules_dryad_multilayer %>% count(module, type) #count how many plants and pollinators are found in each module
modules_count_not_proportion <- modules_count
for(i in (1:nrow(modules_count))){
  if (modules_count$type[i] == "plant"){
    modules_count$n[i] <- modules_count$n[i]/plants
  }
  else{
    modules_count$n[i] <- modules_count$n[i]/pols
  }
}


#distribution by species and proportion
pdf('./graphs/modularity_analysis/proportion_species_in_modules.pdf', 10, 6)
modules_count %>% ggplot(aes(x=as.numeric(module), y=as.numeric(n), fill= type))+
  geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ theme_classic()+
  scale_x_continuous(breaks=seq(1,88,2))+ labs(x= "Module Number", y= "Proportion")+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=9, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()

#distribution by module size
library(plotrix)
module_sizes <- modules_dryad_multilayer %>% count(module)
min(module_sizes$n) #2 species
max(module_sizes$n) #44 species
mean(module_sizes$n) #5.56 sp
std.error(module_sizes$n)#1 sp
sd(module_sizes$n)

###-- Plot modules in x islands and size
modules_in_i <- NULL
for (i in 1:max(modules_dryad_multilayer$module)){
  modules_for_similarity_in_i_num <- modules_dryad_multilayer %>% select(module, layer_id) %>%
    unique() %>% group_by(module) %>% filter(n()==i) %>% select(module) %>% unique() #check which modules occur in i or more layers
  modules_for_similarity_in_i <- modules_dryad_multilayer %>%
    filter(module %in% modules_for_similarity_in_i_num$module) %>% count(module)#only save the modules that are found in i or more layers
  modules_in_i <- rbind(modules_in_i, tibble(module= modules_for_similarity_in_i$module, i=i))
}

num_of_modules_in_i <- modules_in_i %>% count(i)

layers_and_sizes <- right_join(modules_in_i, module_sizes, by= c("module"="module")) #combine number of modules with size of module

#pdf('./graphs/modularity_analysis/species_im_modules_by_layers.pdf', 10, 6)
main_bar_plot <- ggplot(layers_and_sizes, aes(x=i, y=n, fill=factor(module)))+ 
  geom_bar(stat="identity", color="black", show.legend= FALSE)+ theme_classic()+
  scale_x_continuous(breaks=seq(1,14,1))+ 
  labs(x="Number of Islands", y="Number of Nodes")

print(main_bar_plot)
#dev.off()


##-- Plot modules in each island with size proportion


#Calculate the number and proportion of species in each module
modules_dryad_multilayer <- read.csv("./csvs/Islands/modules_in_network_islands_as_layers.csv") 

#Calculate number of nodes in each module per island
N_species_mod<-modules_dryad_multilayer %>% 
  group_by(layer_id, module) %>% 
  summarize(module_size=n_distinct(node_id))

#Create dataframe with total number of nodes per module across islands
Module_size <- modules_dryad_multilayer %>% count(module) #num of nodes in module. 

#merge dataframes and calculate proportion of species in each module
Prop_sp_module <- left_join(N_species_mod,Module_size)

Prop_sp_module_2<-Prop_sp_module %>% 
  mutate(Prop_sp = module_size / n)

#change order accoridng to distances

Prop_sp_module_island<- Prop_sp_module_2 %>% 
  mutate(layer_name= case_when(layer_id == '1' ~ 'WesternSahara',
                               layer_id == '2' ~ 'Fuerteventura',
                               layer_id == '3' ~ 'GranCanaria',
                               layer_id == '4' ~ 'TenerifeSouth',
                               layer_id == '5' ~ 'TenerifeTeno',
                               layer_id == '6' ~ 'Gomera',
                               layer_id == '7' ~ 'Hierro'))


#write.csv(Prop_sp_module_island , "./csvs/Islands/Prop_sp_module_island.csv", row.names = FALSE)

pdf('./graphs/Islands/Prop_sp_modules_per_islands.pdf', 10, 6)
ggplot(Prop_sp_module_island, aes(x = module, y = layer_name, fill=Prop_sp )) +
  geom_tile(color='white') +
  theme_classic() +
  scale_fill_viridis_c(limits = c(0, 1)) +
  scale_x_continuous(breaks=seq(1,88,2)) +
  scale_y_discrete(limits = c("Hierro","Gomera","TenerifeTeno","TenerifeSouth","GranCanaria","Fuerteventura","WesternSahara"))+
  labs(x='Module ID', y="Islands")
dev.off()


##-- Plot modules in each island with size proportion according to island 

#Calculate the number and proportion of species in each module
modules_dryad_multilayer <- read.csv("./csvs/Islands/modules_in_network_islands_as_layers.csv") 

#Calculate number of nodes in each module per island
N_species_mod<-modules_dryad_multilayer %>% 
  group_by(layer_id, module) %>% 
  summarize(module_size=n_distinct(node_id))

#Create dataframe with total number of nodes per module across islands
N_sp_island <- modules_dryad_multilayer %>% count(module) #num of nodes in module. 

Total_sp<- modules_dryad_multilayer %>%
  group_by(layer_id)%>% 
  summarise(Total_sp = n())

#merge dataframes and calculate proportion of species in each module
Prop_sp_island<- left_join(N_species_mod,Total_sp)

Prop_sp_island_2<-Prop_sp_island %>% 
  mutate(Prop_sp = module_size / Total_sp)

#change order accoridng to distances
Prop_sp_in_island<- Prop_sp_island_2 %>% 
  mutate(layer_name= case_when(layer_id == '1' ~ 'WesternSahara',
                               layer_id == '2' ~ 'Fuerteventura',
                               layer_id == '3' ~ 'GranCanaria',
                               layer_id == '4' ~ 'TenerifeSouth',
                               layer_id == '5' ~ 'TenerifeTeno',
                               layer_id == '6' ~ 'Gomera',
                               layer_id == '7' ~ 'Hierro'))


#write.csv(Prop_sp_in_island , "./csvs/Islands/Prop_sp_in_island.csv", row.names = FALSE)

pdf('./graphs/Islands/Prop_sp_islands_in_modules.pdf', 10, 6)
ggplot(Prop_sp_in_island, aes(x = module, y = layer_name, fill=Prop_sp )) +
  geom_tile(color='white') +
  theme_classic() +
  scale_fill_viridis_c(limits = c(0, 0.2)) +
  scale_x_continuous(breaks=seq(1,88,2)) +
  scale_y_discrete(limits = c("Hierro","Gomera","TenerifeTeno","TenerifeSouth","GranCanaria","Fuerteventura","WesternSahara"))+
  labs(x='Module ID', y="Islands")
dev.off()


## ---- Plot of total possibilities of interaction rewiring to test null model 3
interactions_co_occurences <-read.csv("./csvs/Islands/interactions_co_occurences.csv")

interactions_co_occurences_pot <- interactions_co_occurences %>%
  mutate(possible_changes = co_occurrences - interactions) %>% 
  summarise(co_occurrences = sum (co_occurrences), interactions = sum(interactions), possible_changes = sum(possible_changes)) %>% 
              mutate(Id = c(1))
interactions_co_occurences_pot2<- pivot_longer(interactions_co_occurences_pot, names_to = "group", values_to = "number", cols = -Id)

pdf('./graphs/Islands/interactions_co_occurences_pot.pdf', 10, 6)
interactions_co_occurences_pot2 %>%
  ggplot(aes(x=group, y= number, fill=group))+ geom_bar(stat='identity', alpha= 0.6, color= "black")+ theme_bw()+
  scale_y_continuous(limits = c(0, 1200), breaks=seq(0,1200,300)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'))
dev.off()  






