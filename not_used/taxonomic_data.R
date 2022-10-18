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

#this portion of the code comes to organize the taxonomic data known about the species in the network
#using the ncbi data base to try and see whether modules were divided according to taxonomic similarity

#---- adding basic taxonomic data to data frame-----------------------------------------------------------------------------------------------------

#tax_name(sci = "Themnothorax", get = c("family","order"), db = "ncbi") #has taxa information

physical_nodes_with_family <- NULL

nodes_for_taxa <- physical_nodes %>% separate(species, c("genus", "species")) #create version where genus and species are separated

for(i in 1:288){
  print(i)
  current_physical_node <- physical_nodes %>% filter(node_id == i)
  wanted_genus <- nodes_for_taxa %>% filter(node_id == i) %>% select(genus) %>% unlist() #take the genus name in every row
  wanted_family_and_order <- tax_name(sci = wanted_genus, get = c("family", "order"), db = "ncbi") #get family and order of every genus provided
  physical_nodes_with_family <- rbind(physical_nodes_with_family, tibble(node_id = i, type = current_physical_node$type, 
                                                                         species = current_physical_node$species, genus = wanted_genus,
                                                                         family = wanted_family_and_order$family,
                                                                         order = wanted_family_and_order$order))
} #what do i do with missing information?

#write.csv(physical_nodes_with_family, "./csvs/physical_nodes_with_family.csv", row.names = FALSE)

modules_with_taxa_info <- merge(modules_dryad_multilayer$modules, physical_nodes_with_family, by=c("species","node_id","type")) #add taxa info to modules

tot_species_in_module <- modules_with_taxa_info  %>% count(module)

module_taxa_and_num_of_species <- merge(modules_with_taxa_info, tot_species_in_module, by = "module")
names(module_taxa_and_num_of_species)[10] <- "num_of_nodes_in_module"
module_taxa_and_num_of_species$count <- 1

module_taxa_and_num_of_species %>%
  ggplot(aes(x= module, y=count, fill=order))+ geom_bar(stat="identity")+ #stacked
  theme_classic()+ scale_x_continuous(breaks=seq(1,43,2))+ labs(x= "Module", y= "Number of Nodes in Module")

#---- showing taxonomic data on hierarchial modules----------------------------------------------------------------------------------------------------

modules_multi_lvl_with_taxa_info <- merge(modules_dryad_multilayer_multi_level, physical_nodes_with_family, by=c("species","node_id","type"))

tot_species_in_module_multi_lvl <- modules_multi_lvl_with_taxa_info  %>% count(module)

module_multi_lvl_taxa_and_num_of_species <- merge(modules_multi_lvl_with_taxa_info, tot_species_in_module_multi_lvl, by = "module")
names(module_multi_lvl_taxa_and_num_of_species)[10] <- "num_of_nodes_in_module"
module_multi_lvl_taxa_and_num_of_species$count <- 1

module_multi_lvl_taxa_and_num_of_species_count <- module_multi_lvl_taxa_and_num_of_species %>% group_by(module) %>% count(order)

module_1 <- module_multi_lvl_taxa_and_num_of_species %>% filter(module == 1)
module_1_by_module_and_order <- module_1 %>% group_by(order) %>% count(order)
module_2 <- module_multi_lvl_taxa_and_num_of_species %>% filter(module == 2)
module_2_by_module_and_order <- module_2 %>% group_by(order) %>% count(order)

#barplot modules
module_multi_lvl_taxa_and_num_of_species_count %>%
  ggplot(aes(x= module, y=n, fill=order))+ geom_bar(stat="identity")+ #stacked
  theme_classic()+ scale_x_continuous(breaks=seq(1,43,1))+ labs(x= "Module", y= "Number of Nodes in Module")

#pie chart modules
module_1_by_module_and_order %>%
  ggplot(aes(x= "", y=n, fill=order))+ geom_bar(stat="identity", width = 1)+
  theme_void() +coord_polar("y", start = 0)

module_2_by_module_and_order %>%
  ggplot(aes(x= "", y=n, fill=order))+ geom_bar(stat="identity", width = 1)+
  theme_void() +coord_polar("y", start = 0)

#barplot sub modules
module_1 %>%
  ggplot(aes(x= sub_module, y=count, fill=order))+ geom_bar(stat="identity")+ #stacked
  theme_classic()+ scale_x_continuous(breaks=seq(1,41,2))+ labs(x= "Module", y= "Number of Nodes in Module")

module_2 %>%
  ggplot(aes(x= sub_module, y=count, fill=order))+ geom_bar(stat="identity")+ #stacked
  theme_classic()+ scale_x_continuous(breaks=seq(1,41,2))+ labs(x= "Module", y= "Number of Nodes in Module")

#---- only plants-----------------------------------------------------------------------------------------------------

only_plants_with_taxa <- module_taxa_and_num_of_species %>% filter(type == "plant")

only_plants_with_taxa %>%
  ggplot(aes(x= module, y=count, fill=order))+ geom_bar(stat="identity")+ #stacked
  theme_classic()+ scale_x_continuous(breaks=seq(1,37,2))+ labs(x= "Module", y= "Number of Nodes in Module")

#---- only pollinators--------------------------------------------------------------------------------------------------

only_pols_with_taxa <- module_taxa_and_num_of_species %>% filter(type == "pollinator")

only_pols_with_taxa %>%
  ggplot(aes(x= module, y=count, fill=order))+ geom_bar(stat="identity")+ #stacked
  theme_classic()+ scale_x_continuous(breaks=seq(1,37,2))+ labs(x= "Module", y= "Number of Nodes in Module")

#---- built phylogeny tree-----------------------------------------------------------------------------------------

order_names_plants <- only_plants_with_taxa$order %>% unique() %>% na.omit() #create list of plant orders
order_names_pols <- only_pols_with_taxa$order %>% unique() %>% na.omit() #create list of pol orders

order_names_plants_class <- classification(order_names_plants, db = "ncbi")
plants_tree <- class2tree(order_names_plants_class)
plot(plants_tree)

order_names_pols_class <- classification(order_names_pols, db = "ncbi")
pols_tree <- class2tree(order_names_pols_class)
plot(pols_tree)

