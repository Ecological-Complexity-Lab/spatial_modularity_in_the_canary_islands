#### ISLANDS AND JACCARD AS INTEREDGES


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
library(reshape2)
library(extRC)



##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")

dryad_intralayer <- read.csv("./csvs/intralayer_file.csv")
dryad_interlayer <- read.csv("./csvs/interlayer_file.csv") #already has inverted within

##---- layers as islands -------------------------------------------------------------------------------
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

## ----dryad intralayer interlayer both ways-------------------------------------------------------------------------------------
intralayer_inverted <- tibble(values= dryad_intralayer_islands$layer_to, dryad_intralayer_islands$node_to, dryad_intralayer_islands$layer_from, 
                              dryad_intralayer_islands$node_from, dryad_intralayer_islands$weight) #create an inverted copy for directed intralayers
colnames(intralayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight")


## ---- normalize intralayer weights--------------------------------------------------------------------------
#plants in from
tot_plant <- dryad_intralayer_islands %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted <- dryad_intralayer_islands %>% left_join(tot_plant) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight)

#pols in from
tot_pol <- intralayer_inverted %>% 
  group_by(layer_from,node_from) %>% 
  dplyr::summarise(tot=sum(weight))
intralayer_weighted_inverted <- intralayer_inverted %>% left_join(tot_pol) %>% mutate(rel_weight=weight/tot) %>% 
  select(-weight,-tot) %>% dplyr::rename(weight=rel_weight) 

#Create intraedgelist
edgelist_intralayers_both <- bind_rows(intralayer_weighted, intralayer_weighted_inverted) #combine inverted and non inverted versions of intra

##---- create interedges according to jaccard similarity between partners --------------------------------------------------------
# keep only species that occur in 2 or more layers
co_occurrence<- edgelist_intralayers_both%>% 
  group_by(node_from) %>%
  mutate(num_layers_from=n_distinct(layer_from)) %>% 
  filter(num_layers_from>="2")

# a for loop that calculates all the interlayer edges based on jaccard index
interlayers_with_weights_islands <- NULL


for (i in unique(co_occurrence$node_from)) {
  print(i)
  partners_sp <- 
    co_occurrence %>%
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
    mutate(node_from=i, node_to =i) %>%
    select(c(node_from,layer_from=Var1, layer_to=Var2,node_to, weight=value))
  interlayers_with_weights_islands <- rbind(interlayers_with_weights_islands,inter_fid)
}


## ----multilayer_extended_final--------------------------------------------------------------------------------------
dryad_edgelist_complete <- bind_rows(edgelist_intralayers_both, interlayers_with_weights_islands) #combine weighted version of intra with inter

## ----node_metadata--------------------------------------------------------------------------------------------------                                            
pollinators <- sort(unique(intralayer_weighted$node_to)) #adding up only pol who haven't been added yet 
plants <- sort(unique(intralayer_weighted$node_from)) #adding up only plants who haven't been added yet
intersect(pollinators, plants) #making sure I don't have plants in pol or other way around
A <- length(pollinators) # Number of pollinators
P <- length(plants) # Number of plants
S <- A+P

island_names <- c("WesternSahara", #islands as layers
                  "Fuerteventura",
                  "GranCanaria",
                  "TenerifeSouth",
                  "TenerifeTeno",
                  "Gomera",
                  "Hierro")

# Create a table with node metadata
physical_nodes <- tibble(node_id=1:S, #1 till the last species
                         type=c(rep('plant',P),rep('pollinator',A)), #replicate the words P and A times
                         species=c(plants,pollinators)) #add species from plants and pollinators in accordance
layer_metadata <- tibble(layer_id=c(1:7), layer_name=island_names)  #give num to each layer


#write.csv(physical_nodes, "./csvs/Islands/Jac/physical_nodes_islands.csv", row.names = FALSE)
#write.csv(layer_metadata, "./csvs/Islands/Jac/layer_metadata_islands.csv", row.names = FALSE)


# Replace the node names with node_ids
dryad_edgelist_complete_ids <- 
  dryad_edgelist_complete %>% 
  left_join(physical_nodes, by=c('node_from' = 'species')) %>%  # Join for pollinators
  left_join(physical_nodes, by=c('node_to' = 'species')) %>%  # Join for plants
  dplyr::select(-node_from, -node_to) %>% #choose said columns
  dplyr::select(layer_from, node_from=node_id.x, layer_to, node_to=node_id.y, weight) %>% 
  left_join(layer_metadata, by=c('layer_from' = 'layer_name')) %>%  # Join for plants
  left_join(layer_metadata, by=c('layer_to' = 'layer_name')) %>%  # Join for plants
  dplyr::select(-layer_from, -layer_to) %>% 
  dplyr::select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)

#write.csv(dryad_edgelist_complete_ids, "./csvs/Islands/Jac/dryad_edgelist_complete_ids_islands.csv", row.names = FALSE)

#Plot links distribution (directed)
inter_extended<- dryad_edgelist_complete_ids %>% filter(layer_from!=layer_to)
intra_inter_data_for_distibution <- data.frame(values= c(intralayer_weighted$weight, 
                                                         intralayer_weighted_inverted$weight, 
                                                         inter_extended$weight),
                                               group= c(rep("intra plants", nrow(intralayer_weighted)), 
                                                        rep("intra pollinators", nrow(intralayer_weighted_inverted)),
                                                        rep("inter", nrow(inter_extended))))

#write.csv(intra_inter_data_for_distibution, "./csvs/Islands/Jac/intra_inter_data_for_distibution_islands_as_layers.csv",  row.names = FALSE)

pdf('./graphs/Islands/Jac/intra and interlinks_islands_jaccard.pdf', 10, 6)
intra_inter_data_for_distibution %>%
  ggplot(aes(x=values, fill=group))+ geom_histogram(position= "identity", alpha= 0.6, color= "black")+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
dev.off()


#---- basic analysis for richness in island --------------------------------------------------------------------------------
tot_plant_ids <- tot_plant %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>%
  subset(select = -layer_from) %>% count(layer_id)
tot_plant_ids$type <- "plant"

tot_pol_ids <- tot_pol %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>%
  subset(select = -layer_from) %>% count(layer_id)
tot_pol_ids$type <- "pollinator"

richness_in_islands <- rbind(tot_plant_ids, tot_pol_ids)

richness_in_islands %>%
  ggplot(aes(x= layer_id, y=n, fill=type))+ geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ 
  scale_x_continuous(breaks=seq(1,7,1))+ labs(x= "Island ID", y= "Number of Species")+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())

##-- Distribution of species in islands
co_occurrence_species<- dryad_edgelist_complete_ids %>% filter(layer_from == layer_to) %>% 
 mutate(layer_id = layer_from, node_id = node_from) %>% select(layer_id, node_id) %>% 
  mutate(layer_name= case_when(layer_id == '1' ~ 'WesternSahara',
                               layer_id == '2' ~ 'Fuerteventura',
                               layer_id == '3' ~ 'GranCanaria',
                               layer_id == '4' ~ 'TenerifeSouth',
                               layer_id == '5' ~ 'TenerifeTeno',
                               layer_id == '6' ~ 'Gomera',
                               layer_id == '7' ~ 'Hierro'))

co_occurrence_species<- left_join(co_occurrence_species,physical_nodes, by="node_id") %>% select(-species) %>% 
  group_by(layer_id,node_id,type,layer_name) %>% unique() %>% count()

#write.csv(co_occurrence_species, "./Results/Jac/co_occurrence_species.csv", row.names = FALSE)

## -- Comparision between species co-occurrence and number of interedges

#Maximum potential interlayer links according to species co-occurrence
co_occurrence_tot<-co_occurrence_species %>%  group_by(node_id) %>%  
  mutate(tot_island = sum (n)) %>% 
  summarize(Num_pot_interedge= case_when(tot_island == '1' ~ '0',
                               tot_island == '2' ~ '1',
                               tot_island == '3' ~ '3',
                               tot_island == '4' ~ '6',
                               tot_island == '5' ~ '10',
                               tot_island == '6' ~ '15',
                               tot_island == '7' ~ '21'))%>% unique()  #Potential links according to species co-occurence


co_occurrence_tot$Num_pot_interedge<-as.numeric(co_occurrence_tot$Num_pot_interedge)
pot_interedge<-co_occurrence_tot %>% ungroup() %>% summarize(Pot_interedge = sum(Num_pot_interedge))

#Realized interedges links
real_interedge<- dryad_edgelist_complete_ids %>% filter(layer_from!=layer_to) %>% 
                summarize(Real_interedges = n())  

interedges_comp<-cbind(pot_interedge,real_interedge, col=1)
interedges_comp2<- pivot_longer(interedges_comp, names_to = "group", values_to = "Count", cols = -col)

pdf('./graphs/Islands/Jac/realized_interedges.pdf', 10, 6)
interedges_comp2%>%
  ggplot(aes(x=group, y= Count, fill=group))+ geom_bar(stat='identity', alpha= 0.6, color= "black")+ theme_bw()+
  scale_y_continuous(limits = c(0, 500), breaks=seq(0,500,50)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'))
dev.off()  

## Distribution of plants
co_occurrence_plants<- co_occurrence_species %>% filter(type=="plant") %>% group_by(node_id) %>%  
  mutate(tot_island = sum (n)) %>% arrange(desc(tot_island))

# Convert node_id to an ordered factor (change order that plant appears in the x axis)
co_occurrence_plants$node_id<- as.factor(co_occurrence_plants$node_id)
co_occurrence_plants$node_id <- fct_inorder(co_occurrence_plants$node_id)

pdf('./graphs/Islands/co_occurrence_plants.pdf', 10, 6)
co_occurrence_plants%>% 
  ggplot(aes(x=node_id, y= n ,fill= factor(layer_name)))+ geom_bar(stat= "identity")+ theme_classic()+ labs(y="Number of islands", x="Species ID")+
  guides(fill=guide_legend(title="Island"))
dev.off()

## Distribution of pollinators
co_occurrence_pol<- co_occurrence_species %>% filter(type=="pollinator") %>% group_by(node_id) %>%  
  mutate(tot_island = sum (n)) %>% arrange(desc(tot_island))

# Convert node_id to an ordered factor (change order that plant appears in the x axis)
co_occurrence_pol$node_id <- as.factor(co_occurrence_pol$node_id)
co_occurrence_pol$node_id <- fct_inorder(co_occurrence_pol$node_id)

pdf('./graphs/Islands/co_occurrence_pol.pdf', 10, 6)
co_occurrence_pol%>% 
  ggplot(aes(x=node_id, y= n ,fill= factor(layer_name)))+ geom_bar(stat= "identity")+ theme_classic()+ labs(y="Number of islands", x="Species ID")+
  guides(fill=guide_legend(title="Island"))+ theme(axis.text.x=element_blank())
dev.off()



##--- Exploratory graphs modules

##----upload data frame islands  -------
modules_dryad_multilayer <- read.csv("./csvs/Islands/Jac/modules_in_network_islands_as_layers.csv",row.names = 1) 

##-- Check if modules are integrated by the same species across layers
same_species_mod<-modules_dryad_multilayer  %>% select(-flow,-species) %>% 
  group_by(module) %>% distinct(node_id) %>% count() #all modules have two or more different species

df<-data.frame(Group= c("Self-species per module","Multi-species per module"),
               Count = c(0,30))

pdf('./graphs/Islands/Jac/Modules_self-species.pdf', 10, 6)
df%>%
  ggplot(aes(x=Group, y= Count, fill=Group))+ geom_bar(stat='identity', alpha= 0.6, color= "black")+ theme_bw()+
  scale_y_continuous(limits = c(0, 50), breaks=seq(0,50,10)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'))
dev.off()  



#-----------Final Plots -------------

#Distribution of modules in islands
#modules_with_lat_lon <- read.csv("csvs/Islands/Jac/modules_with_lat_lon_islands_as_layers.csv") 

modules_with_lat_lon$layer_id<- as.character(modules_with_lat_lon$layer_id)

modules_with_lat_lon_id <- modules_with_lat_lon %>% 
  mutate(layer_name= case_when(layer_id == '1' ~ 'WesternSahara',
                               layer_id == '2' ~ 'Fuerteventura',
                               layer_id == '3' ~ 'GranCanaria',
                               layer_id == '4' ~ 'TenerifeSouth',
                               layer_id == '5' ~ 'TenerifeTeno',
                               layer_id == '6' ~ 'Gomera',
                               layer_id == '7' ~ 'Hierro'))

modules_with_lat_lon_id <- modules_with_lat_lon_id %>% unique()  %>% group_by(module) %>%  
  mutate(tot_island = sum (count)) %>% arrange(desc(tot_island))

# Convert module to an ordered factor and change the order according to number of islands
modules_with_lat_lon_id$module<- as.factor(modules_with_lat_lon_id$module)
modules_with_lat_lon_id$module <- fct_inorder(modules_with_lat_lon_id$module)

pdf('./graphs/Islands/Jac/Distribution_modules_islands.pdf', 10, 6)
modules_with_lat_lon_id %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_name)))+ geom_bar(stat= "identity")+ theme_classic()+ labs(y="Number of locations", x="Modules")+
  guides(fill=guide_legend(title="Location"))+ theme(axis.text.x = element_blank(),
                                                   axis.text.y=element_text(size=15), axis.title = element_text(size=17),
                                                   legend.text=element_text(size=12.5),legend.title =element_text(size=14))
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
pdf('./graphs/modularity_analysis/proportion_species_in_modules_Jac.pdf', 10, 6)
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
max(module_sizes$n) #109 species
mean(module_sizes$n) #16.33 sp
std.error(module_sizes$n)#4 sp
sd(module_sizes$n)#21


##-- Plot modules in each island with size proportion

#Calculate the number and proportion of species in each module
modules_dryad_multilayer <- read.csv("./csvs/Islands/Jac/modules_in_network_islands_as_layers.csv") 

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


#write.csv(Prop_sp_module_island , "./csvs/Islands/Jac/Prop_sp_module_island.csv", row.names = FALSE)

#pdf('./graphs/Islands/Jac/Prop_sp_modules_per_islands.pdf', 10, 6)
#ggplot(Prop_sp_module_island, aes(x = module, y = layer_name, fill=Prop_sp )) +
#  geom_tile(color='white') +
#  theme_classic() +
#  scale_fill_viridis_c(limits = c(0, 1)) +
#  scale_x_continuous(breaks=seq(1,88,2)) +
#  scale_y_discrete(limits = c("Hierro","Gomera","TenerifeTeno","TenerifeSouth","GranCanaria","Fuerteventura","WesternSahara"))+
#  labs(x='Module ID', y="Islands")
#dev.off()


##-- Plot modules in each island with size proportion according to island 

#Calculate the number and proportion of species in each module
modules_dryad_multilayer <- read.csv("./csvs/Islands/Jac/modules_in_network_islands_as_layers.csv") 

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


#write.csv(Prop_sp_in_island , "./csvs/Islands/Jac/Prop_sp_in_island.csv", row.names = FALSE)

pdf('./graphs/Islands/Jac/Prop_sp_islands_in_modules.pdf', 10, 6)
ggplot(Prop_sp_in_island, aes(x = module, y = layer_name, fill=Prop_sp )) +
  geom_tile(color='white') +
  theme(panel.grid = element_blank(),panel.background = element_blank(), 
    axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y=element_text(size=13, colour = "black"), axis.title = element_text(size=14),
        legend.text=element_text(size=11.5),legend.title =element_text(size=12),
    axis.line = element_line(colour = "black")) +
  scale_fill_viridis(name = "Prop. sp",limits = c(0, 0.50)) +
  scale_x_continuous(breaks=seq(1,88,4)) +
  scale_y_discrete(limits = c("Hierro","Gomera","TenerifeTeno","TenerifeSouth","GranCanaria","Fuerteventura","WesternSahara"))+
  labs(x='Module ID', y="Locations")
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


######## ------ Final figures manuscript (except for the map) ##########

## -- Figure 3
jaccard_similarity_layer_empirical_and_null_km <- read.csv("csvs/Islands/jaccard_similarity_layer_empirical_and_null_km_islands_m1.csv")

# Panel A
Panel_A <- jaccard_similarity_layer_empirical_and_null_km %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+
  labs(x="Distance (Km)", y="Jaccard Similarity")+  
  scale_color_manual(name = "Null Model",  labels = c("E",expression("M"[1]^P),expression("M"[1]^A),
                                                      expression("M"[1]^AP)), values = c("#FB3B1E", "#15B7BC", 
                                                                                         "#72A323", "#BE75FA" )) +
  theme_classic()+ geom_smooth(method= "lm", se=F) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=15, color='black'),
        axis.title = element_text(size=17, color='black'),
        axis.line = element_blank(),
        legend.text.align = 0,
        legend.title =  element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11))
Panel_A

#Panel B
jaccard_similarity_layer_empirical_and_null_km_M2 <- read.csv("./csvs/Islands/jaccard_similarity_layer_empirical_and_null_km_islands.csv") #need to read this to run next part

Panel_B<- jaccard_similarity_layer_empirical_and_null_km_M2%>% 
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
Panel_B


#Panel C
jaccard_similarity_layer_empirical_and_null_km_classic_M3 <- read.csv("./csvs/Islands/jaccard_similarity_layer_empirical_and_null_km_classic_islands.csv") #need to read this to run next part

Panel_C<- jaccard_similarity_layer_empirical_and_null_km_classic_M3 %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  scale_color_manual (name = "Null Model", labels = c("E",expression("M"[3])),
                      values = c("#FB3B1E","#FA86F2"))+
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
Panel_C

# Figure 3 (all panels together)
pdf('./graphs/Islands/Figure_3.pdf', 10, 8)
upper_row<- plot_grid(Panel_A + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")),
                      Panel_B + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
                      ncol = 2, labels = c('(A)', "(B)"))

bottom_row<- plot_grid(NULL,Panel_C + theme(plot.margin = unit(c(0.75,0.25,0.5,0.5), "cm")), NULL,
                       ncol = 3, labels = c("","(C)",""), rel_widths = c(0.48,1,0.48))

plot_grid(upper_row, bottom_row, ncol = 1, rel_heights = c(1, 1))
dev.off()




## -- Figure 4

# Panel A
correlation_empirical_pols <- read.csv("./csvs/Islands/correlation_empirical_pols.csv")
rqsuares_M1_all <- read.csv("./csvs/Islands/rqsuares_M1_all.csv")
rqsuares_M1_all$type <- factor(rqsuares_M1_all$type, levels = c("shuf_plants","shuf_pollinators","shuf_both"))

Panel_A <- rqsuares_M1_all %>% 
  ggplot(aes(x = rsquared, fill = type))+ 
  geom_density(alpha = 0.5)+ 
  geom_vline(xintercept = correlation_empirical_pols$rsquared, linetype = "dashed", color = "#FB3B1E")+
  labs(x= expression("R"^2), y="Density")+  
  scale_fill_manual(name = "Null Model",  labels = c(expression("M"[1]^P),expression("M"[1]^A),
                                                     expression("M"[1]^AP)), values = c("#72A323","#15B7BC", "#A44CD3" ))+
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
Panel_A

#Panel B
correlation_empirical_interactions <- read.csv("./csvs/Islands/correlation_empirical_interactions.csv")
iteration_correlation_interactions <- read.csv("./csvs/Islands/iteration_correlation_interactions_islands.csv")
iteration_correlation_interactions2<-iteration_correlation_interactions %>% mutate(Type = "null_int")

Panel_B<- iteration_correlation_interactions2 %>% ggplot(aes(x = rsquared, fill= Type))+ 
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

#Panel C
correlation_empirical_classic <- read.csv("./csvs/Islands/correlation_empirical_classic_islands.csv")
iteration_correlation_classic <- read.csv("./csvs/Islands/iteration_correlation_classic_islands.csv")
iteration_correlation_classic2<-iteration_correlation_classic %>% mutate(Type = "null_class")

Panel_C<- iteration_correlation_classic2 %>% ggplot(aes(x = rsquared, fill= Type))+ 
  geom_density(alpha = 0.6)+ 
  geom_vline(xintercept = correlation_empirical_classic$rsquared, linetype = "dashed", color = "#FB3B1E") +
  labs(x= expression("R"^2), y="Density")+  
  scale_fill_manual(name = "Null Model",  label = expression("M"^3), values= "#FA86F2")+
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
Panel_C

# Figure 4 (all panels together)
pdf('./graphs/Islands/Figure_4.pdf', 10, 8)
upper_row<- plot_grid(Panel_A + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")),
                      Panel_B + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
                      ncol = 2, labels = c('(A)', "(B)"))

bottom_row<- plot_grid(NULL,Panel_C + theme(plot.margin = unit(c(0.75,0.25,0.5,0.5), "cm")), NULL,
                       ncol = 3, labels = c("","(C)",""), rel_widths = c(0.48,1,0.48))

plot_grid(upper_row, bottom_row, ncol = 1, rel_heights = c(1, 1))
dev.off()
