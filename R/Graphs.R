#### Exploratory analysis and graphs


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

dryad_intralayer <- read.csv("./csvs_nuevo/intralayer_file.csv")

##---- layers as islands -------------------------------------------------------------------------------
# Remove dataframe of mainland
dryad_intralayer_islands <- dryad_intralayer %>% filter (!(layer_from == "WesternSahara1"|
                                                             layer_from == "WesternSahara2"))

# Merge sites belonging to each island
old_names <- c("Fuerteventura1", "Fuerteventura2",
               "GranCanaria1", "GranCanaria2",
               "TenerifeSouth1", "TenerifeSouth2",
               "TenerifeTeno1", "TenerifeTeno2",
               "Gomera1", "Gomera2",
               "Hierro1", "Hierro2")

new_names <- c( "Fuerteventura", "Fuerteventura",
                "GranCanaria", "GranCanaria",
                "Tenerife", "Tenerife",
                "Tenerife", "Tenerife",
                "Gomera", "Gomera",
                "Hierro", "Hierro")

dryad_intralayer_islands$layer_from[dryad_intralayer_islands$layer_from %in% old_names] <- 
  new_names[match(dryad_intralayer_islands$layer_from, old_names)] #change to reflect layer = island

dryad_intralayer_islands$layer_to[dryad_intralayer_islands$layer_to %in% old_names] <- 
  new_names[match(dryad_intralayer_islands$layer_to, old_names)] #change to reflect layer = island

#if node_from, node_to, layer_from, layer_to are all the same need to sum the weight
dryad_intralayer_islands$weight<-as.integer(dryad_intralayer_islands$weight)
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

interlayers_with_weights_islands<-interlayers_with_weights_islands[,c(2,1,3,4,5)]#change order columns

#inverted version
interlayer_inverted <- tibble(values= interlayers_with_weights_islands$layer_to, interlayers_with_weights_islands$node_to, interlayers_with_weights_islands$layer_from, 
                              interlayers_with_weights_islands$node_from, interlayers_with_weights_islands$weight) #create an inverted copy for directed intralayers
colnames(interlayer_inverted) <- c("layer_from", "node_from", "layer_to", "node_to", "weight")

#Create interedgelist
edgelist_interlayers_both <- bind_rows(interlayers_with_weights_islands, interlayer_inverted) #combine inverted and non inverted versions of intra


## ----multilayer_extended_final--------------------------------------------------------------------------------------
dryad_edgelist_complete <- bind_rows(edgelist_intralayers_both, edgelist_interlayers_both) #combine weighted version of intra with inter

## ----node_metadata--------------------------------------------------------------------------------------------------                                            
pollinators <- sort(unique(intralayer_weighted$node_to)) #adding up only pol who haven't been added yet 
plants <- sort(unique(intralayer_weighted$node_from)) #adding up only plants who haven't been added yet
intersect(pollinators, plants) #making sure I don't have plants in pol or other way around
A <- length(pollinators) # Number of pollinators
P <- length(plants) # Number of plants
S <- A+P

island_names <- c("Fuerteventura",#islands as layers
                  "GranCanaria",
                  "Tenerife",
                  "Gomera",
                  "Hierro")

# Create a table with node metadata
physical_nodes <- tibble(node_id=1:S, #1 till the last species
                         type=c(rep('plant',P),rep('pollinator',A)), #replicate the words P and A times
                         species=c(plants,pollinators)) #add species from plants and pollinators in accordance
layer_metadata <- tibble(layer_id=c(1:5), layer_name=island_names)  #give num to each layer

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



#Prop of shared partners plants
prop_parn_plants<- dryad_edgelist_complete_ids %>% filter(layer_from!=layer_to) %>% filter(node_to<=39) %>% 
  summarize(Prop=mean(weight), stand = sd(weight))

#Prop of shared partners pollinators
prop_parn_poll<- dryad_edgelist_complete_ids %>% filter(layer_from!=layer_to) %>% filter(node_to>39) %>% 
  summarize(Prop=mean(weight), stand = sd(weight))


#Plot links distribution (directed)
inter_extended<- dryad_edgelist_complete_ids %>% filter(layer_from!=layer_to)
intra_inter_data_for_distibution <- data.frame(values= c(intralayer_weighted$weight, 
                                                         intralayer_weighted_inverted$weight, 
                                                         inter_extended$weight),
                                               group= c(rep("intra plants", nrow(intralayer_weighted)), 
                                                        rep("intra pollinators", nrow(intralayer_weighted_inverted)),
                                                        rep("inter", nrow(inter_extended))))


#write.csv(intra_inter_data_for_distibution, "./csvs_nuevo/intra_inter_data_for_distibution_islands_as_layers.csv",  row.names = FALSE)

pdf('./graphs/Fig_S1.pdf', 10, 6)
intra_inter_data_for_distibution %>%
  ggplot(aes(x=values, fill=group))+ geom_histogram(position= "identity", alpha= 0.70, color= "black")+ 
  labs(x='Weight', y="Count") +theme_bw()+
  scale_fill_manual(name = "Links",  label = c("Interlayer","Intralayer \nplants", "Intralayer \npollinators"), values = c("#F53416","#008000","#CBADFF"))+
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


  
#---- Exploratory graph of species distribution --------------------------------------------------------------------------------
tot_plant_ids <- tot_plant %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>%
  subset(select = -layer_from) %>% count(layer_id)
tot_plant_ids$type <- "plant"

tot_pol_ids <- tot_pol %>% inner_join(layer_metadata, by= c("layer_from" = "layer_name")) %>%
  subset(select = -layer_from) %>% count(layer_id)
tot_pol_ids$type <- "pollinator"

richness_in_islands <- rbind(tot_plant_ids, tot_pol_ids)

richness_in_islands %>%
  ggplot(aes(x= layer_id, y=n, fill=type))+ geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ 
  scale_x_continuous(breaks=seq(1,5,1))+ labs(x= "Island ID", y= "Number of Species")+ theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
       axis.title = element_text(size=14, color='black'),
       axis.line = element_blank())

## -- Number of interaction per island
Num_int_island<-dryad_edgelist_complete_ids %>% filter(layer_from==layer_to) %>% 
  group_by(layer_from) %>% count() 

## -- Number of interedges involving poll and plants
Num_inter_taxon<-dryad_edgelist_complete_ids %>% filter(layer_from!=layer_to) %>% 
  left_join(physical_nodes,by = c("node_from" ="node_id")) %>% 
  group_by(type) %>% count()

## -- Number of species in the Canary Islands (region)

Number_sps_region<-dryad_edgelist_complete_ids %>% select(node_from) %>% 
  distinct()

Number_pol_plant<-Number_sps_region %>% left_join(physical_nodes,by = c("node_from" ="node_id")) %>% 
  group_by(type) %>% count()

##-- Distribution of species in islands
co_occurrence_species<- dryad_edgelist_complete_ids %>% filter(layer_from == layer_to) %>% 
 mutate(layer_id = layer_from, node_id = node_from) %>% select(layer_id, node_id) %>% 
  mutate(layer_name= case_when(layer_id == '1' ~ 'Fuerteventura',
                               layer_id == '2' ~ 'GranCanaria',
                               layer_id == '3' ~ 'Tenerife',
                               layer_id == '4' ~ 'Gomera',
                               layer_id == '5' ~ 'Hierro'))

co_occurrence_species<- left_join(co_occurrence_species,physical_nodes, by="node_id") %>% select(-species) %>% 
  group_by(layer_id,node_id,type,layer_name) %>% unique() %>% count()


## -- Comparision between species co-occurrence and number of interedges

#Maximum potential interlayer links according to species co-occurrence
co_occurrence_tot<-co_occurrence_species %>%  group_by(node_id) %>%  
  mutate(tot_island = sum (n)) %>% 
  summarize(Num_pot_interedge= case_when(tot_island == '1' ~ '0',
                               tot_island == '2' ~ '1',
                               tot_island == '3' ~ '3',
                               tot_island == '4' ~ '6',
                               tot_island == '5' ~ '10'))%>% unique()  #Potential links according to species co-occurence


co_occurrence_tot$Num_pot_interedge<-as.numeric(co_occurrence_tot$Num_pot_interedge)
pot_interedge<-co_occurrence_tot %>% ungroup() %>% summarize(Pot_interedge = sum(Num_pot_interedge) *2) #because it is directed, so doble op

#Realized interedges links
real_interedge<- dryad_edgelist_complete_ids %>% filter(layer_from!=layer_to) %>% 
                summarize(Real_interedges = n())

interedges_comp<-cbind(pot_interedge,real_interedge, col=1)
interedges_comp2<- pivot_longer(interedges_comp, names_to = "group", values_to = "Count", cols = -col)


interedges_comp2%>%
  ggplot(aes(x=group, y= Count, fill=group))+ geom_bar(stat='identity', alpha= 0.6, color= "black")+ theme_bw()+
  scale_y_continuous(limits = c(0, 900), breaks=seq(0,900,100)) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'))
 

## Distribution of plants
co_occurrence_plants<- co_occurrence_species %>% filter(type=="plant") %>% group_by(node_id) %>%  
  mutate(tot_island = sum (n)) %>% arrange(desc(tot_island))

# Convert node_id to an ordered factor (change order that plant appears in the x axis)
co_occurrence_plants$node_id<- as.factor(co_occurrence_plants$node_id)
co_occurrence_plants$node_id <- fct_inorder(co_occurrence_plants$node_id)

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

co_occurrence_pol%>% 
  ggplot(aes(x=node_id, y= n ,fill= factor(layer_name)))+ geom_bar(stat= "identity")+ theme_classic()+ labs(y="Number of islands", x="Species ID")+
  guides(fill=guide_legend(title="Island"))+ theme(axis.text.x=element_blank())





#---- Exploratory graphs and plots of modules --------------------------------------------------------------------------------

##----upload data frame islands  -------
modules_dryad_multilayer <- read.csv("./csvs_nuevo/module_data_justislands_as_layers.csv") 


#-----------Final Plots -------------

## Number of plant and pollinators per module (Fig. S3) --
Num_sp_module<-modules_dryad_multilayer %>% 
  group_by(module) %>% count(type)

pdf('./graphs/Fig_S3.pdf', 10, 8)
Num_sp_module %>% ggplot(aes(x=as.numeric(module), y=as.numeric(n), fill= type))+
  geom_bar(stat="identity", position= position_dodge2(preserve = "single"))+ theme_classic()+
  scale_x_continuous(breaks=seq(1,33))+ labs(x= "Module ID", y= "Number of species")+
  scale_fill_manual(name = "Trophic group",  label = c("Plant","Pollinator"), values = c("#008000","#DDA0DD"))+
  theme(panel.grid = element_blank(),panel.background = element_blank(), 
        axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y=element_text(size=13, colour = "black"), axis.title = element_text(size=14),
        legend.text=element_text(size=11.5),legend.title =element_text(size=12),
        axis.line = element_line(colour = "black"))
dev.off()


## Distribution of modules in islands (Fig. S5) --

modules_with_lat_lon <- read.csv("./csvs_nuevo/modules_with_lat_lon_justislands_as_layers.csv") 
modules_with_lat_lon$layer_id<- as.character(modules_with_lat_lon$layer_id)

modules_with_lat_lon_id <- modules_with_lat_lon %>% 
  mutate(layer_name= case_when(layer_id == '1' ~ 'Fuerteventura',
                               layer_id == '2' ~ 'GranCanaria',
                               layer_id == '3' ~ 'Tenerife',
                               layer_id == '4' ~ 'Gomera',
                               layer_id == '5' ~ 'Hierro'))

modules_with_lat_lon_id <- modules_with_lat_lon_id %>% unique()  %>% group_by(module) %>%  
  mutate(tot_island = sum (count)) %>% arrange(desc(tot_island))

# Convert module to an ordered factor and change the order according to number of islands
modules_with_lat_lon_id$module<- as.factor(modules_with_lat_lon_id$module)
modules_with_lat_lon_id$module <- fct_inorder(modules_with_lat_lon_id$module)

pdf('./graphs/Fig_S5.pdf', 10, 6)
modules_with_lat_lon_id %>% 
  ggplot(aes(x=module, y= count ,fill= factor(layer_name)))+ geom_bar(stat= "identity")+ theme_classic()+ labs(y="Number of locations", x="Modules")+
  guides(fill=guide_legend(title="Location"))+ theme(axis.text.x = element_blank(),
                                                   axis.text.y=element_text(size=15), axis.title = element_text(size=17),
                                                   legend.text=element_text(size=12.5),legend.title =element_text(size=14))
dev.off()


# Average description of modules --
library(plotrix)

module_sizes <- modules_dryad_multilayer %>% count(module)
min(module_sizes$n)
max(module_sizes$n) 
mean(module_sizes$n) 
std.error(module_sizes$n)
sd(module_sizes$n)

prueba<-emp %>% mutate(mean_turnover = mean(turnover), sd_turnover = sd(turnover))

## -- Plot modules in each island with size proportion according to island (Fig. S4)

#Calculate the number and proportion of species in each module
modules_dryad_multilayer <- read.csv("./csvs_nuevo/module_data_justislands_as_layers.csv") 

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
  mutate(layer_name= case_when(layer_id == '1' ~ 'Fuerteventura',
                               layer_id == '2' ~ 'GranCanaria',
                               layer_id == '3' ~ 'Tenerife',
                               layer_id == '4' ~ 'Gomera',
                               layer_id == '5' ~ 'Hierro'))


#write.csv(Prop_sp_in_island , "./csvs_nuevo/Prop_sp_in_island.csv", row.names = FALSE)

pdf('./graphs/Fig_S4.pdf', 10, 6)
ggplot(Prop_sp_in_island, aes(x = module, y = layer_name, fill=Prop_sp )) +
  geom_tile(color='white') +
  theme(panel.grid = element_blank(),panel.background = element_blank(), 
    axis.text.x = element_text(size=11, colour = "black"),
        axis.text.y=element_text(size=13, colour = "black"), axis.title = element_text(size=14),
        legend.text=element_text(size=11.5),legend.title =element_text(size=12),
    axis.line = element_line(colour = "black")) +
  scale_fill_viridis(name = "Prop. sp",limits = c(0, 0.50)) +
  scale_x_continuous(breaks=seq(1,33)) +
  scale_y_discrete(limits = c("Hierro","Gomera","Tenerife","GranCanaria","Fuerteventura"))+
  labs(x='Module ID', y="Locations")
dev.off()



######## ------ Final figures manuscript (except map and null model figures)  ##########

## -- Figure 3
jaccard_similarity_layer_empirical_and_null_km <- read.csv("./csvs_nuevo/jaccard_similarity_layer_empirical_and_null_km_justislands_m1.csv")

# Panel A
Panel_A <- jaccard_similarity_layer_empirical_and_null_km %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+
  labs(x="Distance (Km)", y="Jaccard Similarity")+  
  scale_y_continuous(breaks=seq(0,0.6,0.1))+
  scale_color_manual(name = "Models",  labels = c("Empirical",expression("M"[1]^AP), expression("M"[1]^P),
                                                      expression("M"[1]^A)),values = c("#FB3B1E", "#15B7BC", 
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
jaccard_similarity_layer_empirical_and_null_km_M2 <- read.csv("./csvs_nuevo/jaccard_similarity_layer_empirical_and_null_km_interactions_justislands.csv") #need to read this to run next part

Panel_B<- jaccard_similarity_layer_empirical_and_null_km_M2%>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ 
  theme_classic()+ geom_smooth(method= "lm", se=F)+
  scale_color_manual (name = "Models", labels = c("Empirical",expression("M"[2])),
                      values = c("#FB3B1E",  "#E6AB02" ))+
  scale_y_continuous(breaks=seq(0,0.6,0.1))+
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
jaccard_similarity_layer_empirical_and_null_km_classic_M3 <- read.csv("./csvs_nuevo/jaccard_similarity_layer_empirical_and_null_km_classic_justislands.csv") #need to read this to run next part

Panel_C<- jaccard_similarity_layer_empirical_and_null_km_classic_M3 %>% 
  ggplot(aes(x= mean_dist_in_km, y= ave, group= type, color= type))+
  geom_point()+ geom_errorbar(aes(ymin= ave-sd, ymax= ave+sd))+ theme_classic()+ geom_smooth(method= "lm", se=F)+
  scale_color_manual (name = "Models", labels = c("Empirical",expression("M"[3])),
                      values = c("#FB3B1E","#FA86F2"))+
  labs(x="Distance (Km)", y="Jaccard Similarity")+ 
  scale_y_continuous(breaks=seq(0,0.6,0.1))+
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
pdf('./graphs/Figure_3.pdf', 10, 8)
upper_row<- plot_grid(Panel_A + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")),
                      Panel_B + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
                      ncol = 2, labels = c('(A)', "(B)"))

bottom_row<- plot_grid(NULL,Panel_C + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
                       ncol = 2, labels = c("(C)","x"))

plot_grid(upper_row, bottom_row, ncol = 1)
dev.off()




## -- Figure 4

# Panel A
# Panel A
NM1 <- read.csv("./csvs_nuevo/slopes_NM1_all_nonscaled.csv")
NM1$type <- factor(NM1$type, levels = c("shuf_plants","shuf_pollinators","shuf_both"))

Panel_A <- NM1 %>% 
  ggplot(aes(x = slope, fill = type))+ 
  geom_density(alpha = 0.5)+ 
  geom_vline(xintercept = -0.00087, linetype = "dashed", color = "#FB3B1E")+ #line R squared empirical
  labs(x= "Slope", y="Density")+  
  scale_fill_manual(name = "Models",  labels = c(expression("M"[1]^P),expression("M"[1]^A),
                                                 expression("M"[1]^AP)), values = c("#72A323","#A44CD3", "#15B7BC"))+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=13, color='black'),
        axis.title = element_text(size=17, color='black'),
        axis.line = element_blank(),
        legend.text.align = 0,
        legend.title =  element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11))
Panel_A

#Panel B
iteration_correlation_interactions <- read.csv("./csvs_nuevo/iteration_correlation_interactions_justislands.csv")
iteration_correlation_interactions2<-iteration_correlation_interactions %>% mutate(Type = "null_int")

Panel_B<- iteration_correlation_interactions2 %>% ggplot(aes(x = slope, fill= Type))+ 
  geom_density(alpha = 0.6)+ 
  geom_vline(xintercept = -0.00087, linetype = "dashed", color = "#FB3B1E") + #line represting rsquared empirical
  labs(x= "Slope", y="Density")+  
  scale_fill_manual(name = "Models",  label = expression("M"[2]), values= "#E6AB02")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=13, color='black'),
        axis.title = element_text(size=17, color='black'),
        axis.line = element_blank(),
        legend.text.align = 0,
        legend.title =  element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11))
Panel_B

#Panel C
iteration_correlation_classic <- read.csv("./csvs_nuevo/iteration_correlation_classic_justislands.csv")
iteration_correlation_classic2<-iteration_correlation_classic %>% mutate(Type = "null_class")

Panel_C<- iteration_correlation_classic2 %>% ggplot(aes(x = slope, fill= Type))+ 
  geom_density(alpha = 0.6)+ 
  geom_vline(xintercept = -0.00087, linetype = "dashed", color = "#FB3B1E") +
  labs(x= "Slope", y="Density")+  
  scale_fill_manual(name = "Models",  label = expression("M"[3]), values= "#FA86F2")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=13, color='black'),
        axis.title = element_text(size=17, color='black'),
        axis.line = element_blank(),
        legend.text.align = 0,
        legend.title =  element_text(size = 13, color = "black"),
        legend.text = element_text(size = 11))
Panel_C

# Figure 4 (all panels together)
pdf('./graphs/Figure_4.pdf', 10, 8)
upper_row<- plot_grid(Panel_A + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")),
                      Panel_B + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
                      ncol = 2, labels = c('(A)', "(B)"))

bottom_row<- plot_grid(NULL,Panel_C + theme(plot.margin = unit(c(0.2,0.25,0.2,0.5), "cm")), 
                       ncol = 2, labels = c("(C)","x"))

plot_grid(upper_row, bottom_row, ncol = 1, rel_heights = c(1, 1))
dev.off()



