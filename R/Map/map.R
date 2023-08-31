



## ---- Map showing structure similarity between islands
library(tmaptools)
library(rgdal)
library(sf)
library(raster)
library(OpenStreetMap)
library(maptiles)
library(sp)
library(glue)
library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(colorspace)
library(sf)
library(maps)
library(tidyselect)

##----get_data--------------------------------------------------------------------------------------------------------
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
source("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands/R/functions.R")



#-- Plots containing 10% biggest modules 
modules_dryad_multilayer <- read.csv("./csvs/Islands/Jac/modules_in_network_islands_as_layers.csv") 
module_data_with_loc <- read.csv("./csvs/Islands/Jac/modules_with_lat_lon_islands_as_layers.csv")
lon_lat_data <- read.csv("./csvs/Islands/Jac/lon_lat_data_islands_as_layers.csv")

#Calculate number of nodes in each module per island
N_species_mod<-modules_dryad_multilayer %>% 
  group_by(layer_id, module) %>% 
  summarize(module_size=n_distinct(node_id))


# Calculate 30% biggest modules 
top_30_percent <- select(module_data_with_loc, c("module", "size_of_module")) %>% distinct %>% arrange(desc(size_of_module)) %>%
  filter(quantile(size_of_module, probs=0.7) < size_of_module) #arrange all modules to see 90 percentile and filter by them
top_30_percent$size_of_module <- top_30_percent$module #don't need the sizes anymore
names(top_30_percent)[names(top_30_percent)=="size_of_module"] <- "group" #change name to group- only the top percentile
pie_chart_data <- top_30_percent %>% right_join(module_data_with_loc, by=c("module"="module")) #join to have coordinates as well
pie_chart_data$group <- as.character(pie_chart_data$group) 
pie_chart_data <- pie_chart_data %>% drop_na()

#Combine dataframes
pie_chart_data2<- left_join(pie_chart_data, N_species_mod, by = c("layer_id","module"))

#Prepare final dataframe
pivot_by_country <- function(data) {
  pie_chart_data2
  s1 = melt(pie_chart_data2, id = c("layer_id", "group"), measure.vars = "module_size")
  s2 = dcast(s1, layer_id ~ group, measure.vars = "module_size", fill = 0)
  s2$Total = rowSums(s2[,2:NCOL(s2)])
  return(s2)
}

pie_chart_data2 <- lon_lat_data %>% right_join(pivot_by_country(pie_chart_data2), by = c("layer_id"="layer_id")) 
pie_chart_data2 <- pie_chart_data2 %>% group_by(lat, Lon) %>% summarise(across(everything(), sum))


#Plot
top_30_scatterpies <- ggplot()+ geom_scatterpie(aes(x=Lon, y=lat, group=layer_id, r=0.2), alpha=0.90, data= pie_chart_data2, 
                                                cols=colnames(pie_chart_data2[,c(4:11)]))+
  theme_void()+ coord_fixed()+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))
top_30_scatterpies


##-- Map of Canary Island region with Jaccard Similarity 

#combine turnover data frame with coordinates
empirical_turnover_for_module_island_shuf_no_self_loop_km <- read.csv("csvs/Islands/Jac/islands_turnover_with_distnace_empirical.csv")
islands_with_lon_lat_dif <- read.csv("./csvs/islands_with_lon_lat_dif.csv")

# Modify islands_with_lon_lat_dif dataframe 
islands_with_lon_lat_dif2<- islands_with_lon_lat_dif %>% 
  mutate(layer_from = recode(layer_from, "1" = "1", "2" = "2", "3"="3", "4east" = "4", "4west" = "5", "6" = "6", "7" = "7"),
         layer_to = recode(layer_to, "1" = "1", "2" = "2", "3"="3", "4east" = "4", "4west" = "5", "5" = "6", "6" = "7"))


#join both data frames
jaccard_similarity_on_map <- merge(empirical_turnover_for_module_island_shuf_no_self_loop_km, islands_with_lon_lat_dif2, 
                                   by = c("layer_from", "layer_to"))


#- Arrange data frame to create lines between dots and pie (DO IT MANUALLY)

pie_coord<- data.frame(xend = c(-14.57,-14.15,-15.72,-16.70,-17.2,-17.48,-18.1),
                       yend = c(26.18,29.42,28.44,28.9,29.08,28.64,28.33))

segm_connect<- cbind(lon_lat_data, pie_coord)

##Create data for world coordinates using map_data() function
world_coordinates <- map_data("world")

Canary<- world_coordinates %>% 
  filter(long > -18.53 & long < -12.9, lat > 24.5 & lat < 32)#Filter Canary Island region

pdf('./graphs/Islands/Jac/Map_without_modules.pdf', 10, 6)
ggplot() +
  geom_map(data = Canary, map = Canary, aes(long, lat, map_id = region), color = "black", 
           fill = "#c2a18b", size = 0.2) +
  geom_segment(aes(x= x, xend= xend, y= y, yend= yend, color = turnover),
                              data= jaccard_similarity_on_map)+ scale_color_viridis()+ #Lines indicating Jaccard similarity across islands

  geom_segment(aes(x = Lon, y = lat, xend = xend, yend = yend),
               color = "black", data = segm_connect) +  # Lines connecting points and pie charts
  
  geom_point(aes(x= Lon, y= lat), shape = 21, fill= 'white', color= "black", data= lon_lat_data )+ #points indicating islands
  
  xlab("Longitude") +xlim(-18.53, -13.5)+ ylab("Latitude")+ ylim(24.9,30)+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "#F0F8FF",colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.4, linetype = "solid", colour = "black"))+
    theme(axis.title=element_text(size=17))+theme(axis.text.x=element_text(size=13))+
    theme(axis.text.y=element_text(size=13))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(color = "Jaccard Similarity")  
dev.off()

pdf('./graphs/Islands/Jac/Map_without_modules_jaccard.pdf', 10, 6)
ggplot() +
  geom_map(data = Canary, map = Canary, aes(long, lat, map_id = region), color = "black", 
           fill = "#c2a18b", size = 0.2) +
  geom_point(aes(x= Lon, y= lat), shape = 21,size = 4,fill= 'white', color= "black", data= lon_lat_data )+ #points indicating islands
  
  xlab("Longitude") +xlim(-18.53, -13.5)+ ylab("Latitude")+ ylim(24.9,30)+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "#F0F8FF",colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.4, linetype = "solid", colour = "black"))+
  theme(axis.title=element_text(size=17))+theme(axis.text.x=element_text(size=13))+
  theme(axis.text.y=element_text(size=13)) 
dev.off()



##-- Final map containing Jaccard Similarity across Islands and the scatter pies 
pdf('./graphs/Islands/Jac/pie_graph.pdf', 10, 6)
top_30_scatterpies + theme(legend.position = "none")
dev.off()

pdf('./graphs/Islands/Map_modules1.pdf', 10, 6)

ggdraw()+ draw_plot(Canary_map)+
  draw_plot(top_30_scatterpies_no_legend, x = -0.088, y = 0.12, scale = 0.70)
dev.off()



## -- Map with interedges and igraph
setwd("/Users/agustin/Desktop/Papers/Canary_Island_Project/spatial_modularity_in_the_canary_islands")
library(ggplot2)
library(viridis)
#interedges
bid <- read_csv('./csvs/Islands/Jac/dryad_edgelist_complete_ids_islands.csv')

bid_inter <- as_tibble(bid %>% filter(layer_from != layer_to))
bid_inter %>% filter(node_from == 8) %>% print(n=Inf)

layer_connectivity <- 
  bid_inter %>% 
  mutate(group=ifelse(node_from<=39,'plant','animal')) %>% 
  #filter(layer_from %in% c(1,4)) %>% 
  arrange(layer_from) %>% 
  arrange(layer_to) %>% 
  group_by(layer_from, layer_to) %>% 
  #group_by(group) %>% 
  summarise(n=length(node_from),w=mean(weight))
  #summarise(n=length(node_from),w=mean(weight), des = sd(weight))

 thprueba<- layer_connectivity %>%  summarize(total=sum(n))

#combine turnover data frame with coordinates
islands_with_lon_lat_dif <- read.csv("./csvs/islands_with_lon_lat_dif.csv")
lon_lat_data <- read.csv("./csvs/Islands/Jac/lon_lat_data_islands_as_layers.csv")

# Modify islands_with_lon_lat_dif dataframe 
islands_with_lon_lat_dif2<- islands_with_lon_lat_dif %>% 
  mutate(layer_from = recode(layer_from, "1" = "1", "2" = "2", "3"="3", "4east" = "4", "4west" = "5", "5" = "6", "6" = "7"),
         layer_to = recode(layer_to, "1" = "1", "2" = "2", "3"="3", "4east" = "4", "4west" = "5", "5" = "6", "6" = "7"))


#join both data frames
islands_with_lon_lat_side1 <- islands_with_lon_lat_dif2[,c(1,2,4:7)]
islands_with_lon_lat_side2 <- islands_with_lon_lat_dif2[,c(2,1,4:7)]
names(islands_with_lon_lat_side2) <- names(islands_with_lon_lat_side1)
islands_with_lon_lat <- rbind(islands_with_lon_lat_side1,islands_with_lon_lat_side2)

islands_with_lon_lat <- 
  as_tibble(islands_with_lon_lat) %>% 
  mutate_at(1:2, as.integer)
connectivity_on_map <- left_join(layer_connectivity, islands_with_lon_lat)

##Create data for world coordinates using map_data() function
world_coordinates <- map_data("world")

Canary<- world_coordinates %>% 
  filter(long > -18.53 & long < -12.9, lat > 24.5 & lat < 32)#Filter Canary Island region

pdf('./graphs/Islands/Jac/map_interedges.pdf', 10, 6)
  ggplot() +
  geom_map(data = Canary, map = Canary, aes(long, lat, map_id = region), color = "black",fill = "#c2a18b", linewidth = 0.2) +
  geom_segment(data= connectivity_on_map, aes(x= x, xend= xend, y= y, yend= yend, color = n))+
  scale_color_viridis()+ #Lines indicating Jaccard similarity across islands
  
  geom_point(aes(x= Lon, y= lat), shape = 21, fill= 'white', color= "black", data= lon_lat_data )+ #points indicating islands
  
  xlab("Longitude") +xlim(-18.53, -13.5)+ ylab("Latitude")+ ylim(24.9,30)+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "#F0F8FF",colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(linewidth = 0.4, linetype = "solid", colour = "black"))+
  theme(axis.title=element_text(size=17))+theme(axis.text.x=element_text(size=13))+
  theme(axis.text.y=element_text(size=13))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(color = "Interlayer links")  

  dev.off()



## ---  Igraph Map with interedges
library(tidyverse)
library(igraph)

# - Create igraph 
bid <- read_csv('./csvs/Islands/Jac/dryad_edgelist_complete_ids_islands.csv')

bid_inter <- as_tibble(bid %>% filter(layer_from != layer_to))
bid_inter %>% filter(node_from == 8) %>% print(n=Inf)

layer_connectivity <- 
  bid_inter %>% 
  mutate(group=ifelse(node_from<=39,'plant','animal')) %>% 
  # filter(layer_from %in% c(1,4)) %>% 
  arrange(layer_from) %>% 
  arrange(layer_to) %>% 
  group_by(layer_from,layer_to,group) %>% 
  summarise(n=length(node_from),w=mean(weight))

layer_connectivity <- as.data.frame(layer_connectivity)
colnames(layer_connectivity) <- c('from', 'to', 'group', 'n', 'w')

g <- graph_from_data_frame(layer_connectivity, directed = F)
list.edge.attributes(g)
E(g)$w
E(g)$color[E(g)$group=='plant'] <- "#008000"
E(g)$color[E(g)$group=='animal'] <- "#BE75FA" 

#coordinates
lon_lat_data <- read.csv("./csvs/Islands/Jac/lon_lat_data_islands_as_layers.csv")
coordinates_igraph <-  lon_lat_data[,c(3,2)]

coordinates_igraph[1,1]<--12.8
coordinates_igraph[1,2]<-27.1
coordinates_igraph[2,1]<--12.4
coordinates_igraph[2,2]<-28.3
coordinates_igraph[3,1]<--14.2
coordinates_igraph[3,2]<-27.9
coordinates_igraph[4,1]<--15.8
coordinates_igraph[4,2]<-28.2
coordinates_igraph[6,1]<--18
coordinates_igraph[7,1]<--20
coordinates_igraph[7,2]<-27.8


pdf('./graphs/Islands/Jac/interedges_igraph.pdf', 10, 6)
plot.igraph(g, edge.width = E(g)$n*0.55,layout=as.matrix(coordinates_igraph))
dev.off()




