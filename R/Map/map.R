## ---- Code to create map showing structure similarity between islands (Fig. 2)


library(dplyr)
library(ggplot2)
library(stringr)
library(colorspace)
library(tidyselect)

##----get_data--------------------------------------------------------------------------------------------------------
setwd("D:/Trabajo/Papers/Canary_Island/spatial_modularity_in_the_canary_islands")
source("D:/Trabajo/Papers/Canary_Island/spatial_modularity_in_the_canary_islands/R/functions.R")

modules_dryad_multilayer <- read.csv("./csvs_nuevo/module_data_justislands_as_layers.csv") 
module_data_with_loc <- read.csv("./csvs_nuevo/modules_with_lat_lon_justislands_as_layers.csv")
lon_lat_data <- read.csv("./csvs_nuevo/lon_lat_data_justislands_as_layers.csv")

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
pie_chart_data2 <- pie_chart_data2 %>% group_by(lat, long) %>% summarise(across(everything(), sum))


#Plot
top_30_scatterpies <- ggplot()+ geom_scatterpie(aes(x=long, y=lat, group=layer_id, r=0.2), alpha=0.90, data= pie_chart_data2, 
                                                cols=colnames(pie_chart_data2[,c(4:11)]))+
  theme_void()+ coord_fixed()+ theme(legend.text = element_text(size = 13), legend.title = element_text(size = 13))
top_30_scatterpies


##-- Map of Canary Island region with Jaccard Similarity (Fig. 4)

#combine turnover data frame with coordinates
empirical_turnover_for_module_island_shuf_no_self_loop_km <- read.csv("./csvs_nuevo/justislands_turnover_with_distnace_empirical.csv")
islands_with_lon_lat_dif <- read.csv("./csvs_nuevo/islands_with_lon_lat_dif.csv")

#Average Jaccard
Jaccard_ave<- empirical_turnover_for_module_island_shuf_no_self_loop_km %>% 
  summarize(Prop = mean(turnover), stan = sd(turnover))

# Modify islands_with_lon_lat_dif dataframe 
islands_with_lon_lat_dif2<- islands_with_lon_lat_dif %>% filter(layer_from != "1") %>% #remove mainland
  mutate(layer_from = recode(layer_from, "2" = "1", "3" = "2", "4east"="3", "4west" = "3", "5" = "4", "6" = "5"),
         layer_to = recode(layer_to, "2" = "1", "3" = "2", "4east"="3", "4west" = "3", "5" = "4", "6" = "5")) #rearrange name islands


#join both data frames
jaccard_similarity_on_map <- merge(empirical_turnover_for_module_island_shuf_no_self_loop_km, islands_with_lon_lat_dif2, 
                                   by = c("layer_from", "layer_to"))


#- Arrange data frame to create lines between dots and pie (DO IT MANUALLY)

pie_coord<- data.frame(xend = c(-14.15,-15.72,-16.70,-17.48,-18.1),
                       yend = c(29.42,28.44,28.9,28.64,28.33))

segm_connect<- cbind(lon_lat_data, pie_coord)

##Create data for world coordinates using map_data() function
world_coordinates <- map_data("world")

Canary<- world_coordinates %>% 
  filter(long > -18.53 & long < -12.9, lat > 27.5 & lat < 29.5)#Filter Canary Island region

Canary_map<-ggplot() +
  geom_map(data = Canary, map = Canary, aes(x=long, y=lat, map_id = region), color = "black", 
           fill = "#c2a18b", size = 0.2) +
  geom_segment(aes(x= x, xend= xend, y= y, yend= yend, color = turnover),
                              data= jaccard_similarity_on_map)+ scale_color_viridis()+ #Lines indicating Jaccard similarity across islands

  geom_segment(aes(x = long, y = lat, xend = xend, yend = yend),
               color = "black", data = segm_connect) +  # Lines connecting points and pie charts
  
  geom_point(aes(x= long, y= lat), shape = 21, fill= 'white', color= "black", data= lon_lat_data )+ #points indicating islands
  
  xlab("Longitude") +xlim(-18.53, -13.5)+ ylab("Latitude")+ ylim(27.5,29.4)+
  theme(panel.border = element_blank(),
        panel.background = element_rect(fill = "#F0F8FF",colour = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.4, linetype = "solid", colour = "black"))+
    theme(axis.title=element_text(size=17))+theme(axis.text.x=element_text(size=13))+
    theme(axis.text.y=element_text(size=13))+ theme(legend.title = element_text(size = 13), legend.text = element_text(size = 13))+
  labs(color = "Jaccard Similarity")  


##-- Final map containing Jaccard Similarity across Islands and the scatter pies 
top_30_scatterpies_no_legend<-top_30_scatterpies + theme(legend.position = "none")

pdf('./graphs/Fig_2.pdf', 10, 8)
ggdraw()+ draw_plot(Canary_map)+
  draw_plot(top_30_scatterpies_no_legend, x = -0.088, y = 0.12, scale = 0.70)
dev.off()
#Then, rearrange manually with adobe to improve the quality of the image



##################### FIGURE S2

## ---  Igraph Map with interedges
library(tidyverse)
library(igraph)

# - Create igraph 
bid <- read_csv('./csvs_nuevo/dryad_edgelist_complete_ids_justislands.csv')

bid_inter <- as_tibble(bid %>% filter(layer_from != layer_to))

layer_connectivity <- 
  bid_inter %>% 
  mutate(group=ifelse(node_from<=39,'plant','animal')) %>% 
  # filter(layer_from %in% c(1,4)) %>% 
  arrange(layer_from) %>% 
  arrange(layer_to) %>% 
  group_by(layer_from,layer_to,group) %>% 
  summarise(n=length(node_from)*2,w=mean(weight)) # *2 because its directed

layer_connectivity <- as.data.frame(layer_connectivity)
colnames(layer_connectivity) <- c('from', 'to', 'group', 'n', 'w')

g <- graph_from_data_frame(layer_connectivity, directed = F)
list.edge.attributes(g)
E(g)$w
E(g)$color[E(g)$group=='plant'] <- "#008000"
E(g)$color[E(g)$group=='animal'] <- "#BE75FA" 


#coordinates
lon_lat_data <- read.csv("./csvs_nuevo/lon_lat_data_justislands_as_layers.csv")
coordinates_igraph <-  lon_lat_data[,c(3,2)]


coordinates_igraph[3,1]<--18
coordinates_igraph[4,1]<--19
coordinates_igraph[5,1]<--20



pdf('./graphs/fig_S2.pdf', 10, 6)
plot.igraph(g, edge.width = E(g)$w*17,layout=as.matrix(coordinates_igraph))
dev.off()
