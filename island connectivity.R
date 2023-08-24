setwd('~/Downloads/')
library(tidyverse)
library(igraph)

bid <- read_csv('dryad_edgelist_complete_ids_islands_directed.csv')
unid <- read_csv('dryad_edgelist_complete_ids_islands_nondirected.csv')

bid_inter <- as_tibble(bid %>% filter(layer_from != layer_to))
bid_inter %>% filter(node_from == 8) %>% print(n=Inf)

unid_inter <- as_tibble(unid %>% filter(layer_from != layer_to))
unid_inter %>% filter(node_from == 8)

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

g <- graph_from_data_frame(layer_connectivity, directed = T)
list.edge.attributes(g)
E(g)$w
E(g)$color[E(g)$group=='plant'] <- 'red'
E(g)$color[E(g)$group=='animal'] <- 'blue'
plot(g, edge.width=E(g)$n, layout=layout_in_circle)



