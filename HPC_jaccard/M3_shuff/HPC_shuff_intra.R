# bash intro
#! /gpfs0/shai/projects/R4/R-4.2.0/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.2.0/lib64/R/library")
print(.libPaths())
print(sessionInfo())
Sys.setlocale(category = "LC_ALL", locale = "")

library(readr)
library(tidyverse)
library(magrittr)

# ------------- parsing arguments -----------
# read args given in command line:
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) { # make sure we have commands
  stop('No arguments were found!') # the script will not run without arguments
} else {
  args <- commandArgs(trailingOnly=TRUE)
  trial_num <- as.numeric(args[1])
}

# ------------ functions ---------------------------
LOG <- function(s, appnd = T) {
  write_lines(s, paste("logs_both/", JOB_ID, "_", trial_num,'_log.txt', sep=''), append = appnd)
}

# ------------- run --------------
LOG('START RUN LOG', appnd = F)
LOG('====================\n')

# this script uses the received arguments to produce a calculation
intralayers_with_ids_for_shuf <- read.csv("./intralayers_with_ids_for_shuf.csv")

##---- classic shuffle within layers ---------------------------------------------
shuf_trial_matrix_classic <- NULL
shuf_null_edge_list_classic <- NULL

for (i in 1:7){
  current_matrix <- intralayers_with_ids_for_shuf %>%
    filter(layer_from == layer_to) %>% filter(layer_from == i) %>% #take 1 layer at a time
    select(node_from, node_to, sum_weight) %>% #select interactions
    dcast(node_from ~ node_to, value.var = "sum_weight", fill = 0) %>%
    column_to_rownames(var="node_from") 
  
  current_matrix <- t(current_matrix) #put pols in rows
  
  print(i) #to keep tab on which layer we're on
  
  for (j in 1:1000){ #1000 iterations
    null <- vegan::nullmodel(current_matrix, method = 'r00_samp') #shuffle within island
    shuffled_matrices <- simulate(null, nsim = 1)
    
    trial_with_shuf_num <- cbind(shuffled_matrices, j) #add trial number 
    shuf_trial_matrix_classic <- rbind(shuf_trial_matrix_classic, trial_with_shuf_num) #create big matrix of all matrices
    edge_list_version_classic <- melt(as.matrix(shuffled_matrices)) %>% filter(value > 0) %>%
      select(node_from=Var1, node_to=Var2, weight=value) #turn the matrix back into a data frame of edge lists
    edge_list_version_classic$trial_number <- j #add trial number to the edge list
    edge_list_version_classic$layer_from <- i #intra so layer from and to are identical
    edge_list_version_classic$layer_to <- i
    shuf_null_edge_list_classic <- rbind(shuf_null_edge_list_classic, edge_list_version_classic) #create mega edge list with all repetitions
  }
}

shuf_null_edge_list_classic <- shuf_null_edge_list_classic[, c(5,1,6,2,3,4)] #change to regular order


write.csv(shuf_null_edge_list_classic, paste("csvs_intra/", trial_num,'_intra_new_shuf_islands_as_layers.csv', sep=''), row.names = FALSE)

LOG('\n====================')
LOG('FINISH RUN LOG')











# this script uses the received arguments to produce a calculation
interlayer_edges_shuf_both <- read.csv("./interlayer_edges_shuf_both_islands_as_layers.csv")

interlayers_shuf_both <- NULL

interlayer_edges_run <- interlayer_edges_shuf_both %>% filter(trial_number == trial_num)
for (i in 1:nrow(interlayer_edges_run)){ #for every species
  current_run <- interlayer_edges_run[i,]
  for (j in 4:ncol(interlayer_edges_run)){ #for all locations
    current_species <- current_run$node_from
    current_location <- current_run$layer_from
    current_location_to <- current_run[j]
   # if (current_species %in% interlayers_shuf$node_from & current_location %in% interlayers_shuf$layer_to &
    #    current_location_to %in% interlayers_shuf$layer_from) next
    #else{interlayers_shuf <- rbind(interlayers_shuf, tibble(trial_number = trial_num, layer_from = current_location,
    #                                                        node_from = current_species, 
    #                                                        layer_to = as.character(current_location_to), 
    #                                                        node_to = current_species))
    interlayers_shuf_both <- rbind(interlayers_shuf_both, tibble(trial_number = trial_num, layer_from = current_location,
                                                            node_from = current_species, 
                                                            layer_to = as.character(current_location_to), 
                                                            node_to = current_species))
   # } 
  }
}

interlayers_new_shuf_both <- interlayers_shuf_both %>% subset(layer_to != "NA")
interlayers_new_shuf_both <- interlayers_new_shuf_both %>% subset(layer_from != layer_to)

write.csv(interlayers_new_shuf_both, paste("csvs_both/", trial_num,'_interlayers_new_shuf_both_islands_as_layers.csv', sep=''), row.names = FALSE)

LOG('\n====================')
LOG('FINISH RUN LOG')
