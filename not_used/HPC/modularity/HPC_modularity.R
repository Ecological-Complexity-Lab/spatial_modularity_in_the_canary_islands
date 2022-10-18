# bash intro
#! /gpfs0/shai/projects/R4/R-4.2.0/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.2.0/lib64/R/library")
print(.libPaths())
print(sessionInfo())
Sys.setlocale(category = "LC_ALL", locale = "")

library(readr)
library(tidyverse)
library(magrittr)
library(infomapecology)

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
  write_lines(s, paste("logs_modularity_pollinators/", JOB_ID, "_", trial_num,'_log.txt', sep=''), append = appnd)
}

# ------------- run --------------
LOG('START RUN LOG', appnd = F)
LOG('====================\n')

##-------------------------------------------------
dryad_edgelist_complete_shuf_pols <- read.csv("./dryad_edgelist_complete_shuf_pols.csv")
dryad_edgelist_complete_shuf_plants <- read.csv("./dryad_edgelist_complete_shuf_plants.csv")
dryad_edgelist_complete_shuf_both <- read.csv("./dryad_edgelist_complete_shuf_both.csv")
physical_nodes <- read.csv("./physical_nodes.csv")
layer_metadata <- read.csv("./layer_metadata.csv")

dryad_multilayer_shuf_1000_pols <- NULL
dryad_multilayer_shuf_1000_plants <- NULL
dryad_multilayer_shuf_1000_both <- NULL

modularity_for_shuf <- function(edge_list, output){
    #current_trial_edgelist <- dryad_edgelist_complete_shuf %>% filter(trial_number == i) #take 1 trial at a time to create multilayer
    current_trial_edgelist <- edge_list %>% filter(trial_number == trial_num) #take 1 trial at a time to create multilayer
    dryad_multilayer_shuf_trial <- create_multilayer_object(extended = current_trial_edgelist, #taking edge list and returning multilayer network
                                                            nodes = physical_nodes,
                                                            layers = layer_metadata,
                                                            intra_output_extended = T)
    
    
    # Input: intra non-extended edge lists and inter extended edge list
    intra_nonextended_shuf_500 <-
      current_trial_edgelist %>% 
      filter(layer_from==layer_to) %>% #only intra
      dplyr::select(layer=layer_from, node_from, node_to, weight)
    inter_extended_shuf_500 <-
      current_trial_edgelist %>% 
      filter(layer_from!=layer_to) #only inter
    
    
    #write.csv(intra_nonextended_shuf_500, "./csvs/intra_nonextended_shuf_500.csv")
    #write.csv(inter_extended_shuf_500, "./csvs/inter_extended_shuf_500.csv")
    
    #create modules for empirical network
    modules_dryad_multilayer_shuf_500 <- run_infomap_multilayer(dryad_multilayer_shuf_trial, 
                                                                infomap_executable = "Infomap",
                                                                flow_model = 'directed',
                                                                relax = F, 
                                                                silent = T, 
                                                                trials = 100,
                                                                seed = 497294, 
                                                                temporal_network = F)
    
    
    #dryad_multilayer_shuf_500 <- rbind(dryad_multilayer_shuf_500, tibble(modules_dryad_multilayer_shuf_500$modules, i)) #save modules for 500 null models
    output <- rbind(output, tibble(modules_dryad_multilayer_shuf_500$modules, trial_num)) #save modules for 500 null models
    
  return(output)
}

#pols
dryad_multilayer_shuf_1000_pols_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_pols, 
                                                              dryad_multilayer_shuf_1000_pols)

write.csv(dryad_multilayer_shuf_1000_pols_output, paste("csvs_modularity_pollinators/", trial_num,'_new_modularity.csv', sep=''), 
          row.names = FALSE)

#plants
dryad_multilayer_shuf_1000_plants_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_plants, 
                                                                dryad_multilayer_shuf_1000_plants)

write.csv(dryad_multilayer_shuf_1000_plants_output, paste("csvs_modularity_plants/", trial_num,'_new_modularity.csv', sep=''), 
          row.names = FALSE)

#both
dryad_multilayer_shuf_1000_both_output <- modularity_for_shuf(dryad_edgelist_complete_shuf_both, 
                                                              dryad_multilayer_shuf_1000_both)

write.csv(dryad_multilayer_shuf_1000_both_output, paste("csvs_modularity_both/", trial_num,'_new_modularity.csv', sep=''), 
          row.names = FALSE)





LOG('\n====================')
LOG('FINISH RUN LOG')