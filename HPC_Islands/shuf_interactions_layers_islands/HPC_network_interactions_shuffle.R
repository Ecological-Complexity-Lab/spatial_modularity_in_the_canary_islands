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
  trial_number <- as.numeric(args[1])
}

# ------------ functions ---------------------------
LOG <- function(s, appnd = T) {
  write_lines(s, paste("logs_interactions/", JOB_ID, "_", trial_number,'_log.txt', sep=''), append = appnd)
}

# ------------- run --------------
LOG('START RUN LOG', appnd = F)
LOG('====================\n')

# this script uses the received arguments to produce a calculation
interlayer_edges_interactions_shuf <- read.csv("./interlayer_edges_interactions_shuf.csv")


interlayers_interactions_shuf <- NULL
interlayer_edges_run <- interlayer_edges_interactions_shuf %>% filter(trial_num == trial_number) #each HPC run gets one trial

for (i in 1:nrow(interlayer_edges_run)){ #for every species
  current_run <- interlayer_edges_run[i,]
  for (j in 4:ncol(interlayer_edges_run)){ #for all locations
    current_species <- current_run$node_from
    current_location <- current_run$layer_from
    current_location_to <- current_run[j]
    interlayers_interactions_shuf <- rbind(interlayers_interactions_shuf, tibble(layer_from = current_location,
                                                                                 node_from = current_species, 
                                                                                 layer_to = as.character(current_location_to), 
                                                                                 node_to = current_species,
                                                                                 trial_number = trial_number))
  }
}

interlayers_interactions_shuf_new <- interlayers_interactions_shuf %>% subset(layer_to != "NA")
interlayers_interactions_shuf_new <- interlayers_interactions_shuf_new %>% subset(layer_from != layer_to)

write.csv(interlayers_interactions_shuf_new, paste("csvs_interactions/", trial_number,'_interlayers_interactions_shuf_new.csv', sep=''), row.names = FALSE)

LOG('\n====================')
LOG('FINISH RUN LOG')
