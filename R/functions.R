

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



##---- modularity_for_shuf: create modules for shuffled null models

# this functions takes edge lists of shuffled matrices and creates a data frame of modularity 
modularity_for_shuf <- function(edge_list, output, range=1:1000){
  for(trial_num in range){
    print(trial_num) #to keep tab on how far along we are
    current_trial_edgelist <- edge_list %>% filter(trial_number == trial_num) #take 1 trial at a time to create multilayer
    dryad_multilayer_shuf_trial <- create_multilayer_object(extended = current_trial_edgelist, #edge list changes every run (shuffled)
                                                            nodes = physical_nodes, #nodes are always the same. we're not deleting nodes.
                                                            layers = layer_metadata, #layers are always the same. we're not deleting layers.
                                                            intra_output_extended = T)
    
    
    #create modules for empirical network
    modules_dryad_multilayer_shuf_1000 <- modified_multi(dryad_multilayer_shuf_trial, 
                                                                infomap_executable = "Infomap",
                                                                flow_model = 'directed',
                                                                relax = F, 
                                                                silent = T, 
                                                                trials = 100,
                                                                seed = 497294, #always the same seed
                                                                temporal_network = F)
    
    
    output <- rbind(output, tibble(modules_dryad_multilayer_shuf_1000$modules, trial_num)) #save modules for 1000 null models
  }  
  return(output)
}


##---- map creation-------------------------------------------------------------------------------- 

##---- pivot by country: this function takes the data and creates a new version of layer_id and num of state nodes in it
pivot_by_country <- function(data) {  
  s1 = melt(data, id = c("layer_id", "group"), measure.vars = "size_of_module")
  s2 = dcast(s1, layer_id ~ group, length)
  s2$Total = rowSums(s2[,2:NCOL(s2)])
  return(s2)
}




#---- infomap modifiled multi if error -i occures----------------------
modified_multi <- function (M, infomap_executable = "Infomap", flow_model = NULL, 
silent = T, trials = 100, seed = NULL, relax = F, multilayer_relax_rate = 0.1, 
multilayer_relax_limit = NULL, multilayer_relax_limit_up = NULL, 
multilayer_relax_limit_down = NULL, temporal_network = F, 
run_standalone = T, remove_auxilary_files = T, ...) 
{
  if (check_infomap(infomap_executable) == F) {
    stop("Error in Infomap stand-alone file.")
  }
  if (class(M) != "multilayer") {
    stop("M must be of class multilayer")
  }
  arguments <- paste("multilayer --tree -2 -N ", trials, 
                     sep = "")
  arguments <- ifelse(!is.null(seed), paste(arguments, "--seed", 
                                            seed), arguments)
  arguments <- ifelse(!is.null(flow_model), paste(arguments, 
                                                  "-f", flow_model), arguments)
  arguments <- ifelse(silent, paste(arguments, "--silent"), 
                      arguments)
  arguments <- paste(arguments, ...)
  if (relax == F) {
    print("Using interlayer edge values to determine flow between layers.")
    write_lines("*Multilayer", "infomap_multilayer.txt")
    write_delim(M$intra, "infomap_multilayer.txt", delim = " ", 
                append = T)
    write_delim(M$inter, "infomap_multilayer.txt", delim = " ", 
                append = T)
  }
  else {
    if (ncol(M$intra) == 5) {
      stop("Cannot use relax rates with extended format of intralayer edges. See function create_multilayer_object.")
    }
    print("Using global relax to determine flow between layers.")
    write_lines("*Intra", "infomap_multilayer.txt")
    write_delim(M$intra, "infomap_multilayer.txt", delim = " ", 
                append = T)
    if (!is.null(M$inter)) {
      if (ncol(M$inter) == 5) {
        stop("Cannot use relax rates with extended format of interlayer edges. See function create_multilayer_object.")
      }
      print("Global relax will be constrained by interlayer edges.")
      write_lines("*Inter", "infomap_multilayer.txt", append = T)
      write_delim(M$inter, "infomap_multilayer.txt", delim = " ", 
                  append = T)
    }
    arguments <- ifelse(!is.null(multilayer_relax_rate), 
                        paste(arguments, "--multilayer-relax-rate", multilayer_relax_rate), 
                        arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit), 
                        paste(arguments, "--multilayer-relax-limit", multilayer_relax_limit), 
                        arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit_up), 
                        paste(arguments, "--multilayer-relax-limit-up", multilayer_relax_limit_up), 
                        arguments)
    arguments <- ifelse(!is.null(multilayer_relax_limit_down), 
                        paste(arguments, "--multilayer-relax-limit-down", 
                              multilayer_relax_limit_down), arguments)
  }
  call <- paste("./", infomap_executable, " infomap_multilayer.txt . ", 
                arguments, sep = "")
  if (run_standalone == T) {
    print(call)
    system(call)
  }
  else {
    print("Please run Infomap online at https://www.mapequation.org/infomap/ using the following arguments (copy-paste):")
    print(arguments)
    invisible(readline(prompt = "After running, download statenodes results and press [ENTER] when done"))
    if (!file.exists("network_states.tree")) {
      stop("Result file network_states.tree was not found. Did you download results?")
    }
    file.rename(from = "network_states.tree", to = "infomap_multilayer_states.tree")
  }
  L_output <- parse_number(read_lines("infomap_multilayer_states.tree")[6])
  modules <- suppressMessages(read_delim("infomap_multilayer_states.tree", 
                                         delim = " ", skip = 11, col_names = c("path", "flow", 
                                                                               "name", "state_id", "node_id", "layer_id")))
  modules %<>% filter(flow > 0) %>% dplyr::select(path, node_id, layer_id, 
                                                  flow) %>% separate(path, into = c("module", "leaf_id"), 
                                                                     sep = ":") %>% mutate_at(.vars = 1:4, as.integer) %>% 
    left_join(M$nodes, "node_id") %>% dplyr::select(node_id, starts_with("module"), 
                                                    everything(), -leaf_id) %>% dplyr::arrange(node_id, layer_id)
  if (temporal_network) {
    print("Reorganizing modules...")
    renamed_moduels <- modules %>% distinct(module, layer_id) %>% 
      arrange(module, layer_id)
    x <- c(1, table(renamed_moduels$module))
    module_birth_layers <- renamed_moduels %>% slice(cumsum(x)) %>% 
      arrange(layer_id, module)
    module_renaming <- data.frame(module = module_birth_layers$module, 
                                  module_renamed = 1:max(module_birth_layers$module))
    modules %<>% left_join(module_renaming, "module") %>% 
      dplyr::select(-module) %>% rename(module = module_renamed)
  }
  if (remove_auxilary_files) {
    print("Removing auxilary files...")
    file.remove("infomap_multilayer_states.tree")
    file.remove("infomap_multilayer.txt")
    file.remove("infomap_multilayer.tree")
  }
  print(paste("Partitioned into ", max(modules$module), " modules.", 
              sep = ""))
  out <- list(call = call, L = L_output, m = max(modules$module), 
              modules = modules)
  class(out) <- "infomap_multilayer"
  return(out)
}

