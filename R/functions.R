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
#library(taxize)
#library(taxizedb)
#library(ggtree)
library(ggpubr)



#---- comparing coefficients between empirical and null model----------------
# hierarchical and non-hierarchical analysis

# this function takes the coefficients and the sd of two data sets and compares them to see how different the slopes are
compare.coeff <- function(b_1,se_1,b_2,se_2){ #b1, b2 are coefficients 
  return((b_1-b_2)/sqrt(se_1^2+se_2^2)) #se1, se2 are standard errors
}


#---- modularity analysis----------------------------------------------------
# run on empirical network and on shuffled networks

#---- pivot by module

# this function flips the data frame to look at it from a layer perspective
pivot_by_module <- function(data){ #creates a data frame with module on the side and layer_id on the top
  s1 = melt(data, id = c("layer_id", "module"))
  s2 = dcast(s1, layer_id ~ module, length)
  s3 = t(s2) 
  s3 <- s3[-1,]
  colnames(s3) <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
  return(s3)
}


##---- edge list per module

# this function creates an edge list of layers found in a module with distances
edge_list_per_module <- function(data,edge_list){
  #gets one row from a data frame and creates an edge list from it
  for (i in (1:13)){
    if (data[i]==0) next #only take layers where the module is present
    else {
      for (j in ((i+1):14)){
        if (data[j]==0) next #only take layers where the module is present
        else {
          edge_list <- rbind(edge_list, tibble(layer_from=i, layer_to=j, module=NA)) #create edge list of all the layer found in a module
        }
      }
    }
  }
  return(edge_list)
}



##---- species distance decay for shuffled null models

# this function calculates Jaccard Similarity between every two layers
#jaccard similarity in species
distnace_decay_shuf <- function(species_all_layers, output){
  for (k in 1:1000){ #1000 iterations
    current_trial <- species_all_layers %>% filter(trial_number == k)
    print(k) #to keep tab on how far along we are
    for (i in (1:13)){
      for (j in ((i+1):14)){
        physical_nodes_in_layer_from <- filter(current_trial, layer_id == i) %>% select(node_id) %>% unlist() #all species in layer from
        physical_nodes_in_layer_to <- filter(current_trial, (layer_id == j)) %>% select(node_id) %>% unlist() #all species in layer to
        int_both <- intersect(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are common to both layers
        uni_both <- union(physical_nodes_in_layer_from, physical_nodes_in_layer_to) #how many nodes are found in both layers total
        turnover <- length(int_both)/length(uni_both)
        output <- rbind(output, tibble(trial_number = k, layer_from = i, layer_to = j, turnover = turnover))
      }
    }
  }
  
  output <- output %>% unique()
  
  return(output)
}


##---- create modules for shuffled null models

# this functions takes edge lists of shuffled matrices and creates a data frame of modularity 
modularity_for_shuf <- function(edge_list, output){
  for(trial_num in 1:1000){
    print(trial_num) #to keep tab on how far along we are
    current_trial_edgelist <- edge_list %>% filter(trial_number == trial_num) #take 1 trial at a time to create multilayer
    dryad_multilayer_shuf_trial <- create_multilayer_object(extended = current_trial_edgelist, #edge list changes every run (shuffled)
                                                            nodes = physical_nodes, #nodes are always the same. we're not deleting nodes.
                                                            layers = layer_metadata, #layers are always the same. we're not deleting layers.
                                                            intra_output_extended = T)
    
    
    #create modules for empirical network
    modules_dryad_multilayer_shuf_1000 <- run_infomap_multilayer(dryad_multilayer_shuf_trial, 
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


##---- module distance decay for shuffled null models between islands


## this function calculates Jaccard Similarity in modules between islands
module_distance_decay_func <- function(multilayer_1000, 
                                       islands_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial) #to keep tab on how far along we are

    modules_for_similarity_shuf <- multilayer_1000 %>% filter(trial_num == trial) #take only 1 trial each run

    #pivot modules
    module_pivoted_shuf <- pivot_by_module(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    #create edge list with distances
    modules_edge_list_shuf <- NULL
    
    for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[k]
      if (is.null(modules_edge_list_shuf)) next
      modules_edge_list_shuf <- modules_edge_list_shuf %>% mutate(module = replace_na(module, current_module)) #add module number
    }
    
    
    edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
    edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA 
    
    #modules similarity pairwise distance between islands
    #change layer number to island identity (aggregate by island)
    edge_list_by_islands_shuf <- edge_list_with_distances_shuf
    old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
    new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
    edge_list_by_islands_shuf$layer_from[edge_list_by_islands_shuf$layer_from %in% old] <- new[match(edge_list_by_islands_shuf$layer_from, old)]
    edge_list_by_islands_shuf$layer_to[edge_list_by_islands_shuf$layer_to %in% old] <- new[match(edge_list_by_islands_shuf$layer_to, old)]
    
    #version with number of modules in layers
    edge_list_by_islands_modules_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to, module) %>%
      summarise(ave_distance= mean(distance_in_meters))
    edge_list_by_islands_modules_shuf$count <- c(1)
    edge_list_by_islands_modules_shuf <- edge_list_by_islands_modules_shuf %>% mutate(number_of_modules= sum(count)) %>%
      select(layer_from, layer_to, module, number_of_modules) 

    ## jaccard similairty shuf
    #total number of modules in each layer
    module_island_turnover_shuf <- NULL
    island_list <- c("1","2","3","4east","4west","5","6")
    
    for (i in island_list){
      for (j in island_list){
        modules_in_island_from_shuf <- filter(edge_list_by_islands_modules_shuf, layer_from == i) %>% select(module) %>% unique() %>% unlist() #all modules in island from
        modules_in_island_to_shuf <- filter(edge_list_by_islands_modules_shuf, layer_from == j) %>% select(module) %>% unique() %>% unlist() #all modules in island to
        int_both <- intersect(modules_in_island_from_shuf, modules_in_island_to_shuf) #how many modules are common in both islands
        uni_both <- union(modules_in_island_from_shuf, modules_in_island_to_shuf) #how many modules are found in both islands in total
        turnover <- length(int_both)/length(uni_both)
        module_island_turnover_shuf <- rbind(module_island_turnover_shuf, tibble(layer_from= i, layer_to= j, turnover= turnover))
      }
    }
    
    edge_list_by_islands_ave_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to) %>%
      summarise(ave_distance= mean(distance_in_meters)) %>% unique() #create an average distance as every two layers in an island are 50-500 meters apart
    
    islands_turnover_with_distnace <- edge_list_by_islands_ave_shuf %>%
      merge(module_island_turnover_shuf, by= c("layer_from", "layer_to")) #merge both versions
    
    islands_turnover_with_distnace_all_trials <- rbind(islands_turnover_with_distnace_all_trials, 
                                                       tibble(islands_turnover_with_distnace, trial_num = trial))
    
  }
  return(islands_turnover_with_distnace_all_trials)
}


##----number of shared modules between islands

# this function calculates the number of shared modules every two islands share
number_of_shared_modules_func <- function(multilayer_1000, 
                                          all_edge_list_island_combine_no_module_shuf){
  for (trial in 1:1000){
    print(trial) #to keep tab on how far along we are
    
    modules_for_similarity_shuf <- multilayer_1000 %>% filter(trial_num == trial) #take only 1 trial
    #pivot modules
    module_pivoted_shuf <- pivot_by_module(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    
    #create edge list with distances
    modules_edge_list_shuf <- NULL

    for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[k]
      if (is.null(modules_edge_list_shuf)) next
      modules_edge_list_shuf <- modules_edge_list_shuf %>% mutate(module = replace_na(module, current_module)) #add module number
    }
    
    
    edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
    edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA 
    
    #modules similarity pairwise distance between islands
    #change layer number to island identity (aggregate by island)
    edge_list_by_islands_shuf <- edge_list_with_distances_shuf
    old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
    new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
    edge_list_by_islands_shuf$layer_from[edge_list_by_islands_shuf$layer_from %in% old] <- new[match(edge_list_by_islands_shuf$layer_from, old)]
    edge_list_by_islands_shuf$layer_to[edge_list_by_islands_shuf$layer_to %in% old] <- new[match(edge_list_by_islands_shuf$layer_to, old)]
    
    #version with number of modules in islands
    edge_list_by_islands_modules_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to, module) %>%
      summarise(ave_distance= mean(distance_in_meters)) 
    edge_list_by_islands_modules_shuf$count <- c(1)
    edge_list_by_islands_modules_shuf <- edge_list_by_islands_modules_shuf %>% mutate(number_of_modules= sum(count)) %>%
      select(layer_from, layer_to, module, number_of_modules) 
    
    
    #version with correct average between islands
    edge_list_by_islands_ave_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to) %>%
      summarise(ave_distance= mean(distance_in_meters)) %>% unique()  #create an average distance as every two layers in an island are 50-500 meters apart
    
    #combine
    edge_list_island_combine_shuf <- edge_list_by_islands_ave_shuf %>%
      merge(edge_list_by_islands_modules_shuf, by= c("layer_from", "layer_to")) #merge both versions 
    
    edge_list_island_combine_no_module_shuf <- edge_list_island_combine_shuf %>% select(-module) %>% unique() #have version where modules aren't present
    
    all_edge_list_island_combine_no_module_shuf <- rbind(all_edge_list_island_combine_no_module_shuf, 
                                                         tibble(layer_from = edge_list_island_combine_no_module_shuf$layer_from,
                                                                layer_to = edge_list_island_combine_no_module_shuf$layer_to,
                                                                ave_distance = edge_list_island_combine_no_module_shuf$ave_distance,
                                                                number_of_modules = edge_list_island_combine_no_module_shuf$number_of_modules,
                                                                trial = trial))
    
    
  }
  return(all_edge_list_island_combine_no_module_shuf)
}


##---- module distance decay shuffled null models between layers

# this function calculates the Jaccard Similarity in modules between layers
module_distance_decay_layer_func <- function(multilayer_1000, 
                                             layers_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial) #to keep tab on how far along we are

    modules_for_similarity_shuf <- multilayer_1000 %>% filter(trial_num == trial) #take only 1 trial

    #pivot modules
    module_pivoted_shuf <- pivot_by_module(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    
    #create edge list with distances
    modules_edge_list_shuf <- NULL

    for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[k]
      if (is.null(modules_edge_list_shuf)) next
      modules_edge_list_shuf <- modules_edge_list_shuf %>% mutate(module = replace_na(module, current_module)) #add module number
    }
    

    edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
    edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA

    #version with # of modules in layers
    edge_list_by_layers_modules_shuf <- edge_list_with_distances_shuf %>% group_by(layer_from, layer_to, module) %>%
      summarise(ave_distance= mean(distance_in_meters)) 
    edge_list_by_layers_modules_shuf$count <- c(1)
    edge_list_by_layers_modules_shuf <- edge_list_by_layers_modules_shuf %>% mutate(number_of_modules= sum(count)) %>%
      select(layer_from, layer_to, module, number_of_modules) 
    
    
    for (i in 1:14){
      for (j in 1:14){
        modules_in_layer_from_shuf <- filter(edge_list_by_layers_modules_shuf, layer_from == i) %>% select(module) %>% unique() %>% unlist() #modules in layer from
        modules_in_layer_to_shuf <- filter(edge_list_by_layers_modules_shuf, layer_from == j) %>% select(module) %>% unique() %>% unlist() #modules in layer to
        int_both <- intersect(modules_in_layer_from_shuf, modules_in_layer_to_shuf) #how many modules are common in both layers
        uni_both <- union(modules_in_layer_from_shuf, modules_in_layer_to_shuf) #how many modules are found in both layers in total
        turnover <- length(int_both)/length(uni_both)
        module_layer_turnover_shuf <- rbind(module_layer_turnover_shuf, tibble(layer_from= i, layer_to= j, turnover= turnover, trial = trial))
      }
    }
    
    edge_list_by_layers_ave_shuf <- edge_list_with_distances_shuf %>% group_by(layer_from, layer_to) %>%
      summarise(ave_distance= mean(distance_in_meters)) %>% unique()
    
    layers_turnover_with_distnace <- edge_list_by_layers_ave_shuf %>%
      merge(module_layer_turnover_shuf, by= c("layer_from", "layer_to")) #merge both versions
    
  }
  return(layers_turnover_with_distnace)
}


##----partner comparison shuffled null models


# this functions calculates the similarity in module partner composition
# module partner composition = which species are found with which other species in the same module
partner_comparison_func <- function(edge_list, partner_comparison_shuf){
  partner_comparison_shuf <- NULL
  for(i in 1:1000){
    print(i) #to keep tab on how far along we are
    partners_null_model_shuf <- NULL
    partners_trial_edge_list <- edge_list %>% filter(trial_num == i) #one trial at a time
    
    max_node <- partners_trial_edge_list %>% slice(which.max(node_id)) %>% select(node_id) #get the max value in the data frame
    max_node <- max_node$node_id #take the value of the max value
    
    module_partners_shuf_modules <- partners_trial_edge_list %>% group_by(trial_num, node_id) %>% #look at chunks of trial number and nodes
      select(module) %>% unique() #select only the unique module for each one. df has trial_num, node_id, module
    
    species_in_module <- module_partners_shuf_modules %>% group_by(trial_num, module) %>% #look from module perspective
      mutate(partners = list(node_id)) #create new column called partners where the value is a list of all the species in the module
    
    for (j in 1:max_node){
      filter_by_node_shuf <- filter(partners_trial_edge_list, partners_trial_edge_list$node_id == j) #only rows where node is
      modules_with_species_shuf <- unique(filter_by_node_shuf$module) #only take unique module numbers
      module_partners_shuf <- filter(partners_trial_edge_list, partners_trial_edge_list$module == modules_with_species_shuf)
      just_partners_shuf <- subset(module_partners_shuf, node_id != j) #remove the current id from the module
      unique_partners_shuf <- unique(just_partners_shuf$node_id) #unique partners in a specific module 
      partners_null_model_shuf <- rbind(partners_null_model_shuf, tibble(node_id = j, partners = unique_partners_shuf)) #all nodes and all their partners
      
      partners_null_model_shuf <- distinct(partners_null_model_shuf)
    }    
    partners_null_model_shuf <- partners_null_model_shuf %>%
      merge(physical_nodes, by="node_id") %>%
      subset(select= -species)
    
    jaccard_val_shuf <- NULL #restart every run
    
    for (k in 1:max(partners_empirical$node_id)){ #for all the nodes
      filter_empirical <- filter(partners_empirical, partners_empirical$node_id == k) #filter only by node number k in empirical
      filter_null_model_shuf <- filter(partners_null_model_shuf, partners_null_model_shuf$node_id == k) #filter only by node number k in null model
      jaccard_similarity <- jaccard(filter_empirical$partners, filter_null_model_shuf$partners) #jaccard on partners of node k in current val
      jaccard_val_shuf <- rbind(jaccard_val_shuf, 
                                tibble(node = k, similarity = jaccard_similarity, type = filter_empirical$type)) #contains node and similarity
    }
    partner_comparison_shuf <- rbind(partner_comparison_shuf, tibble(trial_number=i, similarity=jaccard_val_shuf))
    partner_comparison_shuf <- distinct(partner_comparison_shuf)
  }
  return(partner_comparison_shuf)
}



##---- map creation-------------------------------------------------------------------------------- 

##---- pivot by country

# this function takes the data and creates a new version of layer_id and num of state nodes in it
pivot_by_country <- function(data) {  
  s1 = melt(data, id = c("layer_id", "group"), measure.vars = "size_of_module")
  s2 = dcast(s1, layer_id ~ group, length)
  s2$Total = rowSums(s2[,2:NCOL(s2)])
  return(s2)
}


##---- hierarchical modularity---------------------------------------------------------------------
## hierarchical empirical network
## hierarchical shuf pols plants both
## hierarchical fixed vs empirical


##---- hierarchical infomap

#this function allows for hierarchical modularity
multi_lvl_infomap <- function (M, infomap_executable = "Infomap", flow_model = NULL, 
                               silent = T, trials = 100, seed = NULL, relax = F, multilayer_relax_rate = 0.1, 
                               multilayer_relax_limit = NULL, multilayer_relax_limit_up = NULL, 
                               multilayer_relax_limit_down = NULL, temporal_network = F, 
                               run_standalone = T, remove_auxilary_files = F, ...) 
{
  if (check_infomap(infomap_executable) == F) {
    stop("Error in Infomap stand-alone file.")
  }
  if (class(M) != "multilayer") {
    stop("M must be of class multilayer")
  }
  arguments <- paste("--tree -N ", trials, sep = "")
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
                                         delim = " ", skip = 8, col_names = c("path", "flow", 
                                                                              "name", "state_id", "node_id", "layer_id")))
  modules$path <- sub(":[0-9]+$", "", modules$path) #del last portion of path
  
  modules %<>% filter(flow > 0) %>% select(path, node_id, layer_id, 
                                           flow) %>% separate(path, into = c("module", "sub_module", "sub_sub_module"), 
                                                              sep = ":") %>% mutate_at(.vars = 1:5, as.integer) %>% 
    full_join(M$nodes, "node_id") %>% select(node_id, starts_with("module"), 
                                             everything()) %>% arrange(node_id, layer_id)
  
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
      select(-module) %>% rename(module = module_renamed)
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


##-- pivot by sub module

# this function creates a data frame with sub-module on the side and layer_id on the top
pivot_by_sub_module <- function(data){ 
  s1 = melt(data, id = c("layer_id", "module_sub_module"))
  s2 = dcast(s1, layer_id ~ module_sub_module, length)
  s3 = t(s2)
  colnames(s3) <- s3[1,] #change names of columns to names of layers
  s3 <- s3[-1,]
  return(s3)
}


##---- edge list per sub module

# this function creates an edge list of all the layer found in a module
edge_list_per_sub_module <- function(data,edge_list){

  for (i in (1:13)){
    if (data[i] == 0) next #only take layers where the module is present
    else {
      for (j in (i+1):14){
        if (data[j] == 0) next #only take layers where the module is present
        else {
          edge_list <- rbind(edge_list, tibble(layer_from=i, layer_to=j, module_sub_module=NA)) 
        }
      }
    }
  }
  return(edge_list)
}

##---- hierarchical modularity for shuffled null models

#this function takes edge lists of shuffled matrices and creates a data frame of hierarchical modularity 
modularity_for_shuf_multi_lvl <- function(edge_list, output){
  for(trial_num in 1:1000){
    print(trial_num) #to keep tab on how far along we are
    current_trial_edgelist <- edge_list %>% filter(trial_number == trial_num) #take 1 trial at a time to create multilayer
    dryad_multilayer_shuf_trial <- create_multilayer_object(extended = current_trial_edgelist, #edge list changes every run
                                                            nodes = physical_nodes,
                                                            layers = layer_metadata,
                                                            intra_output_extended = T)
    
    
    #create modules for empirical network
    modules_dryad_multilayer_shuf_1000 <- multi_lvl_infomap(dryad_multilayer_shuf_trial, 
                                                           infomap_executable = "Infomap",
                                                           flow_model = 'directed',
                                                           relax = F, 
                                                           silent = T, 
                                                           trials = 100,
                                                           seed = 497294, 
                                                           temporal_network = F)
    
    
    modules_for_shuf <- modules_dryad_multilayer_shuf_1000$modules
    
    modules_for_shuf$module_sub_module <- paste(modules_for_shuf$module, ".", modules_for_shuf$sub_module) #create new column combining module and sub module
    modules_for_shuf_multi_lvl_sub_module <- modules_for_shuf %>% drop_na(sub_module) #remove every row where there's no sub module
    
    output <- rbind(output, tibble(modules_for_shuf_multi_lvl_sub_module, trial_num)) #save modules for 1000 null models
    
  }  
  return(output)
}


##----hierarchical sub-module distance decay for shuffled and empirical

# this function calculates the Jaccard Similarity in sub-modules between islands
sub_module_distance_decay_func_shuffle <- function(multilayer_1000, 
                                                   islands_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial) #to keep tab on how far along we are

    modules_for_similarity_shuf <- multilayer_1000 %>% filter(trial_num == trial) #take only 1 trial
    
    #pivot modules
    module_pivoted_shuf <- pivot_by_sub_module(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    #create edge list with distances
    modules_edge_list_shuf <- NULL
    
    for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_sub_module(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_sub_module <- rownames(module_pivoted_shuf)[k]
      if (is.null(modules_edge_list_shuf)) next
      modules_edge_list_shuf <- modules_edge_list_shuf %>% mutate(module_sub_module = replace_na(module_sub_module, current_sub_module)) #add module number
    }
    
    edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
    edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA
    
    #modules similarity pairwise distance between islands
    #change layer number to island identity (aggregate by island)
    edge_list_by_islands_shuf <- edge_list_with_distances_shuf
    old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
    new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
    edge_list_by_islands_shuf$layer_from[edge_list_by_islands_shuf$layer_from %in% old] <- new[match(edge_list_by_islands_shuf$layer_from, old)]
    edge_list_by_islands_shuf$layer_to[edge_list_by_islands_shuf$layer_to %in% old] <- new[match(edge_list_by_islands_shuf$layer_to, old)]
    
    #version with # of modules in layers
    edge_list_by_islands_modules_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to, module_sub_module) %>%
      summarise(ave_distance= mean(distance_in_meters)) 
    edge_list_by_islands_modules_shuf$count <- c(1)
    edge_list_by_islands_modules_shuf <- edge_list_by_islands_modules_shuf %>% mutate(number_of_sub_modules= sum(count)) %>%
      select(layer_from, layer_to, module_sub_module, number_of_sub_modules) 
    
    ## jaccard similairty shuf
    #total number of modules in each layer
    module_island_turnover_shuf_multi_lvl <- NULL
    island_list <- c("1","2","3","4east","4west","5","6")
    
    for (i in island_list){
      for (j in island_list){
        sub_modules_in_island_from_shuf <- filter(edge_list_by_islands_modules_shuf, layer_from == i) %>% 
          select(module_sub_module) %>% unique() %>% unlist() #sub-modules in island from
        sub_modules_in_island_to_shuf <- filter(edge_list_by_islands_modules_shuf, layer_from == j) %>% 
          select(module_sub_module) %>% unique() %>% unlist() #sub-modules in island to
        int_both <- intersect(sub_modules_in_island_from_shuf, sub_modules_in_island_to_shuf) #how many nodes are common in both islands
        uni_both <- union(sub_modules_in_island_from_shuf, sub_modules_in_island_to_shuf) #how many nodes are found in both islands in total
        turnover <- length(int_both)/length(uni_both)
        module_island_turnover_shuf_multi_lvl <- rbind(module_island_turnover_shuf_multi_lvl, tibble(layer_from= i, layer_to= j, turnover= turnover))
      }
    }
    
    edge_list_by_islands_ave_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to) %>%
      summarise(ave_distance= mean(distance_in_meters)) %>% unique()
    
    islands_turnover_with_distnace <- edge_list_by_islands_ave_shuf %>%
      merge(module_island_turnover_shuf_multi_lvl, by= c("layer_from", "layer_to")) #merge both versions
    
    islands_turnover_with_distnace_all_trials <- rbind(islands_turnover_with_distnace_all_trials, 
                                                       tibble(islands_turnover_with_distnace, trial_num = trial))
    
  }
  return(islands_turnover_with_distnace_all_trials)
}


##---- pivot sub modules

# this function takes the data and creates a new version of layer_id and num of state nodes in it
pivot_by_island_sub_module <- function(data) { 
  s1 = melt(data, id = c("layer_id", "module_sub_module"), measure.vars = "size_of_sub_module")
  s2 = dcast(s1, layer_id ~ module_sub_module, length)
  print(s2)
  s2$Total = rowSums(s2[,2:NCOL(s2)])
  return(s2)
}

##---- delete all species found in all 14 layers----------------------------------------------------

##---- create modules for species found in less than 14 layers

modularity_for_shuf_del <- function(edge_list, output){
  for(trial_num in 1:1000){
    #current_trial_edgelist <- dryad_edgelist_complete_shuf %>% filter(trial_number == i) #take 1 trial at a time to create multilayer
    current_trial_edgelist <- edge_list %>% filter(trial_number == trial_num) #take 1 trial at a time to create multilayer
    dryad_multilayer_shuf_trial <- create_multilayer_object(extended = current_trial_edgelist, #taking edge list and returning multilayer network
                                                            nodes = physical_nodes_del_shuf,
                                                            layers = layer_metadata,
                                                            intra_output_extended = T)
    
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
  }  
  return(output)
}

##---- distance decay in modules 

module_distance_decay_func_del <- function(multilayer_1000, 
                                           islands_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial)
    
    modules_for_similarity_shuf <- multilayer_1000 %>% filter(trial_num == trial) #%>% #take only 1 trial
    
    #pivot modules
    module_pivoted_shuf <- pivot_by_module(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    #create edge list with distances
    #similarity between state nodes for modules found in 2 or more layers
    modules_edge_list_shuf <- NULL
    
    for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[k]
      if (is.null(modules_edge_list_shuf)) next
      modules_edge_list_shuf <- modules_edge_list_shuf %>% mutate(module = replace_na(module, current_module)) #add module number
    }
    
    #view(modules_edge_list)
    
    
    edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
    edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA and delete layer name
    
    #modules similarity pairwise distance between islands
    edge_list_by_islands_shuf <- edge_list_with_distances_shuf
    old <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
    new <- c("1","1","2","2","3","3","4east","4east","4west","4west","5","5","6","6")
    edge_list_by_islands_shuf$layer_from[edge_list_by_islands_shuf$layer_from %in% old] <- new[match(edge_list_by_islands_shuf$layer_from, old)]
    edge_list_by_islands_shuf$layer_to[edge_list_by_islands_shuf$layer_to %in% old] <- new[match(edge_list_by_islands_shuf$layer_to, old)]
    
    #version with # of modules in layers
    edge_list_by_islands_modules_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to, module) %>%
      summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
    edge_list_by_islands_modules_shuf$count <- c(1)
    edge_list_by_islands_modules_shuf <- edge_list_by_islands_modules_shuf %>% mutate(number_of_modules= sum(count)) %>%
      select(layer_from, layer_to, module, number_of_modules) 
    
    
    #version with correct average between layers
    edge_list_by_islands_ave_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to) %>%
      summarise(ave_distance= mean(distance_in_meters)) %>% unique()
    
    #combine
    edge_list_island_combine_shuf <- edge_list_by_islands_ave_shuf %>%
      merge(edge_list_by_islands_modules_shuf, by= c("layer_from", "layer_to")) #merge both versions 
    
    edge_list_island_combine_no_module_shuf <- edge_list_island_combine_shuf %>% select(-module) %>% unique() #have version where modules aren't present
    

    #total number of modules in each layer
    island_list <- c("1","2","3","4east","4west","5","6")
    
    for (i in island_list){
      for (j in island_list){
        modules_in_island_from_shuf <- filter(edge_list_by_islands_modules_shuf, layer_from == i) %>% select(module) %>% unique() %>% unlist()
        modules_in_island_to_shuf <- filter(edge_list_by_islands_modules_shuf, layer_from == j) %>% select(module) %>% unique() %>% unlist()
        #take all nodes in layer_from and all nodes in layer_to to check turnover
        int_both <- intersect(modules_in_island_from_shuf, modules_in_island_to_shuf) #how many nodes are found in both layers
        uni_both <- union(modules_in_island_from_shuf, modules_in_island_to_shuf)
        turnover <- length(int_both)/length(uni_both)
        module_island_turnover_shuf <- rbind(module_island_turnover_shuf, tibble(layer_from= i, layer_to= j, turnover= turnover, trial = trial))
      }
    }
    
    edge_list_by_islands_ave_shuf <- edge_list_by_islands_shuf %>% group_by(layer_from, layer_to) %>%
      summarise(ave_distance= mean(distance_in_meters)) %>% unique()
    
    islands_turnover_with_distnace <- edge_list_by_islands_ave_shuf %>%
      merge(module_island_turnover_shuf, by= c("layer_from", "layer_to")) #merge both versions
    
  }
  return(islands_turnover_with_distnace)
}

##---- distance decay in modules layers
module_distance_decay_layer_func_del <- function(multilayer_1000, 
                                                 layers_turnover_with_distnace){
  for (trial in 1:1000){
    print(trial)
    
    modules_for_similarity_shuf <- multilayer_1000 %>% filter(trial_num == trial) #%>% #take only 1 trial

    #pivot modules
    module_pivoted_shuf <- pivot_by_module(modules_for_similarity_shuf) #pivot will be done on 1 trial each time
    
    
    #create edge list with distances
    #similarity between state nodes for modules found in 2 or more layers
    modules_edge_list_shuf <- NULL
    
    for (k in (1:nrow(module_pivoted_shuf))){ #run the function for each row in the data frame
      modules_edge_list_shuf <- edge_list_per_module(module_pivoted_shuf[k,], modules_edge_list_shuf) 
      current_module <- rownames(module_pivoted_shuf)[k]
      if (is.null(modules_edge_list_shuf)) next
      modules_edge_list_shuf <- modules_edge_list_shuf %>% mutate(module = replace_na(module, current_module)) #add module number
    }
    
    #view(modules_edge_list)
    
    
    edge_list_with_distances_shuf <- right_join(modules_edge_list_shuf, distances_with_ids, by= c("layer_from", "layer_to")) #combine the edge list with the distances between each two layers
    edge_list_with_distances_shuf <- na.omit(edge_list_with_distances_shuf) #remove NA and delete layer name

    #version with # of modules in layers
    edge_list_by_layers_modules_shuf <- edge_list_with_distances_shuf %>% group_by(layer_from, layer_to, module) %>%
      summarise(ave_distance= mean(distance_in_meters)) #maybe do it differently? should i make all distances within the same island 0?
    edge_list_by_layers_modules_shuf$count <- c(1)
    edge_list_by_layers_modules_shuf <- edge_list_by_layers_modules_shuf %>% mutate(number_of_modules= sum(count)) %>%
      select(layer_from, layer_to, module, number_of_modules) 
    
    
    #version with correct average between layers
    edge_list_by_layers_ave_shuf <- edge_list_with_distances_shuf %>% group_by(layer_from, layer_to) %>%
      summarise(ave_distance= mean(distance_in_meters)) %>% unique()
    
    #combine
    edge_list_layer_combine_shuf <- edge_list_by_layers_ave_shuf %>%
      merge(edge_list_by_layers_modules_shuf, by= c("layer_from", "layer_to")) #merge both versions 
    
    edge_list_layer_combine_no_module_shuf <- edge_list_layer_combine_shuf %>% select(-module) %>% unique() #have version where modules aren't present
    
    #total number of modules in each layer
    
    for (i in 1:14){
      for (j in 1:14){
        modules_in_layer_from_shuf <- filter(edge_list_by_layers_modules_shuf, layer_from == i) %>% select(module) %>% unique() %>% unlist()
        modules_in_layer_to_shuf <- filter(edge_list_by_layers_modules_shuf, layer_from == j) %>% select(module) %>% unique() %>% unlist()
        #take all nodes in layer_from and all nodes in layer_to to check turnover
        int_both <- intersect(modules_in_layer_from_shuf, modules_in_layer_to_shuf) #how many nodes are found in both layers
        uni_both <- union(modules_in_layer_from_shuf, modules_in_layer_to_shuf)
        turnover <- length(int_both)/length(uni_both)
        module_layer_turnover_shuf_del <- rbind(module_layer_turnover_shuf_del, tibble(layer_from= i, layer_to= j, turnover= turnover, trial = trial))
      }
    }
    
    edge_list_by_layers_ave_shuf <- edge_list_with_distances_shuf %>% group_by(layer_from, layer_to) %>%
      summarise(ave_distance= mean(distance_in_meters)) %>% unique()
    
    layers_turnover_with_distnace <- edge_list_by_layers_ave_shuf %>%
      merge(module_layer_turnover_shuf_del, by= c("layer_from", "layer_to")) #merge both versions
    
  }
  return(layers_turnover_with_distnace)
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
    full_join(M$nodes, "node_id") %>% dplyr::select(node_id, starts_with("module"), 
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

#---- proper theme for figures ------------------------------------------------------------------------------------------------
#not an actual function but must be applied to figures for paper
paper_figs_theme <- 
  theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black",fill = NA,size = 1),
        panel.spacing = unit(0.5, "cm", data = NULL),
        axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size=14, color='black'),
        axis.line = element_blank())
paper_figs_theme_no_legend <- 
  paper_figs_theme +
  theme(legend.position = 'none')

