# :wave: About
This repository contains the code and data for the paper: "Local and regional processes drive distance decay in structure in a spatial multilayer plant-pollinator network" - Currently in review

<!--
# :page_facing_up: Paper and citing
Frydman N, Freilikhman S, Talpaz I, Pilosof S. **Practical guidelines and the EMLN R package for handling ecological multilayer networks**. Methods in Ecology and Evolution. 2023. [DOI:10.1111/2041-210X.14225](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.14225).
-->

# Background for project:
Our goal is to test for distance decay in structure between communities using spatial modularity and disentangle the factors that generate this pattern. In multilayer networks, a module can span several local communities. Therefore, distance decay in structure is a pattern in which the farther apart two local communities are, the less they share modules. We use data from local plant-pollinator networks in the Canary Islands, which we represent as a spatial multilayer network. This data set is particularly suitable because a previous study showed distance decay in species and interactions in this system. Distance decay in species composition should strongly affect spatial modularity because if two locations are remote and do not share species, they are unlikely to share modules. We hypothesized that species turnover ($H_1$) and/or interaction rewiring ($H_2$) across locations drive distance decay in structure at the regional scale because both are processes that contribute to shape structure. We further hypothesized that factors occurring within each location (i.e., local factors) could affect structure distance decay ($H_3$) because they influence processes that favor modular structure, such as resource partitioning and coevolution. Finally, we tested whether species' partner similarity across space triggers distance decay in modular structure ($H_4$) because it could increase the chances of sharing modular structures between locations. To test these alternative and not mutually exclusive hypotheses, we developed a set of four null models, which alter different components of the spatial multilayer network. These models allowed us to disentangle local vs regional drivers in spatial modularity.

# :file_folder: Folder breakdown:
* R- All used code.
* csvs_nuevo- Every notable output saved as a csv.
* graphs- Every graph created as a part of this project.


# :computer: R (utility):
## functions
This portion of the code contains functions used in more than one analysis. 
* The functions are organized by their main goal.
* Before every function there is an explanation of the intention of the function or what it calculates.  

# :computer: R (pipeline):
## Empirical
This code turns the data into a multilayer network with 7 layers (locations as layers) and tests distance decay in module composition of the empirical network.
* Gets an interactions csv organized by Noa Fridman based on the interactions recorded and documented in the original research. Also gets a csv of the distance between locations (layers). 
* Creates interlayer and intralayer edges data and normalizes the weights between 0 and 1.
* Interlayers and intralayers edges, nodes, and layers are combined into a multilayer network.
* Modularity analysis is conducted. It also contains basic exploratory analysis regarding modules such as the number of shared modules between locations.
* Tests distance decay in module composition in the empirical network.

***

## Null_M1
This portion of the code runs the null model M<sub>1</sub>, in which we shuffled species between layers in one of 3 versions:
 1. shuffling plants among themselves
 2. shuffling pollinators among themselves
 3. shuffling plants among themselves and then pollinators among themselves

The shuffled networks are then compared to the empirical network to determine whether species turnover is influencing the pattern found.

***

## Null_M2
This portion of the code runs the null model M<sub>2</sub>, in which we shuffled interactions between layers. 

* We randomly shuffled interactions of each pair of species between all the layers in which they co-occur. For example, if a plant and a pollinator co-occurred in layers 1 and 3 but interacted only in layer 1, they would still co-occur in the same layers but may interact in layer 3 after shuffling. 

The shuffled networks are then compared to the empirical network to determine whether interaction rewiring is influencing the pattern found.

***

## Null_M3
This portion of the code runs the null model M<sub>3</sub>, in which we shuffled interactions within layers.

The shuffled networks are then compared to the empirical network to determine whether local factors influence the pattern found.

***

## Null_M4
This portion of the code runs the null model M<sub>4</sub>, in which we set the weight of the already existing interlayer links to a uniform value ranging from 0.1 to 1. This allowed us to test if the extent to which species share partners across locations triggers distance decay in structure.


