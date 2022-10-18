This project focuses on a multilayer network in the Canary Islands made up of pollination interactions dictating
the intralayer structure and signatures of space dictate the interlayer edges
(those may change as the project progresses to encompass more spatial signatures such as ages of the islands, topography of the islands etc,
and not just the geographical distance between the islands).

The network is made up of 288 physical nodes (39 of which are plants and 249 of which are pollinators) 
spanning over 14 spatial layers (sites tested) in 6 islands and a mainland (Western Sahara),
while only two state nodes are found in all 14 sites (Euphorbia balsamifera male and female), 
which were a main criterion for choice of locality. Each island (and the mainland) were sampled in two adjacent sites (~50 to ~500 meters apart) 
to create the 14 layers. There are a total of 969 pollination interactions (weighted based on the abundance of the species in the sites, 
number of visits between each plant-pollinator pair and total number of interactions each state node had in a certain layer. 
For example, if a certain pollinator interacts with two plants in a given layer, 
and has 20 interactions with plant number 1 and 4 interactions with plant number 2. This pollinator has a total of 24 interactions, 
so the weight of the edge connecting into plant number 2 would be 4/24.

There is a total of 1387 signatures of space (interlayer edges) connecting each physical node to itself in a different layer to 
showcase geographical processes. The interlayer edges are also weighted based on a log of the distance between every two layers (base of e), 
normalized using the shortest distance in the network to limit the weights between 0 and 1. This was done to put the interlayer edges 
and intralayer edges on the same scale (0 to 1) as to not prioritize one set of biological and biogeographical processes over the other 
and bias the network that way. 

An explanation on our methods and analysis in a logical story form can be found in the ecomplab drive -> 3rd year projects -> 
Maya -> Spatial modularity in the Canary Islands
