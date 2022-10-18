# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import pandas as pd
import math
import os
from geopy.distance import lonlat, distance  #geo library

#os.chdir("C:/Users/golds/Desktop/uni/shai_lab/dryad_data")

areas=["WesternSahara1", "WesternSahara2", "Fuerteventura1", "Fuerteventura2", "GranCanaria1", "GranCanaria2", "TenerifeSouth1", "TenerifeSouth2", "TenerifeTeno1", "TenerifeTeno2", "Gomera1", "Gomera2", "Hierro1", "Hierro2"]  #all areas that were sampled in the research
pollination_network=pd.read_excel("all_sites.xlsx", index_col=0, sheet_name=areas)

organism_intralayer_edges=[]  #edge list of intalayer edges
for area_name in areas:
    area=pollination_network[area_name]  #acessing the tab in xls
    pollination_data=pd.DataFrame(columns=["layer_from","node_from","layer_to","node_to","weight"])  #df for intralayer. layer from and to should be the same layer
    pollinators_list=list(area) #nodes
    for i in range(0,len(area)):
        pollination=area.iloc[i]  #node_from
        for x in range (0,len(pollinators_list)):
            if pollination[x]!=0:  #node_to if interaction is present
                weight=pollination[x]
                pollination_data=pollination_data.append(pd.DataFrame([[area_name,pollinators_list[x],area_name,pollination.name,weight]], columns=["layer_from","node_from","layer_to","node_to","weight"]), ignore_index=True)
    organism_intralayer_edges.append(pollination_data) #creating interlayer for each layer on its own
organism_intralayer_edges=pd.concat(organism_intralayer_edges) #adding up all interlayers to one edge list

organism_intralayer_edges.to_csv("intralayer_file.csv" , index=False)

print(organism_intralayer_edges)

interlayer_dict={}  #dict where species is key and spatial networks are values
for i in range (0,len(organism_intralayer_edges)):
    edge=organism_intralayer_edges.iloc[i] #every from to interaction in evert spatial network
    if edge[1] not in interlayer_dict.keys(): #if species (from) hasn't been added yet
        interlayer_dict.update({edge[1]:[edge[0]]}) #add the species and its network 
    else:
        if edge[0] not in interlayer_dict[edge[1]]: #if species (from) already exists but network hasn't been added yet
            interlayer_dict[edge[1]].append(edge[0]) #add network to values
        else:
            continue
    if edge[3] not in interlayer_dict.keys(): #if species (to) hasn't been added yet
        interlayer_dict.update({edge[3]:[edge[0]]}) #add the species and its network 
    else:
        if edge[0] not in interlayer_dict[edge[3]]: #if species (to) already exists but network hasn't been added yet
            interlayer_dict[edge[3]].append(edge[0]) #add network to values
        else:
            continue

#print(interlayer_dict)
        
#geo_data=pd.read_excel("gilarranz2014_datos_sierras.xlsx", index_col=0, sheet_name="Coordinates")
#geo_distances=pd.DataFrame(columns=["layer_from","layer_to","distance_in_km"]) #df of all distances between every 2 layers
#for i in range (0,len(geo_data)):
#    for x in range(i+1,len(geo_data)):
#        coords1=(geo_data.iloc[i][3],geo_data.iloc[i][4]) #coordinates of layer from
#        coords2=(geo_data.iloc[x][3],geo_data.iloc[x][4]) #coordinates of layer to
#        geo_distances=geo_distances.append(pd.DataFrame([[geo_data.index[i],geo_data.index[x],distance(coords1,coords2).km]], columns=["layer_from","layer_to","distance_in_km"]), ignore_index=True)
#print(geo_distances)
#geo_distances.to_csv("distances_file.csv", index=False)

#shortest_distance=geo_distances.min()
#print(shortest_distance[2])

####
#no need to create distance matrix caluse it already exists
#turn matrix from meters to km
distances= pd.read_csv("Distance_between_sites_Dryad.csv", index_col=None) #take csv
distances.rename(columns= {"DISTANCE (m)":"distance_in_meters", "From":"layer_from", "To":"layer_to"}, inplace= True) #change column names
distances= distances.replace("_","", regex= True) #remove _ from all the location names
shortest_distance= distances.min()

#distances.to_csv("distances_file.csv", index=False)

def interlayer_edges (d):
    """recieves d and creates interlayer edge calculation"""
    interlayer_edge=(1/(math.log(d)))/(1/(math.log(shortest_distance[2]))) #interlayer weight grows weaker the bigger the distance is between networks
    return interlayer_edge


organism_interlayer_edges=pd.DataFrame(columns=["layer_from","node_from","layer_to","node_to","weight"]) #df of interlayer edges. node from and to should be the same
for key in interlayer_dict.keys(): #for every species
    organism_locations=interlayer_dict.get(key) #list of all networks where species was found
    for i in range(0,len(organism_locations)):
        for x in range(i+1,len(organism_locations)):
            locations_in_geo=distances[(distances["layer_from"]==organism_locations[i]) & (distances["layer_to"]==organism_locations[x])] 
            #locations_in_geo["distance_in_km"].apply(lambda x: float(x))
            interlayer=interlayer_edges(locations_in_geo["distance_in_meters"]) #run funciton on every two locations where species was found
            organism_interlayer_edges=organism_interlayer_edges.append(pd.DataFrame([[organism_locations[i],key,organism_locations[x],key,interlayer]], columns=["layer_from","node_from","layer_to","node_to","weight"]), ignore_index=True)
            
#print(organism_interlayer_edges)

organism_interlayer_edges.to_csv("interlayer_file.csv" , index=False)

import sys
print(sys.version)