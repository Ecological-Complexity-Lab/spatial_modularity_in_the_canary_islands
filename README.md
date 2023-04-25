An explanation on our methods and analysis in a logical story form can be found in the ecomplab drive -> 3rd year projects -> 
Maya -> Spatial modularity in the Canary Islands

# Question

Can we find a signature of distance decay in spatial modules?

# Network building

## Layers as sites
Nodes: species
Intra: The poroportion of visits within a layer of the species, out of all visits of the species in the layer
Inter: Normalized distance

## Layers as islands
Nodes: aggregate within an island all the species.
Intra: First, sum the number of visits in both sites within an island per species pair; then define as before.
Inter: Distance between sites is disregarded. Between two islands: the mean of th distance between the four sites in two islands. Because the distance between sites is negligible in scale comparted to between islands, this is equivalent to a distance between islands.

# Results

We test this querstion when the network layers are islands and when layers are sites.

## Distance decay in species across layers

### Layers as sites

Negative relationship in empirical data. We test this empirical pattern against null model $M_1$:
1. Shuffle plants ($M_1^P$)
2. Shuffle pollinators ($M_1^A$)
3. Shuffle both ($M_1^{AP}$)

**[FIGURE]**

### Layers as islands

Negative relationship in empirical data. We test this empirical pattern against null model $M_1$:
1. Shuffle plants ($M_1^P$)
2. Shuffle pollinators ($M_1^A$)
3. Shuffle both ($M_1^{AP}$)

**[FIGURE]**

## Distance decay in interactions

### Layers as sites
1. Distance decay in interactions **[FIGURE]**
2. A matrix form of the Jaccard used in 1. **[FIGURE]**

### Layers as islands
1. Distance decay in interactions **[FIGURE]**
2. A matrix form of the Jaccard used in 1. **[FIGURE]**

## Modularity analysis

No heirarchical modularity

### Layers as sites

Description of the modules:
- Number of modules- 43
- Size distribution of modules:

![module_sizes_sites](https://github.com/Ecological-Complexity-Lab/spatial_modularity_in_the_canary_islands/blob/main/graphs/modularity_analysis/module_size_distribution.pdf)

- Distribution of number of layers in a module:

![layers_distribution_sites](https://github.com/Ecological-Complexity-Lab/spatial_modularity_in_the_canary_islands/blob/main/graphs/modularity_analysis/number_of_layers_in_module.pdf)

We find distance decay in modules in the empirical network. Compare to null models
- $M_1$
- $M_2$
- $M_3$
- $M_4$

**[4 FIGURES: one per null model]**

Statistical tests by comparing the $R^2$

**[FIGURE]**

### Layers as islands

Description of the modules:
- Number of modules
- Size distribution of modules **[FIGURE]**
- Distribution of number of layers in a module **[FIGURE]**

We find distance decay in modules in the empirical network. Compare to null models
- $M_1$
- $M_2$
- $M_3$
- $M_4$

**[4 FIGURES: one per null model]**

Statistical tests by comparing the $R^2$
