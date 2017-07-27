# GeoNet
### R code for constructing metacommunity networks using geographic coordinates of the sites (local communities within a landscape) and presence-absence of multiple taxa within the sites.

## Usage:
1) Load main code file:
```
      source("path/Main_GeoNet.r")
```
2) Read presence-absence data:
```
      incid <- read_incid(filename, header=T, sep="\t", row_names=T)
```
3) Read coordinates:
```
      coord <- read_coord(filename, header=T, sep="\t", row_names=T)
```
4) Or alternatively, simulate distance matrix:
```
      dist_matrix <- sim_dist(num_sites)
```
5) Build the network:
```
      network <- build_network(incidence, coord, num_perms, pl)
```
6) Plot the network:
```
      plot_network(network, method=NULL)
```

**For more detail see:**

Layeghifard, M, Makarenkov, V, Peres-Neto, PR. 2015. *Spatial and species compositional networks for inferring connectivity patterns in ecological communities*. Global Ecology & Biogeography. 24:6, 718-727. [DOI: 10.1111/geb.12293](http://onlinelibrary.wiley.com/doi/10.1111/geb.12293/full)
