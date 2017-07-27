##############################################################################
#
# R functions for construction of metacommunity networks by taking advantage
# of geographic coordinates of the sites (local communities within a landscape)
# and presence-absence of multiple taxa within the sites.
#
#
# Usage:
# 1) Load this file:
#       source("path/Main_GeoNet.r")
#
# 2) Read presence-absence data:
#       incid <- read_incid(filename, header=T, sep="\t", row_names=T)
# 
# 3) Read coordinates:
#       coord <- read_coord(filename, header=T, sep="\t", row_names=T)
#
# 4) Or alternatively, simulate distance matrix:
#       dist_matrix <- sim_dist(num_sites)
# 
# 5) Build the network:
#       network <- build_network(incidence, coord, num_perms, pl)
# 
# 6) Plot the network:
#       plot_network(network, method=NULL)
# 
# See for more detail:
# Layeghifard, M, Makarenkov, V, Peres-Neto, PR. 2015. Spatial and species
# compositional networks for inferring connectivity patterns in ecological
# communities. Global Ecology & Biogeography. 24:6, 718-727.
##############################################################################

# load required libraries
require(vegan)
# vegan is required for Jaccard dissimilarity function
require(ape)
# ape is required for phylo tree reconstruction
require(igraph)
# igraph is required for converting phylo tree to graph


##############################################################################
#
# R function to read incidence (presence-absence of species within sites)
# into a matrix.
#
# Preferred file format:
# (otherwise, specify details: header, separator & row names)
#
#       Species1 Species2 Species3 Species4 Species5
# Site1        0        1        0        0        0
# Site2        0        0        1        1        1
# Site3        0        0        0        0        1
# Site4        1        0        0        0        0
# Site5        0        0        0        0        1
#
# Usage:
# incid <- read_incid(filename, header=T, sep=",", row_names=T)

read_incid <- function(infile, header=T, sep=",", row_names=T) {
    if (header & row_names) {
        x <- read.table(infile, header = T, sep = sep, row.names = 1)
    } else if (header) {
        x <- read.table(infile, header = T, sep = sep)
    } else {
        x <- read.table(infile, header = F, sep = "\t")
    }
    incid <- data.matrix(x)
}



##############################################################################
#
# R function to read geographic coordinates (longitude & latitudes) of sites
# into a matrix.
#
# Preferred file format:
# (otherwise, specify details: header, separator & row names)
#
#       Longitude Latitude
# Site1    49.908  -81.993
# Site2    49.452  -82.050
# Site3    49.779  -81.924
# Site4    49.953  -81.872
# Site5    49.778  -82.323
#
# Usage:
# coord <- read_coord(filename, header=T, sep=",", row_names=T)

read_coord <- function(infile, header=T, sep=",", row_names=T) {
    if (header & row_names) {
        x <- read.table(infile, header = T, sep = sep, row.names = 1)
    } else if (header) {
        x <- read.table(infile, header = T, sep = sep)
    } else {
        x <- read.table(infile, header = F, sep = "\t")
    }
    coord <- data.matrix(x)
}



##############################################################################
#
# R function to compute distance (as on the surface of earth) between two
# points (sites) from their coordinates (Longitude, Latitude). To be called
# by calc_dist function below.
#
# Usage:
# dist <- geodist(point1, point2)

geodist <- function(point1, point2) {
    radius <- 6371
    rad1 <- point1 * pi/180
    rad2 <- point2 * pi/180
    diam <- sin(rad1[2]) * sin(rad2[2]) + cos(rad1[2]) * cos(rad2[2]) *
        cos(abs(rad1[1]-rad2[1]))	
    dist <- radius * acos(diam)
    return(dist)
}



##############################################################################
#
# R function to compute distance between sites within a metacommunity using
# geodist function above.
#
# Usage:
# dist_matrix <- calc_dist(coordinates)

calc_dist <- function(coord){
    nSites <- dim(coord)[1]
    distMat <- array(0, dim=c(nSites, nSites))
    for (i in 1:(nSites-1)) {
        for (j in (i+1):nSites) {
            distMat[i,j] = geodist(coord[i,],coord[j,])
            distMat[j,i] = distMat[i,j]
        }
    }
    # Normalizing the distance matrix - only for real coordinates
    distMat <- distMat / max(distMat)
    
    # Add names to rows & cols
    row_names <- rownames(coord)
    rownames(distMat) <- row_names
    colnames(distMat) <- row_names
    
    return(distMat)
}



##############################################################################
#
# R function to simulate distance between sites within a metacommunity, where
# output is a simple distance matrix.
#
# Usage:
# dist_matrix <- sim_dist(num_sites)

sim_dist <- function(num_sites){
    sqr_sites <- sqrt(num_sites)
    coord <- array(0, dim=c(num_sites, 2))
    
    count <- 1
    for (i in 1:sqr_sites) {
        for (j in 1:sqr_sites) {
            coord[count, 1] <- i
            coord[count, 2] <- j
            count <- count + 1
        }
    }
    # Distance matrix
    distMat <- vegdist(coord, "euclidean", diag=T, upper=T)
    
    return(as.matrix(distMat))
}



##############################################################################
#
# R function to detect extra links between sites within a metacommunity using
# incidence data a Similarity/Permutation test to add to links found using
# geographic distance.
#
# Usage:
# extra_branches <- find_extra_brs(incidence, num_perms)

find_extra_brs <- function(incid, num_perms) {
    
    taxa <- dim(incid)[1]
    method = "jaccard"
    
    # Empty 3d matrix to hold distance/simmilarity matrices
    allSimMats = array(0, dim=c(taxa,taxa,num_perms+1))
    
    perms = permatfull(incid, fixedmar="row", shuffle="samp", strata=NULL, mtype="count", times=num_perms)
    # Original + permuted matrices
    
    orig_sim <- vegdist(incid, method=method, diag=T, upper=T)
    allSimMats[,,1] = 1 - (as.matrix(orig_sim))
    
    for(p in 1:num_perms) {
        sim <- vegdist(perms$perm[[p]], method=method, diag=T, upper=T)
        allSimMats[,,p+1] <- 1 - (as.matrix(sim))
    }
    
    num.rows = nrow(allSimMats)
    num.cols = ncol(allSimMats)
    num.gens = dim(allSimMats)[3]
    extra_brs = array(0, dim=c(num.rows, num.cols))
    
    for (x in 1:(num.rows-1)) {
        for (y in (x+1):num.rows) {
            my.val = allSimMats[x,y,1]
            all.vals = allSimMats[x,y,2:(num_perms+1)]
            extra_brs[x, y] = (sum(my.val > all.vals) / length(all.vals))
            extra_brs[y, x] = extra_brs[x, y]
        }
    }
    
    extra_brs[extra_brs>0.95] = 1
    extra_brs[extra_brs<=0.95] = 0
    
    # Add names to rows & cols
    row_names <- rownames(incid)
    rownames(extra_brs) <- row_names
    colnames(extra_brs) <- row_names
    
    return(extra_brs)
}



##############################################################################
#
# R function to build the metacommunity network using neighbor-joining tree
# from geographic distance and extra branches detected using incidence data.
#
# Usage:
# network <- build_network(incidence, coord, num_perms, pl)

build_network <- function(incid, coord, num_perms=100, pl=F) {
    # Building distnce matrix from coordinates
    distMat <- calc_dist(coord)
    
    # Building NJ tree -- ape package
    PhyloTree <- nj(distMat)
    
    # Build graph from tree -- igraph package
    PhyloNet <- as.igraph(PhyloTree)
    
    # Find extra branches
    extra_brs <- find_extra_brs(incid, num_perms)
    
    rnames <- rownames(incid)
    num.rows <- dim(incid)[1]
    # Add extra branches to network
    for (x in 1:(num.rows-1)) {
        for (y in (x+1):num.rows) {
            if (extra_brs[x, y] == 1) {
                PhyloNet <- add.edges(PhyloNet, c(as.character(rnames[x]), as.character(rnames[y])))
                #cat(sprintf("%s : %s\n", x, y))
            }
        }
    }
    
    return(PhyloNet)
}



##############################################################################
#
# R function to build the metacommunity network using neighbor-joining tree
# from geographic distance and extra branches detected using incidence data.
#
# Usage:
# plot_network(network, method=NULL)
#
# Available methods: Infomap, Walktrap, Edge betweenness, Fast greedy,
# Leading eigenvector & Label propagation. Default is NULL.

plot_network <- function(network, method=NULL) {
    comms <- NULL
    if (is.null(method)) {
        plot(network)
    } else if (method == "infomap") {
        comms<- infomap.community(network)
    } else if (method == "walktrap") {
        comms<- walktrap.community(network)
    } else if (method == "edge_betweenness") {
        comms<- edge.betweenness.community(network)
    } else if (method == "fastgreedy") {
        comms<- fastgreedy.community(network)
    } else if (method == "leading_eigenvector") {
        comms<- leading.eigenvector.community(network)
    } else if (method == "label_propagation") {
        comms<- label.propagation.community(network)
    }
    # Ploting network & its communities
    if (!is.null(comms)) {
        plot(comms, network)
    }
}
