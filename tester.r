source("C:\\Users\\mlaye\\Dropbox\\Codes\\RScripts\\GeoNet-Renato\\Main_GeoNet.r")

file1 = "data/InTest.csv"
file2 = "data/CoTest.csv"

incid <- read_incid(file1)
coord <- read_coord(file2)

PhyloNet <- build_network(incid, coord, num_perms=100, pl=F)

# YOU COULD ALSO USE ALL THE FUNCTIONS WITHIN SOURCE FILE SEPARATELY.

plot(PhyloNet, layout=layout_with_fr, vertex.size=8,
     vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)
