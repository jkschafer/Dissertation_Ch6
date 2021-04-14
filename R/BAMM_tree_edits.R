library(ape)
library(phytools)


tree <- read.nexus("./Data/consensusTree_10KTrees_Primates_Version3.nex")
is.ultrametric(tree)
is.binary(tree)
tree <- multi2di(tree)
min(tree$edge.length) # is zero
# make non-zero edge length
n <- length(tree$edge.length)
tree$edge.length <- tree$edge.length[1:n-1] + 0.0000001 
min(tree$edge.length) # is not zero
tree <- chronos(tree, lambda = 0) # make tree ultrametric
write.tree(tree, file = "primatetree.tre") # export for BAMM
