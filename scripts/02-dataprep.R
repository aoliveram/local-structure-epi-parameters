library(igraph)
library(network)
library(netplot)
library(intergraph)
library(ergm)

# Read the Simulated data rds file from data/
networks <- readRDS("data/Simulated_1000_networks.rds")
networks <- unclass(networks)

# Computing statistics using ERGM
S_ergm <- parallel::mclapply(networks, \(n) {
  summary_formula(
    n ~ edges + nodematch("gender") + nodematch("grade")
  )
}, mc.cores = 6)


# Computing statistics based on igraph
inetworks <- parallel::mclapply(networks, \(n) {
  asIgraph(n)  
}, mc.cores = 6)

