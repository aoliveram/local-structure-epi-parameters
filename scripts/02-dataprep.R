library(igraph)
library(network)
library(netplot)
library(intergraph)
library(ergm)
library(data.table)

# Read the Simulated data rds file from data/
networks <- readRDS("data/Simulated_1000_networks.rds")
networks <- unclass(networks)

# Taking a sample of 100 networks
set.seed(123)
networks <- networks[sample(1:1000, 100)]

# Computing statistics using ERGM
S_ergm <- parallel::mclapply(networks, \(n) {
  summary_formula(
    n ~ edges + nodematch("gender") + nodematch("grade")
  ) |> as.data.table()
}, mc.cores = 6)


# Computing statistics based on igraph
inetworks <- parallel::mclapply(networks, \(n) {
  asIgraph(n)  
}, mc.cores = 6)

S_igraph <- parallel::mclapply(inetworks, \(n) {
  data.table(
    modularity = modulality(cluster_edge_betweenness(n)),
    transitivity = transitivity(n),
    density = density(n),
    diameter = diameter(n),
    avg_path_length = mean(distance(n)),
    avg_degree = mean(degree(n)),
    avg_betweenness = mean(betweenness(n)),
    avg_closeness = mean(closeness(n)),
    avg_eigenvector = mean(eigen_centrality(n)$vector
  )
}, mc.cores = 6)
