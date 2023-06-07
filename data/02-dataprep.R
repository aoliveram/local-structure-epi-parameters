#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --ntasks=31
#SBATCH --mem=100GB
#SBATCH --job-name=abm-net-02-dataprep
#SBATCH --mail-type=all
#SBATCH --mail-user=george.vegayon@utah.edu

library(igraph)
library(network)
library(netplot)
library(intergraph)
library(ergm)
library(data.table)

ncores <- 30

# Read the Simulated data rds file from data/
networks <- readRDS("data/Simulated_1000_networks.rds")
networks <- unclass(networks)

# # Taking a sample of 100 networks
# set.seed(123)
# networks <- networks[sample(1:1000, 100)]

# Computing statistics using ERGM
S_ergm <- parallel::mclapply(networks, \(n) {
  summary_formula(
    n ~ edges + nodematch("gender") + nodematch("grade") + triangles + balance +
      twopath + gwdegree(decay = .25, fixed = TRUE) + isolates
  ) |> as.list() |> as.data.table()
}, mc.cores = ncores) |> rbindlist()

head(S_ergm)

# Computing statistics based on igraph
S_igraph <- parallel::mclapply(networks, \(n) { 
  n <- asIgraph(n)
  
  data.table(
    modularity = modularity(cluster_fast_greedy(n)),
    transitivity = transitivity(n),
    density = igraph::edge_density(n),
    diameter = diameter(n),
    avg_path_length = igraph::mean_distance(n),
    avg_degree = mean(degree(n)),
    avg_betweenness = mean(betweenness(n)),
    avg_closeness = mean(closeness(n)),
    avg_eigenvector = mean(eigen_centrality(n)$vector),
    components = components(n)$no
  )
}, mc.cores = ncores) |> rbindlist()

setnames(S_ergm, new = paste0("ergm_", names(S_ergm)))
setnames(S_igraph, new = paste0("igraph_", names(S_igraph)))

# Combining the datasets
fwrite(
  cbind(S_igraph, S_ergm),
  file = "data/02-dataprep-network-stats.csv.gz"
  )


