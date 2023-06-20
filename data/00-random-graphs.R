#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --job-name=00-rand-graphs
#SBATCH --output=00-rand-graphs-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=george.vegayon@utah.edu
#SBATCH --ntasks=31
#SBATCH --mem=64G
#SBATCH --time=4:00:00

library(epiworldR)
library(slurmR)
library(intergraph)
library(igraph)

# Set the seed for reproducibility
set.seed(1231)
ncores <- 30

# Load the simulated networks from the RDS file
message("Loading the networks")
networks <- readRDS("data/Simulated_1000_networks.rds")
message("Done loading the networks")

# Simulating scalefree networks with same density and size
# as the original networks
message("Simulating scalefree networks")
random_seeds <- sample.int(1e6, length(networks), replace = FALSE)
networks_sf <- parallel::mclapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-sf.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- networks[[i]]


  # Number of vertices in a "network" object
  n_nodes <- network::network.size(n)
  n_edges <- network::network.edgecount(n)

  m <- igraph::sample_pa(
    n_nodes, power = 1,
    m = rgamma(1, n_edges/n_nodes * 10, 10),
    directed = FALSE
  )

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Copying the vertex attributes of n to m
  for (a in network::list.vertex.attributes(n)) {
    network::set.vertex.attribute(m, a, network::get.vertex.attribute(n, a))
  }

  # Writing the network to a file 
  saveRDS(m, fn, compress = FALSE)

  NULL

}, mc.cores = ncores)
message("Done simulating scalefree networks")

# Simulating smallworkd networks with same density and size
# as the original networks
message("Simulating smallworld networks")
networks_sw <- parallel::mclapply(seq_along(networks), \(i) {

  n <- networks[[i]]

  # Number of vertices in a "network" object
  n_nodes <- network::network.size(n)
  n_edges <- network::network.edgecount(n)

  m <- igraph::sample_smallworld(
    n_nodes, dim = 1,
    nei = rgamma(1, n_edges/n_nodes * 10, 10),
    p = 0.2
    )

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Copying the vertex attributes of n to m
  for (a in network::list.vertex.attributes(n)) {
    network::set.vertex.attribute(m, a, network::get.vertex.attribute(n, a))
  }
  
  saveRDS(m, sprintf("data/graphs/%04i-sw.rds", i), compress = FALSE)

  NULL
    
}, mc.cores = ncores)
message("Done simulating smallworld networks")

# Degree-sequence preserve rewiring of the networks
message("Degree-sequence preserve rewiring of the networks")
networks_r <- parallel::mclapply(seq_along(networks), \(i) {

  n <- networks[[i]]

  # Number of vertices in a "network" object
  n_nodes <- network::network.size(n)
  n_edges <- network::network.edgecount(n)
  
  m <- igraph::rewire(
    intergraph::asIgraph(n),
    keeping_degseq(niter = n_edges*20)
    )

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Copying the vertex attributes of n to m
  for (a in network::list.vertex.attributes(n)) {
    network::set.vertex.attribute(m, a, network::get.vertex.attribute(n, a))
  }

  saveRDS(m, sprintf("data/graphs/%04i-degseq.rds", i), compress = FALSE)

  NULL

}, mc.cores = ncores)
message("Done degree-sequence preserve rewiring of the networks")

message("Done saving the networks")
