#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --job-name=00-rand-graphs
#SBATCH --output=00-rand-graphs-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=george.vegayon@utah.edu
#SBATCH --ntasks=21
#SBATCH --mem=64G
#SBATCH --time=4:00:00

library(epiworldR)
library(slurmR)
library(intergraph)
library(igraph)

# Set the seed for reproducibility
set.seed(1231)
ncores <- 20L

# Load the simulated networks from the RDS file
message("Loading the networks")
networks <- readRDS("data/Simulated_1000_networks.rds")
message("Done loading the networks")

# Saving the networks as individual graphs
for (i in seq_along(networks)) {

  fn <- sprintf("data/graphs/%04i-ergm.rds", i)

  if (file.exists(fn))
    next

  saveRDS(networks[[i]], fn, compress = FALSE)

}

# All networks have the same attributes. So we get them once
v_attrs <- as.data.frame(networks[[1]], unit = "vertices") |>
  as.list()

n_nodes <- network::network.size(networks[[1]])

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
  n_edges <- network::network.edgecount(n)

  m <- igraph::sample_pa(
    n_nodes, power = 1,
    m = n_edges/n_nodes,
    directed = FALSE
  )

  for (a in names(v_attrs))
    m <- set_vertex_attr(m, name = a, value = v_attrs[[a]])

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Writing the network to a file 
  saveRDS(m, fn, compress = FALSE)

  NULL

}, mc.cores = ncores)
message("Done simulating scalefree networks")

# Simulating smallworkd networks with same density and size
# as the original networks
message("Simulating smallworld networks")
networks_sw <- parallel::mclapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-sw.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- networks[[i]]

  # Number of vertices in a "network" object
  n_edges <- network::network.edgecount(n)

  m <- igraph::sample_smallworld(
    n_nodes,
    dim = 1,
    nei = ceiling(mean(sna::degree(n, gmode="graph", cmode = "indegree")))/2,
    p   = 0.1
    )

  for (a in names(v_attrs))
    m <- set_vertex_attr(m, name = a, value = v_attrs[[a]])

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Writing the network to a file 
  saveRDS(m, fn, compress = FALSE)

  NULL
    
}, mc.cores = ncores)
message("Done simulating smallworld networks")

# Degree-sequence preserve rewiring of the networks
message("Degree-sequence preserve rewiring of the networks")
networks_r <- parallel::mclapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-degseq.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- networks[[i]]


  # Number of vertices in a "network" object
  n_edges <- network::network.edgecount(n)
  
  m <- igraph::rewire(
    intergraph::asIgraph(n),
    keeping_degseq(niter = n_edges*20)
    )

  for (a in names(v_attrs))
    m <- set_vertex_attr(m, name = a, value = v_attrs[[a]])

  # Turning m into a network object
  m <- igraph::delete_edge_attr(m, "na") |>
    intergraph::asNetwork()

  # Writing the network to a file 
  saveRDS(m, fn, compress = FALSE)

  NULL

}, mc.cores = ncores)
message("Done degree-sequence preserve rewiring of the networks")

# Erdos-Renyi networks with same density and size
# as the original networks
message("Simulating Erdos-Renyi networks")
networks_er <- parallel::mclapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-er.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- networks[[i]]

  # Number of vertices in a "network" object
  n_edges <- network::network.edgecount(n)

  m <- igraph::sample_gnm(
    n = n_nodes, 
    m = n_edges,
    directed = FALSE,
    loops    = FALSE
  )

  for (a in names(v_attrs))
    m <- set_vertex_attr(m, name = a, value = v_attrs[[a]])

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Writing the network to a file 
  saveRDS(m, fn, compress = FALSE)

  NULL

}, mc.cores = ncores)
message("Done simulating Erdos-Renyi the networks")

message("Done saving the networks")
