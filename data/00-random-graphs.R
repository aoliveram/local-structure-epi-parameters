#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --job-name=00-random-graphs
#SBATCH --output=00-rand-graphs-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=george.vegayon@utah.edu
#SBATCH --mem=20G
#SBATCH --time=4:00:00

library(epiworldR)
library(slurmR)
library(intergraph)
library(igraph)

# Set the seed for reproducibility
set.seed(1231)
ncores <- 30L

# Set the slurmR options
SB_OPTS <- list(
  account    = "vegayon-np",
  partition  = "vegayon-shared-np",
  mem        = "8G"
  )


# Load the simulated networks from the RDS file
message("Loading the networks")

message("Done loading the networks")

# Saving the networks as individual graphs
if (!file.exists(sprintf("data/graphs/%04i-ergm.rds", 1))) {

  networks <- readRDS("data/Simulated_1000_networks.rds")

  for (i in seq_along(networks)) {

    fn <- sprintf("data/graphs/%04i-ergm.rds", i)

    if (file.exists(fn))
      next

    saveRDS(networks[[i]], fn, compress = FALSE)

  }

}

net1 <- readRDS("data/graphs/0001-ergm.rds")
networks <- list.files("data/graphs", pattern = "[0-9]+-ergm\\.rds$", full.names = TRUE)

# All networks have the same attributes. So we get them once
v_attrs <- as.data.frame(net1, unit = "vertices") |>
  as.list()

n_nodes <- network::network.size(net1)

# Simulating scalefree networks with same density and size
# as the original networks
message("Simulating scalefree networks")
random_seeds <- sample.int(1e6, length(networks), replace = FALSE)
networks_sf <- Slurm_lapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-sf.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- readRDS(sprintf("data/graphs/%04i-ergm.rds", i))

  # Number of vertices in a "network" object
  n_edges <- network::network.edgecount(n)

  m <- igraph::sample_pa(
    n_nodes, power = 1,
    m = floor(n_edges/n_nodes),
    directed = FALSE
  )

  for (a in names(v_attrs))
    m <- set_vertex_attr(m, name = a, value = v_attrs[[a]])

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Writing the network to a file 
  saveRDS(m, fn, compress = FALSE)

  NULL

}, sbatch_opt = SB_OPTS, njobs = ncores,
  export = c("v_attrs", "n_nodes", "random_seeds")
  )
message("Done simulating scalefree networks")

# Simulating smallworkd networks with same density and size
# as the original networks
message("Simulating smallworld networks p=0.1")
networks_sw <- Slurm_lapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-swp01.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- readRDS(sprintf("data/graphs/%04i-ergm.rds", i))

  # Number of vertices in a "network" object
  n_edges <- network::network.edgecount(n)

  nneigh <- ceiling(mean(sna::degree(n, gmode="graph", cmode = "indegree")))/2

  # Adding a bit noise
  nneigh <- nneigh + ceiling(runif(1, -.1, .1) * nneigh)
  m <- igraph::sample_smallworld(
    n_nodes,
    dim = 1,
    nei = nneigh,
    p   = 0.1
    )

  for (a in names(v_attrs))
    m <- set_vertex_attr(m, name = a, value = v_attrs[[a]])

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Writing the network to a file 
  saveRDS(m, fn, compress = FALSE)

  NULL
    
}, njobs = ncores, sbatch_opt = SB_OPTS,
  export = c("v_attrs", "n_nodes", "random_seeds")
)
message("Done simulating smallworld p=0.1 networks")

# Small-world with a larger number of rewires ----------------------------------

message("Simulating smallworld networks p=0.2")
networks_sw <- Slurm_lapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-swp02.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- readRDS(sprintf("data/graphs/%04i-ergm.rds", i))

  # Number of vertices in a "network" object
  n_edges <- network::network.edgecount(n)

  nneigh <- ceiling(mean(sna::degree(n, gmode="graph", cmode = "indegree")))/2

  # Adding a bit noise
  nneigh <- nneigh + ceiling(runif(1, -.1, .1) * nneigh)
  m <- igraph::sample_smallworld(
    n_nodes,
    dim = 1,
    nei = nneigh,
    p   = 0.2
    )

  for (a in names(v_attrs))
    m <- set_vertex_attr(m, name = a, value = v_attrs[[a]])

  # Turning m into a network object
  m <- intergraph::asNetwork(m)

  # Writing the network to a file 
  saveRDS(m, fn, compress = FALSE)

  NULL
    
}, njobs = ncores, sbatch_opt = SB_OPTS,
  export = c("v_attrs", "n_nodes", "random_seeds")
  )
  
message("Done simulating smallworld networks p0.2")

# Degree-sequence preserve rewiring of the networks
message("Degree-sequence preserve rewiring of the networks")
networks_r <- Slurm_lapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-degseq.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- readRDS(sprintf("data/graphs/%04i-ergm.rds", i))


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

}, njobs = ncores, sbatch_opt = SB_OPTS,
  export = c("v_attrs", "n_nodes", "random_seeds")
  )
message("Done degree-sequence preserve rewiring of the networks")

# Erdos-Renyi networks with same density and size
# as the original networks
message("Simulating Erdos-Renyi networks")
networks_er <- Slurm_lapply(seq_along(networks), \(i) {

  # Checking out if the file already exists
  fn <- sprintf("data/graphs/%04i-er.rds", i)

  if (file.exists(fn))
    return(NULL)
  
  set.seed(random_seeds[i])

  n <- readRDS(sprintf("data/graphs/%04i-ergm.rds", i))

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

}, njobs = ncores, sbatch_opt = SB_OPTS,
  export = c("v_attrs", "n_nodes", "random_seeds")
)
message("Done simulating Erdos-Renyi the networks")

message("Done saving the networks")
