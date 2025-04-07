#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --job-name=00-random-graphs
#SBATCH --output=00-rand-graphs-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=an.oliveram@udd.cl
#SBATCH --mem=20G
#SBATCH --time=11:00:00

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

# Saving the networks as individual graphs
if (!file.exists(sprintf("data/graphs/%04i-ergm.rds", 1))) {  # Verifica si el primer archivo de red ya existe
  
  networks <- readRDS("data/Simulated_1000_networks.rds")  # Carga las redes simuladas desde un archivo RDS
  
  for (i in seq_along(networks)) {  # Itera sobre cada red en la lista 'networks'
    
    fn <- sprintf("data/graphs/%04i-ergm.rds", i)  # si i=3, fn será data/graphs/0003-ergm.rds
    
    if (file.exists(fn))  # Verifica si el archivo ya existe
      next  # Si el archivo ya existe, pasa al siguiente
    
    saveRDS(networks[[i]], fn, compress = FALSE)  # guarda los datos contenidos en networks[[i]] en el archivo cuyo nombre está almacenado en fn.
    
  }
  
}

networks <- list.files("data/graphs", pattern = "[0-9]+-ergm\\.rds$", full.names = TRUE) #0001-ergm.rds -> coincide, abc-ergm.rds -> no coincide

message("Done loading the networks")

net1 <- readRDS("data/graphs/0001-ergm.rds")

class(net1)
library(network)
network.density(net1)
is.directed(net1)
get.edge.attribute(net1, "weight")
get.edge.attribute(net1, "sign")
list.edge.attributes(net1)
head(as.data.frame(net1, unit = "vertices"))

# All networks have the same attributes. So we get them once
v_attrs <- as.data.frame(net1, unit = "vertices") |>    
  as.list()  #atributos asociados a cada nodo

n_nodes <- network::network.size(net1)

sum(v_attrs$gender == 1)/n_nodes

# Simulating scalefree networks with same density and size --------------------
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

# Simulating smallworkd networks with same density and size -------------------
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

# Degree-sequence preserve rewiring of the networks --------------------------
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

# Erdos-Renyi networks with same density and size -----------------------------
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
