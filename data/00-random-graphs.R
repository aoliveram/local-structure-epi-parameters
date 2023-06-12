#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --job-name=00-rand-graphs
#SBATCH --output=00-rand-graphs-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=george.vegayon@utah.edu
#SBATCH --ntasks=32
#SBATCH --mem=32G

library(epiworldR)
library(slurmR)
library(intergraph)
library(igraph)

# Set the seed for reproducibility
set.seed(1231)
ncores <- 30

# Load the simulated networks from the RDS file
networks <- readRDS("data/Simulated_1000_networks.rds")

# Converting the networks to igraph objects
networks_i <- parallel::mclapply(networks, asIgraph, mc.cores = ncores)

# Simulating scalefree networks with same density and size
# as the original networks
networks_sf <- parallel::mclapply(networks_i, \(n) {
  igraph::sample_pa(igraph::vcount(n), power = 1, m = igraph::ecount(n))
}, mc.cores = ncores)

# Simulating smallworkd networks with same density and size
# as the original networks
networks_sw <- parallel::mclapply(networks_i, \(n) {
  igraph::sample_smallworld(igraph::vcount(n), dim = 1, nei = 2, p = 0.5)
}, mc.cores = ncores)

# Simulating random networks with same density and size
# as the original networks
networks_r <- parallel::mclapply(networks_i, \(n) {
  igraph::sample_gnp(igraph::vcount(n), igraph::edge_density(n))
}, mc.cores = ncores)

# Saving the networks as edgelists
networks_sf <- parallel::mclapply(networks_sf, igraph::as_edgelist, mc.cores = ncores)
networks_sw <- parallel::mclapply(networks_sw, igraph::as_edgelist, mc.cores = ncores)
networks_r <- parallel::mclapply(networks_r, igraph::as_edgelist, mc.cores = ncores)

# Saving the networks as RDS files
saveRDS(networks_sf, "data/Simulated_1000_networks_sf.rds")
saveRDS(networks_sw, "data/Simulated_1000_networks_sw.rds")
saveRDS(networks_r, "data/Simulated_1000_networks_r.rds")

