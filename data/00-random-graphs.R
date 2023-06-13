#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --job-name=00-rand-graphs
#SBATCH --output=00-rand-graphs-%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=george.vegayon@utah.edu
#SBATCH --ntasks=32
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

# Converting the networks to igraph objects
message("Converting the networks to igraph objects")
networks_i <- parallel::mclapply(networks, asIgraph, mc.cores = ncores)
message("Done converting the networks to igraph objects")

# Simulating scalefree networks with same density and size
# as the original networks
message("Simulating scalefree networks")
networks_sf <- parallel::mclapply(networks_i, \(n) {
  m <- igraph::sample_pa(igraph::vcount(n), power = 1,
  m = rgamma(1, ecount(n)/vcount(n) * 10, 10))
  list(net = m, n = vcount(m))
}, mc.cores = ncores)
message("Done simulating scalefree networks")

# Simulating smallworkd networks with same density and size
# as the original networks
message("Simulating smallworld networks")
networks_sw <- parallel::mclapply(networks_i, \(n) {
  m <- igraph::sample_smallworld(
    igraph::vcount(n), dim = 1, nei = rgamma(1, ecount(n)/vcount(n) * 10, 10),
    p = 0.2)
    list(net = m, n = vcount(m))
}, mc.cores = ncores)
message("Done simulating smallworld networks")

# Degree-sequence preserve rewiring of the networks
message("Degree-sequence preserve rewiring of the networks")
networks_r <- parallel::mclapply(networks_i, \(n) {
  m <- igraph::rewire(n, keeping_degseq(niter = ecount(n)*20))
  list(net = m, n = vcount(m))
}, mc.cores = ncores)
message("Done degree-sequence preserve rewiring of the networks")

# Saving the networks as RDS files
message("Saving the networks")
saveRDS(networks_sf, "data/Simulated_1000_networks_sf.rds")
saveRDS(networks_sw, "data/Simulated_1000_networks_sw.rds")
saveRDS(networks_r, "data/Simulated_1000_networks_r.rds")
message("Done saving the networks")
