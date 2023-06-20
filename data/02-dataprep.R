#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --ntasks=51
#SBATCH --mem=128GB
#SBATCH --job-name=02-dataprep
#SBATCH --time=24:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=george.vegayon@utah.edu

library(igraph)
library(network)
# library(netplot)
library(intergraph)
library(ergm)
library(data.table)

ncores <- 50

# Read the Simulated data rds file from data/
message("Loading the networks")
networks <- readRDS("data/Simulated_1000_networks.rds")
networks <- unclass(networks)

# Reading in the edgelists sf, sw, r
networks_sf <- readRDS("data/Simulated_1000_networks_sf.rds")
networks_sw <- readRDS("data/Simulated_1000_networks_sw.rds")
networks_r <- readRDS("data/Simulated_1000_networks_r.rds")
message("Done loading the networks")

# Converting the networks to network objects
message("Converting the networks to network objects")
i2net <- function(net) {
  igraph::as_edgelist(net$net)
}

networks_sf <- parallel::mclapply(networks_sf, i2net, mc.cores = ncores)
networks_sw <- parallel::mclapply(networks_sw, i2net, mc.cores = ncores)
networks_r <- parallel::mclapply(networks_r, i2net, mc.cores = ncores)
message("Done converting the networks to network objects")


# Combining the networks
networks <- c(networks, networks_sf, networks_sw, networks_r)

# Removing the sf, sw, r networks (to save memory)
rm(networks_sf, networks_sw, networks_r)
gc()

# # Taking a sample of 100 networks
# set.seed(123)
# networks <- networks[sample(1:1000, 100)]

# Computing statistics using ERGM
message("Computing statistics using ERGM")
idxs <- 1:length(networks)
pos  <- rep(1:(length(networks)/4), 4)
nettypes <- rep(c("ergm", "sf", "sw", "degseq"), each = 1000)
vattrs <- as.data.frame(networks[[1]], unit = "vertices")
S_ergm <- parallel::mclapply(seq_along(idxs), \(i) {

  n <- networks[[i]]

  if (!inherits(n, "network")) {
    n <- network::network(n, vertex.attr = vattrs, undirected = TRUE)
  }

  summary_formula(
    n ~ edges + nodematch("gender") + nodematch("grade") + triangles + balance +
      twopath + gwdegree(decay = .25, fixed = TRUE) + isolates
  ) |> as.list() |> as.data.table()
}, mc.cores = ncores) |> rbindlist()

head(S_ergm)
message("Done computing statistics using ERGM")

# Computing statistics based on igraph
message("Computing statistics based on igraph")
S_igraph <- parallel::mclapply(networks, \(n) { 

  if (inherits(n, "network")) {
    n <- asIgraph(n)
  } else {
    n <- igraph::graph_from_edgelist(n, directed = FALSE)
  }
  
  data.table(
    modularity      = modularity(cluster_fast_greedy(n)),
    transitivity    = transitivity(n),
    density         = igraph::edge_density(n),
    diameter        = diameter(n),
    avg_path_length = igraph::mean_distance(n),
    avg_degree      = mean(degree(n)),
    avg_betweenness = mean(betweenness(n)),
    avg_closeness   = mean(closeness(n)),
    avg_eigenvector = mean(eigen_centrality(n)$vector),
    components      = components(n)$no
  )
}, mc.cores = ncores) |> rbindlist()

head(S_igraph)

setnames(S_ergm, new = paste0("ergm_", names(S_ergm)))
setnames(S_igraph, new = paste0("igraph_", names(S_igraph)))
message("Done computing statistics based on igraph")

# Combining the datasets
S <- cbind(S_igraph, S_ergm)
S[, netid := seq_len(.N)]
S[, nettype := nettypes]

# Reading simulation results ---------------------------------------------------
message("Processing the simulation results")
simres <- readRDS("data/01-abm-simulation.rds")

simres <- lapply(simres, \(x) {

  if (inherits(x, "error"))
    return(NULL)
  
  res <- tryCatch({
    history <- data.table(x$history)

    # Index of the max prevalence
    peak_idx    <- which.max(x$incidence$Exposed)
    peak_preval <- x$incidence$Exposed[peak_idx]
    peak_time   <- as.integer(rownames(x$incidence)[peak_idx])
    
    # Rt in the peak
    rt_idx     <- with(x$repnum, which.min(abs(peak_time - date)))
    rt         <- x$repnum$avg[rt_idx]
    dispersion <- 1/x$repnum$sd[rt_idx]^2

    # Mean Rt
    r_mean     <- with(x$repnum, sum(avg * n, na.rm = TRUE)/
      sum(n, na.rm = TRUE))

    gentime <- with(x$gentime, sum(avg * n, na.rm = TRUE)/
      sum(n, na.rm = TRUE))

    # Final prevalence
    final_preval <- with(
      x$history, tail(counts[state == "Removed"], 1)
      )

    # Return a data.table
    data.table(
      simid          = x$simid,
      netid          = x$netid,
      peak_time      = peak_time,
      peak_preval    = peak_preval,
      rt             = rt,
      rt_mean        = r_mean,
      dispersion     = dispersion,
      gentime        = gentime,
      final_preval   = final_preval,
      infectiousness = x$param$infectiousness,
      inc_days       = x$param$inc_days,
      recovery_rate  = x$param$recovery_rate
    )
  }, error=function(e) e)

  if (inherits(res, "error"))
    return(NULL)

  res

}) |> rbindlist()
message("Done processing the simulation results")

print(head(S))

print(head(simres))

# Merge the datasets
S <- merge(S, simres, by = "netid", all.x = TRUE)

fwrite(
  S,
  file = "data/02-dataprep-network-stats.csv.gz"
)

message("Done!")
