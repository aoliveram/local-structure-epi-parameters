#!/bin/sh
#SBATCH --account=vegayon-np
#SBATCH --partition=vegayon-shared-np
#SBATCH --ntasks=31
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

ncores <- 30

# Read the Simulated data rds file from data/
message("Loading the networks files")
networks <- list.files(
  "data/graphs", pattern = "[0-9]+-(ergm|sf|sw|degseq)\\.rds$",full.names = TRUE) 

# Computing statistics using ERGM
message("Computing statistics using ERGM")
S_ergm <- parallel::mclapply(seq_along(networks), \(i) {

  n <- networks[[i]]

  # Reading the network
  n <- readRDS(n)

  res <- summary_formula(
    n ~ edges + nodematch("gender") + nodematch("grade") + triangles + balance +
      twopath + gwdegree(decay = .25, fixed = TRUE) + isolates
  ) |> as.list() |> as.data.table()

  res[, netfile := networks[[i]]]

  res
}, mc.cores = ncores) |> rbindlist()

# Pre-appending `ergm_` to all column names, except netfile
setnames(S_ergm, new = paste0("ergm_", names(S_ergm)))
setnames(S_ergm, old = "ergm_netfile", new = "netfile")


head(S_ergm)
message("Done computing statistics using ERGM")

# Computing statistics based on igraph
message("Computing statistics based on igraph")
S_igraph <- parallel::mclapply(networks, \(i) { 

  # Reading the network and converting it into igraph
  n <- readRDS(i)
  n <- intergraph::asIgraph(n)
  
  data.table(
    netfile         = i,
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

setnames(S_igraph, new = paste0("igraph_", names(S_igraph)))
setnames(S_igraph, old = "igraph_netfile", new = "netfile")

message("Done computing statistics based on igraph")

# Combining the datasets
S <- merge(S_igraph, S_ergm, by = "netfile", all = TRUE)

S[, nettype := gsub(".+-([a-z]+)\\.rds", "\\1", netfile)]
S[, netid   := gsub(".+/([0-9]+-[a-z]+)\\.rds", "\\1", netfile)]


# Reading simulation results ---------------------------------------------------
message("Processing the simulation results")
simfiles <- list.files("data/graphs", pattern = "-sim-[0-9]+\\.rds$", full.names = TRUE)

simres <- parallel::mclapply(simfiles, \(fn) {

  x <- readRDS(fn)
  
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
      simfile           = fn,
      simid             = x$simid,
      peak_time         = peak_time,
      peak_preval       = peak_preval,
      rt                = rt,
      rt_mean           = r_mean,
      dispersion        = dispersion,
      gentime           = gentime,
      final_preval      = final_preval,
      transmission_rate = x$param$transmission_rate,
      inc_days          = x$param$inc_days,
      recovery_rate     = x$param$recovery_rate
    )
  }, error=function(e) e)

  if (inherits(res, "error"))
    return(NULL)

  res

}, mc.cores = ncores) |> rbindlist()
message("Done processing the simulation results")

simres[, netid := gsub(".+/([0-9]+-[a-z]+)-sim.+", "\\1", simfile)]

print(head(S))

print(head(simres))

# Merge the datasets
S <- merge(simres, S, by = "netid", all.x = TRUE, all.y = FALSE)

fwrite(
  S,
  file = "data/02-dataprep-network-stats.csv.gz"
)

message("Done!")
